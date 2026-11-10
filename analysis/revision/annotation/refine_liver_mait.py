"""Replay and inspect liver notebook cells 330-336 on a saved neighbor graph.

The input H5AD is opened read-only. By default only the cluster/marker record
is written in audit/. A proposed barcode-to-label table is opt-in and also
stays in audit/: historical cluster numbers do not establish cell identity on
a rerun. This script does not update the publication figures or source H5AD.
No expression matrix is normalized and no neighbor graph is rebuilt.

Clustering calls use the original Scanpy API, not a reimplementation:
https://scanpy.readthedocs.io/en/1.9.x/generated/scanpy.tl.leiden.html
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import sys
import tempfile


GENERIC_MAIT = "TCRab CD8 MAIT"
EFFECTOR = "TCRab CD8 effector"
MAIT_AA = "TCRab CD8aa MAIT"
MAIT_AB = "TCRab CD8ab MAIT"
LABEL_FILENAME = "liver_proposed_cell_labels_after_MAIT_refinement.csv"
RECORDED_VERSIONS = {"scanpy": "1.9.3", "igraph": "0.10.4", "leidenalg": "0.9.0"}


def parse_args():
    repo = Path(__file__).resolve().parents[3]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=repo / "config/paths.local.json")
    parser.add_argument("--input", type=Path, help="Liver H5AD before the documented MAIT refinement.")
    parser.add_argument("--output", type=Path, help="Analysis folder; defaults to OUTPUT_DIR/09_annotation.")
    parser.add_argument("--seed", type=int, default=0, help="The notebook uses Scanpy's default seed 0.")
    parser.add_argument("--graph-key", default=None, help="Override the saved neighbors connectivities key.")
    parser.add_argument("--allow-version-mismatch", action="store_true",
                        help="Explicitly allow graph-clustering versions different from the notebook.")
    parser.add_argument("--stage1-cd8ab-cluster", type=int, default=0,
                        help="Effector subcluster provisionally mapped to CD8ab MAIT; notebook value: 0.")
    parser.add_argument("--stage2-effector-cluster", type=int, default=1,
                        help="Candidate subcluster returned to effector; notebook value: 1.")
    parser.add_argument("--export-proposed-labels", action="store_true",
                        help="Also export proposed labels to audit/ for review; never applied automatically.")
    args = parser.parse_args()
    if args.stage1_cd8ab_cluster < 0 or args.stage2_effector_cluster < 0:
        parser.error("Cluster numbers must be nonnegative.")
    return args


def decode(value):
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def read_column(node, h5py, np):
    if isinstance(node, h5py.Group):
        categories = [decode(value) for value in node["categories"][()]]
        return np.asarray([categories[int(code)] if code >= 0 else "" for code in node["codes"][()]], dtype=object)
    return np.asarray([decode(value) for value in node[()]], dtype=object)


def array_hash(*arrays):
    digest = hashlib.sha256()
    for values in arrays:
        digest.update(str(values.shape).encode("ascii"))
        digest.update(str(values.dtype).encode("ascii"))
        digest.update(values.tobytes())
    return digest.hexdigest()


def string_hash(values):
    return hashlib.sha256("\n".join(values).encode("utf-8")).hexdigest()


def read_markers(handle, selected, h5py, np):
    """Read only four marker columns for the affected cells from stored logcounts."""
    var = handle["var"]
    index = decode(var.attrs.get("_index", "_index"))
    genes = read_column(var[index], h5py, np).tolist()
    names = ["TRAV1-2", "CD4", "CD8A", "CD8B"]
    columns = {genes.index(name): i for i, name in enumerate(names) if name in genes}
    if len(columns) != len(names):
        raise ValueError("The saved expression layer is missing a required MAIT marker.")
    layer = handle["layers/logcounts"]
    if decode(layer.attrs.get("encoding-type", "")) != "csr_matrix":
        raise ValueError("Expected the original sparse CSR layers/logcounts.")
    ptr = layer["indptr"][()]
    values = np.zeros((len(selected), len(names)), dtype=float)
    for position, row in enumerate(selected):
        start, end = int(ptr[row]), int(ptr[row + 1])
        gene_ids = layer["indices"][start:end]
        expression = layer["data"][start:end]
        for gene_id, value in zip(gene_ids, expression):
            column = columns.get(int(gene_id))
            if column is not None:
                values[position, column] = float(value)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("The requested stored marker values must be finite and nonnegative.")
    return names, values


def write_outputs(folder, rows, audit):
    checks = folder / "audit"
    checks.mkdir(parents=True, exist_ok=True)
    # Temporary artifacts are confined to a verified task-specific temp directory.
    temp_root = folder / "temp"
    temp_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="mait_refinement_", dir=temp_root) as temporary:
        temporary = Path(temporary).resolve()
        if temporary.parent != temp_root.resolve():
            raise ValueError("Unexpected temporary directory location.")
        if rows is not None:
            labels = temporary / LABEL_FILENAME
            with labels.open("w", encoding="utf-8-sig", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
        record = temporary / "liver_MAIT_refinement.json"
        record.write_text(json.dumps(audit, indent=2), encoding="utf-8")
        if rows is not None:
            os.replace(labels, checks / LABEL_FILENAME)
        else:
            (checks / LABEL_FILENAME).unlink(missing_ok=True)
        os.replace(record, checks / record.name)
    try:
        temp_root.rmdir()
    except OSError:
        pass  # Other temporary work, if present, belongs to its owning task.


def main():
    args = parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8-sig"))
    source = (args.input or Path(config["liver_h5ad"])).expanduser().resolve()
    folder = (args.output or Path(config["output_dir"]) / "09_annotation").expanduser().resolve()
    versions = {name: metadata.version(name) for name in RECORDED_VERSIONS}
    mismatch = {name: versions[name] for name in versions if versions[name] != RECORDED_VERSIONS[name]}
    if mismatch and not args.allow_version_mismatch:
        raise RuntimeError("Use requirements_mait_refinement.txt in the isolated Python 3.9 runtime; "
                           "graph-clustering versions differ: " + json.dumps(mismatch))

    import h5py
    import numpy as np
    import pandas as pd
    from scipy import sparse
    import anndata as ad
    import leidenalg
    import scanpy as sc

    before_stat = source.stat()
    with h5py.File(source, "r") as handle:
        obs = handle["obs"]
        index = decode(obs.attrs.get("_index", "_index"))
        barcodes = read_column(obs[index], h5py, np)
        original = read_column(obs["celltype_highres"], h5py, np)
        batch = read_column(obs["batch"], h5py, np)
        donor = read_column(obs["Donor ID"], h5py, np)
        if len(set(barcodes)) != len(barcodes):
            raise ValueError("Duplicate cell barcodes are not supported.")
        if not (original == GENERIC_MAIT).any() or not (original == EFFECTOR).any():
            raise ValueError("The input must have the earlier CD8 MAIT and CD8 effector labels.")
        if (original == MAIT_AA).any() or (original == MAIT_AB).any():
            raise ValueError("The input already contains refined MAIT labels; do not replay onto refined input.")
        graph_key = args.graph_key
        if graph_key is None:
            neighbors = handle["uns/neighbors"]
            graph_key = decode(neighbors["connectivities_key"][()]) if "connectivities_key" in neighbors else "connectivities"
        saved = handle["obsp"][graph_key]
        if decode(saved.attrs["encoding-type"]) != "csr_matrix":
            raise ValueError("Expected a CSR saved neighbor graph.")
        data, indices, indptr = saved["data"][()], saved["indices"][()], saved["indptr"][()]
        shape = tuple(int(value) for value in saved.attrs["shape"])
        graph = sparse.csr_matrix((data, indices, indptr), shape=shape)
        if shape != (len(barcodes), len(barcodes)) or not np.isfinite(data).all() or (data < 0).any():
            raise ValueError("Invalid saved neighbor graph.")
        affected = np.flatnonzero(np.isin(original, [GENERIC_MAIT, EFFECTOR]))
        marker_names, marker_values = read_markers(handle, affected, h5py, np)

    # A graph-only AnnData avoids loading the multi-GB expression matrices/uns arrays.
    graph_data = ad.AnnData(X=sparse.csr_matrix((len(barcodes), 0), dtype=np.float32),
                           obs=pd.DataFrame({"celltype_highres": pd.Categorical(original)}, index=barcodes))
    options = dict(random_state=args.seed, directed=True, use_weights=True,
                   n_iterations=-1, partition_type=leidenalg.RBConfigurationVertexPartition)
    sc.tl.leiden(graph_data, resolution=1.0, restrict_to=("celltype_highres", [EFFECTOR]),
                 adjacency=graph, key_added="mait_stage1", **options)
    stage1 = graph_data.obs["mait_stage1"].astype(str).to_numpy()
    stage1_selection = EFFECTOR + "," + str(args.stage1_cd8ab_cluster)
    stage2_selection = MAIT_AB + "," + str(args.stage2_effector_cluster)
    if not (stage1 == stage1_selection).any():
        raise ValueError("The selected stage-1 cluster is absent: " + stage1_selection)
    refined = original.copy()
    refined[stage1 == stage1_selection] = MAIT_AB
    refined[original == GENERIC_MAIT] = MAIT_AA
    graph_data.obs["celltype_highres"] = pd.Categorical(refined)
    sc.tl.leiden(graph_data, resolution=0.3, restrict_to=("celltype_highres", [MAIT_AB]),
                 adjacency=graph, key_added="mait_stage2", **options)
    stage2 = graph_data.obs["mait_stage2"].astype(str).to_numpy()
    # The notebook's conditional assignment is a no-op if no cluster 1 exists.
    refined[stage2 == stage2_selection] = EFFECTOR
    untouched = ~np.isin(original, [GENERIC_MAIT, EFFECTOR])
    if not np.array_equal(original[untouched], refined[untouched]):
        raise AssertionError("A non-MAIT/non-effector cell was changed.")
    if not (refined == MAIT_AB).any() or not np.all(refined[original == GENERIC_MAIT] == MAIT_AA):
        raise AssertionError("The proposed MAIT mapping was not preserved.")
    after_stat = source.stat()
    if (before_stat.st_size, before_stat.st_mtime_ns) != (after_stat.st_size, after_stat.st_mtime_ns):
        raise RuntimeError("The source file changed while it was being read.")

    marker_summary = []
    for stage, labels in [("stage1", stage1), ("stage2", stage2), ("final", refined)]:
        local = labels[affected]
        for group in sorted(set(local)):
            take = local == group
            marker_summary.append({"stage": stage, "group": group, "n_cells": int(take.sum()),
                                   "mean_stored_log_expression": dict(zip(marker_names, marker_values[take].mean(axis=0).tolist())),
                                   "fraction_expressing": dict(zip(marker_names, (marker_values[take] > 0).mean(axis=0).tolist()))})
    counts = [{"batch": value, "CD8aa_MAIT": int(np.sum((batch == value) & (refined == MAIT_AA))),
               "CD8ab_MAIT": int(np.sum((batch == value) & (refined == MAIT_AB)))} for value in sorted(set(batch))]
    rows = None
    if args.export_proposed_labels:
        rows = [{"cell_barcode": str(barcodes[i]), "batch": str(batch[i]), "donor": str(donor[i]),
                 "original_cell_type": str(original[i]), "proposed_cell_type": str(refined[i])} for i in range(len(barcodes))]
    issues = []
    if not (stage2 == stage2_selection).any():
        issues.append("The selected stage-2 effector cluster is absent. The original notebook's saved plot has candidate clusters 0 and 1; the current replay does not reproduce that structure.")
    audit = {"source_h5ad": str(source), "source_bytes": before_stat.st_size,
             "input_h5ad_modified": False, "source_notebook": "01_preprocess_annotate_liver.ipynb",
             "source_cells_one_based": [330, 331, 333, 334, 335, 336],
             "graph_key": graph_key, "graph_sha256": array_hash(data, indices, indptr),
             "cell_barcodes_sha256": string_hash(barcodes), "original_labels_sha256": string_hash(original),
             "parameters": {"resolutions": [1.0, 0.3], "seed": args.seed, "directed": True,
                            "use_weights": True, "n_iterations": -1, "partition": "RBConfigurationVertexPartition"},
             "proposed_mapping": {"stage1_CD8ab_MAIT": stage1_selection,
                                  "stage2_return_to_effector": stage2_selection,
                                  "original_CD8_MAIT": MAIT_AA},
             "versions": {name: metadata.version(name) for name in ["scanpy", "anndata", "igraph", "leidenalg", "numpy", "scipy", "pandas", "h5py"]},
             "python": sys.version, "changed_cells": int(np.sum(original != refined)),
             "stage1_effector_cluster_counts": dict(Counter(stage1[original == EFFECTOR])),
             "stage2_candidate_cluster_counts": dict(Counter(stage2[stage1 == stage1_selection])),
             "stage2_cells_returned_to_effector": int(np.sum(stage2 == stage2_selection)),
             "donor_counts": counts, "marker_summaries": marker_summary,
             "historical_cell_membership_verified": False,
             "validation_status": "Requires cluster-identity review; no labels applied to final results.",
             "replay_structure_differences": issues,
             "proposed_labels_file": "audit/" + LABEL_FILENAME if rows is not None else None,
             "interpretation": "Cluster numbers are specific to a clustering run. Inspect TRAV1-2, CD4, CD8A and CD8B together and compare to the notebook before accepting a mapping. No new marker-positivity threshold is imposed.",
             "result_kind": "Replay of notebook calls on the saved graph; proposed identities, not recovered historical cell assignments."}
    folder.mkdir(parents=True, exist_ok=True)
    write_outputs(folder, rows, audit)
    print(json.dumps({"audit": str(folder / "audit/liver_MAIT_refinement.json"),
                      "proposed_labels": audit["proposed_labels_file"],
                      "validation_status": audit["validation_status"],
                      "replay_structure_differences": issues,
                      "proposed_changed_cells": audit["changed_cells"], "proposed_donor_counts": counts}, indent=2))


if __name__ == "__main__":
    main()
