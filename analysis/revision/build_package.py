"""Build the reviewer package from configured, recovered inputs; never refit DE."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BASELINE = "c21587ad29c7ddb41d30266071435bf5da6f0ccf"

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/paths.local.json")
    parser.add_argument("--steps", default="nichenet,annotation,cellphonedb")
    parser.add_argument("--stage-only", action="store_true")
    parser.add_argument("--install-deps", action="store_true")
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8-sig"))
    steps = set(args.steps.split(","))
    if not steps <= {"nichenet", "annotation", "cellphonedb"}:
        parser.error("Unknown step; choose nichenet,annotation,cellphonedb.")
    def path(key, default=None):
        value = config.get(key) or default
        if value is None:
            raise ValueError(f"Set {key} in {args.config}")
        result = Path(value).expanduser()
        return result if result.is_absolute() else (ROOT / result).resolve()
    work = path("work_dir", ".local/revision")
    output = path("output_dir", "../gut-liver-revision-output")
    work.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    pylib = path("python_package_dir", str(work / "python-packages"))
    env["PYTHONPATH"] = str(pylib) + os.pathsep + env.get("PYTHONPATH", "")
    env["NICHE_R_LIBRARY"] = str(path("r_library_dir", str(work / "runtime/R-library")))
    env["NICHE_VENDOR"] = str(HERE / "nichenet/vendor")
    env["REVISION_OUTPUT_DIR"] = str(output)
    def run(command):
        print("Running", command[1], flush=True)
        subprocess.run([str(value) for value in command], check=True, cwd=ROOT, env=env)
    def stage(source, destination):
        if not source.is_file():
            raise FileNotFoundError(source)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if source.resolve() != destination.resolve():
            temporary = destination.with_name(destination.name + ".tmp")
            shutil.copyfile(source, temporary)
            os.replace(temporary, destination)
    if "nichenet" in steps:
        source = path("deg_dir")
        for receiver in ["gutT", "gutMac", "liverT", "liverMac"]:
            for suffix in ["_DEG.csv", ".gz"]:
                stage(source / (receiver + suffix), work / "data/original_deg" / (receiver + suffix))
        cached = config.get("cached_nichenet_dir")
        if cached:
            for source in path("cached_nichenet_dir").rglob("*.csv"):
                if source.name.endswith(("_nichenet.csv", "_receptors.csv")):
                    stage(source, work / "data/original_nichenet" / source.name)
        env["NICHE_MODEL"] = str(path("ligand_target_matrix"))
        env["NICHE_LR"] = str(path("ligand_receptor_network"))
        for name in ["validate_official_functions.R"]:
            stage(HERE / "nichenet" / name, work / "scripts/nichenet" / name)
        shutil.copytree(HERE / "nichenet/vendor", work / "work/nichenet_environment/nichenetr_v2.0.4", dirs_exist_ok=True)
    baseline_notebooks = work / "data/baseline_notebooks"
    if "cellphonedb" in steps:
        baseline_notebooks.mkdir(parents=True, exist_ok=True)
        for name in ["cellphonedb_gut.ipynb", "cellphonedb_liver.ipynb"]:
            data = subprocess.run(["git", "-C", str(ROOT), "show", f"{BASELINE}:{name}"], check=True, capture_output=True).stdout
            (baseline_notebooks / name).write_bytes(data)
    if args.stage_only:
        print("Inputs staged; no analysis was run.")
        return
    if args.install_deps:
        pylib.mkdir(parents=True, exist_ok=True)
        run([sys.executable, "-m", "pip", "install", "--target", pylib, "numpy", "pandas", "scipy", "matplotlib", "h5py", "Pillow"])
    if "nichenet" in steps:
        rscript = config.get("rscript") or "Rscript"
        if args.install_deps:
            run([rscript, HERE / "nichenet/prepare_environment.R", work])
        run([rscript, HERE / "nichenet/run_nichenet_revision.R", work])
    if "annotation" in steps:
        common = ["--root", work, "--gut", path("gut_h5ad"), "--liver", path("liver_h5ad"), "--output", output / "09_annotation"]
        run([sys.executable, HERE / "annotation/annotation_markers.py", *common])
        run([sys.executable, HERE / "annotation/annotation_texts.py", *common,
             "--gut-notebook", ROOT / "analysis/scrna/01_preprocess_annotate_gut.ipynb",
             "--liver-notebook", ROOT / "analysis/scrna/01_preprocess_annotate_liver.ipynb"])
    if "cellphonedb" in steps:
        common = ["--output-dir", output / "10_cellphonedb"]
        run([sys.executable, HERE / "cellphonedb/cellphone_recover_counts.py", *common,
             "--repo-dir", baseline_notebooks, "--project-data-dir", path("original_result_root")])
        run([sys.executable, HERE / "cellphonedb/cellphone_verify_saved_counts.py", *common,
             "--project-data-dir", path("original_result_root")])
        run([sys.executable, HERE / "cellphonedb/cellphone_plot_counts.py", *common])
    run([sys.executable, HERE / "assemble_outputs.py", "--output", output])
    print("Package saved to", output)

if __name__ == "__main__":
    main()
