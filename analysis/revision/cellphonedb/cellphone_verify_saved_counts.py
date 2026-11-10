import argparse,csv,json
from pathlib import Path
from cellphone_count_data import count_table_path
import numpy as np
BASE=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description='Check recovered counts against available original significant_means tables.')
parser.add_argument('--output-dir',type=Path,default=BASE/'output'/'10_cellphonedb')
parser.add_argument('--project-data-dir',type=Path,default=Path(r'K:\model paper\v2'))
args=parser.parse_args()
OUT=args.output_dir
verification=[]
for key in ['gut_3','liver_4']:
    table=args.project_data_dir/'results'/key/f'statistical_analysis_significant_means__{key}.txt'
    with table.open(newline='',encoding='utf-8-sig') as f:
        reader=csv.DictReader(f,delimiter='\t'); cols=[c for c in reader.fieldnames if '|' in c and c.split('|')[0]!=c.split('|')[1]]
        vals=np.array([[float(row[c]) if row[c] else 0 for c in cols] for row in reader],dtype=float)
    vals=np.nan_to_num(vals,nan=0.0)
    transformed=np.log1p(vals); cutoff=float(np.percentile(transformed[transformed>0],5))
    counts=(transformed>cutoff).sum(axis=0)
    with count_table_path(OUT,key).open(newline='',encoding='utf-8') as f: saved=list(csv.reader(f))
    labels=saved[0][1:]; index={s:i for i,s in enumerate(labels)}; directed=np.zeros((15,15),dtype=int)
    for col,value in zip(cols,counts):
        a,b=col.split('|'); directed[index[a],index[b]]=value
    rebuilt=directed+directed.T; stored=np.array([[int(v) for v in row[1:]] for row in saved[1:]])
    result={'dataset':key,'significant_means_path':str(table),'code_rule':'count entries > 5th percentile of positive log1p significant means, excluding self-pairs; add reverse ordered counts','cutoff_log1p':cutoff,'complete_saved_table_match':bool(np.array_equal(rebuilt,stored)),'different_cells':int(np.count_nonzero(rebuilt!=stored)),'rebuilt_max':int(rebuilt.max()),'stored_max':int(stored.max())}
    verification.append(result); print(json.dumps(result))
audit=OUT/'audit'; audit.mkdir(parents=True,exist_ok=True)
(audit/'interaction_count_checks_against_original_results.json').write_text(json.dumps(verification,indent=2),encoding='utf-8')
if not all(result['complete_saved_table_match'] for result in verification):
    raise ValueError('Recovered counts differ from the available original results; see audit/.')
