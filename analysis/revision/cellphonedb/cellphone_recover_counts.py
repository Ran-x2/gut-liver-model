import argparse, ast, csv, hashlib, json
from html.parser import HTMLParser
from pathlib import Path
from cellphone_count_data import DISPLAY, PRESENCE_FILE, count_table_path

BASE=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description='Recover complete CellPhoneDB count tables from saved notebook HTML.')
parser.add_argument('--repo-dir',type=Path,default=BASE/'repo')
parser.add_argument('--output-dir',type=Path,default=BASE/'output'/'10_cellphonedb')
parser.add_argument('--project-data-dir',type=Path,default=Path(r'K:\model paper\v2'))
parser.add_argument('--export-source-html',action='store_true',help='Optionally retain the source HTML tables in audit/saved_notebook_count_tables/.')
args=parser.parse_args()
REPO=args.repo_dir
OUT=args.output_dir
DATA=OUT/'output'
DATA.mkdir(parents=True, exist_ok=True)
AUDIT=OUT/'audit'; AUDIT.mkdir(parents=True,exist_ok=True)
PLOT_DATA=OUT/'plot_data'; PLOT_DATA.mkdir(parents=True,exist_ok=True)
SPECS=[('gut_3','cellphonedb_gut.ipynb',28,10,'AJD3280'),('gut_4','cellphonedb_gut.ipynb',53,37,'AJG2309'),('liver_3','cellphonedb_liver.ipynb',28,10,'AJD3280'),('liver_4','cellphonedb_liver.ipynb',52,36,'AJG2309')]

class TableParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True); self.rows=[]; self.row=None; self.parts=None
    def handle_starttag(self, tag, attrs):
        if tag=='tr': self.row=[]
        elif tag in ('th','td') and self.row is not None: self.parts=[]
    def handle_data(self, value):
        if self.parts is not None: self.parts.append(value)
    def handle_endtag(self, tag):
        if tag in ('th','td') and self.parts is not None:
            self.row.append(''.join(self.parts).strip()); self.parts=None
        elif tag=='tr' and self.row is not None:
            self.rows.append(self.row); self.row=None

records=[]; recovered={}; display_rows=[]
for key, filename, cell_index, presence_cell, donor in SPECS:
    source_path=REPO/filename; nb=json.loads(source_path.read_text(encoding='utf-8'))
    cell=nb['cells'][cell_index]
    html=''.join(next(o['data']['text/html'] for o in cell['outputs'] if 'text/html' in o.get('data',{})))
    parser=TableParser(); parser.feed(html); rows=parser.rows
    assert len(rows)==16 and all(len(r)==16 for r in rows), (key,[len(r) for r in rows])
    assert not any(v=='...' for r in rows for v in r), key
    labels=rows[0][1:]; rowlabels=[r[0] for r in rows[1:]]
    assert rowlabels==labels and len(set(labels))==15, key
    values=[[float(v) for v in r[1:]] for r in rows[1:]]
    assert all(v>=0 and v.is_integer() for r in values for v in r), key
    values=[[int(v) for v in r] for r in values]
    assert all(values[i][j]==values[j][i] for i in range(15) for j in range(15)), key
    assert all(values[i][i]==0 for i in range(15)), key
    presence_text=''.join(next(o['data']['text/plain'] for o in nb['cells'][presence_cell]['outputs'] if 'text/plain' in o.get('data',{})))
    present=sorted(ast.literal_eval(presence_text))
    absent=sorted(set(labels)-set(present))
    with count_table_path(OUT,key).open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f); writer.writerow(['cell_type']+labels)
        writer.writerows([[label]+row for label,row in zip(labels,values)])
    if args.export_source_html:
        html_folder=AUDIT/'saved_notebook_count_tables'; html_folder.mkdir(parents=True,exist_ok=True)
        (html_folder/(count_table_path(OUT,key).stem+'.html')).write_text(html,encoding='utf-8')
    display_rows.extend([key,label,DISPLAY.get(label,label),int(label in present)] for label in labels)
    manifest={'dataset':key,'donor':donor,'source_notebook':filename,'notebook_sha256':hashlib.sha256(source_path.read_bytes()).hexdigest(),'source_cell_zero_based':cell_index,'source_execution_count':cell.get('execution_count'),'complete_table_shape':[15,15],'all_225_values_explicit':True,'symmetric':True,'diagonal_zero':True,'min_count':min(map(min,values)),'max_count':max(map(max,values)),'present_cell_types':present,'absent_cell_types':absent,'html_sha256':hashlib.sha256(html.encode()).hexdigest(),'numeric_values_from':'saved text/html table, never image pixels','raw_csv':count_table_path(OUT,key).relative_to(OUT).as_posix()}
    records.append(manifest); recovered[key]={'labels':labels,'values':values}
    print(key,'max',manifest['max_count'],'absent',absent,'execution',manifest['source_execution_count'])

csv_path=args.project_data_dir/'gut_3_cpdb_chord.csv'
if csv_path.exists():
    with csv_path.open(newline='',encoding='utf-8-sig') as f: original=list(csv.reader(f))
    observed=recovered['gut_3']; equal=(original[0][1:]==observed['labels'] and [r[0] for r in original[1:]]==observed['labels'] and [[float(v) for v in r[1:]] for r in original[1:]]==observed['values'])
    records[0]['original_chord_csv_path']=str(csv_path); records[0]['original_chord_csv_exact_match']=equal
    print('gut3 original CSV exact match:',equal)
print('gut3 and gut4 stored matrices identical:', recovered['gut_3']==recovered['gut_4'])
with (PLOT_DATA/PRESENCE_FILE).open('w',newline='',encoding='utf-8') as f:
    writer=csv.writer(f); writer.writerow(['dataset','original_cell_type','display_label','present_in_saved_metadata']); writer.writerows(display_rows)
(AUDIT/'saved_notebook_count_table_sources.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
