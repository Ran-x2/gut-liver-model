"""Replot exact historical CellPhoneDB count tables; never rerun inference."""
import argparse, csv, hashlib, json, platform
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
from cellphone_count_data import CASES, COUNT_FIGURE_LEGEND, FIGURE_STEM, PRESENCE_FILE, count_table_path

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir',type=Path,default=Path(__file__).resolve().parents[1]/'output'/'10_cellphonedb')
args=parser.parse_args(); OUT=args.output_dir; FIGS=OUT/'output'; FIGS.mkdir(parents=True,exist_ok=True)
AUDIT=OUT/'audit'; AUDIT.mkdir(parents=True,exist_ok=True)
with (OUT/'plot_data'/PRESENCE_FILE).open(newline='',encoding='utf-8') as f:
    metadata=list(csv.DictReader(f))
presence={(r['dataset'],r['original_cell_type']):bool(int(r['present_in_saved_metadata'])) for r in metadata}
if len(presence)!=len(metadata): raise ValueError('Duplicate cell-population metadata')
DISPLAY={r['original_cell_type']:r['display_label'] for r in metadata}
SPECS=[(key,letter,f'{tissue.title()} | {donor}') for key,letter,tissue,_,donor in CASES]
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.linewidth':0.6,'savefig.facecolor':'white','xtick.major.size':0,'ytick.major.size':0,'figure.max_open_warning':0})
cmap=plt.get_cmap('YlGnBu').copy(); cmap.set_bad('#DEDEDE'); norm=Normalize(0,50)
loaded={}
for key,letter,title in SPECS:
    with count_table_path(OUT,key).open(newline='',encoding='utf-8') as f: rows=list(csv.reader(f))
    labels=rows[0][1:]; values=np.array([[int(v) for v in row[1:]] for row in rows[1:]],dtype=int)
    assert values.shape==(15,15) and np.array_equal(values,values.T)
    mask=np.eye(15,dtype=bool)
    for i,label in enumerate(labels):
        present=presence[key,label]
        if not present: mask[i,:]=True; mask[:,i]=True
    loaded[key]=(labels,values,np.ma.masked_array(values,mask=mask))

def draw(ax,key,title,letter,ticksize=7.2,numbersize=5.9):
    labels,values,masked=loaded[key]
    ax.set_facecolor('#DEDEDE')
    ax.pcolormesh(np.arange(16),np.arange(16),masked,cmap=cmap,norm=norm,shading='flat',edgecolors='white',linewidth=0.32,rasterized=False)
    ax.set_aspect('equal'); ax.set_xlim(0,15); ax.set_ylim(15,0)
    ax.set_xticks(np.arange(15)+0.5,labels=[DISPLAY.get(l,l) for l in labels],rotation=55,ha='right',rotation_mode='anchor',fontsize=ticksize)
    ax.set_yticks(np.arange(15)+0.5,labels=[DISPLAY.get(l,l) for l in labels],fontsize=ticksize)
    ax.tick_params(axis='both',pad=2)
    for spine in ax.spines.values(): spine.set_visible(False)
    for i in range(15):
        for j in range(15):
            if not masked.mask[i,j]:
                ax.text(j+0.5,i+0.5,str(values[i,j]),ha='center',va='center',fontsize=numbersize,color='white' if values[i,j]>=28 else '#101010')
    ax.set_title(f'{letter}   {title}',loc='center',fontsize=10,fontweight='bold',pad=8,fontfamily='Arial')

fig,axes=plt.subplots(2,2,figsize=(8.8,8.9))
fig.subplots_adjust(left=0.115,right=0.985,top=0.965,bottom=0.165,wspace=0.45,hspace=0.42)
for ax,(key,letter,title) in zip(axes.flat,SPECS): draw(ax,key,title,letter)
cax=fig.add_axes([0.27,0.064,0.51,0.014])
cbar=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='horizontal',ticks=np.arange(0,51,10))
cbar.solids.set_rasterized(False)
cbar.outline.set_linewidth(0.5); cbar.ax.tick_params(labelsize=8,length=2,pad=2)
cbar.set_label('Retained significant interactions (both partner assignments)',fontsize=8,labelpad=4)
fig.legend(handles=[Patch(facecolor='#DEDEDE',edgecolor='white',label='Self-pair excluded or cell type absent')],loc='lower center',bbox_to_anchor=(0.5,0.003),frameon=False,fontsize=8)
base=FIGS/FIGURE_STEM
fig.savefig(base.with_suffix('.png'),dpi=600)
plt.close(fig)
base.with_name(base.name+'_legend.txt').write_text(COUNT_FIGURE_LEGEND,encoding='utf-8')
manifest={'figure':FIGURE_STEM,'formats':['png'],'font':'Arial','title_alignment':'center','subtitles':False,'native_width_inches':8.8,'native_height_inches':8.9,'raster_export_dpi':600,'quantitative_common_scale':[0,50],'plotted_quantity':'Exact saved symmetrized significant-interaction counts, after the original 5th-percentile expression filter and before the R chord count<5 suppression','numeric_input':[count_table_path(OUT,key).relative_to(OUT).as_posix() for key,_,_,_,_ in CASES],'new_CellPhoneDB_inference':False,'raw_numeric_tables_modified':False,'display_mask':'Diagonal self-pairs plus absent cell types per saved metadata; masked gray. Stored CSV zeros preserved.','normalization':'None; a shared linear count scale is used for all panels.','versions':{'python':platform.python_version(),'numpy':np.__version__,'matplotlib':matplotlib.__version__},'files':[]}
for path in sorted(FIGS.iterdir()):
    if path.is_file() and path.suffix.lower()=='.png' and not path.stem.endswith('_preview'): manifest['files'].append({'file':path.name,'bytes':path.stat().st_size,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
(AUDIT/'interaction_count_figure_provenance.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(json.dumps({'figure_dir':str(FIGS),'files':[p['file'] for p in manifest['files']],'common_count_scale':[0,50]},indent=2))
