from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).parent
def read(n):return list(csv.DictReader((O/n).open(encoding='utf-8-sig')))
c=[r for r in read('R3_top10_replay_table.csv') if r['actual_score_kcal_mol']]
r=[x for x in read('R3_reference_replay_table.csv') if x['score_kcal_mol']]
fig,axs=plt.subplots(1,2,figsize=(15,7),sharex=True,layout='constrained')
for ax,rows,label,key in [(axs[0],c,'Frozen generated top 10: 9 docked, 0 joint pass','actual_score_kcal_mol'),(axs[1],r,'Reference panel: 8 docked, all technical QC passed','score_kcal_mol')]:
 scores=[float(x[key]) for x in rows];names=[x.get('name',x['compound_id'].replace('GEN_','')) for x in rows]
 colors=['#2a809b' if x['middle_contact_4A']=='True' else '#c99858' for x in rows]
 ax.barh(names,scores,color=colors);ax.invert_yaxis();ax.axvline(-8.079600334,color='#bf4444',ls='--',lw=1.7,label='Frozen score threshold')
 ax.set_xlim(-10.3,0);ax.set_xlabel('Vina first-pose score (kcal/mol; more negative is favorable)');ax.set_title(label,fontsize=12)
 for i,s in enumerate(scores):ax.text(s-0.08,i,f'{s:.3f}',ha='right',va='center',fontsize=9)
 ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
fig.suptitle('R3 same-protocol comparison | M_PUB_ALL, Middle box, e8, seed 20260927',fontsize=15)
fig.supxlabel('Blue: Middle contact <4 A; ochre: no contact. Reference labels do not establish direct TMC1 activity.\n6 dye-reduction drugs docked: 3 joint passes; generated 9: 0. Limited panel, not a biological model benchmark.',fontsize=11)
fig.savefig(O/'R3_same_protocol_score_comparison.png',dpi=180);fig.savefig(O/'R3_same_protocol_score_comparison.pdf')
