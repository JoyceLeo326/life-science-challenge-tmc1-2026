"""Visualize the QC-passed trilaciclib/desmethyl pair and three actual seeds."""
from pathlib import Path
import os
import csv
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FuncFormatter
from PIL import Image, ImageChops

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DATA=ROOT/'04_分子设计'/'TWT_three_seed_parent_product_sensitivity.csv'
IMAGE=Path(r'<PACKAGE_ROOT>\Documents\Codex\TMC1_R2_Scratch_20260929\09_visuals\panel_D_pose_shift.png')
rows=list(csv.DictReader(DATA.open(encoding='utf-8-sig',newline='')))
assert len(rows)==3 and all(r['design_independent_qc_pass']=='True' and r['parent_independent_qc_pass']=='True' for r in rows)
seed=[int(r['seed'])%100 for r in rows]
parent_score=[float(r['parent_vina_kcal_mol']) for r in rows]
design_score=[float(r['design_vina_kcal_mol']) for r in rows]
parent_dist=[float(r['parent_nearest_middle_A']) for r in rows]
design_dist=[float(r['design_nearest_middle_A']) for r in rows]
parent_margin=[float(r['parent_box_margin_A']) for r in rows]
design_margin=[float(r['design_box_margin_A']) for r in rows]

font=FontProperties()
bold=FontProperties()
ink='#17333d'; muted='#526773'; blue='#067da8'; red='#d34f37'; green='#278541'
fig=plt.figure(figsize=(15,7.7),facecolor='white')
fig.text(.04,.956,'更负的对接分数伴随位姿迁移：Trilaciclib 去甲基例',fontproperties=bold,fontsize=22,color=ink,va='top')
fig.text(.04,.904,'同一小鼠六链受体 · 同一中区盒 · Vina e8 · 三个独立种子；6/6 作业通过独立技术质检',fontproperties=font,fontsize=11,color=muted,va='top')

im=Image.open(IMAGE).convert('RGB')
delta=ImageChops.difference(im,Image.new('RGB',im.size,'white'))
bb=delta.point(lambda v:0 if v<20 else 255).getbbox()
pad=40
if bb: im=im.crop((max(0,bb[0]-pad),max(0,bb[1]-pad),min(im.width,bb[2]+pad),min(im.height,bb[3]+pad)))
ax0=fig.add_axes([.055,.245,.50,.59]); ax0.imshow(im,interpolation='lanczos')
ax0.set_xticks([]);ax0.set_yticks([])
for spine in ax0.spines.values(): spine.set_color('#dbe5e9');spine.set_linewidth(1)
fig.text(.055,.846,'A  首姿势叠加（种子 20260929）',fontproperties=bold,fontsize=13,color=ink)
fig.text(.055,.196,'青：亲本 Trilaciclib    赭：去甲基设计物    绿球：预设 Middle 15 残基 Cα',fontproperties=font,fontsize=10,color=ink)
fig.text(.055,.163,'亲本接触 4 个 Middle 残基；设计物 0 个。设计物距盒面仅 0.637 Å。',fontproperties=font,fontsize=10,color=ink)

def style(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_color('#9caeb7')
    ax.grid(axis='y',color='#dfe8eb',linewidth=.7)
    ax.set_xticks(seed,[str(v) for v in seed],fontproperties=font,fontsize=9)
    ax.tick_params(axis='both',labelsize=9,colors=muted)

ax1=fig.add_axes([.665,.535,.285,.29]);style(ax1)
ax1.scatter([s-.045 for s in seed],parent_score,color=blue,s=88,zorder=3,label='亲本')
ax1.scatter([s+.045 for s in seed],design_score,color=red,s=88,zorder=3,label='设计物')
ax1.set_xlim(26.7,29.3);ax1.set_ylim(-10.2,-8.45)
ax1.set_ylabel('Vina 首姿势 / kcal/mol',fontproperties=font,fontsize=10,color=ink)
ax1.set_title('B  计算分数',fontproperties=bold,fontsize=13,color=ink,loc='left',pad=9)

ax2=fig.add_axes([.665,.205,.285,.25]);style(ax2)
ax2.scatter([s-.045 for s in seed],parent_dist,color=blue,s=88,zorder=3)
ax2.scatter([s+.045 for s in seed],design_dist,color=red,s=88,zorder=3)
ax2.axhline(4,color=green,linestyle='--',linewidth=1.5)
ax2.text(29.27,4.08,'4 Å 接触判据',fontproperties=font,fontsize=8,color=green,ha='right')
ax2.set_xlim(26.7,29.3);ax2.set_ylim(2.2,7.7)
ax2.set_ylabel('距 Middle 最近 / Å',fontproperties=font,fontsize=10,color=ink)
ax2.set_title('C  预设区域接触',fontproperties=bold,fontsize=13,color=ink,loc='left',pad=9)
ax2.set_xlabel('种子 202609xx（独立重复）',fontproperties=font,fontsize=9,color=muted)

fig.text(.04,.099,'27/29：设计物分数更负，但距 Middle 6.854/6.863 Å，盒边距仅 0.634/0.637 Å；28：回到 Middle，分数反比同种子亲本高 0.046 kcal/mol。',
         fontproperties=font,fontsize=9.4,color=ink,va='top')
fig.text(.04,.061,'静态对接分数和位姿只描述该计算条件；不能据此宣称同位点亲和力改善、通道激活或功能救援。',
         fontproperties=font,fontsize=9.4,color=muted,va='top')
fig.savefig(HERE/'图09_亲本设计物区域迁移.png',dpi=300,facecolor='white')
fig.savefig(HERE/'图09_亲本设计物区域迁移.svg',facecolor='white')
plt.close(fig)
print('saved second figure')
