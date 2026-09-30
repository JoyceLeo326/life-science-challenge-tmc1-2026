"""Compose PyMOL coordinate renders into an annotated, editable SVG and 300-dpi PNG."""
from pathlib import Path
import os
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from PIL import Image, ImageChops

HERE = Path(__file__).resolve().parent
SCRATCH = Path(r'<PACKAGE_ROOT>\Documents\Codex\TMC1_R2_Scratch_20260929\09_visuals')
FONT = FontProperties()
FONT_B = FontProperties()
INK = '#17333d'
MUTED = '#4f6570'

fig = plt.figure(figsize=(16, 8.4), facecolor='white')
fig.text(.045, .957, '鼠 TMC1–TMIE–CIB2 六链模型：中区与门控环是不同的空间假设',
         fontproperties=FONT_B, fontsize=23, color=INK, va='top')
fig.text(.045, .907, '作者公开 AF2/MD 坐标模型  ·  M_PUB_ALL 对齐重原子  ·  所有面板使用同一坐标系',
         fontproperties=FONT, fontsize=11.5, color=MUTED, va='top')

panels = [
    ('A', '六链复合体总览', 'panel_A_six_chain.png',
     'TMC1 A/B（青） · TMIE C/D（橙） · CIB2 E/F（紫）'),
    ('B', '两个事先定义的搜索区域', 'panel_B_regions.png',
     '蓝框：首轮中区  ·  赭框：门控环探索盒\n黄球：A 链 F234/G235 的 Cα；首轮盒未覆盖这两个残基'),
    ('C', 'F234 附近的 TMC1–TMIE 界面', 'panel_C_gate_interface.png',
     '金色：A 链 F234/G235/Y238/W397\n紫色：TMIE C 链 K47/E48/R56/W58；F234–E48 最近 2.92 Å'),
]
xs = [.045, .365, .685]
for x, (letter, title, file, caption) in zip(xs, panels):
    ax = fig.add_axes([x, .285, .285, .55])
    im = Image.open(SCRATCH / file).convert('RGB')
    # Remove uninformative white margin while retaining every rendered atom/box.
    delta = ImageChops.difference(im, Image.new('RGB', im.size, 'white'))
    bbox = delta.point(lambda v: 0 if v < 20 else 255).getbbox()
    pad = 42
    if bbox:
        bbox = (max(0, bbox[0]-pad), max(0, bbox[1]-pad),
                min(im.width, bbox[2]+pad), min(im.height, bbox[3]+pad))
        im = im.crop(bbox)
    ax.imshow(im, interpolation='lanczos')
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_color('#dce5e8'); spine.set_linewidth(1.2)
    fig.text(x, .858, f'{letter}   {title}', fontproperties=FONT_B, fontsize=13.7,
             color=INK, va='bottom')
    fig.text(x, .255, caption, fontproperties=FONT, fontsize=10.5,
             color=INK, va='top', linespacing=1.7)

fig.text(.045, .105,
         '盒中心与尺寸（Å）：中区 (2.396, −12.964, 0.915), 37.687 × 31.374 × 31.191；探索盒 (2.547, −20.084, 23.319), 22.767 × 20.833 × 27.715。',
         fontproperties=FONT, fontsize=9.2, color=MUTED, va='top')
fig.text(.045, .069,
         '这是静态结构与几何示意。模型未含膜脂、水或离子；界面邻近不等于已证实小分子口袋，也不能推断通道激活或抑制。',
         fontproperties=FONT, fontsize=9.2, color=MUTED, va='top')
fig.savefig(HERE / '图09_六链结构与双区域.png', dpi=300, facecolor='white')
fig.savefig(HERE / '图09_六链结构与双区域.svg', facecolor='white')
plt.close(fig)
print('saved', HERE / '图09_六链结构与双区域.png')
