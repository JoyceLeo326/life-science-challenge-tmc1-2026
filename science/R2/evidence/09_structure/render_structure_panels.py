"""Render coordinates from the frozen mouse six-chain model, using PyMOL 3.1.

Run with the existing isolated PyMOL environment. No docking or coordinate change.
"""
from pathlib import Path
import os
import hashlib
import json
import time

import pymol2
from pymol.cgo import BEGIN, END, VERTEX, COLOR, LINES, LINEWIDTH

HERE = Path(__file__).resolve().parent
SRC = Path(r"<PACKAGE_ROOT>")
SCRATCH = Path(r"<PACKAGE_ROOT>\Documents\Codex\TMC1_R2_Scratch_20260929\09_visuals")
SCRATCH.mkdir(parents=True, exist_ok=True)

OLD_CENTER = (2.3955, -12.9641, 0.9145)
OLD_SIZE = (37.687, 31.3744, 31.1908)
GATE_CENTER = (2.547, -20.084, 23.319)
GATE_SIZE = (22.767, 20.833, 27.715)


def box_cgo(center, size, rgb):
    lo = [c - s / 2 for c, s in zip(center, size)]
    hi = [c + s / 2 for c, s in zip(center, size)]
    vertices = [(x, y, z) for x in (lo[0], hi[0])
                for y in (lo[1], hi[1]) for z in (lo[2], hi[2])]
    edges = [(i, j) for i, a in enumerate(vertices)
             for j, b in enumerate(vertices) if j > i and sum(a[k] != b[k] for k in range(3)) == 1]
    out = [LINEWIDTH, 3.2, BEGIN, LINES, COLOR, *rgb]
    for i, j in edges:
        out.extend([VERTEX, *vertices[i], VERTEX, *vertices[j]])
    out.append(END)
    return out


def set_scene(cmd):
    cmd.set('max_threads', 1)
    cmd.set('orthoscopic', 1)
    cmd.set('antialias', 1)
    cmd.set('ray_shadows', 0)
    cmd.set('depth_cue', 0)
    cmd.set('specular', 0.15)
    cmd.set('ambient', 0.52)
    cmd.set('direct', 0.72)
    cmd.set('cartoon_fancy_helices', 1)
    cmd.set('cartoon_sampling', 10)
    cmd.set('cartoon_loop_radius', 0.20)
    cmd.set('stick_radius', 0.24)
    cmd.bg_color('white')
    cmd.load(str(SRC), 'complex')
    cmd.dss('complex')
    cmd.hide('everything', 'all')
    colors = {
        'A': (0.07, 0.47, 0.54), 'B': (0.38, 0.71, 0.72),
        'C': (0.91, 0.48, 0.25), 'D': (0.97, 0.70, 0.46),
        'E': (0.44, 0.38, 0.65), 'F': (0.71, 0.66, 0.82),
    }
    for ch, rgb in colors.items():
        cmd.set_color('ch_' + ch, rgb)
        cmd.color('ch_' + ch, 'complex and chain ' + ch)
    cmd.show('cartoon', 'complex')
    cmd.set_color('old_box', (0.03, 0.57, 0.72))
    cmd.set_color('gate_box', (0.82, 0.34, 0.20))
    cmd.set_color('gate_res', (0.96, 0.67, 0.20))
    cmd.set_color('tmie_res', (0.61, 0.33, 0.72))
    cmd.orient('complex and chain A+B')
    cmd.turn('y', 22)
    cmd.turn('x', 13)
    return colors


def frame_zoom(cmd, boxes, buffer=4):
    pts = []
    for center, size in boxes:
        lo = [c - s / 2 for c, s in zip(center, size)]
        hi = [c + s / 2 for c, s in zip(center, size)]
        pts.extend((lo, hi))
    for i, pos in enumerate(pts):
        cmd.pseudoatom('_frame', pos=pos)
    cmd.zoom('_frame', buffer=buffer)
    cmd.delete('_frame')


def save(cmd, stem, w=1200, h=920):
    cmd.png(str(SCRATCH / (stem + '.png')), width=w, height=h, dpi=300, ray=1, quiet=1)
    return list(cmd.get_view())


def main():
    start = time.time()
    sha = hashlib.sha256(SRC.read_bytes()).hexdigest()
    with pymol2.PyMOL() as pm:
        cmd = pm.cmd
        colors = set_scene(cmd)
        rotation = list(cmd.get_view()[:9])
        cmd.zoom('complex', buffer=12)
        v1 = save(cmd, 'panel_A_six_chain')

        # Two deliberately separate, documented search regions in same aligned coordinates.
        cmd.hide('everything', 'all')
        cmd.show('cartoon', 'complex and chain A+C')
        cmd.set('cartoon_transparency', 0.14, 'complex and chain A+C')
        cmd.load_cgo(box_cgo(OLD_CENTER, OLD_SIZE, (0.03, 0.57, 0.72)), 'old_box')
        cmd.load_cgo(box_cgo(GATE_CENTER, GATE_SIZE, (0.82, 0.34, 0.20)), 'gate_box')
        cmd.show('spheres', 'complex and chain A and resi 234+235 and name CA')
        cmd.color('gate_res', 'complex and chain A and resi 234+235 and name CA')
        cmd.set('sphere_scale', 0.72, 'complex and chain A and resi 234+235 and name CA')
        cmd.zoom('complex and chain A+C', buffer=10)
        v2 = save(cmd, 'panel_B_regions')

        cmd.delete('old_box'); cmd.delete('gate_box')
        cmd.hide('everything', 'all')
        f = 'complex and chain A and resi 234+235+238+397'
        t = 'complex and chain C and resi 47+48+56+58'
        cmd.select('_local_context', f'byres (complex and chain A+C within 11 of ({f} or {t}))')
        cmd.show('cartoon', '_local_context')
        cmd.set('cartoon_transparency', 0.56, '_local_context')
        cmd.show('sticks', f + ' or ' + t)
        cmd.color('gate_res', f + ' and elem C')
        cmd.color('tmie_res', t + ' and elem C')
        cmd.color('atomic', '(' + f + ' or ' + t + ') and not elem C')
        cmd.set('stick_radius', 0.31, f + ' or ' + t)
        cmd.show('spheres', 'complex and chain A and resi 234+235 and name CA')
        cmd.set('sphere_scale', 0.36, 'complex and chain A and resi 234+235 and name CA')
        cmd.zoom('_local_context', buffer=12)
        v3 = save(cmd, 'panel_C_gate_interface')
        manifest = {
            'created_local': time.strftime('%Y-%m-%d %H:%M:%S %z'),
            'source_coordinate_file': str(SRC), 'source_sha256': sha,
            'source_model': 'Published authors AF2/MD six-protein mouse model, aligned heavy atoms; no lipids, ions, water',
            'chain_mapping': {'A,B': 'mouse TMC1 residues 81-746', 'C,D': 'mouse TMIE residues 44-118', 'E,F': 'mouse CIB2 residues 1-187'},
            'boxes_Angstrom': {'old_middle': {'center': OLD_CENTER, 'size': OLD_SIZE}, 'exploratory_gate': {'center': GATE_CENTER, 'size': GATE_SIZE}},
            'camera_views': {'A': v1, 'B': v2, 'C': v3}, 'initial_rotation': rotation,
            'colors': {ch: list(rgb) for ch, rgb in colors.items()},
            'renderer': 'PyMOL ' + str(cmd.get_version()[0]), 'elapsed_seconds': round(time.time() - start, 2),
        }
    (HERE / 'structure_render_manifest.json').write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps({'elapsed_seconds': manifest['elapsed_seconds'], 'panels': [str(p) for p in SCRATCH.glob('panel_*.png')]}, ensure_ascii=False))


if __name__ == '__main__':
    main()


