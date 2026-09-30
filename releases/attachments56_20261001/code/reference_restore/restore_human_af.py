"""Restore the frozen human H_AF inputs beside the local mouse receptor.

Uses the same frozen alignment functions and Meeko options as the original
scientific pipeline. Coordinate outputs remain outside the code/Git package.
"""
import argparse
import ast
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

EXPECTED_ALIGNED = '93572ba54235e6e338205e6272d81093d5c938ba21594a676b029fa8dd5a8044'
EXPECTED_PDBQT = '35b4acccd44a5c2c5c8aa7e756abe26c081018a2c26d0899190751c23392f2c1'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', required=True)
    parser.add_argument('--reference-dir', type=Path, default=Path(__file__).resolve().parent,
                        help='Directory containing the unchanged frozen mouse restoration method.')
    parser.add_argument('--human-pdb', help='Optional lawful existing fixed source; SHA remains mandatory.')
    args = parser.parse_args()
    reference = args.reference_dir.resolve()
    spec = importlib.util.spec_from_file_location('frozen_mouse_restore', reference/'restore_m_pub_all.py')
    restore = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(restore)
    root = restore.local_output_root(args.output_root)
    raw = root/'raw'; raw.mkdir(exist_ok=True)
    started = time.perf_counter()
    receipt = {'coordinate_payload_publicly_distributed': False,
               'python': platform.python_version(), 'platform': platform.platform(),
               'method': 'Original H_AF self-reference rigid alignment and original Meeko flags; no new modelling, docking or MD.'}
    try:
        installed = {name: importlib.metadata.version(name) for name in restore.VERSIONS}
        receipt['dependencies'] = installed
        if installed != restore.VERSIONS or platform.python_version() != '3.12.10':
            raise RuntimeError('Use Python 3.12.10 and requirements-restore.txt; version mismatch.')
        source = raw/'AF-Q8TDI8-F1-model_v6.pdb'
        method = restore.fetch_verified(restore.HUMAN_URL, source, restore.HUMAN_SHA, args.human_pdb)
        receipt['source'] = {'url': restore.HUMAN_URL, 'sha256': restore.HUMAN_SHA, 'method': method}
        original = reference/'original_prepare_receptors_v2.py'
        restore.require_hash(original, '76a534122dd35c0593ceea26fda721439634ef8d91aa53e671d00edc1b9f512a')
        import numpy as np
        from Bio import Align
        from Bio.Align import substitution_matrices
        from Bio.SeqUtils import seq1
        tree = ast.parse(original.read_text(encoding='utf-8'))
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                     and node.name in ('atoms','ca_map','align')]
        namespace = {'np':np,'Align':Align,'substitution_matrices':substitution_matrices,'seq1':seq1}
        exec(compile(ast.Module(body=functions,type_ignores=[]), original.name, 'exec'),namespace)
        atoms, align = namespace['atoms'],namespace['align']
        frozen = json.loads((reference/'frozen_alignment_inputs.json').read_text('utf-8'))
        tm = {n for a,b in frozen['human_tm_ranges'] for n in range(a,b+1)}
        human = atoms(source)
        pairs,rot,xc,yc,rmsd,fit = align(human,human,tm)
        output = root/'screen_v2/inputs/receptors'; output.mkdir(parents=True,exist_ok=True)
        aligned = output/'H_AF_aligned_heavy.pdb'
        lines=[]
        for atom in human:
            xyz=(atom['xyz']-xc)@rot+yc
            lines.append(atom['line'][:30]+''.join(f'{v:8.3f}' for v in xyz)+atom['line'][54:])
        with aligned.open('w',encoding='utf-8',newline='\r\n') as handle:
            handle.write('\n'.join(lines)+'\nEND\n')
        restore.require_hash(aligned,EXPECTED_ALIGNED)
        receipt['alignment']={'source_heavy_atoms':len(human),'aligned_pairs':len(pairs),
                              'conserved_tm_fit_pairs':len(fit),'tm_ca_rmsd_A':rmsd,
                              'aligned_sha256':EXPECTED_ALIGNED,'matches_historical':True}
        prefix='screen_v2/inputs/receptors/H_AF'
        command=[sys.executable,'-m','meeko.cli.mk_prepare_receptor','--read_pdb',prefix+'_aligned_heavy.pdb',
                 '-o',prefix,'-p','-j','--write_pdb',prefix+'_H.pdb']
        with (root/'human_meeko_prepare.log').open('w',encoding='utf-8') as handle:
            result=subprocess.run(command,cwd=root,stdout=handle,stderr=subprocess.STDOUT)
        receipt['meeko_arguments']=command[2:]
        receipt['meeko_returncode']=result.returncode
        if result.returncode:raise RuntimeError('Human Meeko preparation failed; inspect local log.')
        restore.require_hash(output/'H_AF.pdbqt',EXPECTED_PDBQT)
        receipt['pdbqt']={'sha256':EXPECTED_PDBQT,'matches_historical':True}
        mouse_aligned=output/'M_PUB_ALL_aligned_heavy.pdb'; mouse_pdbqt=output/'M_PUB_ALL.pdbqt'
        if mouse_aligned.exists() and mouse_pdbqt.exists():
            restore.require_hash(mouse_aligned,frozen['expected_aligned_sha256'])
            restore.require_hash(mouse_pdbqt,frozen['expected_pdbqt_sha256'])
            inputs={'mouse_pdbqt':str(mouse_pdbqt.relative_to(root)).replace('\\','/'),
                    'human_pdbqt':str((output/'H_AF.pdbqt').relative_to(root)).replace('\\','/'),
                    'mouse_aligned_pdb':str(mouse_aligned.relative_to(root)).replace('\\','/'),
                    'human_aligned_pdb':str(aligned.relative_to(root)).replace('\\','/')}
            (root/'EXTERNAL_INPUTS.local.json').write_text(json.dumps(inputs,indent=2)+'\n',encoding='utf-8')
            receipt['four_external_inputs_verified']=True
        else:
            receipt['four_external_inputs_verified']=False
        receipt['status']='restored_and_historical_sha256_verified'
    except Exception as error:
        receipt['status']='failed';receipt['error']=str(error)
        raise
    finally:
        receipt['elapsed_seconds']=round(time.perf_counter()-started,3)
        (root/'human_restore_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
