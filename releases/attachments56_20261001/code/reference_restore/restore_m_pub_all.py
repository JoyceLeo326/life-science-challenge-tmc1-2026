"""Fetch fixed reference coordinates locally and restore the historical receptor.

Never distribute this script's downloaded or generated coordinate outputs.
No docking, molecular dynamics, or model inference is performed.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import urllib.request

PACKAGE = Path(__file__).resolve().parent
AUTHOR_URL = ('https://raw.githubusercontent.com/ramirezlab/Drug-design-targeting-TMC1/'
              '74274dad76201258adbd76c7058b705170c4b857/'
              'Mechano_electrical_transducer_channels_complex/MET-channel_complex.pdb')
AUTHOR_SHA = '14b5d4010eb2220ebddca0ef5cd0be346e6ba9bcc6db14311745b684e9807d82'
HUMAN_URL = 'https://alphafold.ebi.ac.uk/files/AF-Q8TDI8-F1-model_v6.pdb'
HUMAN_SHA = 'ba6b438b036b9386f4f04534cfa1dda3b6a35fbe05dd7fb2ec1b415bc6768d3b'
VERSIONS = {'biopython': '1.88', 'gemmi': '0.7.5', 'meeko': '0.8.0',
            'numpy': '2.5.3', 'rdkit': '2025.9.6', 'scipy': '1.18.1'}


def sha256(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def require_hash(path, expected):
    actual = sha256(path)
    if actual != expected:
        raise RuntimeError(f'SHA256 mismatch for {Path(path).name}: {actual} != {expected}')
    return actual


def local_output_root(value):
    root = Path(value).expanduser().resolve()
    # Coordinates may not enter the delivery tree, this package, or a Git checkout.
    if root == PACKAGE or PACKAGE in root.parents or 'completion_staging' in {p.casefold() for p in root.parts}:
        raise ValueError('Choose a separate local output directory outside completion_staging.')
    if any((p / '.git').exists() for p in [root, *root.parents]):
        raise ValueError('Choose an output directory outside any Git checkout.')
    root.mkdir(parents=True, exist_ok=True)
    (root / 'NOT_FOR_UPLOAD.txt').write_text(
        'Downloaded author coordinates and derived coordinates are local research inputs.\n'
        'Do not upload this directory, add it to Git, or put it in public archives.\n',
        encoding='utf-8')
    return root


def fetch_verified(url, target, expected, supplied=None):
    if supplied:
        source = Path(supplied).expanduser().resolve()
        require_hash(source, expected)
        if source != target:
            shutil.copyfile(source, target)
        require_hash(target, expected)
        return 'verified_existing_input'
    if target.exists():
        require_hash(target, expected)
        return 'verified_local_cache'
    temporary = target.with_suffix(target.suffix + '.download')
    request = urllib.request.Request(url, headers={'User-Agent': 'TMC1-reference-restorer/1'})
    try:
        with urllib.request.urlopen(request, timeout=120) as response, temporary.open('wb') as output:
            shutil.copyfileobj(response, output)
        require_hash(temporary, expected)
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)
    return 'fixed_url_download'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', default=str(Path.home() / 'TMC1_reference_restore_local'))
    parser.add_argument('--fetch-only', action='store_true', help='Download and verify source inputs only.')
    parser.add_argument('--source-pdb', help='Optional already downloaded author input; SHA256 is mandatory.')
    parser.add_argument('--human-pdb', help='Optional already downloaded AFDB v6 input; SHA256 is mandatory.')
    args = parser.parse_args()
    root = local_output_root(args.output_root)
    raw = root / 'raw'
    raw.mkdir(exist_ok=True)
    receipt = {'coordinate_payload_publicly_distributed': False,
               'scientific_operations': ['fixed-source retrieval and hash checking'],
               'python': platform.python_version(), 'platform': platform.platform(), 'sources': {}}
    receipt_path = root / 'restore_receipt.json'
    try:
        for label, url, expected, supplied, name in [
            ('author', AUTHOR_URL, AUTHOR_SHA, args.source_pdb, 'public_MET_complex.pdb'),
            ('human_reference', HUMAN_URL, HUMAN_SHA, args.human_pdb, 'AF-Q8TDI8-F1-model_v6.pdb')]:
            method = fetch_verified(url, raw / name, expected, supplied)
            receipt['sources'][label] = {'url': url, 'sha256': expected, 'method': method}
        if args.fetch_only:
            receipt['status'] = 'sources_verified_receptor_not_generated'
            return
        installed = {name: importlib.metadata.version(name) for name in VERSIONS}
        receipt['dependencies'] = installed
        if installed != VERSIONS or platform.python_version() != '3.12.10':
            raise RuntimeError('Use Python 3.12.10 and requirements-restore.txt; version mismatch.')
        import numpy as np
        from Bio import Align
        from Bio.Align import substitution_matrices
        from Bio.SeqUtils import seq1
        # These three function bodies are copied verbatim from the original source.
        original = (PACKAGE / 'original_prepare_receptors_v2.py').read_text(encoding='utf-8')
        require_hash(PACKAGE / 'original_prepare_receptors_v2.py',
                     '76a534122dd35c0593ceea26fda721439634ef8d91aa53e671d00edc1b9f512a')
        import ast
        tree = ast.parse(original)
        functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                     and node.name in ('atoms', 'ca_map', 'align')]
        namespace = {'np': np, 'Align': Align, 'substitution_matrices': substitution_matrices, 'seq1': seq1}
        exec(compile(ast.Module(body=functions, type_ignores=[]), 'original_prepare_receptors_v2.py', 'exec'), namespace)
        atoms, ca_map, align = (namespace[key] for key in ('atoms', 'ca_map', 'align'))
        frozen = json.loads((PACKAGE / 'frozen_alignment_inputs.json').read_text(encoding='utf-8'))
        tm = {n for a, b in frozen['human_tm_ranges'] for n in range(a, b + 1)}
        source = atoms(raw / 'public_MET_complex.pdb')
        human = atoms(raw / 'AF-Q8TDI8-F1-model_v6.pdb')
        pairs, rot, xc, yc, rmsd, fit = align(source, human, tm)
        middle = set(frozen['middle_human_resnums'])
        mapping = dict(pairs)
        present = {mapping[n] for n in ca_map(source) if mapping.get(n) in middle}
        if present != middle:
            raise RuntimeError(f'Missing historical middle-site residues: {middle - present}')
        output = root / 'screen_v2' / 'inputs' / 'receptors'
        output.mkdir(parents=True, exist_ok=True)
        aligned = output / 'M_PUB_ALL_aligned_heavy.pdb'
        lines = []
        for atom in source:
            xyz = (atom['xyz'] - xc) @ rot + yc
            lines.append(atom['line'][:30] + ''.join(f'{v:8.3f}' for v in xyz) + atom['line'][54:])
        # Preserve the historical Windows newline convention on every platform.
        with aligned.open('w', encoding='utf-8', newline='\r\n') as handle:
            handle.write('\n'.join(lines) + '\nEND\n')
        receipt['alignment'] = {'source_heavy_atoms': len(source), 'chainA_ca_count': len(ca_map(source)),
                                'aligned_pairs': len(pairs), 'conserved_tm_fit_pairs': len(fit),
                                'tm_ca_rmsd_A': rmsd, 'middle_residues_present': len(present),
                                'aligned_sha256': sha256(aligned),
                                'aligned_matches_historical': sha256(aligned) == frozen['expected_aligned_sha256']}
        if not receipt['alignment']['aligned_matches_historical']:
            raise RuntimeError('Aligned PDB differs from historical hash; parameterization stopped.')
        # Identical original Meeko flags; execute its module in the same interpreter.
        prefix = 'screen_v2/inputs/receptors/M_PUB_ALL'
        command = [sys.executable, '-m', 'meeko.cli.mk_prepare_receptor', '--read_pdb',
                   'screen_v2/inputs/receptors/M_PUB_ALL_aligned_heavy.pdb', '-o', prefix,
                   '-p', '-j', '--write_pdb', prefix + '_H.pdb']
        receipt['meeko_arguments'] = command[2:]
        with (root / 'meeko_prepare.log').open('w', encoding='utf-8') as log:
            result = subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT)
        receipt['meeko_returncode'] = result.returncode
        if result.returncode:
            raise RuntimeError('Meeko failed; inspect the local meeko_prepare.log.')
        receptor = output / 'M_PUB_ALL.pdbqt'
        actual = sha256(receptor)
        receipt['pdbqt'] = {'sha256': actual, 'expected_sha256': frozen['expected_pdbqt_sha256'],
                            'matches_historical': actual == frozen['expected_pdbqt_sha256']}
        receipt['scientific_operations'] += ['historical rigid alignment', 'Meeko template-based receptor conversion']
        import meeko
        templates = Path(meeko.__file__).resolve().parent / 'data'
        receipt['meeko_template_sha256'] = {p.name: sha256(p) for p in sorted(templates.glob('*.json'))}
        if not receipt['pdbqt']['matches_historical']:
            raise RuntimeError('Generated PDBQT differs from historical SHA256; do not use for frozen replay.')
        receipt['status'] = 'restored_and_historical_sha256_verified'
    except Exception as error:
        receipt['status'] = 'failed'
        receipt['error'] = str(error)
        raise
    finally:
        receipt_path.write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
        print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
