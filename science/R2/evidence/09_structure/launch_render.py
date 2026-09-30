"""Use the first-round PyMOL environment without installing software."""
from pathlib import Path
import os
import subprocess
import sys

here = Path(__file__).resolve().parent
prefix = Path(r'<PACKAGE_ROOT>\.cache\lsc-scientific-render\pymol-env')
env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1',
           QT_QPA_PLATFORM='offscreen', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
env['PATH'] = ';'.join(str(p) for p in [prefix, prefix/'Library/mingw-w64/bin',
    prefix/'Library/usr/bin', prefix/'Library/bin', prefix/'Scripts',
    Path(os.environ.get('SYSTEMROOT', '')), Path(os.environ.get('WINDIR', ''))])
env['CONDA_PREFIX'] = str(prefix)
env.pop('PYTHONPATH', None)
env.pop('PYTHONHOME', None)
script = here / (sys.argv[1] if len(sys.argv) > 1 else 'render_structure_panels.py')
proc = subprocess.run([str(prefix/'python.exe'), '-X', 'utf8', str(script)],
                      env=env, cwd=here, creationflags=subprocess.CREATE_NO_WINDOW,
                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                      text=True, encoding='utf-8', errors='replace')
(here/'render.log').write_text(proc.stdout, encoding='utf-8')
print(proc.stdout)
raise SystemExit(proc.returncode)
