"""Execute the walkthrough using the current scientific Python environment."""
from pathlib import Path
import sys
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager

here = Path(__file__).resolve().parent
source = here / 'TMC1_walkthrough.ipynb'
book = nbformat.read(source, as_version=4)
nbformat.validate(book)
manager = KernelManager(kernel_name='python3')
manager.kernel_spec.argv[0] = sys.executable
NotebookClient(book, timeout=300, kernel_name='python3', km=manager,
               resources={'metadata': {'path': str(here)}}).execute()
output = here.parent.parent / 'notebook_outputs'
output.mkdir(exist_ok=True)
nbformat.write(book, output / 'TMC1_walkthrough.executed.ipynb')
print('PASS: executed Notebook saved outside the code package in notebook_outputs.')
