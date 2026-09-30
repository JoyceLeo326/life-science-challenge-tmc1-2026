"""Single official entry: preserve R2 candidate CSV; recompute R2 and R3 evidence."""
import argparse,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'results.csv');p.add_argument('--progress-delay',type=float,default=0);a=p.parse_args()
    subprocess.run([sys.executable,'-B',str(ROOT/'make_results_r2.py'),'--output',str(a.output),'--progress-delay',str(a.progress_delay)],check=True)
    subprocess.run([sys.executable,'-B',str(ROOT/'verify_r3_evidence.py'),'--output',str(a.output.parent/'r3_summary.csv')],check=True)
    print('Official integrated entry PASS: six frozen candidates unchanged; R3 evidence verified separately.')
if __name__=='__main__':main()
