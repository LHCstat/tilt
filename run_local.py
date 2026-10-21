"""Desktop launcher; synthetic demo and real POSCAR use the same analysis CLI."""
import argparse
from datetime import datetime
from pathlib import Path
import subprocess
import sys


def main():
    parser=argparse.ArgumentParser(description='Select/drag a POSCAR or CONTCAR and generate a new result directory.')
    parser.add_argument('poscar',nargs='?',type=Path)
    parser.add_argument('--azimuth',type=float)
    parser.add_argument('--out',type=Path)
    parser.add_argument('--demo',action='store_true')
    parser.add_argument('--display', choices=('auto','native','lno'), default='auto')
    args=parser.parse_args()
    root=Path(__file__).resolve().parent
    if args.demo and args.poscar:
        parser.error('Choose --demo or a POSCAR, not both')
    source=root/'examples/POSCAR_rx8' if args.demo else args.poscar
    if source is None:
        try:
            import tkinter as tk
            from tkinter import filedialog
            window=tk.Tk();window.withdraw()
            selected=filedialog.askopenfilename(title='Select POSCAR or CONTCAR',parent=window)
            window.destroy()
            if not selected:
                print('No file selected.');return 0
            source=Path(selected)
        except Exception as exc:
            print('File picker unavailable:',exc)
            parser.error('Pass the POSCAR path as an argument instead')
    out=args.out or root/'results'/('tilt_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    cmd=[sys.executable,'-X','utf8',str(root/'ni_octahedra.py'),str(source.resolve()),'--out',str(out.resolve())]
    cmd += ['--display', args.display]
    azimuth=0 if args.demo else args.azimuth
    if azimuth is not None:cmd += ['--azimuth',str(azimuth)]
    else:print('Automatic in-plane reference. For comparing conventional Amam cells, consider --azimuth 45.')
    result=subprocess.run(cmd)
    if result.returncode==0:print('Results:',out.resolve())
    return result.returncode

if __name__=='__main__':sys.exit(main())
