#!/usr/bin/env python3
from pathlib import Path
import shutil
import subprocess
import sys


def main() -> int:
    here = Path(__file__).resolve().parent
    exe = here / "litio_shooting"
    if not exe.exists():
        print(f"No existe ejecutable: {exe}", file=sys.stderr)
        return 1

    p = subprocess.run([str(exe)], cwd=str(here), check=False)
    if p.returncode != 0:
        return p.returncode

    # Guarda una copia etiquetada como salida del wrapper Python.
    src_e = here / "energias_shooting_fortran.dat"
    src_w = here / "funciones_shooting_fortran.dat"
    if src_e.exists():
        shutil.copyfile(src_e, here / "energias_shooting_py.dat")
    if src_w.exists():
        shutil.copyfile(src_w, here / "funciones_shooting_py.dat")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
