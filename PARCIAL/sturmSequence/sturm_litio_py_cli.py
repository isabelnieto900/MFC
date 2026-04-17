#!/usr/bin/env python3
from pathlib import Path
from typing import Optional
import shutil
import subprocess
import sys


def run_channel(exe: Path, cwd: Path, lval: int) -> int:
    p = subprocess.run(
        [str(exe)],
        cwd=str(cwd),
        input=f"{lval}\n",
        text=True,
        check=False,
    )
    return p.returncode


def parse_channel_arg() -> Optional[int]:
    if len(sys.argv) > 1:
        try:
            return int(sys.argv[1])
        except ValueError:
            pass

    if not sys.stdin.isatty():
        data = sys.stdin.read().strip()
        if data:
            token = data.split()[0]
            try:
                return int(token)
            except ValueError:
                return None

    return None


def main() -> int:
    here = Path(__file__).resolve().parent
    exe = here / "sturm_litio_cli"
    if not exe.exists():
        print(f"No existe ejecutable: {exe}", file=sys.stderr)
        return 1

    # Si recibe un canal por argumento/stdin, ejecuta solo ese.
    # Si no, ejecuta los 3 canales para uso interactivo.
    only_l = parse_channel_arg()
    channels = (only_l,) if only_l is not None else (0, 1, 2)

    for l in channels:
        rc = run_channel(exe, here, l)
        if rc != 0:
            return rc
        src = here / "LITIO_STURM.DAT"
        if src.exists() and l in (0, 1):
            shutil.copyfile(src, here / f"LITIO_STURM_PY{l}.DAT")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
