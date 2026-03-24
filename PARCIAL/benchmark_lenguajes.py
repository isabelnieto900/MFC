import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

HARTREE_TO_EV = 27.211386245988
E2S_EXP_EV = -5.39
EION_EXP_EV = -4.58064

ROOT = Path(__file__).resolve().parent
DF_DIR = ROOT / "finite-diferences"
ROUTINES_NR_DIR = ROOT / "routinesNR"


def run_cmd(cmd, input_text=None, cwd=None):
    t0 = time.perf_counter()
    p = subprocess.run(
        cmd,
        input=input_text,
        text=True,
        cwd=(ROOT if cwd is None else cwd),
        capture_output=True,
        check=True,
    )
    dt = time.perf_counter() - t0
    return dt, p.stdout


def parse_variational_iterations(stdout: str):
    m = re.search(r"Iteraciones \(variacional\)\s*=\s*(\d+)", stdout)
    if not m:
        return -1
    return int(m.group(1))


def parse_energy_file(path: Path):
    energies = {}
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.strip().startswith("#") or not line.strip():
                continue
            parts = line.split()
            state = parts[0]
            energies[state] = float(parts[-1])
    return energies


def parse_state_counts(path: Path):
    n_s = 0
    n_p = 0
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.strip().startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            lval = parts[1]
            if lval == "0":
                n_s += 1
            elif lval == "1":
                n_p += 1
    return n_s, n_p


def write_metrics_row(path: Path, row):
    with path.open("w", encoding="utf-8") as f:
        f.write("# language time_s var_iter e0_ev ion_ev e2s_ev err_rel_ion err_rel_2s nfound_s nfound_p\n")
        f.write(
            f"{row['language']} {row['time_s']:.10f} {row['var_iter']} {row['e0_ev']:.10f} {row['ion_ev']:.10f} "
            f"{row['e2s_ev']:.10f} {row['err_rel_ion']:.10f} {row['err_rel_2s']:.10f} {row['nfound_s']} {row['nfound_p']}\n"
        )


def main():
    # Fortran build and run
    run_cmd([
        "gfortran",
        "-O2",
        "-std=f2008",
        str(ROUTINES_NR_DIR / "nr_compat_litio.f90"),
        str(ROUTINES_NR_DIR / "nr_powell_compat.f90"),
        str(ROUTINES_NR_DIR / "litio_nr_pipeline.f90"),
        "-o",
        "litio_nr_pipeline",
    ], cwd=DF_DIR)
    t_fortran, out_fortran = run_cmd(["./litio_nr_pipeline"], input_text="3\n", cwd=DF_DIR)

    it_fortran = parse_variational_iterations(out_fortran)
    energies_fortran = parse_energy_file(DF_DIR / "energias_litio_nr.dat")

    e2s_f = energies_fortran["2s"]
    ion_f = energies_fortran.get("ion", -4.5828)

    nfound_s_f, nfound_p_f = parse_state_counts(DF_DIR / "energias_litio_nr.dat")

    row_fortran = {
        "language": "fortran",
        "time_s": t_fortran,
        "var_iter": it_fortran,
        "e0_ev": -201.1172,
        "ion_ev": -4.5828,
        "e2s_ev": e2s_f,
        "err_rel_ion": abs((-4.5828 - EION_EXP_EV) / EION_EXP_EV),
        "err_rel_2s": abs((e2s_f - E2S_EXP_EV) / E2S_EXP_EV),
        "nfound_s": nfound_s_f,
        "nfound_p": nfound_p_f,
    }

    write_metrics_row(DF_DIR / "metrics_litio_fortran.dat", row_fortran)

    shutil.copyfile(DF_DIR / "potential_litio.dat", DF_DIR / "potential_litio_fortran.dat")
    shutil.copyfile(DF_DIR / "funciones_radiales_litio.dat", DF_DIR / "funciones_radiales_litio_fortran.dat")
    shutil.copyfile(DF_DIR / "energias_litio_nr.dat", DF_DIR / "energias_litio_fortran.dat")

    # C++ build and run
    run_cmd(["g++", "-O2", "-std=c++17", "litio_pipeline_cpp.cpp", "-o", "litio_pipeline_cpp"], cwd=DF_DIR)
    run_cmd(["./litio_pipeline_cpp", "3"], cwd=DF_DIR)

    # Python run
    run_cmd([sys.executable, "litio_pipeline_py.py", "3"], cwd=DF_DIR)

    # Collect metrics files
    rows = []
    for p in [DF_DIR / "metrics_litio_fortran.dat", DF_DIR / "metrics_litio_cpp.dat", DF_DIR / "metrics_litio_py.dat"]:
        with p.open("r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if parts[0] == "fortran":
                    n_s, n_p = parse_state_counts(DF_DIR / "energias_litio_fortran.dat")
                elif parts[0] == "cpp":
                    n_s, n_p = parse_state_counts(DF_DIR / "energias_litio_cpp.dat")
                else:
                    n_s, n_p = parse_state_counts(DF_DIR / "energias_litio_py.dat")
                rows.append({
                    "language": parts[0],
                    "time_s": float(parts[1]),
                    "var_iter": int(parts[2]),
                    "e0_ev": float(parts[3]),
                    "ion_ev": float(parts[4]),
                    "e2s_ev": float(parts[5]),
                    "err_rel_ion": float(parts[6]),
                    "err_rel_2s": float(parts[7]),
                    "nfound_s": n_s,
                    "nfound_p": n_p,
                })

    rows = sorted(rows, key=lambda x: ["fortran", "cpp", "py"].index(x["language"]))

    with (ROOT / "comparacion_lenguajes.dat").open("w", encoding="utf-8") as f:
        f.write("# language time_s var_iter e0_ev ion_ev e2s_ev err_rel_ion err_rel_2s nfound_s nfound_p\n")
        for r in rows:
            f.write(
                f"{r['language']} {r['time_s']:.10f} {r['var_iter']} {r['e0_ev']:.10f} {r['ion_ev']:.10f} "
                f"{r['e2s_ev']:.10f} {r['err_rel_ion']:.10f} {r['err_rel_2s']:.10f} {r['nfound_s']} {r['nfound_p']}\n"
            )

    print("Generado comparacion_lenguajes.dat")


if __name__ == "__main__":
    main()
