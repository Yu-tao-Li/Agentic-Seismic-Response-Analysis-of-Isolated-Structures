import json
import math
from pathlib import Path

import numpy as np


def model_parameters():
    masses = np.array([2.50e5, 2.70e5, 2.70e5, 1.80e5], dtype=float)
    stiffness = np.array([50e6, 245e6, 195e6, 98e6], dtype=float)
    yield_forces = np.array([math.inf, 1225e3, 975e3, 490e3], dtype=float)
    bmat = np.array(
        [[1.0, 0.0, 0.0, 0.0],
         [-1.0, 1.0, 0.0, 0.0],
         [0.0, -1.0, 1.0, 0.0],
         [0.0, 0.0, -1.0, 1.0]],
        dtype=float,
    )
    k0 = bmat.T @ np.diag(stiffness) @ bmat
    return {
        "masses": masses,
        "stiffness": stiffness,
        "yield_forces": yield_forces,
        "post_yield": np.array([math.nan, 0.0, 0.0, 0.0], dtype=float),
        "c_iso": 1000e3,
        "B": bmat,
        "K0": k0,
        "M": np.diag(masses),
    }


def rayleigh_from_initial_modes(mass, stiffness, zeta=0.05):
    vals, _ = np.linalg.eig(np.linalg.solve(mass, stiffness))
    omegas = np.sqrt(np.sort(np.real(vals)))
    w1 = omegas[0]
    w2 = omegas[-1]
    mat = np.array([[1.0 / (2.0 * w1), w1 / 2.0], [1.0 / (2.0 * w2), w2 / 2.0]])
    alpha_m, beta_k = np.linalg.solve(mat, np.array([zeta, zeta]))
    return float(alpha_m), float(beta_k), omegas


def fmt(value):
    return f"{float(value):.12e}"


def append_amplitude(lines, time_vec, ag):
    lines.append("*AMPLITUDE, NAME=AG, TIME=TOTAL TIME")
    pairs = [f"{fmt(t)}, {fmt(a)}" for t, a in zip(time_vec, ag)]
    for i in range(0, len(pairs), 4):
        lines.append(", ".join(pairs[i:i + 4]))


def write_input(script_dir):
    root_dir = script_dir.parent.parent
    input_file = root_dir / "input_data" / "Northridge_01_NO_968.txt"
    raw = np.loadtxt(input_file, dtype=float).reshape(-1)
    if raw.size == 0:
        raise RuntimeError(f"Input ground motion file is empty: {input_file}")
    dt = 0.01
    g = 9.81
    target_pga = 0.40 * g
    scale = target_pga / np.max(np.abs(raw))
    ag = raw * scale
    time_vec = np.arange(raw.size, dtype=float) * dt

    model = model_parameters()
    alpha_m, beta_k, omegas = rayleigh_from_initial_modes(model["M"], model["K0"], 0.05)
    c_story = beta_k * model["stiffness"]
    c_mass = alpha_m * model["masses"]

    lines = []
    lines.extend([
        "*HEADING",
        "Four-story isolated shear building; generated in this run for ABAQUS verification.",
        "*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO, CONTACT=NO",
        "** Nodes are placed one metre apart; truss area is 1 m^2, so E = story stiffness.",
        "*NODE",
        "1, 0.0, 0.0",
        "2, 1.0, 0.0",
        "3, 2.0, 0.0",
        "4, 3.0, 0.0",
        "5, 4.0, 0.0",
        "*NSET, NSET=NRESP",
        "2, 3, 4, 5",
        "*NSET, NSET=NTOP",
        "5",
        "*NSET, NSET=NISO",
        "2",
        "*NSET, NSET=NGROUND",
        "1",
        "*ELEMENT, TYPE=T2D2, ELSET=TRUSS1",
        "1, 1, 2",
        "*ELEMENT, TYPE=T2D2, ELSET=TRUSS2",
        "2, 2, 3",
        "*ELEMENT, TYPE=T2D2, ELSET=TRUSS3",
        "3, 3, 4",
        "*ELEMENT, TYPE=T2D2, ELSET=TRUSS4",
        "4, 4, 5",
        "*ELSET, ELSET=TRUSSES",
        "1, 2, 3, 4",
    ])

    for i in range(4):
        mat = f"MAT{i + 1}"
        lines.extend([
            f"*SOLID SECTION, ELSET=TRUSS{i + 1}, MATERIAL={mat}",
            "1.0",
            f"*MATERIAL, NAME={mat}",
            "*DENSITY",
            "1.0e-12",
            "*ELASTIC",
            f"{fmt(model['stiffness'][i])}, 0.0",
        ])
        if i > 0:
            lines.extend([
                "*PLASTIC",
                f"{fmt(model['yield_forces'][i])}, 0.0",
            ])

    for i, mass in enumerate(model["masses"], start=1):
        node = i + 1
        ele = 100 + i
        lines.extend([
            f"*ELEMENT, TYPE=MASS, ELSET=MASS{i}",
            f"{ele}, {node}",
            f"*MASS, ELSET=MASS{i}",
            fmt(mass),
        ])

    dash_ele = 201
    for story in range(1, 5):
        coeff = c_story[story - 1]
        if story == 1:
            coeff += model["c_iso"]
        n1 = story
        n2 = story + 1
        lines.extend([
            f"*ELEMENT, TYPE=DASHPOT2, ELSET=DASH_STORY{story}",
            f"{dash_ele}, {n1}, {n2}",
            f"*DASHPOT, ELSET=DASH_STORY{story}",
            "1, 1",
            fmt(coeff),
        ])
        dash_ele += 1

    for i, coeff in enumerate(c_mass, start=1):
        node = i + 1
        lines.extend([
            f"*ELEMENT, TYPE=DASHPOT2, ELSET=DASH_MASS{i}",
            f"{dash_ele}, 1, {node}",
            f"*DASHPOT, ELSET=DASH_MASS{i}",
            "1, 1",
            fmt(coeff),
        ])
        dash_ele += 1

    lines.extend([
        "*BOUNDARY",
        "1, 1, 2, 0.0",
        "2, 2, 2, 0.0",
        "3, 2, 2, 0.0",
        "4, 2, 2, 0.0",
        "5, 2, 2, 0.0",
    ])
    append_amplitude(lines, time_vec, ag)
    total_time = time_vec[-1]
    lines.extend([
        "*STEP, NAME=Dynamic, NLGEOM=NO, INC=100000",
        "*DYNAMIC, ALPHA=0.0",
        f"{fmt(dt)}, {fmt(total_time)}, {fmt(1.0e-8)}, {fmt(dt)}",
        "*CLOAD, AMPLITUDE=AG",
    ])
    for i, mass in enumerate(model["masses"], start=1):
        node = i + 1
        lines.append(f"{node}, 1, {fmt(-mass)}")
    lines.extend([
        "*OUTPUT, FIELD, TIME INTERVAL=0.01",
        "*NODE OUTPUT, NSET=NRESP",
        "U, V, A",
        "*ELEMENT OUTPUT, ELSET=TRUSSES",
        "S, E, PE, PEEQ",
        "*END STEP",
        "",
    ])

    inp_path = script_dir / "isolated_shear_building.inp"
    inp_path.write_text("\n".join(lines), encoding="utf-8")

    meta = {
        "input_file": str(input_file),
        "raw_count": int(raw.size),
        "dt": dt,
        "scale_to_mps2": float(scale),
        "normalized_pga_mps2": float(np.max(np.abs(ag))),
        "target_pga_g": 0.40,
        "total_time": float(total_time),
        "masses_kg": model["masses"].tolist(),
        "story_stiffness_N_per_m": model["stiffness"].tolist(),
        "yield_forces_N": [None, 1225e3, 975e3, 490e3],
        "post_yield_ratios": [None, 0.0, 0.0, 0.0],
        "isolation_dashpot_Ns_per_m": model["c_iso"],
        "rayleigh": {
            "zeta": 0.05,
            "mode_pair": [1, 4],
            "alpha_mass": alpha_m,
            "beta_initial_stiffness": beta_k,
            "c_story_dashpots_Ns_per_m": c_story.tolist(),
            "c_mass_dashpots_Ns_per_m": c_mass.tolist(),
            "initial_omega_rad_per_s": omegas.tolist(),
        },
    }
    (script_dir / "abaqus_input_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return inp_path


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    path = write_input(here)
    print(path)
