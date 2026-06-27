import json
import math
from pathlib import Path

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parents[1]
GM_PATH = ROOT_DIR / "input_data" / "GM1.txt"

DT = 0.02
G = 9.81
TARGET_PGA = 0.40 * G
MASS_PER_STORY = 1000.0
STORY_STIFFNESS = 500.0e3
ZETA = 0.05


def rayleigh_coefficients():
    m = MASS_PER_STORY * np.eye(3)
    k = STORY_STIFFNESS
    k_mat = np.array(
        [[2.0 * k, -k, 0.0], [-k, 2.0 * k, -k], [0.0, -k, k]],
        dtype=float,
    )
    eigvals = np.linalg.eigvals(np.linalg.solve(m, k_mat))
    omega = np.sort(np.sqrt(np.real(eigvals)))
    a = np.array(
        [[1.0 / (2.0 * omega[0]), omega[0] / 2.0],
         [1.0 / (2.0 * omega[1]), omega[1] / 2.0]],
        dtype=float,
    )
    alpha_m, beta_k = np.linalg.solve(a, np.array([ZETA, ZETA], dtype=float))
    return float(alpha_m), float(beta_k), omega


def main():
    raw = np.loadtxt(GM_PATH, dtype=float).reshape(-1)
    raw_peak = float(np.max(np.abs(raw)))
    if raw_peak <= 0.0:
        raise ValueError("GM1.txt has zero peak acceleration.")
    ag = raw / raw_peak * TARGET_PGA
    n = int(ag.size)
    time = np.arange(n, dtype=float) * DT
    input_pga = float(np.max(np.abs(ag)))
    eq_force = -MASS_PER_STORY * ag
    total_time = float(time[-1])

    alpha_m, beta_k, omega = rayleigh_coefficients()
    c_story = beta_k * STORY_STIFFNESS
    c_mass = alpha_m * MASS_PER_STORY

    np.savetxt(SCRIPT_DIR / "normalized_ground_acceleration.txt", np.column_stack([time, ag]), fmt="%.12e")
    np.savetxt(SCRIPT_DIR / "abaqus_eq_force_N.txt", np.column_stack([time, eq_force]), fmt="%.12e")

    lines = []
    lines.extend([
        "*Heading",
        "Three-story linear shear building; relative-coordinate earthquake loading",
        "** Units: kg, m, s, N",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
        "1, 0., 0., 0.",
        "2, 0., 1., 0.",
        "3, 0., 2., 0.",
        "999, 0., -1., 0.",
        "*Nset, nset=FLOORS",
        "1, 2, 3",
        "*Nset, nset=TOP",
        "3",
        "*Nset, nset=GROUND",
        "999",
        "*Element, type=MASS, elset=MASS1",
        "101, 1",
        "*Element, type=MASS, elset=MASS2",
        "102, 2",
        "*Element, type=MASS, elset=MASS3",
        "103, 3",
        "*Mass, elset=MASS1",
        "{:.12e}".format(MASS_PER_STORY),
        "*Mass, elset=MASS2",
        "{:.12e}".format(MASS_PER_STORY),
        "*Mass, elset=MASS3",
        "{:.12e}".format(MASS_PER_STORY),
        "*Element, type=SPRING2, elset=SPRING_STORY1",
        "201, 999, 1",
        "*Element, type=SPRING2, elset=SPRING_STORY2",
        "202, 1, 2",
        "*Element, type=SPRING2, elset=SPRING_STORY3",
        "203, 2, 3",
        "*Spring, elset=SPRING_STORY1",
        "1, 1",
        "{:.12e}".format(STORY_STIFFNESS),
        "*Spring, elset=SPRING_STORY2",
        "1, 1",
        "{:.12e}".format(STORY_STIFFNESS),
        "*Spring, elset=SPRING_STORY3",
        "1, 1",
        "{:.12e}".format(STORY_STIFFNESS),
        "*Element, type=DASHPOT2, elset=DASHPOT_STORY1",
        "301, 999, 1",
        "*Element, type=DASHPOT2, elset=DASHPOT_STORY2",
        "302, 1, 2",
        "*Element, type=DASHPOT2, elset=DASHPOT_STORY3",
        "303, 2, 3",
        "*Dashpot, elset=DASHPOT_STORY1",
        "1, 1",
        "{:.12e}".format(c_story),
        "*Dashpot, elset=DASHPOT_STORY2",
        "1, 1",
        "{:.12e}".format(c_story),
        "*Dashpot, elset=DASHPOT_STORY3",
        "1, 1",
        "{:.12e}".format(c_story),
        "*Element, type=DASHPOT2, elset=DASHPOT_MASS1",
        "401, 999, 1",
        "*Element, type=DASHPOT2, elset=DASHPOT_MASS2",
        "402, 999, 2",
        "*Element, type=DASHPOT2, elset=DASHPOT_MASS3",
        "403, 999, 3",
        "*Dashpot, elset=DASHPOT_MASS1",
        "1, 1",
        "{:.12e}".format(c_mass),
        "*Dashpot, elset=DASHPOT_MASS2",
        "1, 1",
        "{:.12e}".format(c_mass),
        "*Dashpot, elset=DASHPOT_MASS3",
        "1, 1",
        "{:.12e}".format(c_mass),
        "*Boundary",
        "999, 1, 1, 0.",
        "*Amplitude, name=EQFORCE, definition=TABULAR, time=TOTAL TIME, smooth=0.",
    ])

    for t, f in zip(time, eq_force):
        lines.append("{:.12e}, {:.12e}".format(float(t), float(f)))

    lines.extend([
        "*Step, name=EQ, nlgeom=NO, inc=3000",
        "*Dynamic, direct, alpha=0.",
        "{:.12e}, {:.12e}".format(DT, total_time),
        "*Cload, amplitude=EQFORCE",
        "1, 1, 1.",
        "2, 1, 1.",
        "3, 1, 1.",
        "*Output, history, frequency=1",
        "*Node Output, nset=TOP",
        "U1, V1, A1",
        "*End Step",
    ])

    (SCRIPT_DIR / "shear_building_abaqus.inp").write_text("\n".join(lines) + "\n", encoding="ascii")

    metadata = {
        "dt": DT,
        "num_samples": n,
        "input_pga_after_normalization": input_pga,
        "mass_per_story_kg": MASS_PER_STORY,
        "story_stiffness_N_per_m": STORY_STIFFNESS,
        "damping_ratio": ZETA,
        "rayleigh_alpha_m": alpha_m,
        "rayleigh_beta_k": beta_k,
        "dashpot_story_coeff_N_s_per_m": c_story,
        "dashpot_mass_coeff_N_s_per_m": c_mass,
        "circular_frequencies_rad_s": omega.tolist(),
        "forcing": "CLOAD amplitude contains -mass_per_story * normalized_ground_acceleration for each floor node",
    }
    (SCRIPT_DIR / "abaqus_model_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    print("Abaqus input created: shear_building_abaqus.inp")
    print("Samples: {}, dt: {:.8f}, total time: {:.8f}, normalized PGA: {:.8f} m/s^2".format(n, DT, total_time, input_pga))
    print("Rayleigh alpha_m: {:.12g}, beta_k: {:.12g}".format(alpha_m, beta_k))
    print("Dashpot coefficients: story {:.12g}, mass {:.12g}".format(c_story, c_mass))


if __name__ == "__main__":
    main()
