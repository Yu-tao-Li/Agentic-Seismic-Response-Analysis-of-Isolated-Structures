import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_PATH = ROOT / "input_data" / "GM1.txt"
INP_PATH = SCRIPT_DIR / "isolated_building.inp"


def build_matrices():
    masses = np.array([500.0, 1000.0, 1000.0, 1000.0])
    stiffnesses = np.array([200e3, 500e3, 500e3, 500e3])
    k1, k2, k3, k4 = stiffnesses
    mass_matrix = np.diag(masses)
    stiffness_matrix = np.array(
        [
            [k1 + k2, -k2, 0.0, 0.0],
            [-k2, k2 + k3, -k3, 0.0],
            [0.0, -k3, k3 + k4, -k4],
            [0.0, 0.0, -k4, k4],
        ],
        dtype=float,
    )
    return masses, stiffnesses, mass_matrix, stiffness_matrix


def modal_rayleigh(mass_matrix, stiffness_matrix, damping_ratio):
    eigvals = np.linalg.eigvals(np.linalg.solve(mass_matrix, stiffness_matrix))
    omegas = np.sqrt(np.sort(np.real(eigvals)))
    system = np.array(
        [[1.0 / (2.0 * omegas[0]), omegas[0] / 2.0], [1.0 / (2.0 * omegas[1]), omegas[1] / 2.0]]
    )
    alpha_m, beta_k = np.linalg.solve(system, np.array([damping_ratio, damping_ratio]))
    return omegas, float(alpha_m), float(beta_k)


def amplitude_lines(time, ag):
    pairs = [f"{t:.12e}, {a:.12e}" for t, a in zip(time, ag)]
    lines = []
    for i in range(0, len(pairs), 4):
        lines.append(", ".join(pairs[i : i + 4]))
    return lines


def write_input_deck(time, ag, masses, stiffnesses, alpha_m, beta_k):
    c_iso = 2e3
    story_dashpots = beta_k * stiffnesses
    mass_dashpots = alpha_m * masses

    lines = [
        "*Heading",
        "Four-story isolated shear-building, SI units",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
        "1, 0., 0., 0.",
        "2, 0., 0., 0.",
        "3, 0., 0., 0.",
        "4, 0., 0., 0.",
        "5, 0., 0., 0.",
        "*Nset, nset=GROUND",
        "1",
        "*Nset, nset=FLOORS",
        "2, 3, 4, 5",
        "*Nset, nset=RESP_NODES",
        "2, 3, 4, 5",
    ]

    for i, mass in enumerate(masses, start=2):
        lines.extend(
            [
                f"*Element, type=MASS, elset=MASS_NODE_{i}",
                f"{100 + i}, {i}",
                f"*Mass, elset=MASS_NODE_{i}",
                f"{mass:.12e}",
            ]
        )

    element_id = 200
    for story, k in enumerate(stiffnesses, start=1):
        lower = story
        upper = story + 1
        lines.extend(
            [
                f"*Element, type=SPRING2, elset=SPRING_{story}",
                f"{element_id}, {lower}, {upper}",
                f"*Spring, elset=SPRING_{story}",
                "1, 1",
                f"{k:.12e}",
            ]
        )
        element_id += 1

    for story, c in enumerate(story_dashpots, start=1):
        lower = story
        upper = story + 1
        if story == 1:
            c += c_iso
        lines.extend(
            [
                f"*Element, type=DASHPOT2, elset=STORY_DASHPOT_{story}",
                f"{element_id}, {lower}, {upper}",
                f"*Dashpot, elset=STORY_DASHPOT_{story}",
                "1, 1",
                f"{c:.12e}",
            ]
        )
        element_id += 1

    for floor, c in enumerate(mass_dashpots, start=2):
        lines.extend(
            [
                f"*Element, type=DASHPOT2, elset=MASS_DASHPOT_{floor}",
                f"{element_id}, 1, {floor}",
                f"*Dashpot, elset=MASS_DASHPOT_{floor}",
                "1, 1",
                f"{c:.12e}",
            ]
        )
        element_id += 1

    lines.extend(
        [
            "*Boundary",
            "1, 1, 3",
            "2, 2, 3",
            "3, 2, 3",
            "4, 2, 3",
            "5, 2, 3",
            "*Amplitude, name=GM_ACCEL, definition=TABULAR, time=TOTAL TIME",
        ]
    )
    lines.extend(amplitude_lines(time, ag))
    lines.extend(
        [
            f"*Step, name=GM_STEP, nlgeom=NO, inc={len(time) + 10}",
            "*Dynamic, direct, alpha=0.",
            f"{time[1] - time[0]:.12e}, {time[-1]:.12e}",
            "*Cload, amplitude=GM_ACCEL",
        ]
    )
    for node, mass in zip(range(2, 6), masses):
        lines.append(f"{node}, 1, {-mass:.12e}")
    lines.extend(
        [
            "*Output, field, frequency=1",
            "*Node Output, nset=RESP_NODES",
            "U, V, A",
            "*End Step",
            "",
        ]
    )
    INP_PATH.write_text("\n".join(lines), encoding="utf-8")


def main():
    dt = 0.02
    g = 9.80665
    target_pga = 0.40 * g
    damping_ratio = 0.05

    raw_gm = np.loadtxt(INPUT_PATH, dtype=float).reshape(-1)
    raw_peak = float(np.max(np.abs(raw_gm)))
    if raw_peak <= 0.0:
        raise ValueError("GM1.txt has zero peak acceleration.")
    ag = raw_gm / raw_peak * target_pga
    time = np.arange(ag.size, dtype=float) * dt

    masses, stiffnesses, mass_matrix, stiffness_matrix = build_matrices()
    omegas, alpha_m, beta_k = modal_rayleigh(mass_matrix, stiffness_matrix, damping_ratio)

    write_input_deck(time, ag, masses, stiffnesses, alpha_m, beta_k)
    np.savetxt(SCRIPT_DIR / "normalized_ground_motion_mps2.txt", np.column_stack([time, ag]), fmt="%.12e", delimiter="\t")

    metadata = {
        "software": "ABAQUS",
        "input_file": str(INP_PATH),
        "source_ground_motion": str(INPUT_PATH),
        "time_step": dt,
        "number_of_samples": int(ag.size),
        "target_pga_mps2": target_pga,
        "pga_after_normalization_mps2": float(np.max(np.abs(ag))),
        "pga_after_normalization_g": float(np.max(np.abs(ag)) / g),
        "masses_kg": masses.tolist(),
        "stiffnesses_N_per_m": stiffnesses.tolist(),
        "isolation_damping_Ns_per_m": 2e3,
        "rayleigh_damping_ratio": damping_ratio,
        "rayleigh_alpha_mass": alpha_m,
        "rayleigh_beta_stiffness": beta_k,
        "modal_frequencies_rad_per_s": omegas.tolist(),
        "damping_realization": (
            "Rayleigh matrix implemented with physical dashpot elements: beta*K as story dashpots "
            "and alpha*M as ground dashpots, plus isolation dashpot."
        ),
    }
    with (SCRIPT_DIR / "model_metadata.json").open("w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    print(f"Wrote {INP_PATH.name}")
    print(f"Samples: {ag.size}, dt: {dt:.6f} s, normalized PGA: {np.max(np.abs(ag)):.12g} m/s^2")
    print(f"Rayleigh alpha={alpha_m:.12g}, beta={beta_k:.12g}")


if __name__ == "__main__":
    main()
