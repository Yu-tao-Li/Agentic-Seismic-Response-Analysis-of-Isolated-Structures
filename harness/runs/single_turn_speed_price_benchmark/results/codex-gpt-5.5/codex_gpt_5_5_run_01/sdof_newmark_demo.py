import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


MASS = 1000.0
STIFFNESS = 200000.0
DAMPING_RATIO = 0.05
GRAVITY = 9.81
TIME_STEP = 0.01
DURATION = 20.0
INITIAL_DISPLACEMENT = 0.0
INITIAL_VELOCITY = 0.0


def make_time_vector(dt, duration):
    if dt <= 0.0:
        raise ValueError("dt must be positive")
    if duration < 0.0:
        raise ValueError("duration must be nonnegative")

    number_of_steps = int(round(duration / dt))
    if not math.isclose(number_of_steps * dt, duration, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("duration must be an integer multiple of dt")

    return np.arange(number_of_steps + 1, dtype=float) * dt


def ground_acceleration(t, g=9.81):
    t = np.asarray(t, dtype=float)
    return (
        0.30 * g * np.sin(2.0 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
        + 0.10 * g * np.sin(2.0 * np.pi * 3.0 * t) * np.exp(-0.20 * t)
    )


def compute_damping(m, k, zeta):
    if m <= 0.0 or k <= 0.0:
        raise ValueError("mass and stiffness must be positive")
    if zeta < 0.0:
        raise ValueError("damping ratio must be nonnegative")

    return 2.0 * zeta * math.sqrt(m * k)


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    if m <= 0.0:
        raise ValueError("mass must be positive")
    if c < 0.0:
        raise ValueError("damping coefficient must be nonnegative")
    if k <= 0.0:
        raise ValueError("stiffness must be positive")
    if dt <= 0.0:
        raise ValueError("dt must be positive")

    ag = np.asarray(ag, dtype=float)
    if ag.ndim != 1 or ag.size == 0:
        raise ValueError("ag must be a nonempty one-dimensional array")

    beta = 0.25
    gamma = 0.5
    n = ag.size

    u = np.zeros(n, dtype=float)
    v = np.zeros(n, dtype=float)
    a_rel = np.zeros(n, dtype=float)

    u[0] = float(u0)
    v[0] = float(v0)
    a_rel[0] = (-m * ag[0] - c * v[0] - k * u[0]) / m

    effective_mass = m + gamma * dt * c + beta * dt * dt * k

    for i in range(n - 1):
        u_predictor = (
            u[i]
            + dt * v[i]
            + dt * dt * (0.5 - beta) * a_rel[i]
        )
        v_predictor = v[i] + dt * (1.0 - gamma) * a_rel[i]

        a_rel[i + 1] = (
            -m * ag[i + 1]
            - c * v_predictor
            - k * u_predictor
        ) / effective_mass

        u[i + 1] = u_predictor + beta * dt * dt * a_rel[i + 1]
        v[i + 1] = v_predictor + gamma * dt * a_rel[i + 1]

    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    arrays = [
        np.asarray(values, dtype=float)
        for values in (t, ag, u, v, a_rel, a_abs)
    ]

    if any(values.ndim != 1 for values in arrays):
        raise ValueError("all response inputs must be one-dimensional")
    if len({values.size for values in arrays}) != 1:
        raise ValueError("all response inputs must have the same length")
    if arrays[0].size == 0:
        raise ValueError("response inputs must not be empty")

    t, ag, u, v, a_rel, a_abs = arrays
    omega_n = math.sqrt(STIFFNESS / MASS)
    natural_period = 2.0 * math.pi / omega_n
    damping_coefficient = compute_damping(
        MASS, STIFFNESS, DAMPING_RATIO
    )
    peak_displacement_index = int(np.argmax(np.abs(u)))

    return {
        "natural_circular_frequency_rad_s": float(omega_n),
        "natural_period_s": float(natural_period),
        "damping_coefficient_N_s_m": float(damping_coefficient),
        "peak_relative_displacement_m": float(np.max(np.abs(u))),
        "peak_relative_velocity_m_s": float(np.max(np.abs(v))),
        "peak_relative_acceleration_m_s2": float(np.max(np.abs(a_rel))),
        "peak_absolute_acceleration_m_s2": float(np.max(np.abs(a_abs))),
        "peak_pseudo_base_shear_N": float(
            np.max(np.abs(STIFFNESS * u))
        ),
        "time_of_peak_displacement_s": float(t[peak_displacement_index]),
    }


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    path = Path(path)
    arrays = [
        np.asarray(values, dtype=float)
        for values in (t, ag, u, v, a_rel, a_abs)
    ]

    if len({values.size for values in arrays}) != 1:
        raise ValueError("all CSV columns must have the same length")

    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "time_s",
                "ground_acceleration_m_s2",
                "relative_displacement_m",
                "relative_velocity_m_s",
                "relative_acceleration_m_s2",
                "absolute_acceleration_m_s2",
            ]
        )
        writer.writerows(zip(*arrays))


def save_summary_json(path, metrics):
    path = Path(path)
    with path.open("w", encoding="utf-8") as file:
        json.dump(metrics, file, indent=2, sort_keys=True)
        file.write("\n")


def plot_response(path, t, ag, u, v, a_abs):
    figure, axes = plt.subplots(
        4, 1, figsize=(10, 11), sharex=True, constrained_layout=True
    )

    axes[0].plot(t, ag, color="tab:gray", linewidth=1.0)
    axes[0].set_ylabel("Ground accel.\n(m/s²)")
    axes[0].set_title("SDOF Newmark Average-Acceleration Response")

    axes[1].plot(t, u, color="tab:blue", linewidth=1.0)
    axes[1].set_ylabel("Relative disp.\n(m)")

    axes[2].plot(t, v, color="tab:green", linewidth=1.0)
    axes[2].set_ylabel("Relative vel.\n(m/s)")

    axes[3].plot(t, a_abs, color="tab:red", linewidth=1.0)
    axes[3].set_ylabel("Absolute accel.\n(m/s²)")
    axes[3].set_xlabel("Time (s)")

    for axis in axes:
        axis.grid(True, alpha=0.3)

    figure.savefig(path, dpi=180)
    plt.close(figure)


def run_analysis(output_dir="sdof_output"):
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    t = make_time_vector(TIME_STEP, DURATION)
    ag = ground_acceleration(t, GRAVITY)
    damping = compute_damping(MASS, STIFFNESS, DAMPING_RATIO)

    u, v, a_rel, a_abs = newmark_average_acceleration(
        MASS,
        damping,
        STIFFNESS,
        ag,
        TIME_STEP,
        INITIAL_DISPLACEMENT,
        INITIAL_VELOCITY,
    )

    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)

    save_csv(
        output_path / "response.csv",
        t,
        ag,
        u,
        v,
        a_rel,
        a_abs,
    )
    save_summary_json(output_path / "summary.json", metrics)
    plot_response(output_path / "response.png", t, ag, u, v, a_abs)

    labels = {
        "natural_circular_frequency_rad_s": "Natural circular frequency (rad/s)",
        "natural_period_s": "Natural period (s)",
        "damping_coefficient_N_s_m": "Damping coefficient (N·s/m)",
        "peak_relative_displacement_m": "Peak relative displacement (m)",
        "peak_relative_velocity_m_s": "Peak relative velocity (m/s)",
        "peak_relative_acceleration_m_s2": (
            "Peak relative acceleration (m/s²)"
        ),
        "peak_absolute_acceleration_m_s2": (
            "Peak absolute acceleration (m/s²)"
        ),
        "peak_pseudo_base_shear_N": "Peak pseudo base shear (N)",
        "time_of_peak_displacement_s": "Time of peak displacement (s)",
    }

    print("SDOF Newmark average-acceleration analysis")
    for key, value in metrics.items():
        print(f"{labels[key]}: {value:.8g}")
    print(f"Output directory: {output_path.resolve()}")

    return metrics


def main():
    run_analysis()


def test_zero_ground_motion():
    ag = np.zeros(101, dtype=float)
    u, v, a_rel, a_abs = newmark_average_acceleration(
        MASS,
        compute_damping(MASS, STIFFNESS, DAMPING_RATIO),
        STIFFNESS,
        ag,
        TIME_STEP,
    )
    assert np.allclose(u, 0.0)
    assert np.allclose(v, 0.0)
    assert np.allclose(a_rel, 0.0)
    assert np.allclose(a_abs, 0.0)


def test_output_shapes():
    t = make_time_vector(TIME_STEP, 1.0)
    ag = ground_acceleration(t)
    u, v, a_rel, a_abs = newmark_average_acceleration(
        MASS,
        compute_damping(MASS, STIFFNESS, DAMPING_RATIO),
        STIFFNESS,
        ag,
        TIME_STEP,
    )
    assert t.shape == ag.shape == u.shape == v.shape
    assert v.shape == a_rel.shape == a_abs.shape


def test_damping_positive():
    damping = compute_damping(MASS, STIFFNESS, DAMPING_RATIO)
    assert damping > 0.0
    assert math.isclose(
        damping,
        2.0 * DAMPING_RATIO * math.sqrt(MASS * STIFFNESS),
    )


if __name__ == "__main__":
    main()
