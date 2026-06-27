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
    if m <= 0.0:
        raise ValueError("mass must be positive")
    if k <= 0.0:
        raise ValueError("stiffness must be positive")
    if zeta < 0.0:
        raise ValueError("damping ratio must be nonnegative")

    natural_circular_frequency = math.sqrt(k / m)
    return 2.0 * zeta * m * natural_circular_frequency


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
    if not np.all(np.isfinite(ag)):
        raise ValueError("ag must contain only finite values")

    beta = 0.25
    gamma = 0.50
    load = -m * ag

    u = np.zeros_like(ag)
    v = np.zeros_like(ag)
    a_rel = np.zeros_like(ag)

    u[0] = u0
    v[0] = v0
    a_rel[0] = (load[0] - c * v0 - k * u0) / m

    effective_mass = m + gamma * dt * c + beta * dt**2 * k

    for i in range(ag.size - 1):
        u_predictor = (
            u[i]
            + dt * v[i]
            + dt**2 * (0.5 - beta) * a_rel[i]
        )
        v_predictor = v[i] + dt * (1.0 - gamma) * a_rel[i]

        a_rel[i + 1] = (
            load[i + 1] - c * v_predictor - k * u_predictor
        ) / effective_mass
        u[i + 1] = u_predictor + beta * dt**2 * a_rel[i + 1]
        v[i + 1] = v_predictor + gamma * dt * a_rel[i + 1]

    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    arrays = [np.asarray(x, dtype=float) for x in (t, ag, u, v, a_rel, a_abs)]
    if any(x.ndim != 1 for x in arrays):
        raise ValueError("all response arrays must be one-dimensional")
    if len({x.size for x in arrays}) != 1 or arrays[0].size == 0:
        raise ValueError("all response arrays must have the same nonzero length")

    t, ag, u, v, a_rel, a_abs = arrays
    peak_displacement_index = int(np.argmax(np.abs(u)))

    return {
        "peak_ground_acceleration_m_per_s2": float(np.max(np.abs(ag))),
        "peak_relative_displacement_m": float(np.max(np.abs(u))),
        "peak_relative_velocity_m_per_s": float(np.max(np.abs(v))),
        "peak_relative_acceleration_m_per_s2": float(np.max(np.abs(a_rel))),
        "peak_absolute_acceleration_m_per_s2": float(np.max(np.abs(a_abs))),
        "peak_pseudo_base_shear_N": float(np.max(np.abs(STIFFNESS * u))),
        "time_of_peak_displacement_s": float(t[peak_displacement_index]),
    }


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    path = Path(path)
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "time_s",
                "ground_acceleration_m_per_s2",
                "relative_displacement_m",
                "relative_velocity_m_per_s",
                "relative_acceleration_m_per_s2",
                "absolute_acceleration_m_per_s2",
            ]
        )
        writer.writerows(zip(t, ag, u, v, a_rel, a_abs))


def save_summary_json(path, metrics):
    path = Path(path)
    with path.open("w", encoding="utf-8") as file:
        json.dump(metrics, file, indent=2, sort_keys=True)
        file.write("\n")


def plot_response(path, t, ag, u, v, a_abs):
    path = Path(path)
    figure, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)

    axes[0].plot(t, ag, color="tab:gray", linewidth=1.0)
    axes[0].set_ylabel("Ground accel.\n(m/s²)")

    axes[1].plot(t, u, color="tab:blue", linewidth=1.0)
    axes[1].set_ylabel("Relative disp.\n(m)")

    axes[2].plot(t, v, color="tab:green", linewidth=1.0)
    axes[2].set_ylabel("Relative vel.\n(m/s)")

    axes[3].plot(t, a_abs, color="tab:red", linewidth=1.0)
    axes[3].set_ylabel("Absolute accel.\n(m/s²)")
    axes[3].set_xlabel("Time (s)")

    for axis in axes:
        axis.grid(True, alpha=0.3)

    figure.suptitle("SDOF Response: Newmark Average Acceleration")
    figure.tight_layout()
    figure.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(figure)


def run_analysis(output_dir="sdof_output"):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

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

    natural_circular_frequency = math.sqrt(STIFFNESS / MASS)
    natural_period = 2.0 * math.pi / natural_circular_frequency

    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)
    metrics.update(
        {
            "mass_kg": MASS,
            "stiffness_N_per_m": STIFFNESS,
            "damping_ratio": DAMPING_RATIO,
            "damping_coefficient_N_s_per_m": damping,
            "natural_circular_frequency_rad_per_s": natural_circular_frequency,
            "natural_period_s": natural_period,
            "time_step_s": TIME_STEP,
            "duration_s": DURATION,
        }
    )

    save_csv(output_dir / "response.csv", t, ag, u, v, a_rel, a_abs)
    save_summary_json(output_dir / "summary.json", metrics)
    plot_response(output_dir / "response.png", t, ag, u, v, a_abs)

    return metrics


def main():
    metrics = run_analysis()

    labels = [
        ("Natural circular frequency", "natural_circular_frequency_rad_per_s", "rad/s"),
        ("Natural period", "natural_period_s", "s"),
        ("Damping coefficient", "damping_coefficient_N_s_per_m", "N·s/m"),
        ("Peak relative displacement", "peak_relative_displacement_m", "m"),
        ("Peak relative velocity", "peak_relative_velocity_m_per_s", "m/s"),
        ("Peak relative acceleration", "peak_relative_acceleration_m_per_s2", "m/s²"),
        ("Peak absolute acceleration", "peak_absolute_acceleration_m_per_s2", "m/s²"),
        ("Peak pseudo base shear", "peak_pseudo_base_shear_N", "N"),
        ("Time of peak displacement", "time_of_peak_displacement_s", "s"),
    ]

    print("SDOF Newmark average-acceleration analysis")
    for label, key, unit in labels:
        print(f"{label}: {metrics[key]:.8g} {unit}")


def test_zero_ground_motion():
    ag = np.zeros(101)
    u, v, a_rel, a_abs = newmark_average_acceleration(
        MASS, compute_damping(MASS, STIFFNESS, DAMPING_RATIO),
        STIFFNESS, ag, TIME_STEP
    )
    assert np.allclose(u, 0.0)
    assert np.allclose(v, 0.0)
    assert np.allclose(a_rel, 0.0)
    assert np.allclose(a_abs, 0.0)


def test_output_shapes():
    t = make_time_vector(TIME_STEP, 1.0)
    ag = ground_acceleration(t)
    outputs = newmark_average_acceleration(
        MASS, compute_damping(MASS, STIFFNESS, DAMPING_RATIO),
        STIFFNESS, ag, TIME_STEP
    )
    assert all(array.shape == t.shape for array in outputs)


def test_damping_positive():
    damping = compute_damping(MASS, STIFFNESS, DAMPING_RATIO)
    assert damping > 0.0


if __name__ == "__main__":
    main()
