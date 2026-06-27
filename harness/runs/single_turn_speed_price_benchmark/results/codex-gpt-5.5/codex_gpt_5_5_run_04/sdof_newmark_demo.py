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
    return np.linspace(0.0, number_of_steps * dt, number_of_steps + 1)


def ground_acceleration(t, g=9.81):
    t = np.asarray(t, dtype=float)
    return (
        0.30 * g * np.sin(2.0 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
        + 0.10 * g * np.sin(2.0 * np.pi * 3.0 * t) * np.exp(-0.20 * t)
    )


def compute_damping(m, k, zeta):
    if m <= 0.0:
        raise ValueError("Mass must be positive")
    if k <= 0.0:
        raise ValueError("Stiffness must be positive")
    if zeta < 0.0:
        raise ValueError("Damping ratio must be nonnegative")

    natural_circular_frequency = math.sqrt(k / m)
    return 2.0 * zeta * m * natural_circular_frequency


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    if m <= 0.0:
        raise ValueError("Mass must be positive")
    if c < 0.0:
        raise ValueError("Damping coefficient must be nonnegative")
    if k <= 0.0:
        raise ValueError("Stiffness must be positive")
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

    u[0] = u0
    v[0] = v0
    a_rel[0] = (-m * ag[0] - c * v0 - k * u0) / m

    a0 = 1.0 / (beta * dt * dt)
    a1 = gamma / (beta * dt)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2.0 * beta) - 1.0
    a4 = gamma / beta - 1.0
    a5 = dt * (gamma / (2.0 * beta) - 1.0)

    effective_stiffness = k + a0 * m + a1 * c
    force = -m * ag

    for i in range(1, n):
        effective_force = (
            force[i]
            + m * (a0 * u[i - 1] + a2 * v[i - 1] + a3 * a_rel[i - 1])
            + c * (a1 * u[i - 1] + a4 * v[i - 1] + a5 * a_rel[i - 1])
        )

        u[i] = effective_force / effective_stiffness
        a_rel[i] = (
            a0 * (u[i] - u[i - 1])
            - a2 * v[i - 1]
            - a3 * a_rel[i - 1]
        )
        v[i] = v[i - 1] + dt * (
            (1.0 - gamma) * a_rel[i - 1] + gamma * a_rel[i]
        )

    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    arrays = [
        np.asarray(values, dtype=float)
        for values in (t, ag, u, v, a_rel, a_abs)
    ]

    if any(values.ndim != 1 for values in arrays):
        raise ValueError("All response arrays must be one-dimensional")
    if len({values.size for values in arrays}) != 1:
        raise ValueError("All response arrays must have equal lengths")
    if arrays[0].size == 0:
        raise ValueError("Response arrays must not be empty")

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
    arrays = [np.asarray(values, dtype=float) for values in (t, ag, u, v, a_rel, a_abs)]
    if len({values.size for values in arrays}) != 1:
        raise ValueError("All CSV arrays must have equal lengths")

    with open(path, "w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
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
        writer.writerows(zip(*arrays))


def save_summary_json(path, metrics):
    with open(path, "w", encoding="utf-8") as json_file:
        json.dump(metrics, json_file, indent=2, sort_keys=True, allow_nan=False)
        json_file.write("\n")


def plot_response(path, t, ag, u, v, a_abs):
    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)

    axes[0].plot(t, ag / GRAVITY, color="black", linewidth=1.0)
    axes[0].set_ylabel("Ground accel. (g)")
    axes[0].set_title("SDOF Response Using Newmark Average Acceleration")

    axes[1].plot(t, 1000.0 * u, color="tab:blue", linewidth=1.0)
    axes[1].set_ylabel("Rel. disp. (mm)")

    axes[2].plot(t, v, color="tab:green", linewidth=1.0)
    axes[2].set_ylabel("Rel. velocity (m/s)")

    axes[3].plot(t, a_abs / GRAVITY, color="tab:red", linewidth=1.0)
    axes[3].set_ylabel("Abs. accel. (g)")
    axes[3].set_xlabel("Time (s)")

    for axis in axes:
        axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.6)

    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


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

    natural_circular_frequency = math.sqrt(STIFFNESS / MASS)
    natural_period = 2.0 * math.pi / natural_circular_frequency

    metrics = {
        "mass_kg": MASS,
        "stiffness_N_per_m": STIFFNESS,
        "damping_ratio": DAMPING_RATIO,
        "damping_coefficient_N_s_per_m": damping,
        "natural_circular_frequency_rad_per_s": natural_circular_frequency,
        "natural_period_s": natural_period,
        "time_step_s": TIME_STEP,
        "duration_s": DURATION,
    }
    metrics.update(compute_response_metrics(t, ag, u, v, a_rel, a_abs))

    save_csv(output_path / "response.csv", t, ag, u, v, a_rel, a_abs)
    save_summary_json(output_path / "summary.json", metrics)
    plot_response(output_path / "response.png", t, ag, u, v, a_abs)

    return metrics


def main():
    metrics = run_analysis()
    for name, value in metrics.items():
        print(f"{name}: {value:.10g}" if isinstance(value, float) else f"{name}: {value}")


def test_zero_ground_motion():
    ag = np.zeros(101)
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
    outputs = newmark_average_acceleration(
        MASS,
        compute_damping(MASS, STIFFNESS, DAMPING_RATIO),
        STIFFNESS,
        ag,
        TIME_STEP,
    )
    assert t.shape == ag.shape
    assert all(result.shape == t.shape for result in outputs)


def test_damping_positive():
    damping = compute_damping(MASS, STIFFNESS, DAMPING_RATIO)
    assert damping > 0.0


if __name__ == "__main__":
    main()
