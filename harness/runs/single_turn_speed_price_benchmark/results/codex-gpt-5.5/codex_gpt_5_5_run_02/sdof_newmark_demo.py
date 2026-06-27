import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def make_time_vector(dt, duration):
    if dt <= 0.0:
        raise ValueError("dt must be positive.")
    if duration < 0.0:
        raise ValueError("duration must be non-negative.")

    steps = int(round(duration / dt))
    if not math.isclose(steps * dt, duration, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("duration must be an integer multiple of dt.")

    return np.arange(steps + 1, dtype=float) * dt


def ground_acceleration(t, g=9.81):
    t = np.asarray(t, dtype=float)
    return (
        0.30 * g * np.sin(2.0 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
        + 0.10 * g * np.sin(2.0 * np.pi * 3.0 * t) * np.exp(-0.20 * t)
    )


def compute_damping(m, k, zeta):
    if m <= 0.0:
        raise ValueError("Mass must be positive.")
    if k <= 0.0:
        raise ValueError("Stiffness must be positive.")
    if zeta < 0.0:
        raise ValueError("Damping ratio cannot be negative.")

    omega_n = math.sqrt(k / m)
    return 2.0 * zeta * m * omega_n


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    if m <= 0.0:
        raise ValueError("Mass must be positive.")
    if c < 0.0:
        raise ValueError("Damping coefficient cannot be negative.")
    if k <= 0.0:
        raise ValueError("Stiffness must be positive.")
    if dt <= 0.0:
        raise ValueError("dt must be positive.")

    ag = np.asarray(ag, dtype=float)
    if ag.ndim != 1 or ag.size == 0:
        raise ValueError("ag must be a non-empty one-dimensional array.")

    beta = 0.25
    gamma = 0.5
    n = ag.size

    u = np.zeros(n, dtype=float)
    v = np.zeros(n, dtype=float)
    a_rel = np.zeros(n, dtype=float)

    u[0] = float(u0)
    v[0] = float(v0)

    force = -m * ag
    a_rel[0] = (force[0] - c * v[0] - k * u[0]) / m

    a0 = 1.0 / (beta * dt**2)
    a1 = gamma / (beta * dt)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2.0 * beta) - 1.0
    a4 = gamma / beta - 1.0
    a5 = dt * (gamma / (2.0 * beta) - 1.0)

    effective_stiffness = k + a0 * m + a1 * c

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

    return u, v, a_rel


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    arrays = [
        np.asarray(values, dtype=float)
        for values in (t, ag, u, v, a_rel, a_abs)
    ]

    if any(values.ndim != 1 for values in arrays):
        raise ValueError("All response inputs must be one-dimensional.")
    if not arrays[0].size:
        raise ValueError("Response arrays cannot be empty.")
    if any(values.shape != arrays[0].shape for values in arrays[1:]):
        raise ValueError("All response arrays must have the same shape.")

    t, ag, u, v, a_rel, a_abs = arrays
    peak_displacement_index = int(np.argmax(np.abs(u)))

    return {
        "peak_ground_acceleration_m_s2": float(np.max(np.abs(ag))),
        "peak_relative_displacement_m": float(np.max(np.abs(u))),
        "peak_relative_velocity_m_s": float(np.max(np.abs(v))),
        "peak_relative_acceleration_m_s2": float(np.max(np.abs(a_rel))),
        "peak_absolute_acceleration_m_s2": float(np.max(np.abs(a_abs))),
        "time_of_peak_displacement_s": float(t[peak_displacement_index]),
    }


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    arrays = [
        np.asarray(values, dtype=float)
        for values in (t, ag, u, v, a_rel, a_abs)
    ]
    if any(values.shape != arrays[0].shape for values in arrays[1:]):
        raise ValueError("All CSV arrays must have the same shape.")

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
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8") as file:
        json.dump(metrics, file, indent=2, sort_keys=True, allow_nan=False)
        file.write("\n")


def plot_response(path, t, ag, u, v, a_abs):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)

    axes[0].plot(t, ag, color="tab:gray", linewidth=1.0)
    axes[0].set_ylabel("Ground accel.\n(m/s²)")

    axes[1].plot(t, u, color="tab:blue", linewidth=1.0)
    axes[1].set_ylabel("Relative disp.\n(m)")

    axes[2].plot(t, v, color="tab:orange", linewidth=1.0)
    axes[2].set_ylabel("Relative vel.\n(m/s)")

    axes[3].plot(t, a_abs, color="tab:red", linewidth=1.0)
    axes[3].set_ylabel("Absolute accel.\n(m/s²)")
    axes[3].set_xlabel("Time (s)")

    for axis in axes:
        axis.grid(True, alpha=0.3)

    fig.suptitle("SDOF Newmark Average-Acceleration Response")
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def run_analysis(output_dir="sdof_output"):
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    g = 9.81
    dt = 0.01
    duration = 20.0
    u0 = 0.0
    v0 = 0.0

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t, g)
    c = compute_damping(m, k, zeta)
    omega_n = math.sqrt(k / m)
    natural_period = 2.0 * math.pi / omega_n

    u, v, a_rel = newmark_average_acceleration(
        m, c, k, ag, dt, u0=u0, v0=v0
    )
    a_abs = a_rel + ag

    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)
    metrics.update(
        {
            "natural_circular_frequency_rad_s": float(omega_n),
            "natural_period_s": float(natural_period),
            "damping_coefficient_n_s_m": float(c),
            "peak_pseudo_base_shear_n": float(np.max(np.abs(k * u))),
        }
    )

    save_csv(output_dir / "response.csv", t, ag, u, v, a_rel, a_abs)
    save_summary_json(output_dir / "summary.json", metrics)
    plot_response(output_dir / "response.png", t, ag, u, v, a_abs)

    print(f"Natural circular frequency: {omega_n:.6f} rad/s")
    print(f"Natural period: {natural_period:.6f} s")
    print(f"Damping coefficient: {c:.6f} N·s/m")
    print(
        "Peak relative displacement: "
        f"{metrics['peak_relative_displacement_m']:.6e} m"
    )
    print(
        "Peak relative velocity: "
        f"{metrics['peak_relative_velocity_m_s']:.6e} m/s"
    )
    print(
        "Peak relative acceleration: "
        f"{metrics['peak_relative_acceleration_m_s2']:.6e} m/s²"
    )
    print(
        "Peak absolute acceleration: "
        f"{metrics['peak_absolute_acceleration_m_s2']:.6e} m/s²"
    )
    print(
        "Peak pseudo base shear: "
        f"{metrics['peak_pseudo_base_shear_n']:.6e} N"
    )
    print(
        "Time of peak displacement: "
        f"{metrics['time_of_peak_displacement_s']:.6f} s"
    )

    return metrics


def main():
    run_analysis()


def test_zero_ground_motion():
    ag = np.zeros(101, dtype=float)
    u, v, a_rel = newmark_average_acceleration(
        m=1000.0,
        c=1000.0,
        k=200000.0,
        ag=ag,
        dt=0.01,
    )
    assert np.allclose(u, 0.0)
    assert np.allclose(v, 0.0)
    assert np.allclose(a_rel, 0.0)


def test_output_shapes():
    t = make_time_vector(0.01, 1.0)
    ag = ground_acceleration(t)
    c = compute_damping(1000.0, 200000.0, 0.05)
    u, v, a_rel = newmark_average_acceleration(
        1000.0, c, 200000.0, ag, 0.01
    )
    assert t.shape == ag.shape == u.shape == v.shape == a_rel.shape
    assert t.size == 101


def test_damping_positive():
    c = compute_damping(1000.0, 200000.0, 0.05)
    assert c > 0.0
    assert math.isclose(c, 2.0 * 0.05 * math.sqrt(1000.0 * 200000.0))


if __name__ == "__main__":
    main()
