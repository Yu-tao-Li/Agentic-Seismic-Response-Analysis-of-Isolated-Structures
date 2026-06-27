```python
import os
import json
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def make_time_vector(dt, duration):
    n_steps = int(round(duration / dt)) + 1
    return np.linspace(0.0, duration, n_steps)


def ground_acceleration(t, g=9.81):
    a1 = 0.30 * g * np.sin(2 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
    a2 = 0.10 * g * np.sin(2 * np.pi * 3.0 * t) * np.exp(-0.20 * t)
    return a1 + a2


def compute_damping(m, k, zeta):
    omega_n = np.sqrt(k / m)
    c = 2.0 * zeta * m * omega_n
    return c, omega_n


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    beta = 1.0 / 4.0
    gamma = 1.0 / 2.0

    n = len(ag)
    u = np.zeros(n)
    v = np.zeros(n)
    a_rel = np.zeros(n)

    p = -m * ag

    a0 = (p[0] - c * v0 - k * u0) / m
    u[0] = u0
    v[0] = v0
    a_rel[0] = a0

    a1 = 1.0 / (beta * dt * dt) * m + gamma / (beta * dt) * c
    a2 = 1.0 / (beta * dt) * m + (gamma / beta - 1.0) * c
    a3 = (1.0 / (2.0 * beta) - 1.0) * m + dt * (gamma / (2.0 * beta) - 1.0) * c
    k_eff = k + a1

    for i in range(n - 1):
        dp = p[i + 1] + a1 * u[i] + a2 * v[i] + a3 * a_rel[i]
        u[i + 1] = dp / k_eff
        v[i + 1] = gamma / (beta * dt) * (u[i + 1] - u[i]) + (1.0 - gamma / beta) * v[i] + dt * (1.0 - gamma / (2.0 * beta)) * a_rel[i]
        a_rel[i + 1] = 1.0 / (beta * dt * dt) * (u[i + 1] - u[i]) - 1.0 / (beta * dt) * v[i] - (1.0 / (2.0 * beta) - 1.0) * a_rel[i]

    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    c, omega_n = compute_damping(m, k, zeta)
    T_n = 2.0 * np.pi / omega_n

    peak_disp = np.max(np.abs(u))
    peak_vel = np.max(np.abs(v))
    peak_arel = np.max(np.abs(a_rel))
    peak_aabs = np.max(np.abs(a_abs))
    peak_shear = k * peak_disp
    idx_peak = np.argmax(np.abs(u))
    time_peak = t[idx_peak]

    metrics = {
        "natural_circular_frequency_rad_s": float(omega_n),
        "natural_period_s": float(T_n),
        "damping_coefficient_Ns_m": float(c),
        "peak_relative_displacement_m": float(peak_disp),
        "peak_relative_velocity_m_s": float(peak_vel),
        "peak_relative_acceleration_m_s2": float(peak_arel),
        "peak_absolute_acceleration_m_s2": float(peak_aabs),
        "peak_pseudo_base_shear_N": float(peak_shear),
        "time_of_peak_displacement_s": float(time_peak),
    }
    return metrics


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    header = "time,ag,u,v,a_rel,a_abs"
    data = np.column_stack((t, ag, u, v, a_rel, a_abs))
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.10e")


def save_summary_json(path, metrics):
    with open(path, "w") as f:
        json.dump(metrics, f, indent=2)


def plot_response(path, t, ag, u, v, a_abs):
    fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

    axes[0].plot(t, ag, linewidth=0.8)
    axes[0].set_ylabel("ag (m/s^2)")
    axes[0].set_title("Ground Acceleration")
    axes[0].grid(True)

    axes[1].plot(t, u, linewidth=0.8, color="tab:blue")
    axes[1].set_ylabel("u (m)")
    axes[1].set_title("Relative Displacement")
    axes[1].grid(True)

    axes[2].plot(t, v, linewidth=0.8, color="tab:orange")
    axes[2].set_ylabel("v (m/s)")
    axes[2].set_title("Relative Velocity")
    axes[2].grid(True)

    axes[3].plot(t, a_abs, linewidth=0.8, color="tab:red")
    axes[3].set_ylabel("a_abs (m/s^2)")
    axes[3].set_title("Absolute Acceleration")
    axes[3].set_xlabel("Time (s)")
    axes[3].grid(True)

    plt.tight_layout()
    plt.savefig(path, dpi=150)
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

    os.makedirs(output_dir, exist_ok=True)

    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t, g)
    c, omega_n = compute_damping(m, k, zeta)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt, u0, v0)
    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)

    csv_path = os.path.join(output_dir, "response.csv")
    json_path = os.path.join(output_dir, "summary.json")
    png_path = os.path.join(output_dir, "response.png")

    save_csv(csv_path, t, ag, u, v, a_rel, a_abs)
    save_summary_json(json_path, metrics)
    plot_response(png_path, t, ag, u, v, a_abs)

    return metrics


def main():
    metrics = run_analysis("sdof_output")
    print("=== SDOF Newmark Analysis Summary ===")
    for key, val in metrics.items():
        print(f"  {key}: {val}")


def test_zero_ground_motion():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    dt = 0.01
    duration = 1.0
    c, omega_n = compute_damping(m, k, zeta)
    n = int(round(duration / dt)) + 1
    ag = np.zeros(n)
    u0, v0 = 0.01, 0.0
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt, u0, v0)
    assert len(u) == n
    assert abs(u[0] - u0) < 1e-12
    assert abs(v[0] - v0) < 1e-12
    max_u = np.max(np.abs(u))
    assert max_u <= abs(u0) + 1e-10, f"Displacement should not grow: max={max_u}, u0={u0}"


def test_output_shapes():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    dt = 0.01
    duration = 5.0
    c, omega_n = compute_damping(m, k, zeta)
    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt)
    n = len(t)
    assert len(u) == n
    assert len(v) == n
    assert len(a_rel) == n
    assert len(a_abs) == n
    assert len(ag) == n
    assert abs(t[-1] - duration) < 1e-10


def test_damping_positive():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    c, omega_n = compute_damping(m, k, zeta)
    assert c > 0
    assert omega_n > 0
    expected_c = 2.0 * zeta * np.sqrt(k * m)
    assert abs(c - expected_c) < 1e-10


if __name__ == "__main__":
    test_zero_ground_motion()
    test_output_shapes()
    test_damping_positive()
    print("All tests passed.")
    main()
```