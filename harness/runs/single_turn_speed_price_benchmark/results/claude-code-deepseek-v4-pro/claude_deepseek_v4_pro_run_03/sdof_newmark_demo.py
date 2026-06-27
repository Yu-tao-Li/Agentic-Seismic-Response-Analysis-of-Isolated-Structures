```python
import os
import json
import numpy as np
import matplotlib.pyplot as plt


def make_time_vector(dt, duration):
    return np.arange(0, duration + dt, dt)


def ground_acceleration(t, g=9.81):
    return (0.30 * g * np.sin(2 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
            + 0.10 * g * np.sin(2 * np.pi * 3.0 * t) * np.exp(-0.20 * t))


def compute_damping(m, k, zeta):
    return 2 * zeta * np.sqrt(k * m)


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    beta = 0.25
    gamma = 0.5

    n = len(ag)
    u = np.zeros(n)
    v = np.zeros(n)
    a_rel = np.zeros(n)

    u[0] = u0
    v[0] = v0
    a_rel[0] = (-m * ag[0] - c * v[0] - k * u[0]) / m

    a1 = m / (beta * dt**2) + gamma * c / (beta * dt)
    a2 = m / (beta * dt) + (gamma / beta - 1) * c
    a3 = (1 / (2 * beta) - 1) * m + dt * (gamma / (2 * beta) - 1) * c

    k_eff = k + a1

    for i in range(n - 1):
        dp = -m * ag[i + 1] + a1 * u[i] + a2 * v[i] + a3 * a_rel[i]
        du = dp / k_eff
        dv = gamma / (beta * dt) * du - gamma / beta * v[i] + dt * (1 - gamma / (2 * beta)) * a_rel[i]
        da = du / (beta * dt**2) - v[i] / (beta * dt) - a_rel[i] / (2 * beta)

        u[i + 1] = u[i] + du
        v[i + 1] = v[i] + dv
        a_rel[i + 1] = a_rel[i] + da

    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    idx_peak_u = np.argmax(np.abs(u))
    omega_n = None
    T_n = None
    c_val = None
    peak_u = np.max(np.abs(u))
    peak_v = np.max(np.abs(v))
    peak_a_rel = np.max(np.abs(a_rel))
    peak_a_abs = np.max(np.abs(a_abs))
    peak_base_shear = None
    t_peak_u = t[idx_peak_u]

    metrics = {
        "natural_circular_frequency_rad_s": omega_n,
        "natural_period_s": T_n,
        "damping_coefficient_Ns_m": c_val,
        "peak_relative_displacement_m": peak_u,
        "peak_relative_velocity_m_s": peak_v,
        "peak_relative_acceleration_m_s2": peak_a_rel,
        "peak_absolute_acceleration_m_s2": peak_a_abs,
        "peak_pseudo_base_shear_N": peak_base_shear,
        "time_of_peak_displacement_s": t_peak_u,
    }
    return metrics


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    header = "t,ag,u,v,a_rel,a_abs"
    data = np.column_stack((t, ag, u, v, a_rel, a_abs))
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.6e")


def save_summary_json(path, metrics):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2, default=str)


def plot_response(path, t, ag, u, v, a_abs):
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    axes[0].plot(t, ag, linewidth=0.5)
    axes[0].set_ylabel("Ground Accel. (m/s^2)")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, u, linewidth=0.5)
    axes[1].set_ylabel("Displacement (m)")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t, v, linewidth=0.5)
    axes[2].set_ylabel("Velocity (m/s)")
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(t, a_abs, linewidth=0.5)
    axes[3].set_ylabel("Abs. Accel. (m/s^2)")
    axes[3].set_xlabel("Time (s)")
    axes[3].grid(True, alpha=0.3)

    fig.suptitle("SDOF Newmark Average Acceleration Response", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_analysis(output_dir="sdof_output"):
    os.makedirs(output_dir, exist_ok=True)

    m = 1000.0
    k = 200000.0
    zeta = 0.05
    g = 9.81
    dt = 0.01
    T = 20.0
    u0 = 0.0
    v0 = 0.0

    c = compute_damping(m, k, zeta)
    omega_n = np.sqrt(k / m)
    T_n = 2 * np.pi / omega_n

    t = make_time_vector(dt, T)
    ag = ground_acceleration(t, g)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt, u0, v0)

    idx_peak_u = np.argmax(np.abs(u))
    peak_u = np.max(np.abs(u))
    peak_v = np.max(np.abs(v))
    peak_a_rel = np.max(np.abs(a_rel))
    peak_a_abs = np.max(np.abs(a_abs))
    peak_base_shear = k * peak_u
    t_peak_u = t[idx_peak_u]

    metrics = {
        "natural_circular_frequency_rad_s": omega_n,
        "natural_period_s": T_n,
        "damping_coefficient_Ns_m": c,
        "peak_relative_displacement_m": peak_u,
        "peak_relative_velocity_m_s": peak_v,
        "peak_relative_acceleration_m_s2": peak_a_rel,
        "peak_absolute_acceleration_m_s2": peak_a_abs,
        "peak_pseudo_base_shear_N": peak_base_shear,
        "time_of_peak_displacement_s": t_peak_u,
    }

    for key, val in metrics.items():
        print(f"{key}: {val}")

    save_csv(os.path.join(output_dir, "response.csv"), t, ag, u, v, a_rel, a_abs)
    save_summary_json(os.path.join(output_dir, "summary.json"), metrics)
    plot_response(os.path.join(output_dir, "response.png"), t, ag, u, v, a_abs)

    print(f"\nOutput written to {output_dir}/")
    return t, ag, u, v, a_rel, a_abs, metrics


def main():
    run_analysis("sdof_output")


def test_zero_ground_motion():
    dt = 0.01
    T = 2.0
    t = make_time_vector(dt, T)
    ag = np.zeros_like(t)
    m, c, k = 1000.0, 2000.0, 200000.0
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt)
    assert np.allclose(u, 0.0, atol=1e-12), "displacement should stay zero"
    assert np.allclose(v, 0.0, atol=1e-12), "velocity should stay zero"
    assert np.allclose(a_rel, 0.0, atol=1e-12), "relative acceleration should stay zero"


def test_output_shapes():
    dt = 0.01
    T = 2.0
    t = make_time_vector(dt, T)
    ag = ground_acceleration(t)
    m, c, k = 1000.0, 2000.0, 200000.0
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt)
    n = len(t)
    assert len(u) == n
    assert len(v) == n
    assert len(a_rel) == n
    assert len(a_abs) == n


def test_damping_positive():
    c = compute_damping(1000.0, 200000.0, 0.05)
    assert c > 0, "damping coefficient must be positive"


if __name__ == "__main__":
    main()
```