```python
"""
SDOF oscillator response analysis using Newmark average acceleration method.
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── Parameters ──────────────────────────────────────────────────────────────

M = 1000.0
K = 200000.0
ZETA = 0.05
G = 9.81
DT = 0.01
DURATION = 20.0
U0 = 0.0
V0 = 0.0


# ── Core functions ──────────────────────────────────────────────────────────

def make_time_vector(dt, duration):
    n_steps = int(round(duration / dt)) + 1
    return np.linspace(0.0, duration, n_steps)


def ground_acceleration(t, g=9.81):
    term1 = 0.30 * g * np.sin(2.0 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
    term2 = 0.10 * g * np.sin(2.0 * np.pi * 3.0 * t) * np.exp(-0.20 * t)
    return term1 + term2


def compute_damping(m, k, zeta):
    omega_n = np.sqrt(k / m)
    c = 2.0 * zeta * m * omega_n
    return c


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    beta = 1.0 / 4.0
    gamma = 1.0 / 2.0

    n = len(ag)
    u = np.zeros(n)
    v = np.zeros(n)
    a = np.zeros(n)

    # Initial acceleration from equation of motion
    a0 = (-m * ag[0] - c * v0 - k * u0) / m
    u[0] = u0
    v[0] = v0
    a[0] = a0

    # Effective stiffness
    k_eff = k + gamma / (beta * dt) * c + 1.0 / (beta * dt ** 2) * m

    for i in range(n - 1):
        dp = (-m * ag[i + 1] + m * (1.0 / (beta * dt ** 2) * u[i]
              + 1.0 / (beta * dt) * v[i]
              + (1.0 / (2.0 * beta) - 1.0) * a[i])
              + c * (gamma / (beta * dt) * u[i]
              + (gamma / beta - 1.0) * v[i]
              + (gamma / (2.0 * beta) - 1.0) * dt * a[i]))

        u[i + 1] = dp / k_eff
        v[i + 1] = (gamma / (beta * dt) * (u[i + 1] - u[i])
                     + (1.0 - gamma / beta) * v[i]
                     + (1.0 - gamma / (2.0 * beta)) * dt * a[i])
        a[i + 1] = (1.0 / (beta * dt ** 2) * (u[i + 1] - u[i])
                     - 1.0 / (beta * dt) * v[i]
                     - (1.0 / (2.0 * beta) - 1.0) * a[i])

    a_rel = a.copy()
    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    m = M
    k = K
    zeta = ZETA
    omega_n = np.sqrt(k / m)
    T_n = 2.0 * np.pi / omega_n
    c = 2.0 * zeta * m * omega_n

    idx_peak_u = int(np.argmax(np.abs(u)))
    metrics = {
        "natural_circular_frequency_rad_s": round(float(omega_n), 6),
        "natural_period_s": round(float(T_n), 6),
        "damping_coefficient_Ns_m": round(float(c), 6),
        "peak_relative_displacement_m": round(float(np.max(np.abs(u))), 8),
        "peak_relative_velocity_m_s": round(float(np.max(np.abs(v))), 8),
        "peak_relative_acceleration_m_s2": round(float(np.max(np.abs(a_rel))), 8),
        "peak_absolute_acceleration_m_s2": round(float(np.max(np.abs(a_abs))), 8),
        "peak_pseudo_base_shear_N": round(float(k * np.max(np.abs(u))), 4),
        "time_of_peak_displacement_s": round(float(t[idx_peak_u]), 4),
    }
    return metrics


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    header = "time_s,ground_acc_m_s2,rel_disp_m,rel_vel_m_s,rel_acc_m_s2,abs_acc_m_s2"
    data = np.column_stack([t, ag, u, v, a_rel, a_abs])
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.8e")


def save_summary_json(path, metrics):
    with open(path, "w") as f:
        json.dump(metrics, f, indent=2)


def plot_response(path, t, ag, u, v, a_abs):
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    axes[0].plot(t, ag, color="tab:red", linewidth=0.8)
    axes[0].set_ylabel("Ground Accel (m/s²)")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, u * 1000.0, color="tab:blue", linewidth=0.8)
    axes[1].set_ylabel("Rel. Disp. (mm)")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t, a_abs, color="tab:green", linewidth=0.8)
    axes[2].set_ylabel("Abs. Accel (m/s²)")
    axes[2].set_xlabel("Time (s)")
    axes[2].grid(True, alpha=0.3)

    fig.suptitle("SDOF Newmark Response", fontsize=13)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_analysis(output_dir="sdof_output"):
    os.makedirs(output_dir, exist_ok=True)

    t = make_time_vector(DT, DURATION)
    ag = ground_acceleration(t, g=G)
    c = compute_damping(M, K, ZETA)
    u, v, a_rel, a_abs = newmark_average_acceleration(M, c, K, ag, DT, U0, V0)
    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)

    csv_path = os.path.join(output_dir, "response.csv")
    json_path = os.path.join(output_dir, "summary.json")
    png_path = os.path.join(output_dir, "response.png")

    save_csv(csv_path, t, ag, u, v, a_rel, a_abs)
    save_summary_json(json_path, metrics)
    plot_response(png_path, t, ag, u, v, a_abs)

    print("Analysis complete. Output written to:", output_dir)
    for k, v in metrics.items():
        print(f"  {k}: {v}")

    return t, ag, u, v, a_rel, a_abs, metrics


# ── Test functions ──────────────────────────────────────────────────────────

def test_zero_ground_motion():
    t = make_time_vector(DT, DURATION)
    ag = np.zeros_like(t)
    c = compute_damping(M, K, ZETA)
    u, v, a_rel, a_abs = newmark_average_acceleration(M, c, K, ag, DT, U0, V0)
    assert np.allclose(u, 0.0, atol=1e-14), "Displacement should be zero"
    assert np.allclose(v, 0.0, atol=1e-14), "Velocity should be zero"
    assert np.allclose(a_abs, 0.0, atol=1e-14), "Acceleration should be zero"
    print("test_zero_ground_motion PASSED")


def test_output_shapes():
    t = make_time_vector(DT, DURATION)
    ag = ground_acceleration(t, g=G)
    c = compute_damping(M, K, ZETA)
    u, v, a_rel, a_abs = newmark_average_acceleration(M, c, K, ag, DT, U0, V0)
    n = len(t)
    assert len(u) == n, f"u length {len(u)} != {n}"
    assert len(v) == n, f"v length {len(v)} != {n}"
    assert len(a_rel) == n, f"a_rel length {len(a_rel)} != {n}"
    assert len(a_abs) == n, f"a_abs length {len(a_abs)} != {n}"
    assert len(ag) == n, f"ag length {len(ag)} != {n}"
    print("test_output_shapes PASSED")


def test_damping_positive():
    c = compute_damping(M, K, ZETA)
    omega_n = np.sqrt(K / M)
    expected_c = 2.0 * ZETA * M * omega_n
    assert c > 0, "Damping coefficient must be positive"
    assert abs(c - expected_c) < 1e-10, "Damping coefficient mismatch"
    print("test_damping_positive PASSED")


# ── Entry point ─────────────────────────────────────────────────────────────

def main():
    test_zero_ground_motion()
    test_output_shapes()
    test_damping_positive()
    print()
    run_analysis()


if __name__ == "__main__":
    main()
```