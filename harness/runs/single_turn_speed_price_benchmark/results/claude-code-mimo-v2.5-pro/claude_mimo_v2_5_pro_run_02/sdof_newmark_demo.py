<tool_call>
<function=write>
<parameter=file_path>harness\runs\single_turn_speed_price_benchmark\results\claude-code-mimo-v2.5-pro\claude_mimo_v2_5_pro_run_02\sdof_newmark_demo.py</parameter>
<parameter=content>import numpy as np
import matplotlib.pyplot as plt
import json
import os


def make_time_vector(dt, duration):
    return np.arange(0.0, duration + dt * 0.5, dt)


def ground_acceleration(t, g=9.81):
    return (
        0.30 * g * np.sin(2 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
        + 0.10 * g * np.sin(2 * np.pi * 3.0 * t) * np.exp(-0.20 * t)
    )


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
    a_rel = np.zeros(n)

    a0 = (-m * ag[0] - c * v0 - k * u0) / m
    u[0] = u0
    v[0] = v0
    a_rel[0] = a0

    k_eff = k + gamma / (beta * dt) * c + 1.0 / (beta * dt**2) * m
    a = 1.0 / (beta * dt) * m + gamma / beta * c
    b = 1.0 / (2.0 * beta) * m + dt * (gamma / (2.0 * beta) - 1.0) * c

    for i in range(n - 1):
        dp = -m * (ag[i + 1] - ag[i]) + a * v[i] + b * a_rel[i]
        du = dp / k_eff
        dv = gamma / (beta * dt) * du - gamma / beta * v[i] + dt * (1.0 - gamma / (2.0 * beta)) * a_rel[i]
        da = 1.0 / (beta * dt**2) * du - 1.0 / (beta * dt) * v[i] - 1.0 / (2.0 * beta) * a_rel[i]

        u[i + 1] = u[i] + du
        v[i + 1] = v[i] + dv
        a_rel[i + 1] = a_rel[i] + da

    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    m = 1000.0
    k = 200000.0
    zeta = 0.05

    omega_n = np.sqrt(k / m)
    T_n = 2.0 * np.pi / omega_n
    c = 2.0 * zeta * m * omega_n

    peak_disp = np.max(np.abs(u))
    peak_vel = np.max(np.abs(v))
    peak_acc_rel = np.max(np.abs(a_rel))
    peak_acc_abs = np.max(np.abs(a_abs))
    peak_base_shear = k * peak_disp
    idx_peak = np.argmax(np.abs(u))
    time_peak_disp = t[idx_peak]

    metrics = {
        "natural_circular_frequency_rad_s": float(omega_n),
        "natural_period_s": float(T_n),
        "damping_coefficient_Ns_m": float(c),
        "peak_relative_displacement_m": float(peak_disp),
        "peak_relative_velocity_m_s": float(peak_vel),
        "peak_relative_acceleration_m_s2": float(peak_acc_rel),
        "peak_absolute_acceleration_m_s2": float(peak_acc_abs),
        "peak_pseudo_base_shear_N": float(peak_base_shear),
        "time_of_peak_displacement_s": float(time_peak_disp),
    }
    return metrics


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    header = "time,ag,u,v,a_rel,a_abs"
    data = np.column_stack((t, ag, u, v, a_rel, a_abs))
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.8e")


def save_summary_json(path, metrics):
    with open(path, "w") as f:
        json.dump(metrics, f, indent=2)


def plot_response(path, t, ag, u, v, a_abs):
    fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)

    axes[0].plot(t, ag, color="tab:blue", linewidth=0.8)
    axes[0].set_ylabel("ag (m/s²)")
    axes[0].set_title("Ground Acceleration")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, u * 1000, color="tab:red", linewidth=0.8)
    axes[1].set_ylabel("u (mm)")
    axes[1].set_title("Relative Displacement")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t, v, color="tab:green", linewidth=0.8)
    axes[2].set_ylabel("v (m/s)")
    axes[2].set_title("Relative Velocity")
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(t, a_abs, color="tab:purple", linewidth=0.8)
    axes[3].set_ylabel("a_abs (m/s²)")
    axes[3].set_xlabel("Time (s)")
    axes[3].set_title("Absolute Acceleration")
    axes[3].grid(True, alpha=0.3)

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
    c = compute_damping(m, k, zeta)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt, u0, v0)

    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)

    csv_path = os.path.join(output_dir, "response.csv")
    json_path = os.path.join(output_dir, "summary.json")
    png_path = os.path.join(output_dir, "response.png")

    save_csv(csv_path, t, ag, u, v, a_rel, a_abs)
    save_summary_json(json_path, metrics)
    plot_response(png_path, t, ag, u, v, a_abs)

    for key, val in metrics.items():
        print(f"  {key}: {val}")

    return metrics


def test_zero_ground_motion():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    dt = 0.01
    duration = 1.0
    t = make_time_vector(dt, duration)
    ag = np.zeros(len(t))
    c = compute_damping(m, k, zeta)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0)
    assert np.allclose(u, 0.0, atol=1e-14), "Displacement should be zero for zero ground motion"
    assert np.allclose(v, 0.0, atol=1e-14), "Velocity should be zero for zero ground motion"
    assert np.allclose(a_rel, 0.0, atol=1e-14), "Acceleration should be zero for zero ground motion"
    print("test_zero_ground_motion passed")


def test_output_shapes():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    g = 9.81
    dt = 0.01
    duration = 20.0
    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t, g)
    c = compute_damping(m, k, zeta)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt)
    n = len(t)
    assert len(u) == n, f"u length {len(u)} != {n}"
    assert len(v) == n, f"v length {len(v)} != {n}"
    assert len(a_rel) == n, f"a_rel length {len(a_rel)} != {n}"
    assert len(a_abs) == n, f"a_abs length {len(a_abs)} != {n}"
    print("test_output_shapes passed")


def test_damping_positive():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    c = compute_damping(m, k, zeta)
    omega_n = np.sqrt(k / m)
    assert c > 0, "Damping coefficient must be positive"
    assert np.isclose(c, 2.0 * zeta * m * omega_n), "Damping coefficient mismatch"
    print("test_damping_positive passed")


if __name__ == "__main__":
    test_zero_ground_motion()
    test_output_shapes()
    test_damping_positive()
    print("\n--- SDOF Newmark Analysis ---")
    run_analysis()
</parameter>
</function>
</tool_call>