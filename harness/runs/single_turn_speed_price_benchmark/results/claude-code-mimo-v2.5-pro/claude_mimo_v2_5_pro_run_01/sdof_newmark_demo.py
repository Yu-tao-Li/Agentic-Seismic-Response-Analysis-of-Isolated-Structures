<tool_call>
<function=Write>
<parameter=file_path>harness\runs\single_turn_speed_price_benchmark\results\claude-code-mimo-v2.5-pro\claude_mimo_v2_5_pro_run_01\sdof_newmark_demo.py</parameter>
<parameter=content>import os
import json
import numpy as np
import matplotlib.pyplot as plt


def make_time_vector(dt, duration):
    n_steps = int(round(duration / dt)) + 1
    return np.linspace(0.0, duration, n_steps)


def ground_acceleration(t, g=9.81):
    ag = (0.30 * g * np.sin(2 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
          + 0.10 * g * np.sin(2 * np.pi * 3.0 * t) * np.exp(-0.20 * t))
    return ag


def compute_damping(m, k, zeta):
    omega_n = np.sqrt(k / m)
    c = 2.0 * zeta * omega_n * m
    return c


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    beta = 1.0 / 4.0
    gamma = 1.0 / 2.0

    n = len(ag)
    u = np.zeros(n)
    v = np.zeros(n)
    a = np.zeros(n)

    a0 = (-m * ag[0] - c * v0 - k * u0) / m
    u[0] = u0
    v[0] = v0
    a[0] = a0

    k_eff = k + gamma / (beta * dt) * c + 1.0 / (beta * dt**2) * m
    a1 = m / (beta * dt) + gamma / beta * c
    a2 = m / (2.0 * beta) + dt * (gamma / (2.0 * beta) - 1.0) * c

    for i in range(n - 1):
        dp = (-m * ag[i + 1] - (-m * ag[i])) + a1 * v[i] + a2 * a[i]
        du = dp / k_eff
        dv = gamma / (beta * dt) * du - gamma / beta * v[i] + dt * (1.0 - gamma / (2.0 * beta)) * a[i]
        da = 1.0 / (beta * dt**2) * du - 1.0 / (beta * dt) * v[i] - 1.0 / (2.0 * beta) * a[i]

        u[i + 1] = u[i] + du
        v[i + 1] = v[i] + dv
        a[i + 1] = a[i] + da

    a_rel = a
    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    m = 1000.0
    k = 200000.0
    zeta = 0.05

    omega_n = np.sqrt(k / m)
    T_n = 2.0 * np.pi / omega_n
    c = compute_damping(m, k, zeta)

    peak_disp = np.max(np.abs(u))
    peak_vel = np.max(np.abs(v))
    peak_rel_acc = np.max(np.abs(a_rel))
    peak_abs_acc = np.max(np.abs(a_abs))
    peak_base_shear = k * peak_disp
    idx_peak_disp = np.argmax(np.abs(u))
    time_peak_disp = t[idx_peak_disp]

    metrics = {
        "natural_circular_frequency_rad_s": float(omega_n),
        "natural_period_s": float(T_n),
        "damping_coefficient_Ns_m": float(c),
        "peak_relative_displacement_m": float(peak_disp),
        "peak_relative_velocity_m_s": float(peak_vel),
        "peak_relative_acceleration_m_s2": float(peak_rel_acc),
        "peak_absolute_acceleration_m_s2": float(peak_abs_acc),
        "peak_pseudo_base_shear_N": float(peak_base_shear),
        "time_of_peak_displacement_s": float(time_peak_disp),
    }
    return metrics


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    with open(path, "w") as f:
        f.write("t_s,ag_m_s2,u_m,v_m_s,a_rel_m_s2,a_abs_m_s2\n")
        for i in range(len(t)):
            f.write(f"{t[i]:.6f},{ag[i]:.6f},{u[i]:.6f},{v[i]:.6f},{a_rel[i]:.6f},{a_abs[i]:.6f}\n")


def save_summary_json(path, metrics):
    with open(path, "w") as f:
        json.dump(metrics, f, indent=2)


def plot_response(path, t, ag, u, v, a_abs):
    fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

    axes[0].plot(t, ag, color="black", linewidth=0.8)
    axes[0].set_ylabel("ag (m/s²)")
    axes[0].set_title("SDOF Newmark Average Acceleration — Response")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, u, color="blue", linewidth=0.8)
    axes[1].set_ylabel("u (m)")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t, v, color="red", linewidth=0.8)
    axes[2].set_ylabel("v (m/s)")
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(t, a_abs, color="green", linewidth=0.8)
    axes[3].set_ylabel("a_abs (m/s²)")
    axes[3].set_xlabel("Time (s)")
    axes[3].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()


def run_analysis(output_dir="sdof_output"):
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    g = 9.81
    dt = 0.01
    duration = 20.0
    u0 = 0.0
    v0 = 0.0

    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t, g)
    c = compute_damping(m, k, zeta)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt, u0, v0)
    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)

    os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(output_dir, "response.csv")
    json_path = os.path.join(output_dir, "summary.json")
    png_path = os.path.join(output_dir, "response.png")

    save_csv(csv_path, t, ag, u, v, a_rel, a_abs)
    save_summary_json(json_path, metrics)
    plot_response(png_path, t, ag, u, v, a_abs)

    print("=== SDOF Newmark Analysis Results ===")
    for key, val in metrics.items():
        print(f"  {key}: {val:.6f}")
    print(f"\nOutput directory: {os.path.abspath(output_dir)}")
    print(f"  - {csv_path}")
    print(f"  - {json_path}")
    print(f"  - {png_path}")

    return t, ag, u, v, a_rel, a_abs, metrics


def main():
    run_analysis("sdof_output")


def test_zero_ground_motion():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    dt = 0.01
    duration = 2.0

    t = make_time_vector(dt, duration)
    ag = np.zeros_like(t)
    c = compute_damping(m, k, zeta)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0)

    assert np.allclose(u, 0.0, atol=1e-12), "Displacement should be zero for zero ground motion"
    assert np.allclose(v, 0.0, atol=1e-12), "Velocity should be zero for zero ground motion"
    assert np.allclose(a_rel, 0.0, atol=1e-12), "Acceleration should be zero for zero ground motion"
    print("test_zero_ground_motion PASSED")


def test_output_shapes():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    g = 9.81
    dt = 0.01
    duration = 1.0

    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t, g)
    c = compute_damping(m, k, zeta)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt)

    n = len(t)
    assert len(u) == n, f"u length {len(u)} != {n}"
    assert len(v) == n, f"v length {len(v)} != {n}"
    assert len(a_rel) == n, f"a_rel length {len(a_rel)} != {n}"
    assert len(a_abs) == n, f"a_abs length {len(a_abs)} != {n}"
    assert len(ag) == n, f"ag length {len(ag)} != {n}"

    t2 = make_time_vector(0.01, 20.0)
    assert len(t2) == 2001, f"Expected 2001 points for 20s at dt=0.01, got {len(t2)}"
    print("test_output_shapes PASSED")


def test_damping_positive():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    c = compute_damping(m, k, zeta)
    assert c > 0, f"Damping coefficient should be positive, got {c}"

    omega_n = np.sqrt(k / m)
    expected_c = 2.0 * zeta * omega_n * m
    assert abs(c - expected_c) < 1e-10, f"Expected c={expected_c}, got {c}"
    print("test_damping_positive PASSED")


if __name__ == "__main__":
    test_zero_ground_motion()
    test_output_shapes()
    test_damping_positive()
    print("\nAll tests passed.\n")
    main()
</parameter>
</function>
</tool_call>