<tool_call>
<function=Write>
<parameter=file_path>harness\runs\single_turn_speed_price_benchmark\results\claude-code-mimo-v2.5-pro\claude_mimo_v2_5_pro_run_05\sdof_newmark_demo.py</parameter>
<parameter=content>import numpy as np
import matplotlib.pyplot as plt
import json
import os
import math


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
    return omega_n, c


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    beta = 1.0 / 4.0
    gamma = 1.0 / 2.0

    n = len(ag)
    u = np.zeros(n)
    v = np.zeros(n)
    a_rel = np.zeros(n)

    a_rel[0] = (-m * ag[0] - c * v0 - k * u0) / m

    k_eff = k + gamma / (beta * dt) * c + 1.0 / (beta * dt ** 2) * m
    a = 1.0 / (beta * dt) * m + gamma / beta * c
    b = 1.0 / (2.0 * beta) * m + dt * (gamma / (2.0 * beta) - 1.0) * c

    u[0] = u0
    v[0] = v0

    for i in range(n - 1):
        dp = -m * (ag[i + 1] - ag[i]) + a * v[i] + b * a_rel[i]
        du = dp / k_eff
        u[i + 1] = u[i] + du
        dv = gamma / (beta * dt) * du - gamma / beta * v[i] + dt * (1.0 - gamma / (2.0 * beta)) * a_rel[i]
        v[i + 1] = v[i] + dv
        da = 1.0 / (beta * dt ** 2) * du - 1.0 / (beta * dt) * v[i] - 1.0 / (2.0 * beta) * a_rel[i]
        a_rel[i + 1] = a_rel[i] + da

    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs, m, k, zeta):
    omega_n, c = compute_damping(m, k, zeta)
    T_n = 2.0 * np.pi / omega_n

    idx_disp = int(np.argmax(np.abs(u)))
    peak_disp = np.abs(u[idx_disp])
    time_peak_disp = t[idx_disp]

    peak_vel = np.max(np.abs(v))
    peak_a_rel = np.max(np.abs(a_rel))
    peak_a_abs = np.max(np.abs(a_abs))
    peak_base_shear = k * peak_disp

    metrics = {
        "natural_circular_frequency_rad_s": float(omega_n),
        "natural_period_s": float(T_n),
        "damping_coefficient_Ns_m": float(c),
        "peak_relative_displacement_m": float(peak_disp),
        "time_of_peak_displacement_s": float(time_peak_disp),
        "peak_relative_velocity_m_s": float(peak_vel),
        "peak_relative_acceleration_m_s2": float(peak_a_rel),
        "peak_absolute_acceleration_m_s2": float(peak_a_abs),
        "peak_pseudo_base_shear_N": float(peak_base_shear),
    }
    return metrics


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    header = "time,ground_acceleration,relative_displacement,relative_velocity,relative_acceleration,absolute_acceleration"
    data = np.column_stack([t, ag, u, v, a_rel, a_abs])
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.10e")


def save_summary_json(path, metrics):
    with open(path, "w") as f:
        json.dump(metrics, f, indent=2)


def plot_response(path, t, ag, u, v, a_abs):
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    axes[0].plot(t, ag, color="black", linewidth=0.7)
    axes[0].set_ylabel("Ag (m/s²)")
    axes[0].set_title("Ground Acceleration")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, u * 1000, color="blue", linewidth=0.7)
    axes[1].set_ylabel("Disp (mm)")
    axes[1].set_title("Relative Displacement")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t, a_abs, color="red", linewidth=0.7)
    axes[2].set_ylabel("A_abs (m/s²)")
    axes[2].set_xlabel("Time (s)")
    axes[2].set_title("Absolute Acceleration")
    axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()


def run_analysis(output_dir="sdof_output"):
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    g = 9.81
    dt = 0.01
    T = 20.0
    u0 = 0.0
    v0 = 0.0

    t = make_time_vector(dt, T)
    ag = ground_acceleration(t, g)
    omega_n, c = compute_damping(m, k, zeta)

    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt, u0, v0)
    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs, m, k, zeta)

    os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(output_dir, "response.csv")
    save_csv(csv_path, t, ag, u, v, a_rel, a_abs)

    json_path = os.path.join(output_dir, "summary.json")
    save_summary_json(json_path, metrics)

    png_path = os.path.join(output_dir, "response.png")
    plot_response(png_path, t, ag, u, v, a_abs)

    print("=== SDOF Newmark Analysis Results ===")
    print(f"  Natural frequency        : {metrics['natural_circular_frequency_rad_s']:.4f} rad/s")
    print(f"  Natural period           : {metrics['natural_period_s']:.4f} s")
    print(f"  Damping coefficient      : {metrics['damping_coefficient_Ns_m']:.4f} Ns/m")
    print(f"  Peak rel. displacement   : {metrics['peak_relative_displacement_m']:.6e} m")
    print(f"    at time                : {metrics['time_of_peak_displacement_s']:.2f} s")
    print(f"  Peak rel. velocity       : {metrics['peak_relative_velocity_m_s']:.6e} m/s")
    print(f"  Peak rel. acceleration   : {metrics['peak_relative_acceleration_m_s2']:.6e} m/s²")
    print(f"  Peak abs. acceleration   : {metrics['peak_absolute_acceleration_m_s2']:.6e} m/s²")
    print(f"  Peak pseudo base shear   : {metrics['peak_pseudo_base_shear_N']:.4f} N")
    print(f"  Outputs written to       : {os.path.abspath(output_dir)}")

    return t, ag, u, v, a_rel, a_abs, metrics


def test_zero_ground_motion():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    dt = 0.01
    duration = 2.0
    omega_n, c = compute_damping(m, k, zeta)
    t = make_time_vector(dt, duration)
    ag = np.zeros_like(t)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0)
    assert np.allclose(u, 0.0, atol=1e-15), "Displacement should be zero for zero ground motion"
    assert np.allclose(v, 0.0, atol=1e-15), "Velocity should be zero for zero ground motion"
    assert np.allclose(a_rel, 0.0, atol=1e-15), "Relative acceleration should be zero for zero ground motion"
    assert np.allclose(a_abs, 0.0, atol=1e-15), "Absolute acceleration should be zero for zero ground motion"
    print("test_zero_ground_motion: PASSED")


def test_output_shapes():
    dt = 0.01
    duration = 20.0
    t = make_time_vector(dt, duration)
    n = len(t)
    ag = ground_acceleration(t)
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    omega_n, c = compute_damping(m, k, zeta)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt)
    assert len(t) == n, f"Time vector length {len(t)} != expected {n}"
    assert len(u) == n, f"Displacement length {len(u)} != expected {n}"
    assert len(v) == n, f"Velocity length {len(v)} != expected {n}"
    assert len(a_rel) == n, f"Rel. acceleration length {len(a_rel)} != expected {n}"
    assert len(a_abs) == n, f"Abs. acceleration length {len(a_abs)} != expected {n}"
    assert len(ag) == n, f"Ground acceleration length {len(ag)} != expected {n}"
    print("test_output_shapes: PASSED")


def test_damping_positive():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    omega_n, c = compute_damping(m, k, zeta)
    assert omega_n > 0, "Natural frequency must be positive"
    assert c > 0, "Damping coefficient must be positive"
    assert zeta < 1.0, "System should be underdamped"
    c_expected = 2.0 * zeta * m * omega_n
    assert abs(c - c_expected) < 1e-10, f"Damping coefficient mismatch: {c} vs {c_expected}"
    print("test_damping_positive: PASSED")


def main():
    t, ag, u, v, a_rel, a_abs, metrics = run_analysis()
    print("\n--- Running tests ---")
    test_zero_ground_motion()
    test_output_shapes()
    test_damping_positive()
    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
</parameter>
</function>
</tool_call>