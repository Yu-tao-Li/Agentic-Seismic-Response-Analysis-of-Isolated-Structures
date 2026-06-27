<tool_call>
<function=Write>
<parameter=file_path>harness\runs\single_turn_speed_price_benchmark\results\claude-code-mimo-v2.5-pro-ultraspeed\claude_mimo_v2_5_pro_ultraspeed_run_01\sdof_newmark_demo.py</parameter>
<parameter=content>import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def make_time_vector(dt, duration):
    return np.arange(0.0, duration + dt * 0.5, dt)


def ground_acceleration(t, g=9.81):
    a = 0.30 * g * np.sin(2 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
    b = 0.10 * g * np.sin(2 * np.pi * 3.0 * t) * np.exp(-0.20 * t)
    return a + b


def compute_damping(m, k, zeta):
    omega_n = np.sqrt(k / m)
    c = 2.0 * zeta * m * omega_n
    return c, omega_n


def newmark_average_accelration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    beta = 1.0 / 4.0
    gamma = 1.0 / 2.0

    n = len(ag)
    u = np.zeros(n)
    v = np.zeros(n)
    a_rel = np.zeros(n)

    u[0] = u0
    v[0] = v0
    a_rel[0] = -ag[0] - (c * v0 + k * u0) / m

    a1 = 1.0 / (beta * dt ** 2) * m + gamma / (beta * dt) * c
    a2 = 1.0 / (beta * dt) * m + (gamma / beta - 1.0) * c
    a3 = (1.0 / (2.0 * beta) - 1.0) * m + dt * (gamma / (2.0 * beta) - 1.0) * c

    k_eff = k + a1

    for i in range(n - 1):
        dp = -m * ag[i + 1] + a2 * v[i] + a3 * a_rel[i]
        u[i + 1] = dp / k_eff
        v[i + 1] = gamma / (beta * dt) * (u[i + 1] - u[i]) + (1.0 - gamma / beta) * v[i] + dt * (1.0 - gamma / (2.0 * beta)) * a_rel[i]
        a_rel[i + 1] = 1.0 / (beta * dt ** 2) * (u[i + 1] - u[i]) - 1.0 / (beta * dt) * v[i] - (1.0 / (2.0 * beta) - 1.0) * a_rel[i]

    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    m = 1000.0
    k = 200000.0
    zeta = 0.05

    c, omega_n = compute_damping(m, k, zeta)
    T_n = 2.0 * np.pi / omega_n

    peak_u = np.max(np.abs(u))
    peak_v = np.max(np.abs(v))
    peak_a_rel = np.max(np.abs(a_rel))
    peak_a_abs = np.max(np.abs(a_abs))
    peak_base_shear = k * peak_u
    idx_peak = np.argmax(np.abs(u))
    t_peak = t[idx_peak]

    metrics = {
        "natural_circular_frequency_rad_s": float(omega_n),
        "natural_period_s": float(T_n),
        "damping_coefficient_Ns_m": float(c),
        "peak_relative_displacement_m": float(peak_u),
        "peak_relative_velocity_m_s": float(peak_v),
        "peak_relative_acceleration_m_s2": float(peak_a_rel),
        "peak_absolute_acceleration_m_s2": float(peak_a_abs),
        "peak_pseudo_base_shear_N": float(peak_base_shear),
        "time_of_peak_displacement_s": float(t_peak),
    }
    return metrics


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    with open(path, "w") as f:
        f.write("t,ag,u,v,a_rel,a_abs\n")
        for i in range(len(t)):
            f.write(f"{t[i]:.6f},{ag[i]:.6e},{u[i]:.6e},{v[i]:.6e},{a_rel[i]:.6e},{a_abs[i]:.6e}\n")


def save_summary_json(path, metrics):
    with open(path, "w") as f:
        json.dump(metrics, f, indent=2)


def plot_response(path, t, ag, u, v, a_abs):
    fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(t, ag, linewidth=0.8)
    axes[0].set_ylabel("ag (m/s²)")
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
    axes[3].set_ylabel("a_abs (m/s²)")
    axes[3].set_xlabel("Time (s)")
    axes[3].set_title("Absolute Acceleration")
    axes[3].grid(True)

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

    c, omega_n = compute_damping(m, k, zeta)
    T_n = 2.0 * np.pi / omega_n

    print(f"Natural circular frequency: {omega_n:.4f} rad/s")
    print(f"Natural period: {T_n:.4f} s")
    print(f"Damping coefficient: {c:.4f} Ns/m")

    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t, g)

    u, v, a_rel, a_abs = newmark_average_accelration(m, c, k, ag, dt, u0, v0)

    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)

    os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(output_dir, "response.csv")
    save_csv(csv_path, t, ag, u, v, a_rel, a_abs)
    print(f"CSV saved to {csv_path}")

    json_path = os.path.join(output_dir, "summary.json")
    save_summary_json(json_path, metrics)
    print(f"JSON saved to {json_path}")

    png_path = os.path.join(output_dir, "response.png")
    plot_response(png_path, t, ag, u, v, a_abs)
    print(f"Plot saved to {png_path}")

    print("\nResponse Metrics:")
    for key, val in metrics.items():
        print(f"  {key}: {val}")

    return t, ag, u, v, a_rel, a_abs, metrics


def main():
    run_analysis("sdof_output")


def test_zero_ground_motion():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    c, omega_n = compute_damping(m, k, zeta)
    dt = 0.01
    t = make_time_vector(dt, 5.0)
    ag = np.zeros(len(t))
    u, v, a_rel, a_abs = newmark_average_accelration(m, c, k, ag, dt, u0=0.0, v0=0.0)
    assert np.allclose(u, 0.0, atol=1e-14), "Displacement should be zero for zero ground motion"
    assert np.allclose(v, 0.0, atol=1e-14), "Velocity should be zero for zero ground motion"
    assert np.allclose(a_rel, 0.0, atol=1e-14), "Acceleration should be zero for zero ground motion"
    print("test_zero_ground_motion: PASSED")


def test_output_shapes():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    c, omega_n = compute_damping(m, k, zeta)
    dt = 0.01
    t = make_time_vector(dt, 10.0)
    ag = ground_acceleration(t)
    u, v, a_rel, a_abs = newmark_average_accelration(m, c, k, ag, dt)
    n = len(t)
    assert len(u) == n, f"u length {len(u)} != {n}"
    assert len(v) == n, f"v length {len(v)} != {n}"
    assert len(a_rel) == n, f"a_rel length {len(a_rel)} != {n}"
    assert len(a_abs) == n, f"a_abs length {len(a_abs)} != {n}"
    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)
    assert "peak_relative_displacement_m" in metrics
    assert "natural_circular_frequency_rad_s" in metrics
    print("test_output_shapes: PASSED")


def test_damping_positive():
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    c, omega_n = compute_damping(m, k, zeta)
    assert c > 0, "Damping coefficient must be positive"
    assert omega_n > 0, "Natural frequency must be positive"
    expected_c = 2.0 * zeta * np.sqrt(k * m)
    assert abs(c - expected_c) < 1e-10, f"c={c}, expected={expected_c}"
    expected_omega = np.sqrt(k / m)
    assert abs(omega_n - expected_omega) < 1e-10, f"omega_n={omega_n}, expected={expected_omega}"
    print("test_damping_positive: PASSED")


if __name__ == "__main__":
    main()
    print("\nRunning tests...")
    test_zero_ground_motion()
    test_output_shapes()
    test_damping_positive()
    print("\nAll tests passed!")
</parameter>
</function>
</tool_call>