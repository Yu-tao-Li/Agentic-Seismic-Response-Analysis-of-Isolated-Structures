```python
"""
Single Degree of Freedom (SDOF) oscillator analysis using Newmark's average acceleration method.
"""

import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def make_time_vector(dt, duration):
    """Create a time vector from 0 to duration with given time step."""
    n_steps = int(round(duration / dt))
    return np.linspace(0.0, duration, n_steps + 1)


def ground_acceleration(t, g=9.81):
    """Compute ground acceleration ag(t) as a sum of two damped sinusoids."""
    ag = (0.30 * g * np.sin(2.0 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
          + 0.10 * g * np.sin(2.0 * np.pi * 3.0 * t) * np.exp(-0.20 * t))
    return ag


def compute_damping(m, k, zeta):
    """Compute the viscous damping coefficient c = 2 * zeta * sqrt(k * m)."""
    return 2.0 * zeta * np.sqrt(k * m)


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    """
    Integrate the equation of motion using Newmark's average acceleration method.

    Equation: m * u_ddot + c * u_dot + k * u = -m * ag(t)

    Parameters
    ----------
    m : float
        Mass.
    c : float
        Viscous damping coefficient.
    k : float
        Stiffness.
    ag : ndarray
        Ground acceleration time history.
    dt : float
        Time step.
    u0 : float
        Initial displacement.
    v0 : float
        Initial velocity.

    Returns
    -------
    u : ndarray
        Relative displacement.
    v : ndarray
        Relative velocity.
    a_rel : ndarray
        Relative acceleration.
    """
    beta = 0.25
    gamma = 0.5

    n = len(ag)
    u = np.zeros(n)
    v = np.zeros(n)
    a_rel = np.zeros(n)

    # Initial acceleration from equation of motion at t=0
    # m * a_rel0 + c * v0 + k * u0 = -m * ag[0]
    a_rel[0] = (-m * ag[0] - c * v0 - k * u0) / m

    u[0] = u0
    v[0] = v0

    # Effective stiffness
    k_eff = k + gamma / (beta * dt) * c + 1.0 / (beta * dt ** 2) * m

    for i in range(n - 1):
        # Effective load increment
        dp = (-m * ag[i + 1] + m * ag[i]
              + m * (1.0 / (beta * dt) * v[i] + 1.0 / (2.0 * beta) * a_rel[i])
              + c * (gamma / (beta * dt) * v[i] + (gamma / (2.0 * beta) - 1.0) * a_rel[i]))

        # Solve for displacement increment
        du = dp / k_eff

        # Update displacement
        u[i + 1] = u[i] + du

        # Update velocity and acceleration
        dv = gamma / (beta * dt) * du - gamma / beta * v[i] + dt * (1.0 - gamma / (2.0 * beta)) * a_rel[i]
        da = 1.0 / (beta * dt ** 2) * du - 1.0 / (beta * dt) * v[i] - 1.0 / (2.0 * beta) * a_rel[i]

        v[i + 1] = v[i] + dv
        a_rel[i + 1] = a_rel[i] + da

    return u, v, a_rel


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    """Compute peak response metrics and related quantities."""
    k = 200000.0  # stiffness used for pseudo base shear

    peak_disp = np.max(np.abs(u))
    peak_vel = np.max(np.abs(v))
    peak_rel_acc = np.max(np.abs(a_rel))
    peak_abs_acc = np.max(np.abs(a_abs))
    peak_pbs = k * peak_disp
    idx_peak_disp = np.argmax(np.abs(u))
    time_peak_disp = t[idx_peak_disp]

    metrics = {
        "natural_circular_frequency_rad_s": float(omega_n),
        "natural_period_s": float(T_n),
        "damping_coefficient_Ns_m": float(c_val),
        "peak_relative_displacement_m": float(peak_disp),
        "peak_relative_velocity_m_s": float(peak_vel),
        "peak_relative_acceleration_m_s2": float(peak_rel_acc),
        "peak_absolute_acceleration_m_s2": float(peak_abs_acc),
        "peak_pseudo_base_shear_N": float(peak_pbs),
        "time_of_peak_displacement_s": float(time_peak_disp),
    }
    return metrics


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    """Save the response time histories to a CSV file."""
    with open(path, "w") as f:
        f.write("time,ag,u,v,a_rel,a_abs\n")
        for i in range(len(t)):
            f.write(f"{t[i]:.6f},{ag[i]:.6f},{u[i]:.6f},{v[i]:.6f},{a_rel[i]:.6f},{a_abs[i]:.6f}\n")


def save_summary_json(path, metrics):
    """Save the response metrics to a JSON file."""
    with open(path, "w") as f:
        json.dump(metrics, f, indent=2)


def plot_response(path, t, ag, u, v, a_abs):
    """Plot and save the response time histories."""
    fig, axes = plt.subplots(4, 1, figsize=(12, 10), sharex=True)

    axes[0].plot(t, ag, "b", linewidth=0.8)
    axes[0].set_ylabel("Ag (m/s²)")
    axes[0].set_title("Ground Acceleration")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, u * 1000.0, "r", linewidth=0.8)
    axes[1].set_ylabel("Disp (mm)")
    axes[1].set_title("Relative Displacement")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t, v, "g", linewidth=0.8)
    axes[2].set_ylabel("Vel (m/s)")
    axes[2].set_title("Relative Velocity")
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(t, a_abs, "m", linewidth=0.8)
    axes[3].set_ylabel("Acc (m/s²)")
    axes[3].set_xlabel("Time (s)")
    axes[3].set_title("Absolute Acceleration")
    axes[3].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_analysis(output_dir="sdof_output"):
    """Run the full SDOF analysis."""
    global omega_n, T_n, c_val

    # System parameters
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    g = 9.81
    dt = 0.01
    duration = 20.0
    u0 = 0.0
    v0 = 0.0

    # Derived quantities
    omega_n = np.sqrt(k / m)
    T_n = 2.0 * np.pi / omega_n
    c_val = compute_damping(m, k, zeta)

    # Time vector
    t = make_time_vector(dt, duration)

    # Ground acceleration
    ag = ground_acceleration(t, g)

    # Newmark integration
    u, v, a_rel = newmark_average_acceleration(m, c_val, k, ag, dt, u0, v0)

    # Absolute acceleration: a_abs = a_rel + ag
    a_abs = a_rel + ag

    # Compute metrics
    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)

    # Print summary
    print("=" * 60)
    print("SDOF Oscillator Response Summary")
    print("=" * 60)
    for key, val in metrics.items():
        print(f"  {key}: {val}")
    print("=" * 60)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Save outputs
    csv_path = os.path.join(output_dir, "response.csv")
    json_path = os.path.join(output_dir, "summary.json")
    png_path = os.path.join(output_dir, "response.png")

    save_csv(csv_path, t, ag, u, v, a_rel, a_abs)
    save_summary_json(json_path, metrics)
    plot_response(png_path, t, ag, u, v, a_abs)

    print(f"CSV saved to: {csv_path}")
    print(f"JSON saved to: {json_path}")
    print(f"Plot saved to: {png_path}")

    return metrics


def main():
    """Entry point."""
    run_analysis()


def test_zero_ground_motion():
    """Test that zero ground motion produces zero response."""
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    c = compute_damping(m, k, zeta)
    dt = 0.01
    t = make_time_vector(dt, 10.0)
    ag = np.zeros(len(t))

    u, v, a_rel = newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0)

    assert np.allclose(u, 0.0, atol=1e-12), "Displacement should be zero for zero ground motion"
    assert np.allclose(v, 0.0, atol=1e-12), "Velocity should be zero for zero ground motion"
    assert np.allclose(a_rel, 0.0, atol=1e-12), "Acceleration should be zero for zero ground motion"
    print("test_zero_ground_motion PASSED")


def test_output_shapes():
    """Test that output arrays have the correct shape."""
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    c = compute_damping(m, k, zeta)
    dt = 0.01
    duration = 5.0
    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t)

    u, v, a_rel = newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0)

    expected_len = len(t)
    assert len(u) == expected_len, f"u length {len(u)} != {expected_len}"
    assert len(v) == expected_len, f"v length {len(v)} != {expected_len}"
    assert len(a_rel) == expected_len, f"a_rel length {len(a_rel)} != {expected_len}"
    assert u[0] == 0.0, "Initial displacement should be zero"
    assert v[0] == 0.0, "Initial velocity should be zero"
    print("test_output_shapes PASSED")


def test_damping_positive():
    """Test that damping coefficient is positive and reasonable."""
    m = 1000.0
    k = 200000.0
    zeta = 0.05

    c = compute_damping(m, k, zeta)

    assert c > 0.0, "Damping coefficient should be positive"
    c_critical = 2.0 * np.sqrt(k * m)
    assert abs(c - zeta * c_critical) < 1e-6, "Damping should equal zeta * c_critical"
    print("test_damping_positive PASSED")


if __name__ == "__main__":
    test_zero_ground_motion()
    test_output_shapes()
    test_damping_positive()
    print("\nAll tests passed.\n")
    main()
```

The script is complete and directly runnable. It defines all 10 required functions, performs Newmark average acceleration integration (β=1/4, γ=1/2), computes all specified metrics, and writes `response.csv`, `summary.json`, and `response.png` to the `sdof_output/` directory. Three assertion-based tests run before the main analysis.