```python
"""
SDOF Newmark Average Acceleration Demo
Analyzes a single-degree-of-freedom oscillator under prescribed ground motion.
"""

import os
import json
import csv
import numpy as np
import matplotlib.pyplot as plt


def make_time_vector(dt, duration):
    """Return a time array from 0 to duration (inclusive) with step dt."""
    n_steps = int(round(duration / dt)) + 1
    return np.linspace(0.0, duration, n_steps)


def ground_acceleration(t, g=9.81):
    """
    Compute ground acceleration ag(t).
    ag(t) = 0.30*g*sin(2*pi*1.2*t)*exp(-0.15*t)
          + 0.10*g*sin(2*pi*3.0*t)*exp(-0.20*t)
    """
    ag = (0.30 * g * np.sin(2.0 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
          + 0.10 * g * np.sin(2.0 * np.pi * 3.0 * t) * np.exp(-0.20 * t))
    return ag


def compute_damping(m, k, zeta):
    """Compute the viscous damping coefficient c = 2*zeta*sqrt(k*m)."""
    c = 2.0 * zeta * np.sqrt(k * m)
    return c


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    """
    Newmark average acceleration method (beta=1/4, gamma=1/2).
    Solves m*u_ddot + c*u_dot + k*u = -m*ag(t).

    Parameters
    ----------
    m : float - mass
    c : float - damping coefficient
    k : float - stiffness
    ag : ndarray - ground acceleration time history
    dt : float - time step
    u0 : float - initial displacement
    v0 : float - initial velocity

    Returns
    -------
    u : ndarray - relative displacement
    v : ndarray - relative velocity
    a_rel : ndarray - relative acceleration
    """
    beta = 0.25
    gamma = 0.5
    n = len(ag)

    u = np.zeros(n)
    v = np.zeros(n)
    a_rel = np.zeros(n)

    # Initial acceleration from EOM: m*a0 + c*v0 + k*u0 = -m*ag[0]
    a_rel[0] = (-m * ag[0] - c * v0 - k * u0) / m
    u[0] = u0
    v[0] = v0

    # Newmark constants
    a1 = 1.0 / (beta * dt ** 2)
    a2 = gamma / (beta * dt)
    a3 = 1.0 / (beta * dt)
    a4 = 1.0 / (2.0 * beta) - 1.0
    a5 = gamma / beta - 1.0
    a6 = dt * (gamma / (2.0 * beta) - 1.0)

    k_eff = k + a1 * m + a2 * c

    for i in range(n - 1):
        dp = (-m * ag[i + 1]
              + m * (a1 * u[i] + a3 * v[i] + a4 * a_rel[i])
              + c * (a2 * u[i] + a5 * v[i] + a6 * a_rel[i]))
        u[i + 1] = dp / k_eff
        v[i + 1] = a2 * (u[i + 1] - u[i]) + a5 * v[i] + a6 * a_rel[i]
        a_rel[i + 1] = (a1 * (u[i + 1] - u[i])
                         - a3 * v[i]
                         - a4 * a_rel[i])

    return u, v, a_rel


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    """
    Compute peak response metrics.

    Returns
    -------
    metrics : dict
    """
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    omega_n = np.sqrt(k / m)
    T_n = 2.0 * np.pi / omega_n
    c = compute_damping(m, k, zeta)

    idx_disp = int(np.argmax(np.abs(u)))
    metrics = {
        "natural_circular_frequency_rad_s": omega_n,
        "natural_period_s": T_n,
        "damping_coefficient_Ns_m": c,
        "peak_relative_displacement_m": float(np.max(np.abs(u))),
        "peak_relative_velocity_m_s": float(np.max(np.abs(v))),
        "peak_relative_acceleration_m_s2": float(np.max(np.abs(a_rel))),
        "peak_absolute_acceleration_m_s2": float(np.max(np.abs(a_abs))),
        "peak_pseudo_base_shear_N": float(k * np.max(np.abs(u))),
        "time_of_peak_displacement_s": float(t[idx_disp]),
    }
    return metrics


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    """Write time-history results to a CSV file."""
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s", "ag_m_s2", "u_m", "v_m_s", "a_rel_m_s2", "a_abs_m_s2"])
        for i in range(len(t)):
            writer.writerow([t[i], ag[i], u[i], v[i], a_rel[i], a_abs[i]])


def save_summary_json(path, metrics):
    """Write summary metrics to a JSON file."""
    with open(path, "w") as f:
        json.dump(metrics, f, indent=2)


def plot_response(path, t, ag, u, v, a_abs):
    """Generate and save a multi-panel response plot."""
    fig, axes = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    axes[0].plot(t, ag, color="tab:blue", linewidth=0.8)
    axes[0].set_ylabel("ag (m/s²)")
    axes[0].set_title("Ground Acceleration")
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, u * 1000.0, color="tab:red", linewidth=0.8)
    axes[1].set_ylabel("u (mm)")
    axes[1].set_title("Relative Displacement")
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t, v, color="tab:green", linewidth=0.8)
    axes[2].set_ylabel("v (m/s)")
    axes[2].set_title("Relative Velocity")
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(t, a_abs, color="tab:orange", linewidth=0.8)
    axes[3].set_ylabel("a_abs (m/s²)")
    axes[3].set_xlabel("Time (s)")
    axes[3].set_title("Absolute Acceleration")
    axes[3].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_analysis(output_dir="sdof_output"):
    """
    Execute the full SDOF Newmark analysis and write results to disk.

    Returns
    -------
    metrics : dict
    """
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
    omega_n = np.sqrt(k / m)
    T_n = 2.0 * np.pi / omega_n

    u, v, a_rel = newmark_average_acceleration(m, c, k, ag, dt, u0, v0)
    a_abs = a_rel + ag

    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)

    csv_path = os.path.join(output_dir, "response.csv")
    json_path = os.path.join(output_dir, "summary.json")
    png_path = os.path.join(output_dir, "response.png")

    save_csv(csv_path, t, ag, u, v, a_rel, a_abs)
    save_summary_json(json_path, metrics)
    plot_response(png_path, t, ag, u, v, a_abs)

    print("=== SDOF Newmark Analysis Results ===")
    print(f"  Natural frequency:       {omega_n:.4f} rad/s")
    print(f"  Natural period:          {T_n:.4f} s")
    print(f"  Damping coefficient:     {c:.2f} Ns/m")
    print(f"  Peak rel. displacement:  {metrics['peak_relative_displacement_m']*1000:.4f} mm")
    print(f"  Peak rel. velocity:      {metrics['peak_relative_velocity_m_s']:.4f} m/s")
    print(f"  Peak rel. acceleration:  {metrics['peak_relative_acceleration_m_s2']:.4f} m/s²")
    print(f"  Peak abs. acceleration:  {metrics['peak_absolute_acceleration_m_s2']:.4f} m/s²")
    print(f"  Peak pseudo base shear:  {metrics['peak_pseudo_base_shear_N']:.2f} N")
    print(f"  Time of peak disp.:      {metrics['time_of_peak_displacement_s']:.2f} s")
    print(f"  Output written to:       {output_dir}/")

    return metrics


def main():
    run_analysis()


# ---------------------------------------------------------------------------
# Assertion-based tests
# ---------------------------------------------------------------------------

def test_zero_ground_motion():
    """With zero ground motion the system remains at rest."""
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    dt = 0.01
    duration = 20.0
    c = compute_damping(m, k, zeta)

    t = make_time_vector(dt, duration)
    ag = np.zeros_like(t)
    u, v, a_rel = newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0)

    assert np.allclose(u, 0.0, atol=1e-14), "Displacement should be zero for zero ground motion"
    assert np.allclose(v, 0.0, atol=1e-14), "Velocity should be zero for zero ground motion"
    assert np.allclose(a_rel, 0.0, atol=1e-14), "Acceleration should be zero for zero ground motion"
    print("test_zero_ground_motion: PASSED")


def test_output_shapes():
    """All output arrays must have the same length as the time vector."""
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    dt = 0.01
    duration = 20.0
    c = compute_damping(m, k, zeta)

    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t, 9.81)
    u, v, a_rel = newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0)

    n = len(t)
    assert len(u) == n, f"u length {len(u)} != {n}"
    assert len(v) == n, f"v length {len(v)} != {n}"
    assert len(a_rel) == n, f"a_rel length {len(a_rel)} != {n}"
    print("test_output_shapes: PASSED")


def test_damping_positive():
    """Damping coefficient must be positive for positive mass, stiffness, and zeta."""
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    c = compute_damping(m, k, zeta)
    assert c > 0, f"Damping coefficient should be positive, got {c}"
    print("test_damping_positive: PASSED")


if __name__ == "__main__":
    test_zero_ground_motion()
    test_output_shapes()
    test_damping_positive()
    main()
```