I'm writing the complete `sdof_newmark_demo.py` file with all required functions, parameters, and output generation.

```python
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import json
import os


def make_time_vector(dt, duration):
    return np.arange(0, duration + dt, dt)


def ground_acceleration(t, g=9.81):
    ag = (0.30 * g * np.sin(2 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
          + 0.10 * g * np.sin(2 * np.pi * 3.0 * t) * np.exp(-0.20 * t))
    return ag


def compute_damping(m, k, zeta):
    omega_n = np.sqrt(k / m)
    c = 2 * zeta * m * omega_n
    return c, omega_n


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    beta = 1.0 / 4.0
    gamma = 1.0 / 2.0

    n = len(ag)
    u = np.zeros(n)
    v = np.zeros(n)
    a_rel = np.zeros(n)
    a_abs = np.zeros(n)

    u[0] = u0
    v[0] = v0
    a_rel[0] = (-c * v0 - k * u0 - m * ag[0]) / m
    a_abs[0] = a_rel[0] + ag[0]

    k_eff = m / (beta * dt * dt) + gamma * c / (beta * dt) + k

    for i in range(1, n):
        a1 = m / (beta * dt) + gamma * c / beta
        a2 = m / (2 * beta) + dt * c * (gamma / (2 * beta) - 1)

        u_pred = u[i - 1] + dt * v[i - 1] + (dt * dt) * (0.5 - beta) * a_rel[i - 1]
        u[i] = u_pred

        dp_eff = -m * ag[i]
        du = 0.0

        resid = dp_eff - (m * a_rel[i - 1] + c * v[i - 1] + k * u[i - 1])
        du += resid / k_eff

        dv = gamma / (beta * dt) * du - gamma / beta * v[i - 1] + dt * (1 - gamma / (2 * beta)) * a_rel[i - 1]
        da = 1.0 / (beta * dt * dt) * du - 1.0 / (beta * dt) * v[i - 1] - 1.0 / (2 * beta) * a_rel[i - 1]

        u[i] = u[i - 1] + du
        v[i] = v[i - 1] + dv
        a_rel[i] = a_rel[i - 1] + da

    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    omega_n = np.sqrt(200000.0 / 1000.0)
    Tn = 2 * np.pi / omega_n
    c_val = 2 * 0.05 * 1000.0 * omega_n

    idx_u_peak = np.argmax(np.abs(u))
    u_peak = np.abs(u[idx_u_peak])
    t_u_peak = t[idx_u_peak]

    v_peak = np.max(np.abs(v))
    a_rel_peak = np.max(np.abs(a_rel))
    a_abs_peak = np.max(np.abs(a_abs))

    k = 200000.0
    pseudo_base_shear = k * u_peak

    metrics = {
        "omega_n_rad_per_s": float(omega_n),
        "Tn_s": float(Tn),
        "c_Ns_per_m": float(c_val),
        "u_max_m": float(u_peak),
        "v_max_m_per_s": float(v_peak),
        "a_rel_max_m_per_s2": float(a_rel_peak),
        "a_abs_max_m_per_s2": float(a_abs_peak),
        "k_u_max_N": float(pseudo_base_shear),
        "t_u_peak_s": float(t_u_peak),
    }
    return metrics


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    header = "t,ag,u,v,a_rel,a_abs"
    data = np.column_stack((t, ag, u, v, a_rel, a_abs))
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.8e")


def save_summary_json(path, metrics):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)


def plot_response(path, t, ag, u, v, a_abs):
    fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(t, ag, linewidth=0.8, color="black")
    axes[0].set_ylabel("ag (m/s²)")
    axes[0].grid(True, linestyle=":", alpha=0.6)

    axes[1].plot(t, u, linewidth=0.8, color="steelblue")
    axes[1].set_ylabel("u (m)")
    axes[1].grid(True, linestyle=":", alpha=0.6)

    axes[2].plot(t, v, linewidth=0.8, color="darkgreen")
    axes[2].set_ylabel("v (m/s)")
    axes[2].grid(True, linestyle=":", alpha=0.6)

    axes[3].plot(t, a_abs, linewidth=0.8, color="darkred")
    axes[3].set_ylabel("a_abs (m/s²)")
    axes[3].set_xlabel("Time (s)")
    axes[3].grid(True, linestyle=":", alpha=0.6)

    fig.suptitle("SDOF Newmark Average Acceleration Response", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_analysis(output_dir="sdof_output"):
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    g = 9.81
    dt = 0.01
    duration = 20.0

    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t, g)
    c, omega_n = compute_damping(m, k, zeta)

    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt)

    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)

    os.makedirs(output_dir, exist_ok=True)

    save_csv(os.path.join(output_dir, "response.csv"), t, ag, u, v, a_rel, a_abs)
    save_summary_json(os.path.join(output_dir, "summary.json"), metrics)
    plot_response(os.path.join(output_dir, "response.png"), t, ag, u, v, a_abs)

    print("--- SDOF Newmark Average Acceleration Analysis ---")
    print(f"Natural circular frequency : {metrics['omega_n_rad_per_s']:.4f} rad/s")
    print(f"Natural period              : {metrics['Tn_s']:.4f} s")
    print(f"Damping coefficient         : {metrics['c_Ns_per_m']:.4f} N·s/m")
    print(f"Peak relative displacement  : {metrics['u_max_m']:.6e} m")
    print(f"Peak relative velocity      : {metrics['v_max_m_per_s']:.6e} m/s")
    print(f"Peak relative acceleration  : {metrics['a_rel_max_m_per_s2']:.6e} m/s²")
    print(f"Peak absolute acceleration  : {metrics['a_abs_max_m_per_s2']:.6e} m/s²")
    print(f"Peak pseudo base shear k*u  : {metrics['k_u_max_N']:.6e} N")
    print(f"Time of peak displacement   : {metrics['t_u_peak_s']:.4f} s")
    print(f"\nOutput written to: {output_dir}/")

    return metrics


def main():
    run_analysis("sdof_output")


def test_zero_ground_motion():
    m, c, k = 1000.0, 2000.0, 200000.0
    dt = 0.01
    t = np.arange(0, 5.0 + dt, dt)
    ag = np.zeros_like(t)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt)
    assert np.allclose(u, 0.0, atol=1e-12), "u should stay zero with zero input"
    assert np.allclose(v, 0.0, atol=1e-12), "v should stay zero with zero input"
    assert np.allclose(a_rel, 0.0, atol=1e-12), "a_rel should stay zero with zero input"


def test_output_shapes():
    m, c, k = 1000.0, 2000.0, 200000.0
    dt = 0.01
    duration = 5.0
    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt)
    assert len(u) == len(ag), "displacement length mismatch"
    assert len(v) == len(ag), "velocity length mismatch"
    assert len(a_rel) == len(ag), "rel accel length mismatch"
    assert len(a_abs) == len(ag), "abs accel length mismatch"
    assert (a_abs[0] - a_rel[0]) - ag[0] == 0.0, "a_abs = a_rel + ag mismatch"


def test_damping_positive():
    m, k, zeta = 1000.0, 200000.0, 0.05
    c, _ = compute_damping(m, k, zeta)
    assert c > 0.0, "damping coefficient must be positive"


if __name__ == "__main__":
    main()
```