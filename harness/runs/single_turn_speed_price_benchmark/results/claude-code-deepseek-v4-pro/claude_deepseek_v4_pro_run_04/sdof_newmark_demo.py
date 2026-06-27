I'll write the complete Python file with all 10 required functions plus 3 test functions.

```python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
import os
import csv


def make_time_vector(dt, duration):
    n_steps = int(duration / dt) + 1
    return np.linspace(0, duration, n_steps)


def ground_acceleration(t, g=9.81):
    ag = (0.30 * g * np.sin(2 * np.pi * 1.2 * t) * np.exp(-0.15 * t)
          + 0.10 * g * np.sin(2 * np.pi * 3.0 * t) * np.exp(-0.20 * t))
    return ag


def compute_damping(m, k, zeta):
    return 2.0 * zeta * np.sqrt(m * k)


def newmark_average_acceleration(m, c, k, ag, dt, u0=0.0, v0=0.0):
    beta = 0.25
    gamma = 0.5

    n = len(ag)
    u = np.zeros(n)
    v = np.zeros(n)
    a_rel = np.zeros(n)

    u[0] = u0
    v[0] = v0
    a_rel[0] = (-m * ag[0] - c * v0 - k * u0) / m

    m_eff = m + gamma * dt * c + beta * dt * dt * k

    for i in range(n - 1):
        u_tilde = u[i] + dt * v[i] + (0.5 - beta) * dt * dt * a_rel[i]
        v_tilde = v[i] + (1.0 - gamma) * dt * a_rel[i]

        a_rel[i + 1] = (-m * ag[i + 1] - c * v_tilde - k * u_tilde) / m_eff
        u[i + 1] = u_tilde + beta * dt * dt * a_rel[i + 1]
        v[i + 1] = v_tilde + gamma * dt * a_rel[i + 1]

    a_abs = a_rel + ag
    return u, v, a_rel, a_abs


def compute_response_metrics(t, ag, u, v, a_rel, a_abs):
    idx_peak_u = np.argmax(np.abs(u))
    peak_u = np.abs(u[idx_peak_u])
    peak_v = np.max(np.abs(v))
    peak_a_rel = np.max(np.abs(a_rel))
    peak_a_abs = np.max(np.abs(a_abs))

    omega_n = None
    T_n = None
    c = None

    idx_max_disp = np.argmax(np.abs(u))

    return {
        'peak_relative_displacement_m': float(peak_u),
        'peak_relative_velocity_m_s': float(peak_v),
        'peak_relative_acceleration_m_s2': float(peak_a_rel),
        'peak_absolute_acceleration_m_s2': float(peak_a_abs),
        'time_of_peak_displacement_s': float(t[idx_max_disp]),
    }


def save_csv(path, t, ag, u, v, a_rel, a_abs):
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time_s', 'ground_accel_m_s2', 'disp_m',
                         'velocity_m_s', 'rel_accel_m_s2', 'abs_accel_m_s2'])
        for i in range(len(t)):
            writer.writerow([
                f'{t[i]:.6f}',
                f'{ag[i]:.8e}',
                f'{u[i]:.8e}',
                f'{v[i]:.8e}',
                f'{a_rel[i]:.8e}',
                f'{a_abs[i]:.8e}',
            ])


def save_summary_json(path, metrics):
    with open(path, 'w') as f:
        json.dump(metrics, f, indent=2)


def plot_response(path, t, ag, u, v, a_abs):
    fig, axes = plt.subplots(4, 1, figsize=(10, 10), sharex=True)

    axes[0].plot(t, ag, linewidth=0.8)
    axes[0].set_ylabel('Ground Accel [m/s$^2$]')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(t, u * 1000, linewidth=0.8)
    axes[1].set_ylabel('Displacement [mm]')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(t, v, linewidth=0.8)
    axes[2].set_ylabel('Velocity [m/s]')
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(t, a_abs, linewidth=0.8)
    axes[3].set_ylabel('Abs. Accel [m/s$^2$]')
    axes[3].set_xlabel('Time [s]')
    axes[3].grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_analysis(output_dir='sdof_output'):
    m = 1000.0
    k = 200000.0
    zeta = 0.05
    g = 9.81
    dt = 0.01
    duration = 20.0
    u0 = 0.0
    v0 = 0.0

    omega_n = np.sqrt(k / m)
    T_n = 2.0 * np.pi / omega_n
    c = compute_damping(m, k, zeta)

    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t, g)

    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt, u0, v0)

    metrics = compute_response_metrics(t, ag, u, v, a_rel, a_abs)
    metrics['mass_kg'] = m
    metrics['stiffness_N_m'] = k
    metrics['damping_ratio'] = zeta
    metrics['damping_coefficient_N_s_m'] = c
    metrics['natural_circular_frequency_rad_s'] = omega_n
    metrics['natural_period_s'] = T_n
    metrics['time_step_s'] = dt
    metrics['duration_s'] = duration

    peak_pseudo_base_shear = k * metrics['peak_relative_displacement_m']
    metrics['peak_pseudo_base_shear_N'] = peak_pseudo_base_shear

    os.makedirs(output_dir, exist_ok=True)

    save_csv(os.path.join(output_dir, 'response.csv'), t, ag, u, v, a_rel, a_abs)
    save_summary_json(os.path.join(output_dir, 'summary.json'), metrics)
    plot_response(os.path.join(output_dir, 'response.png'), t, ag, u, v, a_abs)

    print(f'Natural circular frequency: {omega_n:.4f} rad/s')
    print(f'Natural period: {T_n:.4f} s')
    print(f'Damping coefficient: {c:.4f} N·s/m')
    print(f'Peak relative displacement: {metrics["peak_relative_displacement_m"]:.6f} m')
    print(f'Peak relative velocity: {metrics["peak_relative_velocity_m_s"]:.6f} m/s')
    print(f'Peak relative acceleration: {metrics["peak_relative_acceleration_m_s2"]:.6f} m/s^2')
    print(f'Peak absolute acceleration: {metrics["peak_absolute_acceleration_m_s2"]:.6f} m/s^2')
    print(f'Peak pseudo base shear (k*u): {peak_pseudo_base_shear:.2f} N')
    print(f'Time of peak displacement: {metrics["time_of_peak_displacement_s"]:.4f} s')
    print(f'Output written to: {output_dir}/')

    return metrics


def main():
    run_analysis(output_dir='sdof_output')


def test_zero_ground_motion():
    m, k, zeta = 1000.0, 200000.0, 0.05
    dt, duration = 0.01, 5.0
    c = compute_damping(m, k, zeta)
    t = make_time_vector(dt, duration)
    ag = np.zeros_like(t)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt, u0=0.01, v0=0.0)
    assert u[-1] < 0.01, 'Free vibration should decay from initial 0.01 m'
    print('test_zero_ground_motion: PASSED')


def test_output_shapes():
    m, k, zeta = 1000.0, 200000.0, 0.05
    dt, duration = 0.01, 5.0
    c = compute_damping(m, k, zeta)
    t = make_time_vector(dt, duration)
    ag = ground_acceleration(t, g=9.81)
    u, v, a_rel, a_abs = newmark_average_acceleration(m, c, k, ag, dt)
    n = len(t)
    assert len(u) == n, 'Displacement length mismatch'
    assert len(v) == n, 'Velocity length mismatch'
    assert len(a_rel) == n, 'Rel acceleration length mismatch'
    assert len(a_abs) == n, 'Abs acceleration length mismatch'
    print('test_output_shapes: PASSED')


def test_damping_positive():
    c = compute_damping(1000.0, 200000.0, 0.05)
    assert c > 0, 'Damping coefficient must be positive'
    print('test_damping_positive: PASSED')


if __name__ == '__main__':
    main()
```