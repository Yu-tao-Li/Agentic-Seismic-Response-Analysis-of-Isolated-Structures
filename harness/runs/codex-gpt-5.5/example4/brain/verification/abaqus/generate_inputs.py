#!/usr/bin/env python
"""Generate ABAQUS input - direct earthquake force method.
Apply F_i = -m_i * ag(t) to each structural mass. Base fixed.
Equation: M*a_rel + C*v_rel + K*u_rel = -M*ag
"""
import os, sys
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import numpy as np

# Parameters
mass_isolation = 250e3
mass_stories = [270e3, 270e3, 180e3]
k_stories = [245e6, 195e6, 98e6]
fy_stories = [1225e3, 975e3, 490e3]
keq = 50e6
c_iso = 1000e3
dt = 0.01; g = 9.81; PGA = 0.4 * g

# Load ground acceleration
acc = np.loadtxt('../../input_data/Northridge_01_NO_968.txt')
acc = acc / np.max(np.abs(acc)) * PGA
n_steps = len(acc)

# Amplitude file: normalized to 1.0 at peak, applied with -m_i factor
with open('eq_ampl.txt', 'w') as f:
    for i, a in enumerate(acc):
        f.write(f'{i*dt:.10f}, {a:.10f}\n')

# Rayleigh coefficients
alpha_m = 0.505929
beta_k = 0.003477
c_m = [alpha_m * m for m in [mass_isolation] + mass_stories]
c_k_vals = [beta_k * keq] + [beta_k * k for k in k_stories]
dy = [fy_stories[i] / k_stories[i] for i in range(3)]

lines = []
def add(s): lines.append(s)

add("*HEADING")
add("4-story structure - direct force seismic analysis")
add("*NODE")
for i in range(1, 6):
    add(f"{i}, {i-1}.0, 0.0")
add("*NSET, NSET=BASE"); add("1")
add("*NSET, NSET=ISO"); add("2")
add("*NSET, NSET=S2"); add("3")
add("*NSET, NSET=S3"); add("4")
add("*NSET, NSET=S4"); add("5")
add("*NSET, NSET=ALLFLOORS"); add("2,3,4,5")

# Springs
add("*ELEMENT, TYPE=SPRING2, ELSET=EL-ISO")
add("101, 1, 2")
add("*SPRING, ELSET=EL-ISO")
add("1, 1")
add(f" {keq:.1f}")

for idx in range(3):
    fy = fy_stories[idx]; d_y = dy[idx]
    n_a = idx + 2; n_b = idx + 3
    add(f"*ELEMENT, TYPE=SPRING2, ELSET=EL-S{idx+2}")
    add(f"{201+idx}, {n_a}, {n_b}")
    add(f"*SPRING, ELSET=EL-S{idx+2}, NONLINEAR")
    add(f"1, 1")
    add(f" {-fy:.1f}, -100.0")
    add(f" {-fy:.1f}, {-d_y:.10f}")
    add(f" 0.0, 0.0")
    add(f" {fy:.1f}, {d_y:.10f}")
    add(f" {fy:.1f}, 100.0")

# Dashpots
add("*ELEMENT, TYPE=DASHPOT2, ELSET=EL-DASH-ISO")
add("301, 1, 2")
add("*DASHPOT, ELSET=EL-DASH-ISO")
add("1, 1")
add(f" {c_iso:.1f}")

for idx, ck in enumerate(c_k_vals):
    n_a = idx + 1; n_b = idx + 2
    add(f"*ELEMENT, TYPE=DASHPOT2, ELSET=EL-DK-{idx}")
    add(f"{401+idx}, {n_a}, {n_b}")
    add(f"*DASHPOT, ELSET=EL-DK-{idx}")
    add(f"1, 1")
    add(f" {ck:.1f}")

for idx, cm in enumerate(c_m):
    n_id = idx + 2
    add(f"*ELEMENT, TYPE=DASHPOT1, ELSET=EL-DM-{idx}")
    add(f"{501+idx}, {n_id}")
    add(f"*DASHPOT, ELSET=EL-DM-{idx}")
    add(f"1")
    add(f" {cm:.1f}")

# Mass elements (structure only, no large mass at base)
mass_list = [mass_isolation] + mass_stories
for idx, m in enumerate(mass_list):
    n_id = idx + 2
    add(f"*ELEMENT, TYPE=MASS, ELSET=EL-M-{idx}")
    add(f"{601+idx}, {n_id}")
    add(f"*MASS, ELSET=EL-M-{idx}")
    add(f" {m:.1f}")

# Amplitude
add("*AMPLITUDE, NAME=EQ, INPUT=eq_ampl.txt")

# Fix base completely, fix Y for all floors
add("*BOUNDARY")
add("BASE, 1, 2")
add("ALLFLOORS, 2, 2")

# Direct earthquake forces: F_i = -m_i * ag(t) on each floor
add("*STEP, NAME=SEISMIC, NLGEOM=NO, INC=2000")
add("*DYNAMIC, ALPHA=0.0, HAFTOL=1.0E6")
add(f"{dt:.6f}, {n_steps*dt:.6f}, 1e-15, {dt:.6f}")

add("*CLOAD, AMPLITUDE=EQ")
add(f"ISO, 1, {-mass_isolation:.1f}")
add(f"S2, 1, {-mass_stories[0]:.1f}")
add(f"S3, 1, {-mass_stories[1]:.1f}")
add(f"S4, 1, {-mass_stories[2]:.1f}")

add("*RESTART, WRITE, FREQUENCY=0")
add("*OUTPUT, FIELD, FREQUENCY=10")
add("*NODE OUTPUT, NSET=ALLFLOORS")
add("U, V, A")
add("*NODE OUTPUT, NSET=BASE")
add("RF")
add("*OUTPUT, HISTORY, FREQUENCY=1")
add("*NODE OUTPUT, NSET=S4")
add("U1, A1")
add("*NODE OUTPUT, NSET=ISO")
add("U1")
add("*END STEP")

with open('model.inp', 'w') as f:
    f.write('\n'.join(lines))

print(f"Generated model.inp ({n_steps} steps)")
print("Direct earthquake forces: F_i = -m_i * ag(t) on each structural mass")
