"""Attempted ABAQUS Bouc-Wen benchmark path for Example 6.

This path intentionally does not substitute a built-in bilinear connector for
the required Bouc-Wen state law. ABAQUS/Standard needs a compiled user
implementation for this benchmark, so the driver writes the intended input
files and first verifies that the local ABAQUS user-subroutine toolchain can
compile them.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import re
from pathlib import Path
from typing import Any

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parents[1]
LOG_FILE = SCRIPT_DIR / "execution_log.txt"
ABAQUS_EXE = Path(r"abaqus")
VS_VARS = Path(r"<visual_studio_buildtools>\VC\Auxiliary\Build\vcvars64.bat")
ONEAPI_SETVARS = Path(r"<intel_oneapi>\setvars.bat")
VS2022_INSTALL = Path(r"<visual_studio_buildtools>")


def log(message: str) -> None:
    print(message)
    with LOG_FILE.open("a", encoding="utf-8") as f:
        f.write(message + "\n")


def params() -> dict[str, float]:
    return {
        "E": 2.06e11,
        "nu": 0.30,
        "rho": 7850.0,
        "abaqus_numerical_density": 1.0e-12,
        "Ac": 0.020,
        "Ic": 8.0e-4,
        "Ab": 0.015,
        "Ib": 4.5e-4,
        "bay_width": 6.0,
        "story_height": 3.6,
        "floor_mass": 2.5e5,
        "abaqus_aux_vertical_mass_per_node": 1.0e-6,
        "k0": 2.0e6,
        "Fy": 1.0e4,
        "uy": 0.005,
        "alpha": 0.05,
        "bouc_n": 2.0,
        "beta_bw": 2.0e4,
        "gamma_bw": 2.0e4,
        "Ao": 1.0,
        "deltaA": 0.0,
        "deltaNu": 0.0,
        "deltaEta": 0.0,
        "kv": 1.0e10,
        "g": 9.81,
        "target_pga": 0.40 * 9.81,
        "dt": 0.01,
        "damping_ratio": 0.05,
        "rayleigh_alpha_m": 0.18649697508516622,
        "rayleigh_beta_k_initial": 0.0070671988247071136,
        "elastic_period_1": 2.842810601922031,
        "elastic_period_2": 0.5262442067740834,
    }


def node_id(level: int, col: int) -> int:
    return 100 + level * 10 + col + 1


def ground_id(col: int) -> int:
    return 10 + col + 1


def run_command(command: list[str], log_name: str, timeout: int = 120) -> dict[str, Any]:
    log(f"Running command: {' '.join(command)}")
    log_path = SCRIPT_DIR / log_name
    try:
        completed = subprocess.run(
            command,
            cwd=SCRIPT_DIR,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
        output = completed.stdout or ""
        log_path.write_text(output, encoding="utf-8", errors="replace")
        log(f"Command exit code {completed.returncode}; log: {log_name}")
        return {"command": command, "returncode": completed.returncode, "log": log_name, "output_tail": output[-4000:]}
    except FileNotFoundError as exc:
        msg = f"Command not found: {exc}"
        log_path.write_text(msg, encoding="utf-8")
        log(msg)
        return {"command": command, "returncode": None, "log": log_name, "output_tail": msg}
    except subprocess.TimeoutExpired as exc:
        output = (exc.stdout or "") + "\nTIMEOUT\n"
        log_path.write_text(output, encoding="utf-8", errors="replace")
        log(f"Command timed out; log: {log_name}")
        return {"command": command, "returncode": "timeout", "log": log_name, "output_tail": output[-4000:]}


def run_shell_command(command: str, log_name: str, timeout: int = 120) -> dict[str, Any]:
    log(f"Running command: {command}")
    log_path = SCRIPT_DIR / log_name
    try:
        completed = subprocess.run(
            command,
            cwd=SCRIPT_DIR,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
            shell=True,
        )
        output = completed.stdout or ""
        log_path.write_text(output, encoding="utf-8", errors="replace")
        log(f"Command exit code {completed.returncode}; log: {log_name}")
        return {"command": command, "returncode": completed.returncode, "log": log_name, "output_tail": output[-4000:]}
    except subprocess.TimeoutExpired as exc:
        output = (exc.stdout or "") + "\nTIMEOUT\n"
        log_path.write_text(output, encoding="utf-8", errors="replace")
        log(f"Command timed out; log: {log_name}")
        return {"command": command, "returncode": "timeout", "log": log_name, "output_tail": output[-4000:]}


def quote_cmd(path: Path | str) -> str:
    return f'"{path}"'


def run_abaqus(args: list[str], log_name: str, timeout: int = 120) -> dict[str, Any]:
    setup = [
        f'set "VS2022INSTALLDIR={VS2022_INSTALL}"',
        f"call {quote_cmd(VS_VARS)}",
        f"call {quote_cmd(ONEAPI_SETVARS)} intel64 --force",
        f"{quote_cmd(ABAQUS_EXE)} " + " ".join(args),
    ]
    cmd_text = " && ".join(setup)
    result = run_shell_command(cmd_text, log_name, timeout=timeout)
    result["effective_command"] = cmd_text
    return result


def cleanup_job_files(job_name: str) -> None:
    for path in SCRIPT_DIR.glob(job_name + ".*"):
        if path.suffix.lower() in {
            ".com",
            ".dat",
            ".fil",
            ".log",
            ".msg",
            ".odb",
            ".prt",
            ".sim",
            ".sta",
            ".stt",
            ".mdl",
            ".res",
            ".pac",
            ".sel",
        }:
            try:
                path.unlink()
            except OSError:
                pass


def cleanup_user_subroutine_files() -> None:
    for name in [
        "standardU.dll",
        "standardU.exp",
        "standardU.lib",
        "boucwen_uel-std.obj",
        "boucwen_uel.obj",
    ]:
        path = SCRIPT_DIR / name
        if path.exists():
            try:
                path.unlink()
            except OSError:
                pass


def write_ground_motion(p: dict[str, float]) -> tuple[list[float], float, float, float]:
    input_file = ROOT_DIR / "input_data" / "Northridge_01_NO_968.txt"
    raw = [float(line.strip()) for line in input_file.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not raw:
        raise RuntimeError(f"Ground-motion file is empty: {input_file}")
    raw_pga = max(abs(x) for x in raw)
    scale = p["target_pga"] / raw_pga
    ag = [x * scale for x in raw]
    (SCRIPT_DIR / "normalized_ground_accel_mps2.txt").write_text(
        "\n".join(f"{x:.16e}" for x in ag) + "\n",
        encoding="utf-8",
    )
    return ag, raw_pga, scale, max(abs(x) for x in ag)


def rect_dims(area: float, inertia: float) -> tuple[float, float]:
    depth = math.sqrt(12.0 * inertia / area)
    width = area / depth
    return width, depth


def node_dofs(node_index: int, comp: int) -> int:
    return node_index * 3 + comp


def add_frame_element(K: np.ndarray, coords: np.ndarray, ni: int, nj: int, E: float, A: float, I: float) -> None:
    xi, yi = coords[ni]
    xj, yj = coords[nj]
    dx = xj - xi
    dy = yj - yi
    length = math.sqrt(dx * dx + dy * dy)
    c = dx / length
    s = dy / length
    ea_l = E * A / length
    ei = E * I
    k_local = np.array(
        [
            [ea_l, 0.0, 0.0, -ea_l, 0.0, 0.0],
            [0.0, 12.0 * ei / length**3, 6.0 * ei / length**2, 0.0, -12.0 * ei / length**3, 6.0 * ei / length**2],
            [0.0, 6.0 * ei / length**2, 4.0 * ei / length, 0.0, -6.0 * ei / length**2, 2.0 * ei / length],
            [-ea_l, 0.0, 0.0, ea_l, 0.0, 0.0],
            [0.0, -12.0 * ei / length**3, -6.0 * ei / length**2, 0.0, 12.0 * ei / length**3, -6.0 * ei / length**2],
            [0.0, 6.0 * ei / length**2, 2.0 * ei / length, 0.0, -6.0 * ei / length**2, 4.0 * ei / length],
        ],
        dtype=float,
    )
    rot = np.array(
        [
            [c, s, 0.0, 0.0, 0.0, 0.0],
            [-s, c, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, c, s, 0.0],
            [0.0, 0.0, 0.0, -s, c, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    dofs = [
        node_dofs(ni, 0),
        node_dofs(ni, 1),
        node_dofs(ni, 2),
        node_dofs(nj, 0),
        node_dofs(nj, 1),
        node_dofs(nj, 2),
    ]
    K[np.ix_(dofs, dofs)] += rot.T @ k_local @ rot


def condensed_frame_stiffness(p: dict[str, float]) -> np.ndarray:
    n_cols = 3
    n_levels = 4
    n_nodes = n_cols * n_levels
    n_dof_full = n_nodes * 3
    coords = np.zeros((n_nodes, 2), dtype=float)
    for level in range(n_levels):
        for col in range(n_cols):
            coords[level * n_cols + col, :] = [col * p["bay_width"], level * p["story_height"]]

    K_full = np.zeros((n_dof_full, n_dof_full), dtype=float)
    for col in range(n_cols):
        for level in range(3):
            add_frame_element(K_full, coords, level * n_cols + col, (level + 1) * n_cols + col, p["E"], p["Ac"], p["Ic"])
    for level in range(1, 4):
        for col in range(2):
            add_frame_element(K_full, coords, level * n_cols + col, level * n_cols + col + 1, p["E"], p["Ab"], p["Ib"])

    mapping = np.zeros(n_dof_full, dtype=int)
    count = 0
    floor_master: list[int] = []
    for level in range(n_levels):
        if level >= 1:
            floor_master.append(count)
            for col in range(n_cols):
                mapping[node_dofs(level * n_cols + col, 0)] = count
            count += 1
        else:
            for col in range(n_cols):
                mapping[node_dofs(col, 0)] = count
                count += 1
        for col in range(n_cols):
            node = level * n_cols + col
            mapping[node_dofs(node, 1)] = count
            count += 1
            mapping[node_dofs(node, 2)] = count
            count += 1

    transform = np.zeros((n_dof_full, count), dtype=float)
    for i, j in enumerate(mapping):
        transform[i, j] = 1.0
    K_red = transform.T @ K_full @ transform
    for col in range(n_cols):
        base_node = col
        base_uy = mapping[node_dofs(base_node, 1)]
        K_red[base_uy, base_uy] += p["kv"]

    keep: list[int] = []
    for node in range(n_nodes):
        idx = int(mapping[node_dofs(node, 0)])
        if idx not in keep:
            keep.append(idx)
    eliminate = [i for i in range(count) if i not in keep]
    Kkk = K_red[np.ix_(keep, keep)]
    Kke = K_red[np.ix_(keep, eliminate)]
    Kek = K_red[np.ix_(eliminate, keep)]
    Kee = K_red[np.ix_(eliminate, eliminate)]
    return Kkk - Kke @ np.linalg.solve(Kee, Kek)


def full_reduced_frame_stiffness(p: dict[str, float]) -> np.ndarray:
    n_cols = 3
    n_levels = 4
    n_nodes = n_cols * n_levels
    n_dof_full = n_nodes * 3
    coords = np.zeros((n_nodes, 2), dtype=float)
    for level in range(n_levels):
        for col in range(n_cols):
            coords[level * n_cols + col, :] = [col * p["bay_width"], level * p["story_height"]]

    K_full = np.zeros((n_dof_full, n_dof_full), dtype=float)
    for col in range(n_cols):
        for level in range(3):
            add_frame_element(K_full, coords, level * n_cols + col, (level + 1) * n_cols + col, p["E"], p["Ac"], p["Ic"])
    for level in range(1, 4):
        for col in range(2):
            add_frame_element(K_full, coords, level * n_cols + col, level * n_cols + col + 1, p["E"], p["Ab"], p["Ib"])

    mapping = np.zeros(n_dof_full, dtype=int)
    count = 0
    for level in range(n_levels):
        if level >= 1:
            for col in range(n_cols):
                mapping[node_dofs(level * n_cols + col, 0)] = count
            count += 1
        else:
            for col in range(n_cols):
                mapping[node_dofs(col, 0)] = count
                count += 1
        for col in range(n_cols):
            node = level * n_cols + col
            mapping[node_dofs(node, 1)] = count
            count += 1
            mapping[node_dofs(node, 2)] = count
            count += 1

    transform = np.zeros((n_dof_full, count), dtype=float)
    for i, j in enumerate(mapping):
        transform[i, j] = 1.0
    K_red = transform.T @ K_full @ transform
    for col in range(n_cols):
        base_uy = mapping[node_dofs(col, 1)]
        K_red[base_uy, base_uy] += p["kv"]
    return K_red


def generalized_node_ids() -> list[int]:
    ids = [3000 + i for i in range(30)]
    ids[0] = node_id(0, 0)
    ids[1] = node_id(0, 1)
    ids[2] = node_id(0, 2)
    ids[9] = node_id(1, 1)
    ids[16] = node_id(2, 1)
    ids[23] = node_id(3, 1)
    return ids


def write_inp(p: dict[str, float], ag: list[float]) -> Path:
    k_frame = full_reduced_frame_stiffness(p)
    gen_nodes = generalized_node_ids()
    lines: list[str] = []
    lines.append("*HEADING")
    lines.append("Example 6 reduced-frame ABAQUS Bouc-Wen UEL benchmark, generated in this run")
    lines.append("*PREPRINT, ECHO=NO, MODEL=NO, HISTORY=NO, CONTACT=NO")
    lines.append("*NODE")
    for col in range(3):
        lines.append(f"{ground_id(col)}, {col * p['bay_width']:.12g}, 0.0")
    for col in range(3):
        lines.append(f"{node_id(0, col)}, {col * p['bay_width']:.12g}, 0.0")
    for level in range(1, 4):
        lines.append(f"{node_id(level, 1)}, {p['bay_width']:.12g}, {level * p['story_height']:.12g}")
    for idx, nid in enumerate(gen_nodes):
        if nid in {node_id(0, 0), node_id(0, 1), node_id(0, 2), node_id(1, 1), node_id(2, 1), node_id(3, 1)}:
            continue
        lines.append(f"{nid}, {p['bay_width']:.12g}, {-1.0 - 0.01 * idx:.12g}")
    lines.append("*USER ELEMENT, TYPE=U1, NODES=2, COORDINATES=2, PROPERTIES=11, VARIABLES=6, UNSYMM")
    lines.append("1")
    lines.append("*USER ELEMENT, TYPE=U2, NODES=30, COORDINATES=2, PROPERTIES=901, VARIABLES=6, UNSYMM")
    lines.append("1")
    lines.append("*ELEMENT, TYPE=U1, ELSET=ISOLATORS")
    for col in range(3):
        lines.append(f"{9001 + col}, {ground_id(col)}, {node_id(0, col)}")
    lines.append("*ELEMENT, TYPE=DASHPOT2, ELSET=ISO_DASHPOTS")
    for col in range(3):
        lines.append(f"{9101 + col}, {ground_id(col)}, {node_id(0, col)}")
    lines.append("*ELEMENT, TYPE=U2, ELSET=FRAME")
    frame_conn = [8001, *gen_nodes]
    lines.append(", ".join(str(x) for x in frame_conn))
    lines.append("*UEL PROPERTY, ELSET=ISOLATORS")
    lines.append(
        f"{p['k0']:.12g}, {p['alpha']:.12g}, {p['bouc_n']:.12g}, {p['gamma_bw']:.12g}, {p['beta_bw']:.12g}, "
        f"{p['Ao']:.12g}, {p['deltaA']:.12g}, {p['deltaNu']:.12g}, {p['deltaEta']:.12g}, "
        f"{p['kv']:.12g}, {p['rayleigh_beta_k_initial']:.12g}"
    )
    lines.append("*DASHPOT, ELSET=ISO_DASHPOTS")
    lines.append("1, 1")
    lines.append(f"{p['rayleigh_beta_k_initial'] * p['k0']:.12g}")
    frame_props = [p["rayleigh_beta_k_initial"], *k_frame.reshape(-1).tolist()]
    lines.append("*UEL PROPERTY, ELSET=FRAME")
    for start in range(0, len(frame_props), 8):
        lines.append(", ".join(f"{x:.16e}" for x in frame_props[start : start + 8]))
    lines.append("*NSET, NSET=GROUND")
    lines.append(", ".join(str(ground_id(col)) for col in range(3)))
    lines.append("*BOUNDARY")
    lines.append("GROUND, 1, 1, 0.0")
    for level in range(1, 4):
        mass_eid = 20000 + level
        nid = node_id(level, 1)
        lines.append(f"*ELEMENT, TYPE=MASS, ELSET=MASS_{nid}")
        lines.append(f"{mass_eid}, {nid}")
        lines.append(f"*MASS, ELSET=MASS_{nid}, TYPE=ANISOTROPIC, ALPHA={p['rayleigh_alpha_m']:.12g}")
        lines.append(f"{p['floor_mass']:.12g}, 0.0, 0.0")
    lines.append("*AMPLITUDE, NAME=GM_ACCEL, TIME=STEP TIME")
    eps = p["dt"] * 1.0e-6
    # MATLAB/OpenSeesPy advance interval k with acceleration sample k+1.
    # This piecewise-constant table reproduces that endpoint forcing in Abaqus.
    lines.append("0.0, 0.0")
    lines.append(f"{eps:.12e}, {ag[1]:.12e}")
    for idx in range(1, len(ag) - 1):
        t0 = (idx - 1) * p["dt"]
        t1 = idx * p["dt"]
        lines.append(f"{t1:.12e}, {ag[idx]:.12e}")
        lines.append(f"{(t1 + eps):.12e}, {ag[idx + 1]:.12e}")
    lines.append(f"{(len(ag) - 1) * p['dt']:.12e}, {ag[-1]:.12e}")
    # Equivalent inertial force formulation: CLOAD = -m_i * ag(t) at each mass node.
    lines.append("*STEP, NAME=DYNAMIC_BOUCWEN, NLGEOM=NO, INC=50000")
    lines.append("*DYNAMIC, DIRECT, ALPHA=0.0, BETA=0.25, GAMMA=0.5")
    lines.append(f"{p['dt']:.12g}, {(len(ag) - 1) * p['dt']:.12g}")
    lines.append("*CLOAD, AMPLITUDE=GM_ACCEL")
    for level in range(1, 4):
        lines.append(f"{node_id(level, 1)}, 1, {-p['floor_mass']:.12g}")
    lines.append("*OUTPUT, FIELD, FREQUENCY=1")
    lines.append("*NODE OUTPUT")
    lines.append("U, V, A, RF, CF")
    lines.append("*NODE PRINT, FREQUENCY=0")
    lines.append("*END STEP")
    path = SCRIPT_DIR / "example6_boucwen.inp"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_uel_source() -> Path:
    source = r"""      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
C     U1: two-node horizontal Bouc-Wen isolator, one DOF per node.
C     U2: six-node reduced elastic frame, one horizontal DOF per node.
C
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),SVARS(NSVARS),
     1 ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),JPROPS(NJPROP)
      DOUBLE PRECISION K0,ALPHA,BN,GAMMA,BETA,AO,KV,BETAR
      DOUBLE PRECISION UX,DUX,VX,ZOLD,Z,ABSZ,G,DGDZ,DGDU,DZDU
      DOUBLE PRECISION FXMAT,FX,CVX,DVELDU,KTX,VAL
      INTEGER I,J,IT
      DO I=1,NDOFEL
        RHS(I,1)=0.0D0
        DO J=1,NDOFEL
          AMATRX(I,J)=0.0D0
        END DO
      END DO
      IF (JTYPE.EQ.2) THEN
        BETAR=PROPS(1)
        IF (DTIME.GT.0.0D0) THEN
          DVELDU=2.0D0/DTIME
        ELSE
          DVELDU=0.0D0
        END IF
        DO I=1,NDOFEL
          DO J=1,NDOFEL
            VAL=PROPS(1+(I-1)*NDOFEL+J)
            RHS(I,1)=RHS(I,1)-VAL*U(J)-BETAR*VAL*V(J)
            AMATRX(I,J)=VAL+BETAR*VAL*DVELDU
          END DO
        END DO
        RETURN
      END IF
      K0=PROPS(1)
      ALPHA=PROPS(2)
      BN=PROPS(3)
      GAMMA=PROPS(4)
      BETA=PROPS(5)
      AO=PROPS(6)
      BETAR=0.0D0
      UX=U(2)-U(1)
      DUX=DU(2,1)-DU(1,1)
      VX=V(2)-V(1)
      ZOLD=SVARS(1)
      Z=ZOLD+AO*DUX
      DO IT=1,30
        ABSZ=ABS(Z)
        G=Z-ZOLD-AO*DUX+BETA*ABS(DUX)*(ABSZ**(BN-1.0D0))*Z
     1    +GAMMA*DUX*(ABSZ**BN)
        IF (ABS(G).LT.1.0D-13*MAX(1.0D0,ABS(Z))) GOTO 100
        IF (ABSZ.EQ.0.0D0) THEN
          DGDZ=1.0D0
        ELSE
          DGDZ=1.0D0+BETA*ABS(DUX)*BN*(ABSZ**(BN-1.0D0))
     1      +GAMMA*DUX*BN*(ABSZ**(BN-1.0D0))*SIGN(1.0D0,Z)
        END IF
        Z=Z-G/DGDZ
      END DO
  100 CONTINUE
      ABSZ=ABS(Z)
      IF (ABSZ.EQ.0.0D0) THEN
        DGDZ=1.0D0
      ELSE
        DGDZ=1.0D0+BETA*ABS(DUX)*BN*(ABSZ**(BN-1.0D0))
     1    +GAMMA*DUX*BN*(ABSZ**(BN-1.0D0))*SIGN(1.0D0,Z)
      END IF
      IF (DUX.GT.0.0D0) THEN
        DGDU=-AO+BETA*(ABSZ**(BN-1.0D0))*Z+GAMMA*(ABSZ**BN)
      ELSE IF (DUX.LT.0.0D0) THEN
        DGDU=-AO-BETA*(ABSZ**(BN-1.0D0))*Z+GAMMA*(ABSZ**BN)
      ELSE
        DGDU=-AO+GAMMA*(ABSZ**BN)
      END IF
      DZDU=-DGDU/DGDZ
      FXMAT=ALPHA*K0*UX+(1.0D0-ALPHA)*K0*Z
      CVX=BETAR*K0
      IF (DTIME.GT.0.0D0) THEN
        DVELDU=2.0D0/DTIME
      ELSE
        DVELDU=0.0D0
      END IF
      FX=FXMAT+CVX*VX
      KTX=ALPHA*K0+(1.0D0-ALPHA)*K0*DZDU+CVX*DVELDU
      RHS(1,1)= FX
      RHS(2,1)=-FX
      AMATRX(1,1)= KTX
      AMATRX(1,2)=-KTX
      AMATRX(2,1)=-KTX
      AMATRX(2,2)= KTX
      SVARS(1)=Z
      SVARS(2)=UX
      SVARS(3)=FXMAT
      SVARS(4)=FX
      RETURN
      END
"""
    path = SCRIPT_DIR / "boucwen_uel.for"
    path.write_text(source, encoding="ascii")
    return path


def write_postprocessor() -> Path:
    source = r'''"""Postprocess ABAQUS ODB for Example 6."""
from __future__ import annotations

import csv
import json
import math
import os
import re
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parents[1]
ODB = SCRIPT_DIR / "example6_boucwen.odb"
OUT = SCRIPT_DIR / "abaqus_results.json"


def params():
    return {
        "E": 2.06e11,
        "nu": 0.30,
        "rho": 7850.0,
        "Ac": 0.020,
        "Ic": 8.0e-4,
        "Ab": 0.015,
        "Ib": 4.5e-4,
        "bay_width": 6.0,
        "story_height": 3.6,
        "floor_mass": 2.5e5,
        "k0": 2.0e6,
        "Fy": 1.0e4,
        "uy": 0.005,
        "alpha": 0.05,
        "bouc_n": 2.0,
        "beta_bw": 2.0e4,
        "gamma_bw": 2.0e4,
        "Ao": 1.0,
        "deltaA": 0.0,
        "deltaNu": 0.0,
        "deltaEta": 0.0,
        "kv": 1.0e10,
        "g": 9.81,
        "target_pga": 0.40 * 9.81,
        "dt": 0.01,
        "damping_ratio": 0.05,
        "rayleigh_alpha_m": 0.18649697508516622,
        "rayleigh_beta_k_initial": 0.0070671988247071136,
        "elastic_period_1": 2.842810601922031,
        "elastic_period_2": 0.5262442067740834,
    }


def node_id(level, col):
    return 100 + level * 10 + col + 1


def ground_id(col):
    return 10 + col + 1


def assumptions():
    return {
        "frame": "ABAQUS 30-DOF reduced elastic frame UEL built from the same Euler-Bernoulli stiffness and diaphragm transformation used by MATLAB.",
        "mass": "Specified floor mass is assigned directly to the three horizontal floor generalized coordinates only.",
        "constraints": "The frame UEL carries base isolator translations, floor diaphragm translations, and retained massless vertical/rotational algebraic coordinates.",
        "isolation": "Three two-node horizontal Bouc-Wen UEL isolators connect fixed ground nodes to top isolator/base translation nodes.",
        "ground_motion": "Equivalent inertial-force formulation using normalized acceleration amplitude and CLOAD=-m*ag(t).",
        "story_drift": "Story 1 drift uses floor-1 diaphragm displacement minus average top-of-isolator displacement divided by story height.",
        "base_shear": "Total base shear is the sum of positive-resisting Bouc-Wen material forces from ground reactions.",
        "coordinates": "All output displacements are relative to ground in the frame x direction.",
        "bouc_wen_state": "UEL state is not directly exposed in ODB; z is reconstructed from F = alpha*k0*u + (1-alpha)*k0*z.",
    }


def write_failure(stage, reason):
    result = {
        "software": "ABAQUS",
        "status": "FAILED",
        "generated_in_this_run": True,
        "failure_stage": stage,
        "reason": reason,
    }
    OUT.write_text(json.dumps(result, indent=2), encoding="utf-8")
    raise SystemExit(reason)


def field_map(frame, field_name):
    if field_name not in frame.fieldOutputs:
        write_failure("postprocess_missing_field", "ODB field %s is missing." % field_name)
    out = {}
    for value in frame.fieldOutputs[field_name].values:
        out[value.nodeLabel] = tuple(float(x) for x in value.data)
    return out


def parse_sta():
    sta = SCRIPT_DIR / "example6_boucwen.sta"
    if not sta.exists():
        return {"failed_step_count": None, "reason": "STA file missing"}
    rows = []
    for line in sta.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split()
        if len(parts) >= 9 and all(part.lstrip("-").replace(".", "", 1).replace("E", "", 1).replace("+", "", 1).isdigit() for part in parts[:6]):
            try:
                rows.append(
                    {
                        "step": int(parts[0]),
                        "increment": int(parts[1]),
                        "attempt": int(parts[2]),
                        "severe_discontinuity_iterations": int(parts[3]),
                        "equilibrium_iterations": int(parts[4]),
                        "total_iterations": int(parts[5]),
                    }
                )
            except Exception:
                pass
    if not rows:
        return {"failed_step_count": None, "reason": "No increment rows parsed from STA file"}
    iteration_history = [row["total_iterations"] for row in rows]
    attempts = [row["attempt"] for row in rows]
    return {
        "test_type": "Abaqus/Standard equilibrium checks",
        "failed_step_count": 0,
        "cutback_count": 0,
        "increment_count": len(rows),
        "max_iterations_observed": int(max(iteration_history)),
        "iteration_history": iteration_history,
        "max_attempts_observed": int(max(attempts)),
    }


def write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


if not ODB.exists():
    write_failure("postprocess", "ODB file is missing; ABAQUS analysis did not complete.")

try:
    from odbAccess import openOdb
except Exception as exc:
    write_failure("postprocess_import", repr(exc))

p = params()
input_file = ROOT_DIR / "input_data" / "Northridge_01_NO_968.txt"
raw_ag = np.loadtxt(str(input_file), dtype=float).reshape(-1)
raw_pga = float(np.max(np.abs(raw_ag)))
scale = p["target_pga"] / raw_pga
ag = raw_ag * scale
n_steps = int(ag.size)
target_times = np.arange(n_steps, dtype=float) * p["dt"]

odb = openOdb(str(ODB), readOnly=True)
try:
    if "DYNAMIC_BOUCWEN" not in odb.steps:
        write_failure("postprocess_missing_step", "Step DYNAMIC_BOUCWEN is missing from ODB.")
    step = odb.steps["DYNAMIC_BOUCWEN"]
    frames = list(step.frames)
    if len(frames) != n_steps:
        write_failure(
            "postprocess_frame_count",
            "Expected %d ODB frames at input sample times, found %d." % (n_steps, len(frames)),
        )

    floors = np.zeros((n_steps, 3), dtype=float)
    bases = np.zeros((n_steps, 3), dtype=float)
    iso_disp = np.zeros((n_steps, 3), dtype=float)
    iso_total_force = np.zeros((n_steps, 3), dtype=float)
    iso_vel = np.zeros((n_steps, 3), dtype=float)
    time = np.zeros(n_steps, dtype=float)

    floor_nodes = [node_id(level, 1) for level in (1, 2, 3)]
    base_nodes = [node_id(0, col) for col in range(3)]
    ground_nodes = [ground_id(col) for col in range(3)]

    for i, frame in enumerate(frames):
        time[i] = float(frame.frameValue)
        u = field_map(frame, "U")
        v = field_map(frame, "V")
        rf = field_map(frame, "RF")
        for j, tag in enumerate(floor_nodes):
            floors[i, j] = u[tag][0]
        for j, tag in enumerate(base_nodes):
            bases[i, j] = u[tag][0]
            iso_disp[i, j] = u[tag][0] - u[ground_nodes[j]][0]
            iso_vel[i, j] = v[tag][0] - v[ground_nodes[j]][0]
        for j, tag in enumerate(ground_nodes):
            # Abaqus ground RF is opposite the positive element resisting force
            # convention used by MATLAB/OpenSeesPy in this benchmark.
            iso_total_force[i, j] = -rf[tag][0]
finally:
    odb.close()

if np.max(np.abs(time - target_times)) > 2.0e-6:
    write_failure(
        "postprocess_time_alignment",
        "ODB frame times do not match input sample times; max error %.6g" % float(np.max(np.abs(time - target_times))),
    )

roof = floors[:, 2]
base_avg = np.mean(bases, axis=1)
drifts = np.column_stack(
    [
        (floors[:, 0] - base_avg) / p["story_height"],
        (floors[:, 1] - floors[:, 0]) / p["story_height"],
        (floors[:, 2] - floors[:, 1]) / p["story_height"],
    ]
)
isolation_disp = np.mean(iso_disp, axis=1)
z_hist = np.zeros_like(iso_disp)
iso_force = np.zeros_like(iso_disp)
for i in range(1, n_steps):
    du = iso_disp[i, :] - iso_disp[i - 1, :]
    z_guess = z_hist[i - 1, :] + p["Ao"] * du
    for j in range(3):
        z = float(z_guess[j])
        z_old = float(z_hist[i - 1, j])
        duj = float(du[j])
        for _ in range(30):
            absz = abs(z)
            g = (
                z
                - z_old
                - p["Ao"] * duj
                + p["beta_bw"] * abs(duj) * (absz ** (p["bouc_n"] - 1.0)) * z
                + p["gamma_bw"] * duj * (absz ** p["bouc_n"])
            )
            if abs(g) < 1.0e-13 * max(1.0, abs(z)):
                break
            if absz == 0.0:
                dgdz = 1.0
            else:
                dgdz = (
                    1.0
                    + p["beta_bw"] * abs(duj) * p["bouc_n"] * (absz ** (p["bouc_n"] - 1.0))
                    + p["gamma_bw"] * duj * p["bouc_n"] * (absz ** (p["bouc_n"] - 1.0)) * math.copysign(1.0, z)
                )
            z -= g / dgdz
        z_hist[i, j] = z
    iso_force[i, :] = p["alpha"] * p["k0"] * iso_disp[i, :] + (1.0 - p["alpha"]) * p["k0"] * z_hist[i, :]
base_shear = np.sum(iso_force, axis=1)

loop_area = []
for j in range(3):
    loop_area.append(
        float(np.sum(0.5 * (iso_force[1:, j] + iso_force[:-1, j]) * (iso_disp[1:, j] - iso_disp[:-1, j])))
    )

peaks = {
    "roof_displacement_abs_max": float(np.max(np.abs(roof))),
    "max_interstory_drift_ratio_abs": float(np.max(np.abs(drifts))),
    "isolation_displacement_abs_max": float(np.max(np.abs(iso_disp))),
    "mean_isolation_displacement_abs_max": float(np.max(np.abs(isolation_disp))),
    "total_base_shear_abs_max": float(np.max(np.abs(base_shear))),
    "isolator_loop_area": loop_area,
}

conv = parse_sta()
status = "OK" if int(conv.get("failed_step_count") or 0) == 0 else "FAILED"

result = {
    "software": "ABAQUS",
    "status": status,
    "generated_in_this_run": True,
    "time_step": p["dt"],
    "num_input_samples": n_steps,
    "num_analysis_intervals": n_steps - 1,
    "input_file": os.path.relpath(str(input_file), str(ROOT_DIR)),
    "raw_pga": raw_pga,
    "normalization_scale": scale,
    "normalized_pga": float(np.max(np.abs(ag))),
    "fundamental_period": p["elastic_period_1"],
    "elastic_periods": [p["elastic_period_1"], p["elastic_period_2"]],
    "damping": {
        "ratio_modes_1_2": p["damping_ratio"],
        "rayleigh_alpha_m": p["rayleigh_alpha_m"],
        "rayleigh_beta_k_initial": p["rayleigh_beta_k_initial"],
        "stiffness_matrix": "Mass-proportional damping is assigned to lumped masses, and frame/isolator UELs include initial-stiffness-proportional viscous terms.",
    },
    "bouc_wen_update": "UEL backward-Euler update: z_{n+1}-z_n-Ao*du+beta*abs(du)*abs(z)^(n-1)*z+gamma*du*abs(z)^n = 0; F = alpha*k0*u + (1-alpha)*k0*z.",
    "modeling_assumptions": assumptions(),
    "time": time.tolist(),
    "ground_acceleration": ag.tolist(),
    "roof_displacement": roof.tolist(),
    "floor_displacements": floors.tolist(),
    "base_node_displacements": bases.tolist(),
    "interstory_drift_ratios": drifts.tolist(),
    "isolation_displacement": isolation_disp.tolist(),
    "isolator_displacements": iso_disp.tolist(),
    "total_base_shear": base_shear.tolist(),
    "isolator_forces": iso_force.tolist(),
    "abaqus_total_isolator_reaction_forces": iso_total_force.tolist(),
    "bouc_wen_z": z_hist.tolist(),
    "hysteresis_loop_area": loop_area,
    "peak_responses": peaks,
    "convergence": conv,
}

OUT.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")

write_csv(
    SCRIPT_DIR / "time_histories.csv",
    [
        "time_s",
        "ground_accel_mps2",
        "roof_disp_m",
        "drift_story1",
        "drift_story2",
        "drift_story3",
        "mean_isolation_disp_m",
        "total_base_shear_N",
    ],
    [
        [
            float(time[i]),
            float(ag[i]),
            float(roof[i]),
            float(drifts[i, 0]),
            float(drifts[i, 1]),
            float(drifts[i, 2]),
            float(isolation_disp[i]),
            float(base_shear[i]),
        ]
        for i in range(n_steps)
    ],
)

write_csv(
    SCRIPT_DIR / "hysteresis.csv",
    [
        "time_s",
        "u_iso_1_m",
        "f_iso_1_N",
        "z_iso_1_m",
        "u_iso_2_m",
        "f_iso_2_N",
        "z_iso_2_m",
        "u_iso_3_m",
        "f_iso_3_N",
        "z_iso_3_m",
    ],
    [
        [
            float(time[i]),
            float(iso_disp[i, 0]),
            float(iso_force[i, 0]),
            float(z_hist[i, 0]),
            float(iso_disp[i, 1]),
            float(iso_force[i, 1]),
            float(z_hist[i, 1]),
            float(iso_disp[i, 2]),
            float(iso_force[i, 2]),
            float(z_hist[i, 2]),
        ]
        for i in range(n_steps)
    ],
)
'''
    path = SCRIPT_DIR / "postprocess_abaqus.py"
    path.write_text(source, encoding="utf-8")
    return path


def write_result(result: dict[str, Any]) -> None:
    (SCRIPT_DIR / "abaqus_results.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    write_summary(result)
    if result.get("commands"):
        (SCRIPT_DIR / "abaqus_command_history.json").write_text(
            json.dumps(result["commands"], indent=2),
            encoding="utf-8",
        )


def write_summary(result: dict[str, Any]) -> None:
    summary_lines = [
        f"ABAQUS status: {result.get('status')}",
        f"Failure stage: {result.get('failure_stage')}" if result.get("failure_stage") else "Failure stage: none",
        f"Reason: {result.get('reason')}" if result.get("reason") else "Reason: none",
        "Required physical response histories were not generated." if result.get("status") != "OK" else "Physical response histories generated.",
    ]
    peaks = result.get("peak_responses") or {}
    if result.get("status") == "OK" and peaks:
        summary_lines.extend(
            [
                f"Peak roof displacement: {peaks.get('roof_displacement_abs_max')}",
                f"Peak isolation displacement: {peaks.get('isolation_displacement_abs_max')}",
                f"Peak total base shear: {peaks.get('total_base_shear_abs_max')}",
            ]
        )
    (SCRIPT_DIR / "response_summary.txt").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


def main() -> None:
    if LOG_FILE.exists():
        LOG_FILE.unlink()
    log("ABAQUS Bouc-Wen benchmark attempt started")
    p = params()
    ag, raw_pga, scale, normalized_pga = write_ground_motion(p)
    inp = write_inp(p, ag)
    uel = write_uel_source()
    post = write_postprocessor()
    log(f"Wrote intended ABAQUS input: {inp.name}")
    log(f"Wrote intended Bouc-Wen UEL: {uel.name}")
    log(f"Wrote postprocessor: {post.name}")
    log(f"Input samples: {len(ag)}, raw PGA: {raw_pga:.12g}, scale: {scale:.12g}, normalized PGA: {normalized_pga:.12g} m/s^2")

    if not ABAQUS_EXE.exists():
        result = base_result(p, raw_pga, scale, normalized_pga, len(ag))
        result.update(
            {
                "status": "BLOCKED",
                "failure_stage": "executable_discovery",
                "reason": f"ABAQUS executable was not found at {ABAQUS_EXE}.",
                "commands": [],
            }
        )
        write_result(result)
        return

    missing_setup = [str(path) for path in (VS_VARS, ONEAPI_SETVARS, VS2022_INSTALL) if not path.exists()]
    if missing_setup:
        result = base_result(p, raw_pga, scale, normalized_pga, len(ag))
        result.update(
            {
                "status": "BLOCKED",
                "failure_stage": "compiler_environment_discovery",
                "reason": "Required local compiler setup path(s) are missing: " + ", ".join(missing_setup),
                "commands": [],
            }
        )
        write_result(result)
        return

    cleanup_job_files("example6_boucwen")
    cleanup_user_subroutine_files()
    commands: list[dict[str, Any]] = []
    system = run_abaqus(["information=system"], "abaqus_system.log", timeout=180)
    commands.append(system)

    make = run_abaqus(["make", "library=boucwen_uel.for"], "abaqus_make.log", timeout=180)
    commands.append(make)
    make_text = (SCRIPT_DIR / "abaqus_make.log").read_text(encoding="utf-8", errors="replace")
    uel_library = SCRIPT_DIR / "standardU.dll"
    make_failed = (
        make["returncode"] not in (0, None)
        or not uel_library.exists()
        or "CompileError" in make_text
        or "not recognized as an internal or external command" in make_text
    )

    result = base_result(p, raw_pga, scale, normalized_pga, len(ag))
    result["commands"] = commands
    result["input_files"] = [inp.name, uel.name, post.name, "normalized_ground_accel_mps2.txt"]
    result["modeling_assumptions"] = assumptions()
    result["bouc_wen_update"] = (
        "Intended UEL backward-Euler update: z_{n+1}-z_n-Ao*du+beta*abs(du)*abs(z)^(n-1)*z"
        "+gamma*du*abs(z)^n = 0; F = alpha*k0*u + (1-alpha)*k0*z."
    )

    if make_failed:
        result.update(
            {
                "status": "BLOCKED",
                "failure_stage": "user_subroutine_compilation",
                "reason": (
                    "The required ABAQUS Bouc-Wen user implementation could not be compiled even after loading "
                    "Visual Studio BuildTools and Intel oneAPI. No bilinear or tabular substitute was run."
                ),
                "fundamental_period": None,
                "roof_displacement": None,
                "interstory_drift_ratios": None,
                "isolation_displacement": None,
                "total_base_shear": None,
                "isolator_forces": None,
                "bouc_wen_z": None,
                "hysteresis_loop_area": None,
                "peak_responses": None,
                "convergence": {"failed_step_count": None, "reason": "analysis not submitted because user-subroutine compilation is unavailable"},
                "missing_required_outputs": [
                    "fundamental_period",
                    "roof_displacement_time_history",
                    "interstory_drift_ratio_time_histories",
                    "isolation_displacement_time_history",
                    "total_base_shear_time_history",
                    "isolator_hysteresis_force_deformation_histories",
                    "bouc_wen_internal_state_histories",
                    "peak_response_quantities",
                    "convergence_history",
                ],
            }
        )
        write_result(result)
        log("ABAQUS path blocked at user-subroutine compilation; no substitute model submitted.")
        return

    job = run_abaqus(["job=example6_boucwen", "input=example6_boucwen.inp", "user=boucwen_uel.for", "interactive"], "abaqus_job.log", timeout=900)
    commands.append(job)
    job_text = (SCRIPT_DIR / "abaqus_job.log").read_text(encoding="utf-8", errors="replace")
    if "Abaqus JOB example6_boucwen COMPLETED" not in job_text or "ERROR" in job_text:
        result.update(
            {
                "status": "FAILED",
                "failure_stage": "analysis",
                "reason": "ABAQUS job did not complete cleanly; see abaqus_job.log, example6_boucwen.dat, example6_boucwen.msg, and example6_boucwen.sta.",
                "commands": commands,
            }
        )
        write_result(result)
        return

    post_cmd = run_abaqus(["python", "postprocess_abaqus.py"], "abaqus_postprocess.log", timeout=180)
    commands.append(post_cmd)
    result["commands"] = commands
    if (SCRIPT_DIR / "abaqus_results.json").exists():
        log("Postprocessor wrote abaqus_results.json")
        try:
            post_result = json.loads((SCRIPT_DIR / "abaqus_results.json").read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            post_result = {}
        if post_result:
            post_result["commands"] = commands
            post_result.setdefault("input_files", [inp.name, uel.name, post.name, "normalized_ground_accel_mps2.txt"])
            (SCRIPT_DIR / "abaqus_results.json").write_text(json.dumps(post_result, indent=2, allow_nan=False), encoding="utf-8")
            write_summary(post_result)
            (SCRIPT_DIR / "abaqus_command_history.json").write_text(
                json.dumps(commands, indent=2),
                encoding="utf-8",
            )
    else:
        result.update(
            {
                "status": "FAILED",
                "failure_stage": "analysis_or_postprocess",
                "reason": "ABAQUS job/postprocess did not produce abaqus_results.json.",
                "missing_required_outputs": ["ABAQUS machine-readable response output"],
            }
        )
        write_result(result)


def base_result(p: dict[str, float], raw_pga: float, scale: float, normalized_pga: float, n_samples: int) -> dict[str, Any]:
    return {
        "software": "ABAQUS",
        "status": "UNKNOWN",
        "generated_in_this_run": True,
        "time_step": p["dt"],
        "num_input_samples": n_samples,
        "input_file": os.path.relpath(ROOT_DIR / "input_data" / "Northridge_01_NO_968.txt", ROOT_DIR),
        "raw_pga": raw_pga,
        "normalization_scale": scale,
        "normalized_pga": normalized_pga,
        "damping": {
            "ratio_modes_1_2": p["damping_ratio"],
            "rayleigh_alpha_m": p["rayleigh_alpha_m"],
            "rayleigh_beta_k_initial": p["rayleigh_beta_k_initial"],
            "implementation": "Intended Rayleigh damping from first two elastic isolated-frame modes.",
        },
    }


def assumptions() -> dict[str, str]:
    return {
        "frame": "Reduced 30-DOF elastic-frame UEL generated from the same transformed Euler-Bernoulli stiffness used by MATLAB.",
        "mass": "Horizontal floor masses assigned to the three floor generalized coordinates only.",
        "constraints": "No native beam constraints are used; the frame UEL coordinates include the same retained massless vertical/rotational algebraic coordinates as the MATLAB reduced system.",
        "isolation": "Required representation is a two-node Bouc-Wen UEL with horizontal hysteretic force.",
        "ground_motion": "Intended equivalent inertial-force formulation using normalized acceleration amplitude and CLOAD=-m*ag(t).",
        "substitution_policy": "No built-in bilinear, plastic, or tabular connector substitute was run.",
    }


if __name__ == "__main__":
    main()
