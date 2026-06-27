C
C Bouc-Wen Isolator UEL for ABAQUS/Standard
C 2-node element with 2 DOFs per node (UX, UY)
C Horizontal (DOF 1): Bouc-Wen hysteretic
C Vertical (DOF 2): Linear elastic
C
C Properties:
C   PROPS(1) = k0    (initial horizontal stiffness, N/m)
C   PROPS(2) = kv    (vertical stiffness, N/m)
C   PROPS(3) = alpha (post-yield stiffness ratio)
C   PROPS(4) = n_bw  (Bouc-Wen smoothness exponent)
C   PROPS(5) = beta  (Bouc-Wen parameter)
C   PROPS(6) = gamma (Bouc-Wen parameter)
C   PROPS(7) = Ao    (Bouc-Wen reference amplitude)
C
C State variables (per element, saved after convergence):
C   SVARS(1) = z     (Bouc-Wen hysteretic state, length units)
C   SVARS(2) = u_conv (last converged isolator displacement)
C   SVARS(3) = F_h   (horizontal isolator force, N)
C   SVARS(4) = F_v   (vertical isolator force, N)
C
C CRITICAL: During Newton iterations within an increment, we must always
C compute the Bouc-Wen state from the LAST CONVERGED state, not from
C intermediate Newton states. SAVE variables track the converged state
C at the start of each increment.
C
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1  PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2  KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3  LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1  SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2  DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3  JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),
     4  LFLAGS(*),JPROPS(*)
C
C     --- Local variables ---
      REAL*8 k0, kv, alpha, n_bw, beta_bw, gamma_bw, Ao
      REAL*8 u_cur, du_iso, z_conv, u_conv, z_new, F_h, F_v
      REAL*8 k_tan_h, k_tan_v
      REAL*8 abs_z, abs_du_iso
      REAL*8 R_bw, J_bw, dz_corr
      INTEGER bw_iter
C
C     --- SAVE: persist across UEL calls, per-element ---
      INTEGER, SAVE :: last_kinc(1000)
      REAL*8,  SAVE :: u_conv_save(1000), z_conv_save(1000)
      INTEGER ielem, ikinc
      LOGICAL is_new_increment
C
C     Extract properties
      k0       = PROPS(1)
      kv       = PROPS(2)
      alpha    = PROPS(3)
      n_bw     = PROPS(4)
      beta_bw  = PROPS(5)
      gamma_bw = PROPS(6)
      Ao       = PROPS(7)
C
C     Extract current displacements
C     Node 1: U(1)=ux1, U(2)=uy1  (bottom / ground)
C     Node 2: U(3)=ux2, U(4)=uy2  (top / isolator top)
      u_cur = U(3) - U(1)
C
C     Use element number for save array indexing
      ielem = JELEM
C
      IF (LFLAGS(3).EQ.4) THEN
C       --- Stiffness matrix request (eigenvalue analysis) ---
        k_tan_h = k0
        k_tan_v = kv
C
        DO I=1,NDOFEL
          DO J=1,NDOFEL
            AMATRX(I,J) = 0.0D0
          END DO
        END DO
C
        AMATRX(1,1) = AMATRX(1,1) + k_tan_h
        AMATRX(1,3) = AMATRX(1,3) - k_tan_h
        AMATRX(3,1) = AMATRX(3,1) - k_tan_h
        AMATRX(3,3) = AMATRX(3,3) + k_tan_h
C
        AMATRX(2,2) = AMATRX(2,2) + k_tan_v
        AMATRX(2,4) = AMATRX(2,4) - k_tan_v
        AMATRX(4,2) = AMATRX(4,2) - k_tan_v
        AMATRX(4,4) = AMATRX(4,4) + k_tan_v
C
        RETURN
      END IF
C
      IF (LFLAGS(3).EQ.5) THEN
C       --- Mass matrix request ---
        DO I=1,NDOFEL
          DO J=1,NDOFEL
            AMATRX(I,J) = 0.0D0
          END DO
        END DO
        RETURN
      END IF
C
C     --- Normal increment: compute RHS and AMATRX ---
C
C     Detect new increment: when KINC changes, save the converged state
      IF (ielem .LE. 1000) THEN
        ikinc = KINC
        is_new_increment = (ikinc .NE. last_kinc(ielem))
        IF (is_new_increment) THEN
          last_kinc(ielem) = ikinc
          u_conv_save(ielem) = SVARS(2)
          z_conv_save(ielem) = SVARS(1)
        END IF
        u_conv = u_conv_save(ielem)
        z_conv = z_conv_save(ielem)
      ELSE
C       Fallback: use SVARS directly (may corrupt state during iterations)
        u_conv = SVARS(2)
        z_conv = SVARS(1)
      END IF
C
C     Compute incremental displacement from LAST CONVERGED state
      du_iso = u_cur - u_conv
      abs_du_iso = ABS(du_iso)
C
C     Bouc-Wen state update: backward Euler, Newton sub-iteration
C     Always start from the last converged z (z_conv)
      z_new = z_conv
      DO bw_iter = 1, 50
        abs_z = ABS(z_new)
        R_bw = z_new - z_conv - Ao*du_iso
     1      + beta_bw * abs_du_iso * abs_z * z_new
     2      + gamma_bw * du_iso * z_new * z_new
        J_bw = 1.0D0
     1       + 2.0D0 * beta_bw * abs_du_iso * abs_z
     2       + 2.0D0 * gamma_bw * du_iso * z_new
        dz_corr = -R_bw / J_bw
        z_new = z_new + dz_corr
        IF (ABS(dz_corr) .LT. 1.0D-14) EXIT
      END DO
C
C     Compute horizontal force: F_h = alpha*k0*u + (1-alpha)*k0*z
      F_h = alpha * k0 * u_cur + (1.0D0 - alpha) * k0 * z_new
C
C     Tangent stiffness: use elastic k0 for robustness (matching MATLAB)
      k_tan_h = k0
C
C     Vertical force and stiffness
      F_v = kv * (U(4) - U(2))
      k_tan_v = kv
C
C     Update state variables (will be saved on convergence)
      SVARS(1) = z_new
      SVARS(2) = u_cur
      SVARS(3) = F_h
      SVARS(4) = F_v
C
C     --- Assemble RHS (internal force vector) ---
C     Node 1 (bottom): forces applied in negative direction
C     Node 2 (top): forces applied in positive direction
C
      RHS(1,1) = RHS(1,1) - F_h
      RHS(3,1) = RHS(3,1) + F_h
C
      RHS(2,1) = RHS(2,1) - F_v
      RHS(4,1) = RHS(4,1) + F_v
C
C     --- Assemble AMATRX (tangent stiffness matrix) ---
      DO I=1,NDOFEL
        DO J=1,NDOFEL
          AMATRX(I,J) = 0.0D0
        END DO
      END DO
C
      AMATRX(1,1) = AMATRX(1,1) + k_tan_h
      AMATRX(1,3) = AMATRX(1,3) - k_tan_h
      AMATRX(3,1) = AMATRX(3,1) - k_tan_h
      AMATRX(3,3) = AMATRX(3,3) + k_tan_h
C
      AMATRX(2,2) = AMATRX(2,2) + k_tan_v
      AMATRX(2,4) = AMATRX(2,4) - k_tan_v
      AMATRX(4,2) = AMATRX(4,2) - k_tan_v
      AMATRX(4,4) = AMATRX(4,4) + k_tan_v
C
      RETURN
      END
