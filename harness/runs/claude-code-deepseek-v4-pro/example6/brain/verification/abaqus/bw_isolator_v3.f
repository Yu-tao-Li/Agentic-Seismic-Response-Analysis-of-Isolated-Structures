C
C Bouc-Wen Isolator UEL for ABAQUS/Standard (v3)
C 2-node element with 2 DOFs per node (UX, UY)
C Horizontal (DOF 1): Bouc-Wen hysteretic with Rayleigh damping
C Vertical (DOF 2): Linear elastic
C
C v3: Damping force in RHS only (no damping tangent in AMATRX)
C     This avoids stiff damping tangent causing convergence failure
C
C Properties:
C   PROPS(1) = k0    (initial horizontal stiffness, N/m)
C   PROPS(2) = kv    (vertical stiffness, N/m)
C   PROPS(3) = alpha (post-yield stiffness ratio)
C   PROPS(4) = n_bw  (Bouc-Wen smoothness exponent)
C   PROPS(5) = beta  (Bouc-Wen parameter)
C   PROPS(6) = gamma (Bouc-Wen parameter)
C   PROPS(7) = Ao    (Bouc-Wen reference amplitude)
C   PROPS(8) = C_damp (horizontal damping coefficient, N*s/m)
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
      REAL*8 k0, kv, alpha, n_bw, beta_bw, gamma_bw, Ao, C_damp
      REAL*8 u_cur, du_iso, z_conv, u_conv, z_new, F_h, F_v
      REAL*8 k_tan_h, k_tan_v
      REAL*8 abs_z, abs_du_iso
      REAL*8 R_bw, J_bw, dz_corr
      REAL*8 v_rel, F_damp_h
      INTEGER bw_iter
C
      INTEGER, SAVE :: last_kinc(1000)
      REAL*8,  SAVE :: u_conv_save(1000), z_conv_save(1000)
      INTEGER ielem, ikinc
      LOGICAL is_new_increment
C
      k0       = PROPS(1)
      kv       = PROPS(2)
      alpha    = PROPS(3)
      n_bw     = PROPS(4)
      beta_bw  = PROPS(5)
      gamma_bw = PROPS(6)
      Ao       = PROPS(7)
      C_damp   = PROPS(8)
C
      u_cur = U(3) - U(1)
      v_rel = V(3) - V(1)
      ielem = JELEM
C
C     --- Eigenvalue: stiffness only ---
      IF (LFLAGS(3).EQ.4) THEN
        k_tan_h = k0
        k_tan_v = kv
        DO I=1,NDOFEL
          DO J=1,NDOFEL
            AMATRX(I,J) = 0.0D0
          END DO
        END DO
        AMATRX(1,1) = k_tan_h
        AMATRX(1,3) = -k_tan_h
        AMATRX(3,1) = -k_tan_h
        AMATRX(3,3) = k_tan_h
        AMATRX(2,2) = k_tan_v
        AMATRX(2,4) = -k_tan_v
        AMATRX(4,2) = -k_tan_v
        AMATRX(4,4) = k_tan_v
        RETURN
      END IF
C
C     --- Mass matrix: zero ---
      IF (LFLAGS(3).EQ.5) THEN
        DO I=1,NDOFEL
          DO J=1,NDOFEL
            AMATRX(I,J) = 0.0D0
          END DO
        END DO
        RETURN
      END IF
C
C     --- Normal increment ---
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
        u_conv = SVARS(2)
        z_conv = SVARS(1)
      END IF
C
      du_iso = u_cur - u_conv
      abs_du_iso = ABS(du_iso)
C
C     Bouc-Wen state update: backward Euler, Newton sub-iteration
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
      F_h = alpha * k0 * u_cur + (1.0D0 - alpha) * k0 * z_new
      F_v = kv * (U(4) - U(2))
C
C     Stiffness-proportional Rayleigh damping force
      F_damp_h = C_damp * v_rel
C
      SVARS(1) = z_new
      SVARS(2) = u_cur
      SVARS(3) = F_h
      SVARS(4) = F_v
C
C     Assemble RHS with damping
      RHS(1,1) = RHS(1,1) - F_h - F_damp_h
      RHS(3,1) = RHS(3,1) + F_h + F_damp_h
      RHS(2,1) = RHS(2,1) - F_v
      RHS(4,1) = RHS(4,1) + F_v
C
C     Assemble AMATRX (elastic stiffness only, no damping tangent)
      k_tan_h = k0
      k_tan_v = kv
C
      DO I=1,NDOFEL
        DO J=1,NDOFEL
          AMATRX(I,J) = 0.0D0
        END DO
      END DO
C
      AMATRX(1,1) = k_tan_h
      AMATRX(1,3) = -k_tan_h
      AMATRX(3,1) = -k_tan_h
      AMATRX(3,3) = k_tan_h
      AMATRX(2,2) = k_tan_v
      AMATRX(2,4) = -k_tan_v
      AMATRX(4,2) = -k_tan_v
      AMATRX(4,4) = k_tan_v
C
      RETURN
      END
