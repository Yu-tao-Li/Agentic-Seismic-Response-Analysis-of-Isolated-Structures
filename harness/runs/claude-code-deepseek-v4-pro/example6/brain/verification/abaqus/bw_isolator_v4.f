C
C Bouc-Wen Isolator UEL v4: Rayleigh damping with limited tangent
C Damping force in RHS + limited damping contribution in AMATRX
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
      REAL*8 k_tan_h, k_tan_v, K_damp_add
      REAL*8 abs_z, abs_du_iso
      REAL*8 R_bw, J_bw, dz_corr
      REAL*8 v_rel, F_damp_h
      REAL*8 hht_gamma, hht_beta, dt_min, dt_use
      INTEGER bw_iter
      PARAMETER (hht_gamma = 0.55D0, hht_beta = 0.275625D0)
      PARAMETER (dt_min = 0.001D0)
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
      IF (LFLAGS(3).EQ.5) THEN
        DO I=1,NDOFEL
          DO J=1,NDOFEL
            AMATRX(I,J) = 0.0D0
          END DO
        END DO
        RETURN
      END IF
C
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
      F_damp_h = C_damp * v_rel
C
      SVARS(1) = z_new
      SVARS(2) = u_cur
      SVARS(3) = F_h
      SVARS(4) = F_v
C
      RHS(1,1) = RHS(1,1) - F_h - F_damp_h
      RHS(3,1) = RHS(3,1) + F_h + F_damp_h
      RHS(2,1) = RHS(2,1) - F_v
      RHS(4,1) = RHS(4,1) + F_v
C
C     Tangent: elastic + limited damping contribution
C     Limit dt for damping tangent to avoid blow-up at small time steps
      k_tan_h = k0
      k_tan_v = kv
      dt_use = MAX(DTIME, dt_min)
      K_damp_add = C_damp * hht_gamma / (hht_beta * dt_use)
C     Cap the damping tangent to 5x elastic stiffness
      IF (K_damp_add .GT. 5.0D0 * k0) K_damp_add = 5.0D0 * k0
C
      DO I=1,NDOFEL
        DO J=1,NDOFEL
          AMATRX(I,J) = 0.0D0
        END DO
      END DO
C
      AMATRX(1,1) = k_tan_h + K_damp_add
      AMATRX(1,3) = -k_tan_h - K_damp_add
      AMATRX(3,1) = -k_tan_h - K_damp_add
      AMATRX(3,3) = k_tan_h + K_damp_add
      AMATRX(2,2) = k_tan_v
      AMATRX(2,4) = -k_tan_v
      AMATRX(4,2) = -k_tan_v
      AMATRX(4,4) = k_tan_v
C
      RETURN
      END
