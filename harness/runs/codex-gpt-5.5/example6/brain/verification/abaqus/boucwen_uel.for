      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
C     Example 6 Bouc-Wen isolation element for Abaqus/Standard.
C     Two-node 2D user element, active DOFs [u1,u2] at each node:
C       U(1)=ground ux, U(2)=ground uy, U(3)=top ux, U(4)=top uy.
C
C     Horizontal force:
C       F = alpha*k0*u + (1-alpha)*k0*z + betaK*k0*v
C     with z updated by the same backward-Euler OpenSees-form equation used
C     by the MATLAB benchmark:
C       z_{n+1}-z_n-du*(Ao-|z_{n+1}|^n*(gamma+beta*sign(du*z))) = 0.
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),SVARS(*),
     1 ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
C
      DOUBLE PRECISION K0,ALPHA,NBW,GAMMA,BETA,A0,KV,BETAK
      DOUBLE PRECISION UREL,VREL,DUH,ZOLD,UOLD,ZNEW,DZDU,KT
      DOUBLE PRECISION FH,FV,CH,CV,DVDU,OFFDIAG,ABSZ,PSI,PHI
      DOUBLE PRECISION RES,DPHIDZ,DGDZ,DZ,SIGNDUZ,SIGNZ,TOL
      DOUBLE PRECISION UCOMM(20000),ZCOMM(20000)
      DOUBLE PRECISION UTRIAL(20000),ZTRIAL(20000)
      DOUBLE PRECISION ONE,ZERO
      INTEGER I,J,ITER,IDX,INIT,KARR(20000)
      SAVE UCOMM,ZCOMM,UTRIAL,ZTRIAL,KARR,INIT
C
      DATA INIT /0/
C
      PARAMETER (ZERO=0.0D0, ONE=1.0D0)
C
      DO I=1,NDOFEL
         DO J=1,NRHS
            RHS(I,J)=ZERO
         END DO
         DO J=1,NDOFEL
            AMATRX(I,J)=ZERO
         END DO
      END DO
C
      K0    = PROPS(1)
      ALPHA = PROPS(2)
      NBW   = PROPS(3)
      GAMMA = PROPS(4)
      BETA  = PROPS(5)
      A0    = PROPS(6)
      KV    = PROPS(7)
      BETAK = PROPS(8)
C
      UREL = U(3) - U(1)
      VREL = V(3) - V(1)
C
      IF (INIT .EQ. 0) THEN
         DO I=1,20000
            UCOMM(I) = ZERO
            ZCOMM(I) = ZERO
            UTRIAL(I) = ZERO
            ZTRIAL(I) = ZERO
            KARR(I) = -1
         END DO
         INIT = 1
      END IF
C
      IDX = JELEM
      IF (IDX .LT. 1 .OR. IDX .GT. 20000) IDX = 1
      IF (KARR(IDX) .EQ. -1) THEN
         UCOMM(IDX) = ZERO
         ZCOMM(IDX) = ZERO
         UTRIAL(IDX) = ZERO
         ZTRIAL(IDX) = ZERO
         KARR(IDX) = KINC
      ELSE IF (KARR(IDX) .NE. KINC) THEN
         UCOMM(IDX) = UTRIAL(IDX)
         ZCOMM(IDX) = ZTRIAL(IDX)
         KARR(IDX) = KINC
      END IF
C
      UOLD = UCOMM(IDX)
      ZOLD = ZCOMM(IDX)
      DUH = UREL - UOLD
      ZNEW = ZOLD
      TOL = 1.0D-12
C
C     Local Newton solve for the implicit Bouc-Wen state.
      DO ITER=1,40
         ABSZ = ABS(ZNEW)
         IF (ABSZ .LT. 1.0D-14) ABSZ = 1.0D-14
         IF (DUH*ZNEW .GT. ZERO) THEN
            SIGNDUZ = ONE
         ELSE
            SIGNDUZ = -ONE
         END IF
         IF (ZNEW .GT. ZERO) THEN
            SIGNZ = ONE
         ELSE
            SIGNZ = -ONE
         END IF
         PSI = GAMMA + BETA*SIGNDUZ
         PHI = A0 - ABSZ**NBW*PSI
         RES = ZNEW - ZOLD - DUH*PHI
         DPHIDZ = -NBW*ABSZ**(NBW-ONE)*SIGNZ*PSI
         DGDZ = ONE - DUH*DPHIDZ
         IF (ABS(DGDZ) .LT. 1.0D-18) DGDZ = SIGN(1.0D-18,DGDZ)
         DZ = -RES/DGDZ
         ZNEW = ZNEW + DZ
         IF (ABS(DZ) .LE. TOL*MAX(ONE,ABS(ZNEW))) GOTO 20
      END DO
   20 CONTINUE
C
      ABSZ = ABS(ZNEW)
      IF (ABSZ .LT. 1.0D-14) ABSZ = 1.0D-14
      IF (DUH*ZNEW .GT. ZERO) THEN
         SIGNDUZ = ONE
      ELSE
         SIGNDUZ = -ONE
      END IF
      IF (ZNEW .GT. ZERO) THEN
         SIGNZ = ONE
      ELSE
         SIGNZ = -ONE
      END IF
      PSI = GAMMA + BETA*SIGNDUZ
      PHI = A0 - ABSZ**NBW*PSI
      DPHIDZ = -NBW*ABSZ**(NBW-ONE)*SIGNZ*PSI
      DGDZ = ONE - DUH*DPHIDZ
      IF (ABS(DGDZ) .LT. 1.0D-18) DGDZ = SIGN(1.0D-18,DGDZ)
      DZDU = PHI/DGDZ
      KT = ALPHA*K0 + (ONE-ALPHA)*K0*DZDU
      CH = BETAK*K0
      FH = ALPHA*K0*UREL + (ONE-ALPHA)*K0*ZNEW + CH*VREL
      FV = KV*(U(4) - U(2))
      CV = ZERO
C
C     Abaqus/Standard dynamic direct supplies PARAMS(2)=Newmark beta.
C     Include d(vrel)/d(urel) for the velocity-proportional tangent.
      DVDU = ZERO
      IF (DTIME .GT. ZERO .AND. LFLAGS(1) .EQ. 11) THEN
         IF (PARAMS(2) .GT. ZERO) DVDU = PARAMS(3)/(PARAMS(2)*DTIME)
      END IF
      OFFDIAG = KT + CH*DVDU
C
C     Residual contribution is external minus internal force.
      RHS(1,1) =  FH
      RHS(3,1) = -FH
      RHS(2,1) =  FV
      RHS(4,1) = -FV
C
      AMATRX(1,1) =  OFFDIAG
      AMATRX(1,3) = -OFFDIAG
      AMATRX(3,1) = -OFFDIAG
      AMATRX(3,3) =  OFFDIAG
      AMATRX(2,2) =  KV + CV*DVDU
      AMATRX(2,4) = -KV - CV*DVDU
      AMATRX(4,2) = -KV - CV*DVDU
      AMATRX(4,4) =  KV + CV*DVDU
C
C     Store committed trial state and audit quantities.
      UTRIAL(IDX) = UREL
      ZTRIAL(IDX) = ZNEW
      SVARS(1) = UCOMM(IDX)
      SVARS(2) = ZCOMM(IDX)
      SVARS(3) = FH
      SVARS(4) = KT
      SVARS(5) = VREL
      SVARS(6) = DUH
C
      ENERGY(2) = ENERGY(2) + 0.5D0*FH*DUH
      RETURN
      END
