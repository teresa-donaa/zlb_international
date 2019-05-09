MODULE gatsm
!
USE constants
USE observations
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE theta_to_param ( theta, delta0, delta1_1, delta1_2, PhiQ_1, PhiQ_2, &
        muQ_1, muQ_2, Gamm, Sigma, Sigma_1, Sigma_2, mu, mu_1, mu_2, lambda_1, lambda_2, &
        Phi, Phi_1, Phi_2, LLambda_1, LLambda_2, Omega, omega_dex, r_inf )
    !
    ! Transforms theta to real parameters 
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)
    REAL(8), INTENT(OUT) :: delta0(num_C)
    REAL(8), DIMENSION(num_XC(1)), INTENT(OUT) :: delta1_1, mu_1, muQ_1, lambda_1
    REAL(8), DIMENSION(num_XC(2)), INTENT(OUT) :: delta1_2, mu_2, muQ_2, lambda_2
    REAL(8), INTENT(OUT) :: mu(num_X)
    REAL(8), DIMENSION(num_XC(1),num_XC(1)), INTENT(OUT) :: PhiQ_1, Sigma_1, Phi_1, LLambda_1
    REAL(8), DIMENSION(num_XC(2),num_XC(2)), INTENT(OUT) :: PhiQ_2, Sigma_2, Phi_2, LLambda_2
    REAL(8), DIMENSION(num_X,num_X), INTENT(OUT) :: Phi, Gamm, Sigma
    REAL(8), INTENT(OUT) :: Omega(num_C*num_K,num_C*num_K)
    REAL(8), INTENT(OUT) :: omega_dex
    REAL(8), INTENT(OUT) :: r_inf
    !
    ! Declaring local variables
    !
    INTEGER :: il, iu, h, cc, ir, ic, j
    REAL(8) :: tmpPhiQ(num_PhiQ), tmpOmega(num_C)
    REAL(8), DIMENSION(num_X,num_X) :: S1, S, A, K
    CHARACTER(len = 1) :: key
    ! 
    ! Beginning execution
    !
    ! delta0
    !
    il = 0
    IF (switch_delta01 .EQ. 1) THEN
        delta0 = theta(il+1:il+num_C)
    END IF
    !
    IF (switch_delta02 .EQ. 1) THEN
        delta0 = 0.d0
    END IF
    iu = num_delta0
    !
    ! delta1_1, delta1_2
    !
    il = iu
    delta1_1 = EXP(theta(il+1:il+num_XC(1)))
    delta1_2 = EXP(theta(il+num_XC(1)+1:il+num_XC(1)+num_XC(2)))
    iu = il+num_delta1
    !
    ! PhiQ_1, PhiQ_2
    !
    IF (switch_PhiQ1 .EQ. 1) THEN
        !
        il = iu
        PhiQ_1 = 0.d0
        h = 0
        DO ir = 1, num_XC(1)
            !
            DO ic = 1, ir
                !
                h = h+1
                PhiQ_1(ir,ic) = theta(il+h)
                IF (ir .EQ. ic) THEN
                    !
                    PhiQ_1(ir,ir) = (dueinvpi*ATAN(PhiQ_1(ir,ir))+1.d0)/2.d0
                    IF ((ir .GT. 1) .AND. (ir .LE. num_XG)) PhiQ_1(ir,ir) = PhiQ_1(ir-1,ir-1)*PhiQ_1(ir,ir)
                    IF ((ir .GT. (num_XG+1)) .AND. (ir .LE. num_XC(1))) &
                        PhiQ_1(ir,ir) = PhiQ_1(ir-1,ir-1)*PhiQ_1(ir,ir)
                    !
                END IF
                !
            END DO
            !
        END DO
        iu = iu+num_XC(1)*(num_XC(1)+1)/2
        !
        il = iu
        PhiQ_2 = 0.d0
        h = 0
        DO ir = 1, num_XC(2)
            !
            DO ic = 1, ir
                !
                h = h+1
                PhiQ_2(ir,ic) = theta(il+h)
                IF (ir .EQ. ic) THEN
                    !
                    PhiQ_2(ir,ir) = (dueinvpi*ATAN(PhiQ_2(ir,ir))+1.d0)/2.d0
                    IF ((ir .GT. 1) .AND. (ir .LE. num_XG)) PhiQ_2(ir,ir) = PhiQ_2(ir-1,ir-1)*PhiQ_2(ir,ir)
                    IF ((ir .GT. (num_XG+1)) .AND. (ir .LE. num_XC(2))) &
                        PhiQ_2(ir,ir) = PhiQ_2(ir-1,ir-1)*PhiQ_2(ir,ir)
                    !
                END IF
                !
            END DO
            !
        END DO
        iu = iu+num_XC(2)*(num_XC(2)+1)/2
    END IF
    !
    IF (switch_PhiQ2 .EQ. 1) THEN
        !
        il = iu
        PhiQ_1 = posdef_matrix(num_XC(1),theta(il+1:il+num_XC(1)**2))
        iu = iu+num_XC(1)**2
        !
        il = iu
        PhiQ_2 = posdef_matrix(num_XC(2),theta(il+1:il+num_XC(2)**2))
        iu = iu+num_XC(2)**2
    END IF
    !
    ! muQ
    !
    il = iu
    IF (switch_muQ1 .EQ. 1) THEN
        !
        muQ_1 = theta(il+1:il+num_XC(1))
        muQ_2 = theta(il+num_XC(1)+1:il+num_XC(1)+num_XC(2))
        !
    END IF
    !
    IF (switch_muQ2 .EQ. 1) THEN
        !
        muQ_1 = 0.d0
        muQ_2 = 0.d0
        !
    END IF
    iu = iu+num_muQ
    !
    ! Gamma
    !
    il = iu
    Gamm = 0.d0
    DO ir = 1, num_X
        !
        Gamm(ir,ir) = 1.d0
        !
    END DO
    iu = iu+num_Gamm
    !
    ! Sigma, Sigma_1, Sigma_2
    !
    Sigma = 0.d0
    DO ir = 1, num_X
        !
        DO ic = 1, ir
            !
            DO h = 1, num_X
                !
                Sigma(ir,ic) = Sigma(ir,ic)+Gamm(ir,h)*Gamm(ic,h)
                !
            END DO
            IF (ir .NE. ic) Sigma(ic,ir) = Sigma(ir,ic)
            !
        END DO
        !
    END DO
    Sigma_1 = Sigma(sel_X1,sel_X1)
    Sigma_2 = Sigma(sel_X2,sel_X2)
    !
    ! mu, mu_1, mu_2, lambda_1, lambda_2
    !
    il = iu
    IF (switch_mu1 .EQ. 1) mu = theta(il+1:il+num_X)
    IF (switch_mu2 .EQ. 1) mu = 0.d0
    iu = iu+num_mu
    mu_1 = mu(sel_X1)
    mu_2 = mu(sel_X2)
    lambda_1 = mu_1-muQ_1
    lambda_2 = mu_2-muQ_2
    !
    ! Phi, Phi_1, Phi_2, LLambda_1, LLambda_2
    !
    il = iu
    IF (switch_Phi1 .EQ. 1) Phi = posdef_matrix(num_X,theta(il+1:il+num_X**2))
    !
    IF (switch_Phi2 .EQ. 1) THEN
        !
        Phi = 0.d0
        h = 0
        DO ir = 1, num_X
            !
            DO ic = 1, ir
                ! Cycle if (ir,ic) corresponds to local factors
                ! Only works for models with (XG,1,1) factors
                IF ((ir .GT. num_XC(1)) .AND. (ic .GT. num_XG) .AND. (ic .LE. num_XC(1))) CYCLE
                h = h+1
                Phi(ir,ic) = theta(il+h)
                IF (ir .EQ. ic) THEN
                    !
                    Phi(ir,ir) = (dueinvpi*ATAN(Phi(ir,ir))+1.d0)/2.d0
                    IF ((ir .GT. 1) .AND. (ir .LE. num_XG)) Phi(ir,ir) = Phi(ir-1,ir-1)*Phi(ir,ir)
                    IF ((ir .GT. (num_XG+1)) .AND. (ir .LE. num_XC(1))) &
                        Phi(ir,ir) = Phi(ir-1,ir-1)*Phi(ir,ir)
                    IF ((ir .GT. (num_XC(1)+1)) .AND. (ir .LE. num_X)) &
                        Phi(ir,ir) = Phi(ir-1,ir-1)*Phi(ir,ir)
                    !
                END IF
                !
            END DO
            !
        END DO
        !
    END IF
    iu = iu+num_Phi
    Phi_1 = Phi(sel_X1,sel_X1)
    Phi_2 = Phi(sel_X2,sel_X2)
    LLambda_1 = Phi_1-PhiQ_1
    LLambda_2 = Phi_2-PhiQ_2
    !
    ! Omega
    !
    il = iu
    Omega = 0.d0
    tmpOmega = EXP(theta(il+1:il+num_C))
    DO cc = 1, num_C
        !
        DO ir = 1, num_K
            !
            j = (cc-1)*num_K+ir
            Omega(j,j) = tmpOmega(cc)**2
            !
        END DO
        !
    END DO
    iu = iu+num_Omega
    !
    ! omega_dex
    !
    il = iu
    omega_dex = EXP(theta(il+1))**2
    iu = iu+num_omega_dex
    !
    ! r_inf
    !
    il = iu
    IF (switch_r_inf0 .EQ. 1) r_inf = -99.d0
    !
    IF (switch_r_inf1 .EQ. 1) r_inf = r_inf_fixed
    !
    IF (switch_r_inf2 .EQ. 1) r_inf = theta(il+1)
    iu = iu+num_r_inf_true
    !
    ! End execution and returning control
    !
    END SUBROUTINE theta_to_param
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE Kalman_params ( delta0, delta1_1, delta1_2, PhiQ_1, PhiQ_2, muQ_1, muQ_2, &
        Sigma_1, Sigma_2, mu, Phi, Sigma, &
        avec_1, avec_2, Bmat_1, Bmat_2, sigmaQvec_1, sigmaQvec_2, muXI, PhiXI, SigmaXI )
    !
    ! Transforms real parameters to Kalman filter parameters 
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: delta0(num_C)
    REAL(8), DIMENSION(num_XC(1)), INTENT(IN) :: delta1_1, muQ_1
    REAL(8), DIMENSION(num_XC(2)), INTENT(IN) :: delta1_2, muQ_2
    REAL(8), DIMENSION(num_XC(1),num_XC(1)), INTENT(IN) :: PhiQ_1, Sigma_1
    REAL(8), DIMENSION(num_XC(2),num_XC(2)), INTENT(IN) :: PhiQ_2, Sigma_2
    REAL(8), INTENT(IN) :: mu(num_X)
    REAL(8), DIMENSION(num_X,num_X), INTENT(IN) :: Phi, Sigma
    REAL(8), INTENT(OUT) :: avec_1(num_K), avec_2(num_K)
    REAL(8), INTENT(OUT) :: Bmat_1(num_XC(1),num_K), Bmat_2(num_XC(2),num_K)
    REAL(8), INTENT(OUT) :: sigmaQvec_1(num_K), sigmaQvec_2(num_K)
    REAL(8), INTENT(OUT) :: muXI(num_XI)
    REAL(8), DIMENSION(num_XI,num_XI), INTENT(OUT) :: PhiXI, SigmaXI
    ! 
    ! Beginning execution
    !
    ! Measurement equation parameters
    !
    CALL Kalman_params_absigmaQvec(num_XC(1),delta0(1),delta1_1, &
        PhiQ_1,muQ_1,Sigma_1,avec_1,Bmat_1,sigmaQvec_1)
    CALL Kalman_params_absigmaQvec(num_XC(2),delta0(2),delta1_2, &
        PhiQ_2,muQ_2,Sigma_2,avec_2,Bmat_2,sigmaQvec_2)
    !
    ! Transition equation parameters
    !
    muXI(:num_X) = mu
    muXI(num_X+1:) = 0.d0
    !
    PhiXI(:num_X,:num_X) = Phi
    PhiXI(num_X+1:,:num_X) = eyeX
    PhiXI(:,num_X+1:) = 0.d0
    !
    SigmaXI = 0.d0
    SigmaXI(:num_X,:num_X) = Sigma    
    !
    ! End execution and returning control
    !
    END SUBROUTINE Kalman_params
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE Kalman_params_absigmaQvec ( n, delta0cc, delta1cc, PhiQcc, muQcc, Sigmacc, &
        aveccc, Bmatcc, sigmaQveccc )
    !
    ! Transforms real parameters to Kalman filter parameters 
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: delta0cc
    REAL(8), INTENT(IN) :: delta1cc(n)
    REAL(8), INTENT(IN) :: muQcc(n)
    REAL(8), INTENT(IN) :: PhiQcc(n,n)
    REAL(8), INTENT(IN) :: Sigmacc(n,n)
    REAL(8), INTENT(OUT) :: aveccc(num_K)
    REAL(8), INTENT(OUT) :: Bmatcc(n,num_K)
    REAL(8), INTENT(OUT) :: sigmaQveccc(num_K)
    !
    ! Declaring local variables
    !
    INTEGER :: iTau, ik, ir, ic
    REAL(8), DIMENSION(n) :: b, sum_b, prev_b
    REAL(8) :: sigmaQ2   
    !
    ! Beginning execution
    !
    ! Initialization: Tau = 0
    !
    b = delta1cc
    sum_b = 0.d0
    sigmaQ2 = 0.d0
    !
    ! Beginning loop over maturities
    !
    iTau = 0
    DO ik = 1, Tau_vec(num_K)
        !
        ! \sum_{j=0}^{ik-1} b_j
        !
        sum_b = sum_b+b
        !
        ! \sum_{j=0}^{ik-1} b_j' Sigma b_j
        !
        sigmaQ2 = 0.d0
        DO ir = 1, n
            DO ic = 1, n
                sigmaQ2 = sigmaQ2+b(ir)*b(ic)*Sigmacc(ir,ic)
            END DO
        END DO
        !
        ! b_ik
        !
        prev_b = b
        b = 0.d0
        DO ir = 1, n
            DO ic = 1, n
                b(ir) = b(ir)+PHIQcc(ic,ir)*prev_b(ic)
            END DO
        END DO
        !
        ! Check if an observed maturity has been reached
        !
        IF (ik .EQ. Tau_vec(iTau+1)) THEN
            iTau = iTau+1
            Bmatcc(:,iTau) = b
            aveccc(iTau) = delta0cc
            DO ir = 1, n
                aveccc(iTau) = aveccc(iTau)+sum_b(ir)*muQcc(ir)
                DO ic = 1, n
                    aveccc(iTau) = aveccc(iTau)-0.5d0*sum_b(ir)*Sigmacc(ir,ic)*sum_b(ic)
                END DO
            END DO
            sigmaQveccc(iTau) = sigmaQ2
        END IF
        !
    END DO
    !    
    ! Computing sigmaQvec
    !
    sigmaQveccc = SQRT(sigmaQveccc)
    !
    ! End execution and returning control
    !
    END SUBROUTINE Kalman_params_absigmaQvec
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE Kalman_fhat ( r_inf, avec_1, avec_2, Bmat_1, Bmat_2, sigmaQvec_1, sigmaQvec_2, &
        xitt, hvec, Hmat )
    !
    ! Computes fitted yields and Jacobian of fitted yields w.r.t. states
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: r_inf
    REAL(8), DIMENSION(num_K), INTENT(IN) :: avec_1, avec_2, sigmaQvec_1, sigmaQvec_2
    REAL(8), INTENT(IN) :: Bmat_1(num_XC(1),num_K), Bmat_2(num_XC(2),num_K)
    REAL(8), INTENT(IN) :: xitt(num_XI)                 ! Factors
    REAL(8), INTENT(OUT) :: hvec(num_C*num_K)           ! Fitted values
    REAL(8), INTENT(OUT) :: Hmat(num_XI,num_C*num_K)     ! Jacobian of fitted values w.r.t. states
    !
    ! Declaring local variables
    !
    INTEGER :: ic
    REAL(8) :: hvec_c(num_K), Hmat_1(num_XC(1),num_K), Hmat_2(num_XC(2),num_K)
    ! 
    ! Beginning execution
    !
    hvec = 0.d0
    Hmat = 0.d0
    !
    ! GATSM
    !
    IF (switch_gatsm .EQ. 1) THEN
        !
        CALL Kalman_fhat_gatsm(num_XC(1),avec_1,Bmat_1,xitt(sel_X1),hvec_c,Hmat_1)
        hvec(:num_K) = hvec_c
        Hmat(sel_X1,:num_K) = Hmat_1
        !
        CALL Kalman_fhat_gatsm(num_XC(2),avec_2,Bmat_2,xitt(sel_X2),hvec_c,Hmat_2)
        hvec(num_K+1:2*num_K) = hvec_c
        Hmat(sel_X2,num_K+1:2*num_K) = Hmat_2
        !
    END IF
    !
    ! SR-GATSM
    !
    IF (switch_sr .EQ. 1) THEN
        !
        CALL Kalman_fhat_sr(num_XC(1),r_inf,avec_1,Bmat_1,sigmaQvec_1,xitt(sel_X1),hvec_c,Hmat_1)
        hvec(:num_K) = hvec_c
        Hmat(sel_X1,:num_K) = Hmat_1
        !
        CALL Kalman_fhat_sr(num_XC(2),r_inf,avec_2,Bmat_2,sigmaQvec_2,xitt(sel_X2),hvec_c,Hmat_2)
        hvec(num_K+1:2*num_K) = hvec_c
        Hmat(sel_X2,num_K+1:2*num_K) = Hmat_2
        !
    END IF
    !
    ! End execution and returning control
    !
    END SUBROUTINE Kalman_fhat
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE Kalman_fhat_gatsm ( n, avec, Bmat, xtt, hvec, Hmat )
    !
    ! Computes fitted yields and Jacobian of fitted yields w.r.t. states in the GATSM model
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: avec(num_K)
    REAL(8), INTENT(IN) :: Bmat(n,num_K)
    REAL(8), INTENT(IN) :: xtt(n)
    REAL(8), INTENT(OUT) :: hvec(num_K)                 ! Fitted yields
    REAL(8), INTENT(OUT) :: Hmat(n,num_K)               ! Jacobian of fitted yields w.r.t. states
    !
    ! Declaring local variables
    !
    INTEGER :: ir, ik
    ! 
    ! Beginning execution
    !
    ! Beginning loop over maturities
    !
    DO ik = 1, num_K
        !
        hvec(ik) = avec(ik)
        DO ir = 1, n
            !
            hvec(ik) = hvec(ik)+Bmat(ir,ik)*xtt(ir)
            !
        END DO
        Hmat(:,ik) = Bmat(:,ik)
        !
    END DO
    !
    ! End execution and returning control
    !
    END SUBROUTINE Kalman_fhat_gatsm
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE Kalman_fhat_sr ( n, r_inf, avec, Bmat, sigmaQvec, xtt, hvec, Hmat )
    !
    ! Computes fitted yields and Jacobian of fitted yields w.r.t. states in the SR-GATSM model
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: r_inf                        ! Short rate lower bound
    REAL(8), INTENT(IN) :: avec(num_K)
    REAL(8), INTENT(IN) :: Bmat(n,num_K)
    REAL(8), INTENT(IN) :: sigmaQvec(num_K)
    REAL(8), INTENT(IN) :: xtt(n)
    REAL(8), INTENT(OUT) :: hvec(num_K)                 ! Fitted yields
    REAL(8), INTENT(OUT) :: Hmat(n,num_K)               ! Jacobian of fitted yields w.r.t. states
    !
    ! Declaring local variables
    !
    INTEGER :: ir, ik
    REAL(8) :: z, cdfN01_z, pdfN01_z
    ! 
    ! Beginning execution
    !
    ! Beginning loop over maturities
    !
    DO ik = 1, num_K
        !
        ! z is the argument of g and g'
        !
        z = avec(ik)
        DO ir = 1, n
            z = z+Bmat(ir,ik)*xtt(ir)
        END DO
        z = (z-r_inf)/sigmaQvec(ik)
        !
        ! Phi(z) and phi(z)
        !
        cdfN01_z = cdf_N01(z,.FALSE.)
        pdfN01_z = pdf_N01(z)
        !
        ! Adding to hvec and Hmat
        !
        hvec(ik) = r_inf+sigmaQvec(ik)*(z*cdfN01_z+pdfN01_z)
        Hmat(:,ik) = cdfN01_z*Bmat(:,ik)
        !
    END DO
    !
    ! End execution and returning control
    !
    END SUBROUTINE Kalman_fhat_sr
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE Kalman_ehat ( r_inf, delta0, delta1_1, delta1_2, lambda_1, lambda_2, &
        LLambda_1, LLambda_2, mu_1, mu_2, Phi_1, Phi_2, &
        xitt, kvec, Kmat, lambdattlag_1, eps_1, lambdattlag_2, eps_2 )
    !
    ! Computes fitted depreciation and Jacobian of fitted depreciation w.r.t. states
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: r_inf
    REAL(8), INTENT(IN) :: delta0(num_C)
    REAL(8), DIMENSION(num_XC(1)), INTENT(IN) :: delta1_1, mu_1, lambda_1
    REAL(8), DIMENSION(num_XC(2)), INTENT(IN) :: delta1_2, mu_2, lambda_2
    REAL(8), DIMENSION(num_XC(1),num_XC(1)), INTENT(IN) :: Phi_1, LLambda_1
    REAL(8), DIMENSION(num_XC(2),num_XC(2)), INTENT(IN) :: Phi_2, LLambda_2
    REAL(8), INTENT(IN) :: xitt(num_XI)                         ! Factors
    REAL(8), INTENT(OUT) :: kvec                                ! Fitted depreciation
    REAL(8), INTENT(OUT) :: Kmat(num_XI)                        ! Jacobian of fitted depreciation w.r.t. states
    REAL(8), DIMENSION(num_XC(1)), INTENT(OUT) ::  lambdattlag_1, eps_1
    REAL(8), DIMENSION(num_XC(2)), INTENT(OUT) ::  lambdattlag_2, eps_2
    !
    ! Declaring local variables
    !
    REAL(8), DIMENSION(num_X) :: xtt, xttlag
    REAL(8), DIMENSION(num_XC(1)) :: xtt_1, xttlag_1, dr_1
    REAL(8), DIMENSION(num_XC(2)) :: xtt_2, xttlag_2, dr_2
    REAL(8) :: r_1, r_2                 ! REMEMBER: 1 IS US, 2 IS UK !!!
    INTEGER :: ir, ic
    REAL(8) :: Kmat_tmp(num_X,num_C)
    ! 
    ! Beginning execution
    !
    kvec = 0.d0
    Kmat = 0.d0
    !
    ! Factors
    !
    ! tt
    xtt = xitt(:num_X)
    xtt_1 = xtt(sel_X1)
    xtt_2 = xtt(sel_X2)
    !
    ! tt-1
    xttlag = xitt(num_X+1:)
    xttlag_1 = xttlag(sel_X1)
    xttlag_2 = xttlag(sel_X2)
    !
    ! Short rates in a GATSM model
    !
    r_1 = delta0(1)+SUM(delta1_1*xttlag_1)
    r_2 = delta0(2)+SUM(delta1_2*xttlag_2)
    dr_1 = delta1_1
    dr_2 = delta1_2
    !
    ! Short rates in a SR-GATSM model
    !
    IF (switch_sr .EQ. 1) THEN
        !
        IF (r_1 .LT. r_inf) THEN
            !
            r_1 = r_inf
            dr_1 = 0.d0
            !
        END IF
        IF (r_2 .LT. r_inf) THEN
            !
            r_2 = r_inf
            dr_2 = 0.d0
            !
        END IF
        !
    END IF
    !
    ! Risk premia and factors' innovations
    ! BEWARE: THE INNOVATIONS ARE CORRECT ONLY IF GAMMA = IDENTITY MATRIX !!!
    !
    DO ir = 1, num_XC(1)
        !
        lambdattlag_1(ir) = lambda_1(ir)
        eps_1(ir) = xtt_1(ir)-mu_1(ir)
        DO ic = 1, num_XC(1)
            !
            lambdattlag_1(ir) = lambdattlag_1(ir)+LLambda_1(ir,ic)*xttlag_1(ic)
            eps_1(ir) = eps_1(ir)-Phi_1(ir,ic)*xttlag_1(ic)
            !
        END DO
        !
    END DO
    DO ir = 1, num_XC(2)
        !
        lambdattlag_2(ir) = lambda_2(ir)
        eps_2(ir) = xtt_2(ir)
        DO ic = 1, num_XC(2)
            !
            lambdattlag_2(ir) = lambdattlag_2(ir)+LLambda_2(ir,ic)*xttlag_2(ic)
            eps_2(ir) = eps_2(ir)-Phi_2(ir,ic)*xttlag_2(ic)
            !
        END DO
        !
    END DO
    !
    ! kvec
    !
    kvec = r_1-r_2+ &
        0.5d0*(SUM(lambdattlag_1**2)-SUM(lambdattlag_2**2))+ &
        (SUM(lambdattlag_1*eps_1)-SUM(lambdattlag_2*eps_2))
    !
    ! Kmat
    !
    ! xtt
    Kmat_tmp = 0.d0
    DO ir = 1, num_XC(1)
        !
        Kmat_tmp(sel_X1(ir),1) = lambdattlag_1(ir)
        !
    END DO
    DO ir = 1, num_XC(2)
        !
        Kmat_tmp(sel_X2(ir),2) = lambdattlag_2(ir)
        !
    END DO
    Kmat(:num_X) = Kmat_tmp(:,1)-Kmat_tmp(:,2)
    !
    ! xttlag
    Kmat_tmp = 0.d0
    DO ir = 1, num_XC(1)
        !
        Kmat_tmp(sel_X1(ir),1) = dr_1(ir)
        DO ic = 1, num_XC(1)
            !
            Kmat_tmp(sel_X1(ir),1) = Kmat_tmp(sel_X1(ir),1)+ &
                LLambda_1(ic,ir)*(lambdattlag_1(ic)+eps_1(ic))- &
                Phi_1(ic,ir)*lambdattlag_1(ic)
            !
        END DO
        !
    END DO
    DO ir = 1, num_XC(2)
        !
        Kmat_tmp(sel_X2(ir),2) = dr_2(ir)
        DO ic = 1, num_XC(2)
            !
            Kmat_tmp(sel_X2(ir),2) = Kmat_tmp(sel_X2(ir),2)+ &
                LLambda_2(ic,ir)*(lambdattlag_2(ic)+eps_2(ic))- &
                Phi_2(ic,ir)*lambdattlag_2(ic)
            !
        END DO
        !
    END DO
    Kmat(num_X+1:) = Kmat_tmp(:,1)-Kmat_tmp(:,2)
    !
    ! End execution and returning control
    !
    END SUBROUTINE Kalman_ehat
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE Kalman_vec ( theta, objf_vec, xif, xiu, Pf, Pu )
    !
    ! Modified Kalman filter - loglikelihood contributions
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)
    REAL(8), INTENT(OUT) :: objf_vec(num_T)
    REAL(8), INTENT(OUT) :: xif(num_XI,num_T)
    REAL(8), INTENT(OUT) :: xiu(num_XI,num_T)
    REAL(8), INTENT(OUT) :: Pf(num_XI,num_XI,num_T)
    REAL(8), INTENT(OUT) :: Pu(num_XI,num_XI,num_T)
    !
    ! Declaring local variables
    !
    INTEGER :: tt, ktt, info, ir, ic, h, k, il, iu, jl, ju, i, j
    INTEGER, DIMENSION(num_C*num_K) :: Tautt
    INTEGER, DIMENSION(num_Y) :: ipiv2K
    INTEGER :: ipivN(num_XI), ipivNq(num_XI**2)
    REAL(8) :: delta0(num_C), r_inf, tmp, log_det_Ltt, logLtt, logL(num_T), omega_dex
    REAL(8), DIMENSION(num_X) :: mu
    REAL(8), DIMENSION(num_X,num_X) :: Gamm, Sigma, Phi
    REAL(8), DIMENSION(num_XI,num_C*num_K) :: Hmat
    REAL(8), DIMENSION(num_XC(1)) :: mu_1, muQ_1, delta1_1, lambda_1, lambdattlag_1, eps_1
    REAL(8), DIMENSION(num_XC(2)) :: mu_2, muQ_2, delta1_2, lambda_2, lambdattlag_2, eps_2
    REAL(8), DIMENSION(num_XC(1),num_XC(1)) :: PhiQ_1, Sigma_1, Phi_1, LLambda_1
    REAL(8), DIMENSION(num_XC(2),num_XC(2)) :: PhiQ_2, Sigma_2, Phi_2, LLambda_2
    REAL(8), DIMENSION(num_XI) :: muXI, XI0
    REAL(8), DIMENSION(num_XI,num_XI) :: PhiXI, SigmaXI, P0, loc_A
    REAL(8), DIMENSION(num_C*num_K) :: hvec
    REAL(8), DIMENSION(num_C*num_K,num_C*num_K) :: Omega
    REAL(8), DIMENSION(num_Y) :: ustar, ystar, mvecstar
    REAL(8), DIMENSION(num_XI,num_Y) :: Mmatstar, Qtt, Gtt
    REAL(8), DIMENSION(num_Y,num_Y) :: Omegastar, Ltt, triLtt
    REAL(8) :: Qttprime(num_Y,num_XI), z_i(num_Y,1)
    REAL(8) :: loc_A2(num_XI**2,num_XI**2), loc_B(num_XI,1), loc_B2(num_XI**2,1)
    INTEGER :: lwork
    REAL(8) :: wr(num_X), wi(num_X), Phi_tmp(num_X,num_X), vl(num_X,num_X), vr(num_X,num_X), work(10*num_X)
    REAL(8), DIMENSION(num_K) :: avec_1, avec_2, sigmaQvec_1, sigmaQvec_2
    REAL(8), DIMENSION(num_XC(1),num_K) :: Bmat_1
    REAL(8), DIMENSION(num_XC(2),num_K) :: Bmat_2
    REAL(8) :: kvec                                ! Fitted depreciation
    REAL(8) :: Kmat(num_XI)                        ! Jacobian of fitted depreciation w.r.t. states
    !
    ! Beginning execution
    !
    ! Extracting model parameters
    !
    CALL theta_to_param(theta,delta0,delta1_1,delta1_2,PhiQ_1,PhiQ_2,muQ_1,muQ_2, &
        Gamm,Sigma,Sigma_1,Sigma_2,mu,mu_1,mu_2,lambda_1,lambda_2, &
        Phi,Phi_1,Phi_2,LLambda_1,LLambda_2,Omega,omega_dex,r_inf)
    !
    ! Computing Kalman filter parameters
    !
    CALL Kalman_params(delta0,delta1_1,delta1_2,PhiQ_1,PhiQ_2,muQ_1,muQ_2, &
        Sigma_1,Sigma_2,mu,Phi,Sigma, &
        avec_1,avec_2,Bmat_1,Bmat_2,sigmaQvec_1,sigmaQvec_2,muXI,PhiXI,SigmaXI)
    !
    ! Initial conditions
    !
    IF (switch_init_Kalman .EQ. 1) THEN
        !
        ! X0 = 0, P0 = 100*Id, as in Wu-Xia (2014)
        !
        XI0 = 0.d0
        P0 = 1.d2*eyeXI
        !
    ELSE IF (switch_init_Kalman .EQ. 2) THEN
        !
        ! X0 = long run mean
        ! P0 = P0 = 100*Id
        !
        loc_A = eyeXI-PhiXI
        loc_B(:,1) = muXI
        CALL DGESV(num_XI,1,loc_A,num_XI,ipivN,loc_B,num_XI,info)
        XI0 = loc_B(:,1)
        P0 = 1.d2*eyeXI
        !
    ELSE IF (switch_init_Kalman .EQ. 3) THEN
        !
        ! X0 = long run mean
        ! P0 = long run variance
        ! (this does not work if Phi has eigenvalues larger than 1 in modulus)
        !
        loc_A = eyeXI-PhiXI
        loc_B(:,1) = muXI
        CALL DGESV(num_XI,1,loc_A,num_XI,ipivN,loc_B,num_XI,info)
        XI0 = loc_B(:,1)
        !
        loc_A2 = eyeXI2
        loc_B2 = 0.d0
        h = 0
        DO i = 1, num_XI
            DO j = 1, num_XI
                il = (i-1)*num_XI+1
                iu = i*num_XI
                jl = (j-1)*num_XI+1
                ju = j*num_XI
                loc_A2(il:iu,jl:ju) = eyeXI2(il:iu,jl:ju)-PhiXI(i,j)*PhiXI
                h = h+1
                IF (i .EQ. j) loc_B2(h,1) = SigmaXI(i,i)
            END DO
        END DO
        CALL DGESV(num_XI**2,1,loc_A2,num_XI**2,ipivNq,loc_B2,num_XI**2,info)
        h = 0
        DO j = 1, num_XI
            DO i = 1, num_XI
                h = h+1
                P0(i,j) = loc_B2(h,1)
            END DO
        END DO
        !
    END IF
    !
    ! Starting the recursion
    ! Unconditional mean and variance of state variable vector
    ! NB: 
    !       xif = \hat \xi_{t|t-1},    Pf = P_{t|t-1}
    !       xiu = \hat \xi_{t|t},      Pu = P_{t|t}
    !
    xif(:,1) = XI0
    Pf(:,:,1) = P0
    !
    ! Kalman filter
    !
    DO tt = 1, num_T
        !
        ! Computing the loglikelihood contribution
        !
        ! Computing vector of predicted yields (hvec) and slope matrix (Hmat)
        !
        CALL Kalman_fhat(r_inf,avec_1,avec_2,Bmat_1,Bmat_2,sigmaQvec_1,sigmaQvec_2, &
            xif(:,tt),hvec,Hmat)
        CALL Kalman_ehat(r_inf,delta0,delta1_1,delta1_2,lambda_1,lambda_2, &
            LLambda_1,LLambda_2,mu_1,mu_2,Phi_1,Phi_2, &
            xif(:,tt),kvec,Kmat,lambdattlag_1,eps_1,lambdattlag_2,eps_2)    
        !
        ! Observations for time t
        !
        ktt = K_vec(tt)
        Tautt = Tau_mat(:,tt)
        !
        ! Selected vectors and matrices
        !
        ystar = 0.d0
        ystar(:ktt) = F(Tautt(:ktt),tt)
        ystar(ktt+1) = dexrate(tt)
        !
        mvecstar = 0.d0
        mvecstar(:ktt) = hvec(Tautt(:ktt))
        mvecstar(ktt+1) = kvec
        !
        Mmatstar = 0.d0
        Mmatstar(:num_X,:ktt) = Hmat(:,Tautt(:ktt))
        Mmatstar(:,ktt+1) = Kmat
        !
        Omegastar = 0.d0
        Omegastar(:ktt,:ktt) = Omega(Tautt(:ktt),Tautt(:ktt))
        Omegastar(ktt+1,ktt+1) = omega_dex
        !
        ! u = [F - hvec, dexrate - kvec]
        !
        ustar = 0.d0
        ustar(:ktt+1) = ystar(:ktt)-mvecstar(:ktt)
        !
        ! Ltt = Mmatstar'*Pf*Mmatstar + Omegastar
        !
        Ltt(:ktt+1,:ktt+1) = Omegastar(:ktt+1,:ktt+1)
        DO ir = 1, ktt+1
            DO ic = 1, ktt+1
                DO h = 1, num_XI
                    DO k = 1, num_XI
                        Ltt(ir,ic) = Ltt(ir,ic)+Mmatstar(h,ir)*Mmatstar(k,ic)*Pf(h,k,tt)
                    END DO
                END DO
            END DO
        END DO
        !
        ! Choleski decomposition of Ltt(:ktt+1,:ktt+1)
        !
        triLtt = Ltt
        DO ic = 1, ktt+1
            tmp = 0.d0
            DO h = 1, ic-1
                tmp = tmp+triLtt(ic,h)**2
            END DO
            triLtt(ic,ic) = SQRT(triLtt(ic,ic)-tmp)
            DO ir = ic+1, ktt+1
                tmp = 0.d0
                DO h = 1, ic-1
                    tmp = tmp+triLtt(ir,h)*triLtt(ic,h)
                END DO
                triLtt(ir,ic) = (triLtt(ir,ic)-tmp)/triLtt(ic,ic)
            END DO
        END DO
        DO ic = 2, ktt+1
            triLtt(1:ic-1,ic) = 0.d0
        END DO
        !
        ! Loglikelihood contribution
        !
        log_det_Ltt = 0.d0
        DO h = 1, ktt+1
            log_det_Ltt = log_det_Ltt+2.d0*LOG(triLtt(h,h))
        END DO
        !
        z_i(:,1) = ustar
        CALL DGESV(ktt+1,1,triLtt(:ktt+1,:ktt+1),ktt+1,ipiv2K(:ktt+1),z_i(:ktt+1,1),ktt+1,info)
        logLtt = -0.5d0*((ktt+1)*LOG(2.d0*pi)+log_det_Ltt+SUM(z_i**2))
        logL(tt) = logLtt
        !
        ! Update step
        !
        ! Qtt = Pf*Mmatstar
        !
        Qtt = 0.d0
        DO ir = 1, num_XI
            DO ic = 1, ktt+1
                DO h = 1, num_XI
                    Qtt(ir,ic) = Qtt(ir,ic)+Pf(ir,h,tt)*Mmatstar(h,ic)
                END DO
            END DO
        END DO
        !
        ! Qttprime = Qtt'
        !
        Qttprime = 0.d0
        DO ir = 1, num_XI
            DO ic = 1, ktt+1
                Qttprime(ic,ir) = Qtt(ir,ic)
            END DO
        END DO
        !
        ! Solve Ltt(:ktt+1,:ktt+1) Gtt(:,:ktt+1)' = Qttprime(:ktt+1,:) w.r.t Gtt(:,:ktt+1)'
        !
        CALL DGESV(ktt+1,num_XI,Ltt(:ktt+1,:ktt+1),ktt+1,ipiv2K(:ktt+1),Qttprime(:ktt+1,:),ktt+1,info)
        !
        Gtt = 0.d0
        DO ir = 1, ktt+1
            DO ic = 1, num_XI
                Gtt(ic,ir) = Qttprime(ir,ic)
            END DO
        END DO
        !
        ! Update step
        !
        xiu(:,tt) = xif(:,tt)
        DO ir = 1, num_XI
            DO h = 1, ktt+1
                xiu(ir,tt) = xiu(ir,tt)+Gtt(ir,h)*ustar(h)
            END DO
        END DO
        !
        Pu(:,:,tt) = Pf(:,:,tt)
        DO ir = 1, num_XI
            DO ic = 1, ir
                DO h = 1, ktt+1
                    DO k = 1, num_XI
                        Pu(ir,ic,tt) = Pu(ir,ic,tt)-Gtt(ir,h)*Mmatstar(k,h)*Pf(k,ic,tt)
                    END DO
                END DO
            END DO
        END DO
        DO ir = 1, num_XI-1
            DO ic = ir+1, num_XI
                Pu(ir,ic,tt) = Pu(ic,ir,tt)
            END DO
        END DO
        !
        ! Forecast step
        !
        IF (tt .LT. num_T) THEN
            !
            xif(:,tt+1) = muXI
            DO ir = 1, num_XI
                DO h = 1, num_XI
                    xif(ir,tt+1) = xif(ir,tt+1)+PhiXI(ir,h)*xiu(h,tt)
                END DO
            END DO
            !
            Pf(:,:,tt+1) = SigmaXI
            DO ir = 1, num_XI
                DO ic = 1, ir
                    DO h = 1, num_XI
                        DO k = 1, num_XI
                            Pf(ir,ic,tt+1) = Pf(ir,ic,tt+1)+PhiXI(ir,h)*Pu(h,k,tt)*PhiXI(ic,k)
                        END DO
                    END DO
                END DO
            END DO
            DO ir = 1, num_XI-1
                DO ic = ir+1, num_XI
                    Pf(ir,ic,tt+1) = Pf(ic,ir,tt+1)
                END DO
            END DO
        END IF
        !
    END DO
    !
    ! Minus average loglikelihood
    !
    objf_vec = logL
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE Kalman_vec
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE ll ( theta, objf, xif, xiu, Pf, Pu )
    !
    ! Modified Kalman filter - loglikelihood contributions
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)
    REAL(8), INTENT(OUT) :: objf
    REAL(8), INTENT(OUT) :: xif(num_XI,num_T)
    REAL(8), INTENT(OUT) :: xiu(num_XI,num_T)
    REAL(8), INTENT(OUT) :: Pf(num_XI,num_XI,num_T)
    REAL(8), INTENT(OUT) :: Pu(num_XI,num_XI,num_T)
    !
    ! Declaring local variables
    !
    REAL(8) :: objf_vec(num_T)
    !
    ! Beginning execution
    !
    CALL Kalman_vec(theta,objf_vec,xif,xiu,Pf,Pu)
    !
    ! Minus average loglikelihood
    !
    objf = SUM(objf_vec)/num_T
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE ll
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION cdf_N01 ( x, upper )
    !
    ! cdf_N01 computes the cumulative density of the standard normal distribution.
    ! It can be differentiated automatically by TAPENADE
    !
    ! cdf_N01 is based on the ALNORM function
    ! Original FORTRAN77 version by David Hill.
    ! FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, is one endpoint of the semi-infinite interval
    !    over which the integration takes place.
    !
    !    Input, logical UPPER, determines whether the upper or lower
    !    interval is to be integrated:
    !    .TRUE.  => integrate from X to + Infinity;
    !    .FALSE. => integrate from - Infinity to X.
    !
    !    Output, real ( kind = 8 ) ALNORM, the integral of the standard normal
    !    distribution over the desired interval.
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: x
    LOGICAL, INTENT(IN) :: upper
    !
    ! Declaring local variables
    !
    LOGICAL :: up
    REAL(8) :: y, z
    !
    ! Declaring parameters
    !
    REAL(8), PARAMETER :: a1 = 5.75885480458D+00
    REAL(8), PARAMETER :: a2 = 2.62433121679D+00
    REAL(8), PARAMETER :: a3 = 5.92885724438D+00
    REAL(8), PARAMETER :: b1 = -29.8213557807D+00
    REAL(8), PARAMETER :: b2 = 48.6959930692D+00
    REAL(8), PARAMETER :: c1 = -0.000000038052D+00
    REAL(8), PARAMETER :: c2 = 0.000398064794D+00
    REAL(8), PARAMETER :: c3 = -0.151679116635D+00
    REAL(8), PARAMETER :: c4 = 4.8385912808D+00
    REAL(8), PARAMETER :: c5 = 0.742380924027D+00
    REAL(8), PARAMETER :: c6 = 3.99019417011D+00
    REAL(8), PARAMETER :: con = 1.28D+00
    REAL(8), PARAMETER :: d1 = 1.00000615302D+00
    REAL(8), PARAMETER :: d2 = 1.98615381364D+00
    REAL(8), PARAMETER :: d3 = 5.29330324926D+00
    REAL(8), PARAMETER :: d4 = -15.1508972451D+00
    REAL(8), PARAMETER :: d5 = 30.789933034D+00
    REAL(8), PARAMETER :: ltone = 7.0D+00
    REAL(8), PARAMETER :: p = 0.398942280444D+00
    REAL(8), PARAMETER :: q = 0.39990348504D+00
    REAL(8), PARAMETER :: r = 0.398942280385D+00
    REAL(8), PARAMETER :: utzero = 18.66D+00
    !
    ! Declaring function's type
    !
    REAL(8) :: cdf_N01
    !
    ! Beginning execution
    !
    up = upper
    z = x
    !
    IF (z < 0.0D+00) THEN
        up = .not. up
        z = - z
    END IF
    !
    IF ((ltone < z) .AND. ((.NOT. up) .OR. (utzero < z))) THEN
        !
        IF (up) THEN
            cdf_N01 = 0.0D+00
        ELSE
            cdf_N01 = 1.0D+00
        END IF
        RETURN
        !
    END IF
    !
    y = 0.5D+00*z*z
    !
    IF (z <= con) THEN
        cdf_N01 = 0.5D+00-z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))))
    ELSE
        cdf_N01 = r*EXP(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))))
    END IF
    IF (.NOT. up) cdf_N01 = 1.0D+00-cdf_N01
    RETURN
    !
    ! Ending execution
    !
    END FUNCTION cdf_N01
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION pdf_N01 ( x )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: x
    !
    ! Declaring function's type
    !
    REAL(8) :: pdf_N01
    !
    ! Beginning execution
    !
    pdf_N01 = 1.d0/SQRT(2.d0*pi)*EXP(-x**2/2.d0)
    !
    ! Ending execution and returning control
    !
    END FUNCTION pdf_N01
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION posdef_matrix ( n, x )
    !
    ! Given a vector with n^2 elements, generates a positive definite matrix
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: n
    REAL(8), INTENT(IN) :: x(n**2)
    !
    ! Declaring function's type
    !
    REAL(8) :: posdef_matrix(n,n)
    !
    ! Declaring local variables
    !
    INTEGER :: h, ir, ic
    REAL(8), DIMENSION(n,n) :: S1, S, A, K
    !
    ! Beginning execution
    !
    ! S1 = lower triangular (N x N) => N(N+1)/2 parameters
    !
    h = 0
    S1 = 0.d0
    DO ir = 1, n
        !
        DO ic = 1, ir
            !
            h = h+1
            S1(ir,ic) = x(h)
            !
        END DO
        !
    END DO
    !
    ! S = S1 S1' = spd (N x N)
    !
    S = 0.d0
    DO ir = 1, n
        !
        DO ic = 1, n
            !
            DO h = 1, n
                !
                S(ir,ic) = S(ir,ic)+S1(ir,h)*S1(ic,h)
                !
            END DO
            !
        END DO
        !
    END DO
    !
    ! A = antisymmetric matrix (A = -A') (N x N) => N(N-1)/2 parameters
    !
    h = 0
    A = 0.d0
    DO ir = 2, n
        !
        DO ic = 1, ir-1
            !
            h = h+1
            A(ir,ic) = x(n*(n+1)/2+h)
            A(ic,ir) = -A(ir,ic)
            !
        END DO
        !
    END DO
    !
    ! K = S + A = general real positive definite matrix (N x N)
    !
    K = S+A
    !
    ! Phi = exp(-K)
    !
    CALL r8mat_expm1(n,-K,posdef_matrix)
    !
    ! Ending execution and returning control
    !
    END FUNCTION posdef_matrix
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE gatsm