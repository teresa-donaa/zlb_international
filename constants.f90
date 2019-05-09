MODULE constants
!
IMPLICIT NONE
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Selecting optimization method and other switches
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
INTEGER, PARAMETER :: to0 = 0                   ! Test mode
INTEGER, PARAMETER :: to1 = 0                   ! First estimation round
INTEGER, PARAMETER :: to2 = 0                   ! Second estimation round
INTEGER, PARAMETER :: to3 = 0                   ! Final: compute skewnesses amd kurtoses of simulated spreads
INTEGER, PARAMETER :: to4 = 1                   ! Final: factor trajectories and diagnostics
INTEGER, PARAMETER :: to5 = 0                   ! Final: effects of a change in the lower bound
INTEGER, PARAMETER :: to6 = 0                   ! Final: asymptotic variance
!
INTEGER, PARAMETER :: num_stime = to0*1+to1*100+to2*1+to4*0+to6*0
                                                ! Total number of completed estimation trials
!
! Selection of optimization algorithm
!
INTEGER, PARAMETER :: switch_lbfgs = 1          ! = 1: L-BFGS optimization ON; = 0 L-BFGS optimization OFF
INTEGER, PARAMETER :: switch_politope = 1       ! = 1: POLITOPE optimization ON; = 0 POLITOPE optimization OFF
!
! States' initialization in Kalman filtering
!
INTEGER, PARAMETER :: switch_init_Kalman = 2    ! = 1: x0 = 0, P0 = 100*eye(N) (as in Wu-Xia 2014)
                                                ! = 2: x0 = long run mean, P0 = 100*eye(N)
                                                ! = 3: x0 = long run mean, P0 = long run variance
                                                ! (this does not work with a nonstationary Phi)
!
! Initial values
!
INTEGER, PARAMETER :: switch_inival1 = 0        ! = 1: IniVal randomly chosen
INTEGER, PARAMETER :: switch_inival2 = 1        ! = 1: IniVal chosen as in Hamilton-Wu
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring sample characteristics 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Period: 2001m1 - 2016m10 (same for every country)
!
INTEGER, PARAMETER :: num_T_tot = 320           ! Total number of observations in data file
INTEGER, PARAMETER :: ind_obs1 = 1              ! Index of first observation to be used in estimation
INTEGER, PARAMETER :: ind_obsT = 320            ! Index of last observation to be used in estimation
INTEGER, PARAMETER :: num_T = &                 ! Total number of observations to be used in estimation
    ind_obsT-ind_obs1+1
!
! Maturities (same for every country)
!
INTEGER, PARAMETER :: num_K = 6                 ! Number of maturities
INTEGER, PARAMETER :: Tau_vec(num_K) = &        ! Vector of observed maturities
!     3m, 6m, 1y, 3y, 5y, 10y
    (/ 3,  6, 12, 36, 60, 120 /)
!
! Countries
!
INTEGER, PARAMETER :: num_K_US = 6              ! Number of available US maturities
INTEGER, PARAMETER :: data_US_cols(num_K) = &   ! Selected US maturities
!     3m 6m 1y 3y 5y 10y
!      1  2  3  4  5   6
    (/ 1, 2, 3, 4, 5,  6 /)
!
INTEGER, PARAMETER :: num_K_UK = 6              ! Number of available UK maturities
INTEGER, PARAMETER :: data_UK_cols(num_K) = &   ! Selected UK maturities
!     3m 6m 1y 3y 5y 10y
!      1  2  3  4  5   6
    (/ 1, 2, 3, 4, 5,  6 /)
!
INTEGER, PARAMETER :: num_C = 2                 ! Total number of countries in the joint model
CHARACTER(len=2), PARAMETER :: lab_c(num_C) = & ! Country labels
    (/ "US", "UK" /)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Factors specification and selection 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Factors
INTEGER, PARAMETER :: num_X = 4                 ! Total number of factors in the joint model
INTEGER, PARAMETER :: num_XG = 2                ! Number of global factors
INTEGER, PARAMETER :: num_XL = 2                ! Number of local factors
INTEGER, PARAMETER :: num_XC(num_C) = &         ! Number of country factors
    (/ 3, 3 /)
INTEGER, PARAMETER :: sel_X(num_X,num_C) = &    ! (N x C) selection matrix of factors entering in each country's TS
    RESHAPE( (/ &                               ! List indexes of factors used for each country, pad 0s to num_X
!       1   2   3   4   5   6   7   8           ! NB: This is written transposed here! 
        1,  2,  3,  0,  &
        1,  2,  4,  0   &
    /), (/ num_X, num_C /) )
INTEGER, PARAMETER :: sel_X1(num_XC(1)) = sel_X(:num_XC(1),1)
INTEGER, PARAMETER :: sel_X2(num_XC(2)) = sel_X(:num_XC(2),2)
INTEGER, PARAMETER :: num_XI = 2*num_X          ! Number of factors in the companion form model
INTEGER, PARAMETER :: num_Y = num_C*num_K+1     ! (Maximum) Number of observed variables
!
! Model specification
INTEGER, PARAMETER :: switch_gatsm       = 0    ! = 1: GATSM model
INTEGER, PARAMETER :: switch_sr_rinf_fix = 1    ! = 1: SR-GATSM, r_inf fixed (same for all c)
INTEGER, PARAMETER :: switch_sr_rinf_est = 0    ! = 1: SR-GATSM, r_inf estimated (same for all c)
INTEGER, PARAMETER :: switch_sr = &             ! = 1: Shadow Rate model
    switch_sr_rinf_fix+switch_sr_rinf_est
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Parameterization
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Specification of delta0
INTEGER, PARAMETER :: switch_delta01 = 1        ! = 1: delta0 unconstrained
INTEGER, PARAMETER :: switch_delta02 = 0        ! = 1: delta0 = 0
!
! Specification of PhiQ
INTEGER, PARAMETER :: switch_PhiQ1 = 0          ! = 1: PhiQ lower triangular
INTEGER, PARAMETER :: switch_PhiQ2 = 1          ! = 1: PhiQ full
!
! Specification of muQ
INTEGER, PARAMETER :: switch_muQ1 = 1           ! = 1: muQ unconstrained
INTEGER, PARAMETER :: switch_muQ2 = 0           ! = 1: muQ = 0
!
! Specification of mu
INTEGER, PARAMETER :: switch_mu1 = 0            ! = 1: mu unconstrained
INTEGER, PARAMETER :: switch_mu2 = 1            ! = 1: mu = 0
!
! Specification of Phi
INTEGER, PARAMETER :: switch_Phi1 = 0           ! = 1: Phi full
INTEGER, PARAMETER :: switch_Phi2 = 1           ! = 1: Phi l.tr., with zero elements btw local factors
!
! Specification of r_inf:
INTEGER, PARAMETER :: switch_r_inf0 = &         ! = 1: GATSM
    1*switch_gatsm+ &
    0*(switch_sr_rinf_fix+switch_sr_rinf_est)
INTEGER, PARAMETER :: switch_r_inf1 = &         ! = 1: r_inf fixed at r_inf_fixed (same for all c)
    1*switch_sr_rinf_fix+ &
    0*(switch_gatsm+switch_sr_rinf_est)
REAL(8), PARAMETER :: r_inf_fixed = 0.d0        ! Fixed value for the lower bound
INTEGER, PARAMETER :: switch_r_inf2 = &         ! = 1: r_inf to be estimated (same for all c)
    1*switch_sr_rinf_est+ &
    0*(switch_gatsm+switch_sr_rinf_fix)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Defining the number of parameters
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
INTEGER, PARAMETER :: num_delta0 = &
    switch_delta01*num_C+ &
    switch_delta02*0
INTEGER, PARAMETER :: num_delta1 = num_XC(1)+num_XC(2)
INTEGER, PARAMETER :: num_PhiQ = &              ! Only works for models with (Xg,1,1) factors
    switch_PhiQ1*(num_XC(1)*(num_XC(1)+1)/2+num_XC(2)*(num_XC(2)+1)/2)+ &
    switch_PhiQ2*(num_XC(1)**2+num_XC(2)**2)
INTEGER, PARAMETER :: num_muQ = &
    switch_muQ1*(num_XC(1)+num_XC(2))+ &
    switch_muQ2*0
INTEGER, PARAMETER :: num_Gamm = 0
INTEGER, PARAMETER :: num_mu = &
    switch_mu1*num_X+ &
    switch_mu2*0
INTEGER, PARAMETER :: num_Phi = &               ! Only works for models with (Xg,1,1) factors
    switch_Phi1*num_X**2+ &
    switch_Phi2*(num_X*(num_X+1)/2-(num_XC(1)-num_XG)*(num_XC(2)-num_XG))
INTEGER, PARAMETER :: num_Omega = num_C
INTEGER, PARAMETER :: num_omega_dex = 1
INTEGER, PARAMETER :: num_r_inf_true = &
    0*switch_r_inf0+ &    
    0*switch_r_inf1+ &
    1*switch_r_inf2
INTEGER, PARAMETER :: num_r_inf = MAX(1,num_r_inf_true)
!
! Total number of model parameters
INTEGER, PARAMETER :: num_theta = &
    num_delta0+num_delta1+num_PhiQ+num_muQ+num_Gamm+num_mu+num_Phi+ &
    num_Omega+num_omega_dex+num_r_inf_true
!
! Total number of parameters used to compute the asymptotic variance
INTEGER, PARAMETER :: num_param = &
     num_C &                        ! delta0
    +num_XC(1)+num_XC(2) &          ! delta1
    +num_XC(1)**2+num_XC(2)**2 &    ! PhiQ
    +num_XC(1)+num_XC(2) &          ! muQ
    +num_X*(num_X+1)/2 &            ! Gamma
    +num_X &                        ! mu
    +num_X*num_X &                  ! Phi
    +num_C &                        ! Omega
    +1 &                            ! Omega_dex
    +num_r_inf_true                 ! r_inf
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters for starting point selection
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
REAL(8), PARAMETER :: theta_inf(num_theta) = -5.d-1
REAL(8), PARAMETER :: theta_sup(num_theta) = 5.d-1
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters for L-BFGS
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
REAL(8), PARAMETER :: factr = (to0+to1)*1.d+7+(to2+to4)*1.d+1
REAL(8), PARAMETER :: pgtol = (to0+to1)*1.0d-4+(to2+to4)*1.d-5
REAL(8), PARAMETER :: max_iters_lbfgs = (100*num_theta)*to1+(300*num_theta)*to2
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters for POLITOPE
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
INTEGER, PARAMETER :: max_iters_politope = (100*num_theta)*(to0+to1)+(600*num_theta)*to2
                                                ! Maximum number of iterations in politope
REAL(8), PARAMETER :: pert_theta(num_theta) = 0.5d0     
                                                ! Relative perturbation of the parameters
INTEGER, PARAMETER :: max_repeat_politope = num_theta*(to0+to1)+8*num_theta*to2
                                                ! Maximum number of politope restarts in second round optimization
REAL(8), PARAMETER :: tol = 1.d-7*(to0+to1)+1.d-11*to2
                                                ! Tolerance in second round optimization
REAL(8), PARAMETER :: tol_conv = 1.d-6*(to0+to1)+1.d-10*to2
                                                ! Convergence tolerance in second round optimization
REAL(8), PARAMETER :: tol_politope_p = tol
REAL(8), PARAMETER :: crit_politope_conv_p = tol_conv
REAL(8), PARAMETER :: tol_politope_y = tol
REAL(8), PARAMETER :: crit_politope_conv_y = tol_conv
INTEGER, PARAMETER :: rtol_formula = 2          ! Chooses the formula used to compute rtol
                                                ! See lines 150-190 in simplex_M.f90
INTEGER, PARAMETER :: crit_conv_formula = 1     ! Politope minimizations are restarted looking 
                                                ! at improvements in: 
                                                ! 1 = y
                                                ! 2 = p
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters for skewness and kurtosis in simulated yields distribution
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
INTEGER, PARAMETER :: SKttYIELDS = 1            ! From 2001 January onwards
INTEGER, PARAMETER :: num_SKFY = num_K          ! Number of yields to simulate
INTEGER, PARAMETER :: vec_SKFY(num_SKFY) = &    ! Vector of maturities of yields to simulate
    (/ 1, 2, 3, 4, 5, 6 /)
INTEGER, PARAMETER :: num_SKHY = 16             ! Number of horizons at which to simulate yields
INTEGER, PARAMETER :: vec_SKHY(num_SKHY) = &    ! Vector of horizons at which to simulate yields
    (/ 1, 2, 3, 6, 9, 12, 15, 18, 21, 24, 30, 36, 42, 48, 54, 60 /)
INTEGER, PARAMETER :: num_CC = 3                ! Number of spreads to be simulated
CHARACTER(len=4), PARAMETER :: SKlab_c(num_CC) = & ! Spread labels
    (/ "ITDE", "ITOI", "OIDE" /)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters for ETZ, LOH, LOP
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
REAL(8), PARAMETER :: sr_cutoff = 0.d0          ! Short rate cutoff
INTEGER, PARAMETER :: sr_months = 12            ! Number of months sr must stay above sr_cutoff
INTEGER, PARAMETER :: num_H = 600               ! Longest simulation horizon
INTEGER, PARAMETER :: num_S = 20000             ! Number of replications
!
INTEGER, PARAMETER :: ttYIELDS = 171            ! Date at which to forecast yields
INTEGER, PARAMETER :: num_FY = num_K            ! Number of yields to forecast
INTEGER, PARAMETER :: vec_FY(num_FY) = &        ! Vector of maturities of yields to forecast
    (/ 1, 2, 3, 4, 5, 6 /)
INTEGER, PARAMETER :: num_HY = 10               ! Number of horizons at which to forecast yields
INTEGER, PARAMETER :: vec_HY(num_HY) = &        ! Vector of horizons at which to forecast yields
    (/ 6, 12, 18, 24, 30, 36, 42, 48, 54, 60 /)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters to compute the effect of a change in the lower bound
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
INTEGER, PARAMETER :: LBttYIELDS = 1            ! From 2001 January onwards
REAL(8), PARAMETER :: deltaLB = -0.2d0          ! Assumed change in the lower bound
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters about the output files
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
CHARACTER(len=30), PARAMETER :: file_data_us = 'us.txt'
CHARACTER(len=30), PARAMETER :: file_data_uk = 'uk.txt'
CHARACTER(len=30), PARAMETER :: file_data_e = 'e.txt'
CHARACTER(len=30), PARAMETER :: file_res = 'res.txt'
CHARACTER(len=30), PARAMETER :: file_ll_params = 'll_params.txt'
CHARACTER(len=30), PARAMETER :: file_loglik = 'loglik.txt'
CHARACTER(len=30), PARAMETER :: file_theta = 'theta.txt'
CHARACTER(len=30), PARAMETER :: file_res_to1 = 'res_to1.txt'
CHARACTER(len=30), PARAMETER :: file_res_to2 = 'res_to2.txt'
CHARACTER(len=30), PARAMETER :: file_states = 'states.txt'
CHARACTER(len=30), PARAMETER :: file_kc = 'kalman_coeff.txt'
CHARACTER(len=30), PARAMETER :: file_sr = 'short_rates.txt'
CHARACTER(len=30), PARAMETER :: file_fhat = 'fhat.txt'
CHARACTER(len=30), PARAMETER :: file_ehat = 'ehat.txt'
CHARACTER(len=30), PARAMETER :: file_politope = 'politope.txt'
CHARACTER(len=30), PARAMETER :: file_datesP = 'datesP.txt'
CHARACTER(len=30), PARAMETER :: file_datesQ = 'datesQ.txt'
CHARACTER(len=30), PARAMETER :: file_resdatesP = 'resdatesP.txt'
CHARACTER(len=30), PARAMETER :: file_resdatesQ = 'resdatesQ.txt'
CHARACTER(len=30), PARAMETER :: file_resSRP = 'resSRP.txt'
CHARACTER(len=30), PARAMETER :: file_resSRQ = 'resSRQ.txt'
CHARACTER(len=30), PARAMETER :: file_resY10yLOP = 'resY10yLOP.txt'
CHARACTER(len=30), PARAMETER :: file_resY10yLOQ = 'resY10yLOQ.txt'
CHARACTER(len=30), PARAMETER :: file_resY6mLOP = 'resY6mLOP.txt'
CHARACTER(len=30), PARAMETER :: file_resY6mLOQ = 'resY6mLOQ.txt'
CHARACTER(len=30), PARAMETER :: file_resYIELDP = 'resYIELDP.txt'
CHARACTER(len=30), PARAMETER :: file_resYIELDQ = 'resYIELDQ.txt'
CHARACTER(len=30), PARAMETER :: file_resSKEWP = 'resSKEWP.txt'
CHARACTER(len=30), PARAMETER :: file_resSKEWQ = 'resSKEWQ.txt'
CHARACTER(len=30), PARAMETER :: file_resKURTP = 'resKURTP.txt'
CHARACTER(len=30), PARAMETER :: file_resKURTQ = 'resKURTQ.txt'
CHARACTER(len=30), PARAMETER :: file_resLB = 'resLB.txt'
CHARACTER(len=30), PARAMETER :: file_jmat = 'jmat.txt'
CHARACTER(len=30), PARAMETER :: file_invjmat = 'invjmat.txt'
CHARACTER(len=30), PARAMETER :: file_dl = 'dl.txt'
CHARACTER(len=30), PARAMETER :: file_imat = 'imat.txt'
CHARACTER(len=30), PARAMETER :: file_dparam_dtheta = 'dparam_dtheta.txt'
CHARACTER(len=30), PARAMETER :: file_param = 'param.txt'
CHARACTER(len=30), PARAMETER :: file_var_param = 'var_param.txt'
CHARACTER(len=30), PARAMETER :: file_fin_res = 'fin_res.txt'
!
INTEGER, PARAMETER :: unit_data_us = 51
INTEGER, PARAMETER :: unit_data_uk = 52
INTEGER, PARAMETER :: unit_data_e = 53
INTEGER, PARAMETER :: unit_res = 21
INTEGER, PARAMETER :: unit_ll_params = 24
INTEGER, PARAMETER :: unit_loglik = 30
INTEGER, PARAMETER :: unit_theta = 5
INTEGER, PARAMETER :: unit_res_to1 = 22
INTEGER, PARAMETER :: unit_res_to2 = 23
INTEGER, PARAMETER :: unit_states = 25
INTEGER, PARAMETER :: unit_kc = 26
INTEGER, PARAMETER :: unit_sr = 27
INTEGER, PARAMETER :: unit_fhat = 28
INTEGER, PARAMETER :: unit_ehat = 288
INTEGER, PARAMETER :: unit_politope = 4
INTEGER, PARAMETER :: unit_datesP = 291
INTEGER, PARAMETER :: unit_datesQ = 292
INTEGER, PARAMETER :: unit_resdatesP = 281
INTEGER, PARAMETER :: unit_resdatesQ = 282
INTEGER, PARAMETER :: unit_resSRP = 271
INTEGER, PARAMETER :: unit_resSRQ = 272
INTEGER, PARAMETER :: unit_resY10yLOP = 261
INTEGER, PARAMETER :: unit_resY10yLOQ = 262
INTEGER, PARAMETER :: unit_resY6mLOP = 251
INTEGER, PARAMETER :: unit_resY6mLOQ = 252
INTEGER, PARAMETER :: unit_resYIELDP = 241
INTEGER, PARAMETER :: unit_resYIELDQ = 242
INTEGER, PARAMETER :: unit_resSKEWP = 341
INTEGER, PARAMETER :: unit_resSKEWQ = 342
INTEGER, PARAMETER :: unit_resKURTP = 441
INTEGER, PARAMETER :: unit_resKURTQ = 442
INTEGER, PARAMETER :: unit_resLB = 541
INTEGER, PARAMETER :: unit_jmat = 127
INTEGER, PARAMETER :: unit_invjmat = 128
INTEGER, PARAMETER :: unit_dl = 129
INTEGER, PARAMETER :: unit_imat = 126
INTEGER, PARAMETER :: unit_dparam_dtheta = 123
INTEGER, PARAMETER :: unit_param = 121
INTEGER, PARAMETER :: unit_var_param = 125
INTEGER, PARAMETER :: unit_fin_res = 124
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring constants
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
REAL(8), PARAMETER :: minimum_p = 1.d-10        ! Minimum probability value
REAL(8), PARAMETER :: pi = 3.14159265358979323846264338328d0    ! Pi
REAL(8), PARAMETER :: invsqrtpi = 0.564189583547756286948079451561d0    ! 1/Sqrt[Pi]
REAL(8), PARAMETER :: sqrt2 = 1.41421356237309504880168872421d0 ! Sqrt[2]
REAL(8), PARAMETER :: invsqrt2pi = 0.398942280401432677939946059934d0   ! 1/Sqrt[2*Pi]
REAL(8), PARAMETER :: dueinvpi = 0.636619772367581343075535053490d0 ! 2/Pi
REAL(8), PARAMETER :: pisudue = 1.570796326794896619231321691639d0  ! Pi/2
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE constants
