MODULE observations
!
! Global data 
!
USE constants
!
IMPLICIT NONE
!
INTEGER, SAVE :: date(num_T)
REAL(8), SAVE :: F(num_C*num_K,num_T)               ! forward rates
REAL(8), SAVE :: exrate(num_T)                      ! exchange rate
INTEGER, SAVE :: Tau_mat(num_C*num_K,num_T)         ! Tau_mat = (TauS_mat, TauR_mat)
INTEGER, SAVE :: K_vec(num_T)                       ! K_vec = KS_vec+KR_vec
INTEGER, SAVE :: Tau(num_C*num_K)                   ! Observed maturities
REAL(8), SAVE :: eye_theta(num_theta,num_theta)
REAL(8), SAVE :: eyeN(num_X,num_X)
REAL(8), SAVE :: eyeN2(num_X**2,num_X**2)
!
! Ending module
!
END MODULE observations
