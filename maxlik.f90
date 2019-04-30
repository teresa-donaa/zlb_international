MODULE maxlik
!
USE constants
USE observations
USE printing_routines
USE gatsm
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    !SUBROUTINE estimate_lbfgs ( i_stime, theta, objf, grad, task )
    !!
    !IMPLICIT NONE
    !! 
    !! Declaring dummy variables
    !!
    !INTEGER, INTENT(IN) :: i_stime
    !REAL(8), INTENT(INOUT) :: theta(num_theta)
    !REAL(8), INTENT(OUT) :: objf
    !REAL(8), INTENT(OUT) :: grad(num_theta)
    !CHARACTER(len=60), INTENT(OUT) :: task
    !!
    !! Declaring local variables
    !!
    !INTEGER :: isave(44), iter, i
    !INTEGER, PARAMETER :: n = num_theta, m = 20, iprint = -1
    !INTEGER, ALLOCATABLE :: nbd(:), iwa(:)
    !REAL(8) :: dsave(29)
    !REAL(8), ALLOCATABLE :: wa(:), l(:), u(:)
    !CHARACTER(len=60) :: csave
    !LOGICAL :: lsave(4)
    !!
    !! Beginning execution
    !!
    !ALLOCATE ( nbd(n), iwa(3*n) )
    !ALLOCATE ( wa(2*m*n + 5*n + 11*m*m + 8*m), l(n), u(n) )
    !nbd = 0
    !task = 'START'
    !iter = 0
    !DO WHILE ((task(1:2) .EQ. 'FG') .OR. (task .EQ. 'NEW_X') .OR. (task .EQ. 'START')) 
    !    !
    !    CALL setulb ( n, m, theta, l, u, nbd, objf, grad, factr, pgtol, &
    !        wa, iwa, task, iprint, csave, lsave, isave, dsave )
    !    ! 
    !    IF (task(1:5) .EQ. 'NEW_X') THEN
    !        CALL print_screen(i_stime,iter,objf,theta,grad)
    !        IF (iter .GT. max_iters_lbfgs) EXIT
    !        iter = iter+1
    !    END IF
    !    IF (task(1:2) .eq. 'FG') CALL loglik_lbfgs(n,theta,objf,grad)
    !    !
    !END DO
    !DEALLOCATE ( nbd, iwa )
    !DEALLOCATE ( wa, l, u )
    !!
    !! Ending execution and returning control
    !!
    !END SUBROUTINE estimate_lbfgs
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    SUBROUTINE loglik_lbfgs ( nn, theta, fc, gc ) 
!    !
!    ! Declaring dummy variables
!    !
!    INTEGER, INTENT(IN) :: nn
!    REAL(8), INTENT(IN) :: theta(nn)        
!    REAL(8), INTENT(OUT) :: fc
!    REAL(8), INTENT(OUT) :: gc(nn)
!    !
!    ! Declaring local variables
!    !
!    REAL(8) :: objf, objfd(num_theta)
!    REAL(8), DIMENSION(num_X,num_T) :: xu, xf
!    REAL(8), DIMENSION(num_X,num_X,num_T) :: Pu, Pf
!    REAL(8) :: delta0(num_C), r_inf, Omega(num_C*num_K,num_C*num_K), work(10*num_X)
!    REAL(8), DIMENSION(num_X) ::  muQ, mu, wr, wi
!    REAL(8), DIMENSION(1,num_X) :: vl, vr
!    REAL(8), DIMENSION(num_X,num_X) :: PhiQ, Phi, Gamm, Sigma
!    REAL(8), DIMENSION(num_X,num_C) :: delta1
!    INTEGER :: info
!    !
!    ! Beginning execution
!    !
!    ! Call to loglikelihood subroutine
!    !
!    CALL ll(theta,objf,xf,xu,Pf,Pu)
!!    CALL ll_dv(theta,eye_theta,objf,objfd,xf,xu,Pf,Pu,num_theta)
!    !
!    fc = -objf
!    gc = -objfd
!    IF (to0 .EQ. 1) THEN
!        CALL open_write_file(unit_loglik,file_loglik)
!        WRITE(unit_loglik,1) fc
!1       FORMAT(ES15.8)
!        CLOSE(UNIT=unit_loglik)
!    END IF
!    !
!    ! Ending execution and returning control
!    !
!    END SUBROUTINE loglik_lbfgs
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION loglik_politope ( theta ) 
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)   
    !
    ! Declaring function's type
    !
    REAL(8) :: loglik_politope
    !
    ! Declaring local variables
    !
    REAL(8) :: objf
    REAL(8), DIMENSION(num_X,num_T) :: xu, xf
    REAL(8), DIMENSION(num_X,num_X,num_T) :: Pu, Pf
    !
    ! Beginning execution
    !
    ! Call to loglikelihood subroutine
    !
    CALL ll(theta,objf,xf,xu,Pf,Pu)
    IF(ISNAN(objf)) objf = -1.d25
    IF(objf .GT. 1.d3) objf = -1.d25
    loglik_politope = -objf
    !
    ! Ending execution and returning control
    !
    END FUNCTION loglik_politope
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE maxlik
