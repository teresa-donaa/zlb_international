MODULE printing_routines
!
USE constants
USE observations
USE gatsm
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE open_read_file ( unit_number, file_name )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unit_number
    CHARACTER(len=30), INTENT(IN) :: file_name
    !
    ! Declaring local variables
    !
    INTEGER :: open_err                 ! Open file error code
    !
    ! Beginning execution
    !
    OPEN ( UNIT=unit_number, FILE=file_name, STATUS='old', IOSTAT=open_err )
    IF (open_err .NE. 0) THEN
        WRITE(*,1234) file_name
        STOP
    END IF
    1234 FORMAT ('Unable to open input file ', A20)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE open_read_file
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE open_write_file ( unit_number, file_name )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unit_number
    CHARACTER(len=30), INTENT(IN) :: file_name
    !
    ! Declaring local variables
    !
    INTEGER :: open_err                 ! Open file error code
    !
    ! Beginning execution
    !
    OPEN ( UNIT=unit_number, FILE=file_name, STATUS='replace', IOSTAT=open_err )
    IF (open_err .NE. 0) THEN
        WRITE(*,1234) file_name
        STOP
    END IF
    1234 FORMAT ('Unable to open input file ', A20)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE open_write_file
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE print_screen_politope ( i_stime, iter, objf )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: i_stime      ! Estimation trial number
    INTEGER, INTENT(IN) :: iter         ! Iteration number
    REAL(8), INTENT(IN) :: objf         ! Latent criterion function at the optimum
    !
    ! Beginning execution
    !
    WRITE(*,19) i_stime, iter, objf
19  FORMAT('Stima ', I3, ', Iter: ', I5, ' ; -LL/N: ', ES15.8)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE print_screen_politope 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    SUBROUTINE print_screen ( i_stime, iter, objf, theta, grad )
!    !
!    IMPLICIT NONE
!    !
!    ! Declaring dummy variables
!    !
!    INTEGER, INTENT(IN) :: i_stime      ! Estimation trial number
!    INTEGER, INTENT(IN) :: iter         ! Iterationnumber
!    REAL(8), INTENT(IN) :: objf         ! Latent criterion function at the optimum
!    REAL(8), INTENT(IN) :: theta(num_theta)
!    REAL(8), INTENT(IN) :: grad(num_theta)
!    !
!    ! Beginning execution
!    !
!    ! Printing to screen
!    !
!    WRITE(*,19) i_stime, iter, objf, SQRT(SUM(grad**2))
!19  FORMAT('Stima ', I2, ', Iter: ', I5, ' ; -LL/N: ', ES15.8,' ; ngrad: ', ES12.5)
!    !
!    ! Ending execution and returning control
!    !
!    END SUBROUTINE print_screen
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE print_res ( i_stime, objf, theta, task, grad )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: i_stime      ! Estimation trial number
    REAL(8), INTENT(IN) :: objf         ! Latent criterion function at the optimum
    REAL(8), INTENT(IN) :: theta(num_theta)
    CHARACTER(len=60), INTENT(IN) :: task
    REAL(8), INTENT(IN) :: grad(num_theta)
    !
    ! Declaring local variables
    !
    INTEGER :: i, ic
    REAL(8) :: delta0(num_C), r_inf, omega_dex
    REAL(8), DIMENSION(num_X) :: mu
    REAL(8), DIMENSION(num_XC(1)) :: mu_1, muQ_1, delta1_1, lambda_1
    REAL(8), DIMENSION(num_XC(2)) :: mu_2, muQ_2, delta1_2, lambda_2
    REAL(8), DIMENSION(num_X,num_X) :: Gamm, Sigma, Phi
    REAL(8), DIMENSION(num_XC(1),num_XC(1)) :: PhiQ_1, Sigma_1, Phi_1, LLambda_1
    REAL(8), DIMENSION(num_XC(2),num_XC(2)) :: PhiQ_2, Sigma_2, Phi_2, LLambda_2
    REAL(8), DIMENSION(num_C*num_K,num_C*num_K) :: Omega
    !
    ! Beginning execution
    !
    ! Print res file
    !
    WRITE (unit_res,35) i_stime, objf, &
        theta, task, grad
35  FORMAT ( I3, " # ", ES25.16E3, &
        " #", <num_theta>(1X, ES25.16E3), " # ", A60, " # ", <num_theta>(ES25.16E3,1X))
    !
    ! Print ll_params file
    !
    CALL theta_to_param(theta,delta0,delta1_1,delta1_2,PhiQ_1,PhiQ_2, &
        muQ_1,muQ_2,Gamm,Sigma,Sigma_1,Sigma_2,mu,mu_1,mu_2,lambda_1,lambda_2, &
        Phi,Phi_1,Phi_2,LLambda_1,LLambda_2,Omega,omega_dex,r_inf)
    WRITE (unit_ll_params,36) i_stime, objf, &
        r_inf, delta0, delta1_1, delta1_2, &
        (PhiQ_1(i,:), i = 1, num_XC(1)), (PhiQ_2(i,:), i = 1, num_XC(2)), &
        muQ_1, muQ_2, (Gamm(i,:), i = 1, num_X), &
        mu, (Phi(i,:), i = 1, num_X), (SQRT(Omega(i,i)), i = 1, num_C*num_K), SQRT(omega_dex), task
36  FORMAT(I3, " # ", ES25.16E3, &
        " # ", <num_r_inf+num_C+num_XC(1)+num_XC(2)+ &
            num_XC(1)**2+num_XC(2)**2+ &
            num_XC(1)+num_XC(2)+num_X**2+ &
            num_X+num_X**2+num_C*num_K+1>(ES25.16E3, 1X), " # ", A60)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE print_res
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE head_ll_params ( )
    !
    IMPLICIT NONE
    !
    ! Declaring local variables
    !
    INTEGER :: i, j, ic
    !
    ! Beginning execution
    !
    WRITE(unit_ll_params,11) (i, i = 1, num_C), &
        (i, i = 1, num_XC(1)), (i, i = 1, num_XC(2)), &
        ((i, j, j = 1, num_XC(1)), i = 1, num_XC(1)), &
        ((i, j, j = 1, num_XC(2)), i = 1, num_XC(2)), &
        (i, i = 1, num_XC(1)), (i, i = 1, num_XC(2)), &
        ((i, j, j = 1, num_X), i = 1, num_X), (i, i = 1, num_X), &
        ((i, j, j = 1, num_X), i = 1, num_X), &
        ((lab_c(ic), i, i, i = 1, num_K), ic = 1, num_C)
11  FORMAT("  i #                    '-ll/T # ", &
        "                    r_inf ", &
        <num_C>("                delta0(", I1, ") "), &
        <num_XC(1)>("              delta1_1(", I1, ") "), &
        <num_XC(2)>("              delta1_2(", I1, ") "), &
        <num_XC(1)**2>("              PhiQ_1(", I1, ",", I1, ") "), &
        <num_XC(2)**2>("              PhiQ_2(", I1, ",", I1, ") "), &
        <num_XC(1)>("                 muQ_1(", I1, ") "), &
        <num_XC(2)>("                 muQ_2(", I1, ") "), &
        <num_X**2>("               Gamma(", I1, ",", I1, ") "), &    
        <num_X>("                    mu(", I1, ") "), &
        <num_X**2>("                 Phi(", I1, ",", I1, ") "), &   
        <num_C>(<num_K>("            Omega", A3, "(", I1, ",", I1, ") ")), &
        "                omega_dex ", &
        " # task")
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE head_ll_params
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    SUBROUTINE print_final_results ( param, stderr, objf, grad )
!    !
!    ! Computes final results and writes them on file
!    !
!    IMPLICIT NONE
!    !
!    ! Declaring dummy variables
!    !
!    REAL(8), INTENT(IN) :: param(num_param)             ! Estimates of parameters
!    REAL(8), INTENT(IN) :: stderr(num_param)            ! Estimated asymptotic standard errors
!    REAL(8), INTENT(IN) :: objf                         ! Optimized latent criterion
!    REAL(8), INTENT(IN) :: grad(num_param)              ! Estimated asymptotic standard errors
!    !
!    ! Declaring local variables
!    !
!    INTEGER :: ind, ir, ic                              ! Loop indexes
!    !
!    ! Beginning execution
!    !
!    ! Computing statistics for the parameter estimates
!    !
!    ! delta0
!    !
!    WRITE (unit_fin_res,3525)
!3525 FORMAT ( /, &
!        '-----------------------------------------------------------------', /, &
!        '    Parameter     Estimate        AsySE   AsyT-Ratio     Gradient', /, &
!        '-----------------------------------------------------------------' )
!    !
!    ind = 0
!    DO ir = 1, num_C
!        ind = ind+1
!        WRITE (unit_fin_res, 9000) ir, param(ind), stderr(ind), &
!            param(ind)/stderr(ind), grad(ind)
!        9000 FORMAT ( 1X, 'delta0(', I1, ')   ', 4(1X, ES12.5) )
!    END DO
!    !
!    ! delta1
!    !
!    DO ic = 1, num_C
!        DO ir = 1, num_X
!            ind = ind+1
!            WRITE (unit_fin_res, 9001) ir, ic, param(ind), stderr(ind), &
!                param(ind)/stderr(ind), grad(ind)
!        9001 FORMAT ( 1X, 'delta1(', I1, ',', I1, ') ', 4(1X, ES12.5) )
!        END DO
!    END DO
!    !
!    ! PhiQ
!    !
!    DO ir = 1, num_X
!        DO ic = 1, ir
!            ind = ind+1
!            WRITE (unit_fin_res, 9002) ir, ic, param(ind), stderr(ind), &
!                param(ind)/stderr(ind), grad(ind)
!        9002 FORMAT ( 1X, 'PhiQ(', I1, ',', I1, ')   ', 4(1X, ES12.5) )
!        END DO
!    END DO
!    !
!    ! muQ
!    !
!    DO ir = 1, num_X
!        ind = ind+1
!        WRITE (unit_fin_res, 9003) ir, param(ind), stderr(ind), &
!            param(ind)/stderr(ind), grad(ind)
!        9003 FORMAT ( 1X, 'muQ(', I1, ')      ', 4(1X, ES12.5) )
!    END DO
!    !
!    ! Gamma
!    !
!    DO ir = 1, num_X
!        DO ic = 1, ir
!            ind = ind+1
!            WRITE (unit_fin_res, 9004) ir, ic, param(ind), stderr(ind), &
!                param(ind)/stderr(ind), grad(ind)
!        9004 FORMAT ( 1X, 'Gamma(', I1, ',', I1, ')  ', 4(1X, ES12.5) )
!        END DO
!    END DO
!    !
!    ! mu
!    !
!    DO ir = 1, num_X
!        ind = ind+1
!        WRITE (unit_fin_res, 9005) ir, param(ind), stderr(ind), &
!            param(ind)/stderr(ind), grad(ind)
!        9005 FORMAT ( 1X, 'mu(', I1, ')       ', 4(1X, ES12.5) )
!    END DO
!    !
!    ! Phi
!    !
!    DO ir = 1, num_X
!        DO ic = 1, num_X
!            ind = ind+1
!            WRITE (unit_fin_res, 9006) ir, ic, param(ind), stderr(ind), &
!                param(ind)/stderr(ind), grad(ind)
!        9006 FORMAT ( 1X, 'Phi(', I1, ',', I1, ')    ', 4(1X, ES12.5) )
!        END DO
!    END DO
!    !
!    ! Omega
!    !
!    DO ir = 1, num_C
!        ind = ind+1
!        WRITE (unit_fin_res, 9007) ir, param(ind), stderr(ind), &
!            param(ind)/stderr(ind), grad(ind)
!        9007 FORMAT ( 1X, 'Omega(', I1, ')    ', 4(1X, ES12.5) )
!    END DO
!    !
!    ! r_inf
!    !
!    DO ir = 1, num_r_inf
!        ind = ind+1
!        WRITE (unit_fin_res, 9008) ir, param(ind), stderr(ind), &
!            param(ind)/stderr(ind), grad(ind)
!        9008 FORMAT ( 1X, 'r_inf(', I1, ')    ', 4(1X, ES12.5) )
!    END DO
!    !
!    ! Ending execution and returning control
!    !
!    END SUBROUTINE print_final_results
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE printing_routines
