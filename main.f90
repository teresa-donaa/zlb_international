PROGRAM main
!
USE constants
USE observations
USE printing_routines
USE load_data
USE starting_points
USE maxlik
USE simplex
USE diagnostics
!USE skewness_kurtosis_yields
!USE lb_change_effects
!USE asymptotic_variance
!
IMPLICIT NONE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring variables
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
INTEGER :: i_stime, i, j, iter, errcode, nlines
CHARACTER(len=2) :: ichar
INTEGER :: seed(2)
REAL(8) :: theta(num_theta), objf, min_objf, grad(num_theta), llcheck
!REAL(8), DIMENSION(num_X,num_T) :: xf, xu
!REAL(8), DIMENSION(num_X,num_X,num_T) :: Pf, Pu
CHARACTER(len=60) :: task
!REAL(8), DIMENSION(num_param) :: param, stderr, grad_param
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Beginning execution
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
IF ((to0 .EQ. 1) .OR. (to1 .EQ. 1) .OR. (to2 .EQ. 1)) THEN
    !
    ! Loading data
    !
    CALL input_data()
    !
    ! Initialising random number generator
    !
    IF (to1 .EQ. 1) THEN
        !
        CALL GETARG(1,ichar)
        seed = INUM(ichar)
        IF (ALL(seed .EQ. 0)) seed = 3
        !
    END IF    
    !
    ! If to1, checking if a res.txt file already exists
    !
    IF (to1 .EQ. 1) THEN
    !    !
    !    OPEN(UNIT = unit_res,FILE = file_res,READONLY,IOSTAT = errcode)
    !    nlines = 0 
    !    IF (errcode .GT. 0) THEN
    !        !
            CALL open_write_file(unit_res,file_res)
            CALL open_write_file(unit_ll_params,file_ll_params)
            CALL head_ll_params()
    !        !
    !    ELSE
    !        !
    !        DO
    !            !
    !            READ(unit_res,*,IOSTAT = errcode)
    !            IF (errcode .NE. 0) EXIT
    !            nlines = nlines+1
    !            !
    !        END DO
    !        CLOSE(UNIT = unit_res)
    !        OPEN(UNIT = unit_res,FILE = file_res,POSITION = 'APPEND')
    !        OPEN(UNIT = unit_ll_params,FILE = file_ll_params,POSITION = 'APPEND')
    !        !
    !    END IF
    !    !
    END IF
    !
    ! If to2, open output files
    !
    IF (to2 .EQ. 1) THEN
        !
        CALL open_write_file(unit_res,file_res)
        CALL open_write_file(unit_ll_params,file_ll_params)
        CALL head_ll_params()
        !
    END IF
    !
    ! Starting loop 
    !
    min_objf = 1.d3
    DO i_stime = 1, num_stime
        !
        ! Creating random starting values of the parameters 
        !
        CALL admissible_starting_point(seed,theta,llcheck)
        IF ((to1 .EQ. 1) .AND. (i_stime .LE. nlines)) CYCLE 
        objf = 1.d3
        !
        ! POLITOPE Estimation 
        !
        IF (switch_politope .EQ. 1) THEN
            !
!            objf = loglik_politope(theta)
            CALL open_write_file(unit_politope,file_politope)
            CALL politope(i_stime,loglik_politope,theta,pert_theta,1,250,objf,iter)
            CLOSE(UNIT=unit_politope)
            min_objf = MIN(objf,min_objf)
            !
        END IF
        !
        ! L-BFGS Estimation 
        !
!        IF ((switch_lbfgs .EQ. 1) .AND. (objf .LE. min_objf+2.d0)) THEN
!            !
!            CALL estimate_lbfgs(i_stime,theta,objf,grad,task)
!            min_objf = MIN(objf,min_objf)
!            !
!        END IF
        !
        ! Printing intermediate trace output 
        !
        CALL print_res(i_stime,objf,theta,task,grad)
        ! 
    END DO 
    !
    ! Closes optimization stage
    !
    CLOSE(UNIT=unit_res)
    CLOSE(UNIT=unit_ll_params)
    !
END IF
!!
!! Computing skewness matrix
!!
!IF (to3 .EQ. 1) THEN
!    !
!    ! Reading parameters estimates and data
!    !
!    CALL admissible_starting_point(seed,theta,llcheck)
!    CALL input_data()
!    !
!    ! Computing and writing diagnostics
!    !
!    CALL compute_skewness_kurtosis_yields(theta,llcheck,objf,grad)
!    !
!END IF
!
! Computing diagnostics
!
IF (to4 .EQ. 1) THEN
    !
    ! Reading parameters estimates and data
    !
    CALL admissible_starting_point(seed,theta,llcheck)
    CALL input_data()
    !
    ! Computing and writing diagnostics
    !
    CALL compute_diagnostics(theta,llcheck,objf,grad)
!    CALL compute_liftoff_yields(theta,llcheck,objf,grad)
    !
END IF
!!
!! Computing effects of a change in the lower bound
!!
!IF (to5 .EQ. 1) THEN
!    !
!    ! Reading parameters estimates and data
!    !
!    CALL admissible_starting_point(seed,theta,llcheck)
!    CALL input_data()
!    !
!    ! Computing and writing diagnostics
!    !
!    CALL compute_lb_change_effects(theta,llcheck,objf,grad)
!    !
!END IF
!!
!! Computing diagnostics
!!
!IF (to6 .EQ. 1) THEN
!    !
!    ! Reading parameters estimates and data
!    !
!    CALL admissible_starting_point(seed,theta,llcheck)
!    CALL input_data()
!    !
!    ! Computing and writing asymptotic standard errors
!    !
!    CALL compute_asymptotic_variance(theta,llcheck,objf,grad_param,param,stderr)
!    CALL open_write_file(unit_fin_res,file_fin_res)
!    CALL print_final_results(param,stderr,objf,grad_param)
!    CLOSE(UNIT=unit_fin_res)
!    !
!END IF
!
END PROGRAM main
