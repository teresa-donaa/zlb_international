MODULE load_data
!
USE constants
USE observations
USE printing_routines
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE input_data ( )
	!
	IMPLICIT NONE
	! 
	! Declaring local variables
	!
	INTEGER :: io_status, datett  				                    ! Open error code
	INTEGER :: tt, ik, im, ic, h                                    ! Loop index
	REAL(8) :: fc(num_K,num_T), e_t, de_t
	!
	! Beginning execution
	!
	! Reading forward rates
	!
    h = 1
    CALL input_yields_country(unit_data_us,file_data_us,num_K_US,data_US_cols,fc)
    F((h-1)*num_K+1:h*num_K,:) = fc
    !
    h = 2
    CALL input_yields_country(unit_data_uk,file_data_uk,num_K_UK,data_UK_cols,fc)
    F((h-1)*num_K+1:h*num_K,:) = fc
    !
    ! Reading exchange rates
    !
    CALL open_read_file(unit_data_e,file_data_e)
	!
    ! Skip the first row of the data files (containing names and maturities)
	READ(unit_data_e,101)
101 FORMAT()
	!
	! Reading remaining lines
    tt = 0
	DO 
        READ(unit_data_e,*,IOSTAT = io_status) datett, e_t, de_t
        IF (io_status .LT. 0) EXIT
        tt = tt+1
        IF ((tt .GE. ind_obs1) .AND. (tt .LE. ind_obsT)) THEN
	        PRINT*, 'tt = ', tt, ' ; e_t = ', e_t
            exrate(tt) = e_t
            dexrate(tt) = de_t
        END IF
    END DO
	CLOSE(UNIT = unit_data_e)
	!
    ! Creating maturities
    !
    DO im = 1, num_C
        Tau((im-1)*num_K+1:im*num_K) = Tau_vec
    END DO
	!
    ! Generating:
    !   Tau_mat = matrix of indexes of available observations
    !   K_vec   = vector of numbers of available observations
    !
    Tau_mat = 0
    K_vec = 0
    DO tt = 1, num_T
        im = 0
        DO ik = 1, num_C*num_K
            IF (F(ik,tt) .GT. -90.d0) THEN
                im = im+1
                Tau_mat(im,tt) = ik
            END IF
        END DO
        K_vec(tt) = im
    END DO
    !
    ! Constructing identity matrices of orders num_theta, num_X and num_X^2
    !
    eye_theta = 0.d0
    DO im = 1, num_theta
        eye_theta(im,im) = 1.d0
    END DO
    eyeX = 0.d0
    DO im = 1, num_X
        eyeX(im,im) = 1.d0
    END DO
    eyeXI = 0.d0
    DO im = 1, num_XI
        eyeXI(im,im) = 1.d0
    END DO
    eyeXI2 = 0.d0
    DO im = 1, num_XI**2
        eyeXI2(im,im) = 1.d0
    END DO
    !
    ! Ending execution and returning control
    !
	END SUBROUTINE input_data
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE input_yields_country ( unit, file, num_mat, ind_cols, yc )
	!
	IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unit
    CHARACTER(len=30), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: num_mat
    INTEGER, INTENT(IN) :: ind_cols(num_K)
    REAL(8), INTENT(OUT) :: yc(num_K,num_T)
	! 
	! Declaring local variables
	!
	INTEGER :: io_status  				                        ! Open error code
	INTEGER :: tt, datett
	REAL(8) :: ytt(num_mat)
	!
	! Beginning execution
	!
	! Reading yields
	!
    CALL open_read_file(unit,file)
	!
    ! Skip the first row of the data files (containing names and maturities)
	READ(unit,101)
101 FORMAT()
	!
	! Reading remaining lines
    tt = 0
	DO 
        READ(unit,*,IOSTAT=io_status) datett, ytt
        IF (io_status .LT. 0) EXIT
        tt = tt+1
	    PRINT*, 'it = ', tt
        IF ((tt .GE. ind_obs1) .AND. (tt .LE. ind_obsT)) THEN
            yc(:,tt) = ytt(ind_cols)
        END IF
    END DO
	CLOSE(UNIT=unit)
    !
    ! Ending execution and returning control
    !
	END SUBROUTINE input_yields_country
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE load_data