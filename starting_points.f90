MODULE starting_points
!
USE constants
USE observations
USE printing_routines
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	SUBROUTINE admissible_starting_point ( seed, theta, llcheck )
	!
	! Computes an admissible starting parameters vector
	!
	IMPLICIT NONE
	!
	! Declaring dummy variables
	!
	INTEGER, INTENT(INOUT) :: seed(2)				! Seed for r.n. generation
	REAL(8), INTENT(OUT) :: theta(num_theta)		! Admissible value of the parameters vector
	REAL(8), INTENT(OUT) :: llcheck
    !
    ! Declaring external routine
    !
    REAL(8), EXTERNAL :: r8_uniform_01
    !
	! Declaring local variables
	!
	INTEGER :: i									! Index
    CHARACTER(len=2) :: ichar
	!
	! Beginning execution
	!
    ! Modalità di test
    !
    IF (to0 .EQ. 1) THEN
        !
        CALL open_read_file(unit_theta,file_theta)
        READ(unit_theta,*) (theta(i), i = 1, num_theta)
        CLOSE(UNIT=unit_theta)
        !
    END IF
    !
	! Step 1: Starting points chosen randomly
	!
    IF (to1 .EQ. 1) THEN
		!
        CALL initialize()
        CALL set_initial_seed(seed(1),seed(2))
        !
        ! Totally random choice
        !
        IF (switch_inival1 .EQ. 1) THEN
            !
		    DO i = 1, num_theta
                theta(i) = theta_inf(i)+(theta_sup(i)-theta_inf(i))*r8_uniform_01()
            END DO 
            !
        END IF
        IF (switch_inival2 .EQ. 1) THEN
            !
            CALL HWinival(seed,theta)
            !
        END IF
		!
    END IF
	!
	! Step 2: Starting point read from the res_to1.txt file
	!
	IF (to2 .EQ. 1) THEN
        !
        CALL GETARG(1,ichar)
        i = INUM(ichar)
        IF (i .EQ. 0) i = 1
        !
        CALL open_read_file(unit_res_to1,file_res_to1)
        IF (i .GT. 1) THEN
            READ(unit_res_to1, 2534) 
2534        FORMAT(<i-1>(/))
            BACKSPACE(UNIT=unit_res_to1)
        END IF
        READ(unit_res_to1,*) llcheck, theta
        CLOSE(UNIT=unit_res_to1)
        !
	END IF 
	!
	! Step 3 or 4: Starting point read from the res_to2.txt file
	!
	IF ((to3 .EQ. 1) .OR. (to4 .EQ. 1)) THEN
        !
        CALL GETARG(1,ichar)
        i = INUM(ichar)
        !
        CALL open_read_file(unit_res_to2,file_res_to2)
        IF (i .GT. 1) THEN
            READ(unit_res_to2, 2535) 
2535        FORMAT(<i-1>(/))
            BACKSPACE(UNIT=unit_res_to2)
        END IF
        READ(unit_res_to2,11) llcheck, theta
11      FORMAT(6X,ES25.18,3X,<num_theta>(ES25.18,1X))
        CLOSE(UNIT=unit_res_to2)
        CALL open_write_file(unit_theta,file_theta)
        WRITE(unit_theta,117) theta
117     FORMAT(<num_theta>( ES25.18,1X))
        CLOSE(unit_theta)
        !
    END IF 
    !
    END SUBROUTINE admissible_starting_point
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE HWinival ( seed, theta )
	!
	! Computes an admissible starting parameters vector following Hamilton and Wu (2012), p. 319 col. 1
	!
	IMPLICIT NONE
	!
	! Declaring dummy variables
	!
	INTEGER, INTENT(INOUT) :: seed(2)				! Seed for r.n. generation
	REAL(8), INTENT(OUT) :: theta(num_theta)		! Admissible value of the parameters vector
    !
    ! Declaring external routine
    !
    REAL(8), EXTERNAL :: r8_uniform_01
    !
	! Declaring local variables
	!
    REAL(8) :: PhiQ_1(num_XC(1),num_XC(2)), PhiQ_2(num_XC(2),num_XC(2))
    REAL(8) :: Phi(num_X,num_X), subF(num_T)
    INTEGER :: cc, il, iu, ir, ic, h, ncc
    !
	! Beginning execution
	!
    ! delta0
    !
    il = 0
    IF (switch_delta01 .EQ. 1) THEN                 ! delta0 unconstrained
        DO cc = 1, num_C
            subF = F((cc-1)*num_K+1,:)
            theta(il+cc) = SUM(subF,MASK = subF .GT. -90.d0)/COUNT(MASK = subF .GT. -90.d0)
        END DO
    END IF
    iu = num_delta0
    !
    ! delta1
    !
    il = iu
    theta(il+1:il+num_delta1) = LOG(1.d-4)
    iu = il+num_delta1
    !
    ! PhiQ_1
    !
    ! First, create PhiQ_1 and PhiQ_2
    !
    PhiQ_1 = 0.d0
    DO ir = 1, num_XC(1)
        PhiQ_1(ir,ir) = 0.5d0+0.5d0*r8_uniform_01()
    END DO
    PhiQ_2 = 0.d0
    DO ir = 1, num_XC(2)
        PhiQ_2(ir,ir) = 0.5d0+0.5d0*r8_uniform_01()
    END DO
    !
    ! Second, extract theta from PhiQ_1
    !
    il = iu
    h = 0
    IF (switch_PhiQ1 .EQ. 1) THEN           ! PhiQ_1 lower triangular
        DO ir = 1, num_XC(1)
            DO ic = 1, ir
                h = h+1
                IF (ir .EQ. ic) THEN
                    theta(il+h) = TAN(pisudue*(2.d0*PhiQ_1(ir,ir)-1.d0))
                ELSE IF (ir .GT. ic) THEN
                    theta(il+h) = PhiQ_1(ir,ic)
                END IF
            END DO
        END DO
        iu = iu+num_XC(1)*(num_XC(1)+1)/2
    END IF
    IF (switch_PhiQ2 .EQ. 1) THEN           ! PhiQ_1 full
        !
        ! S1 = lower triangular; N(N+1)/2 parameters
        !
        h = 0
        DO ir = 1, num_XC(1)
            DO ic = 1, ir
                h = h+1
                IF (ir .EQ. ic) THEN
                    theta(il+h) = SQRT(-LOG(PhiQ_1(ir,ic)))
                ELSE IF (ir .GT. ic) THEN
                    theta(il+h) = 0.d0
                END IF
            END DO
        END DO
        !
        ! A = antisymmetric matrix (A = -A'); N(N-1)/2 parameters
        !
        il = il+num_XC(1)*(num_XC(1)+1)/2
        h = 0
        DO ir = 2, num_XC(1)
            DO ic = 1, ir-1
                h = h+1
                theta(il+h) = 0.d0
            END DO
        END DO
        iu = iu+num_XC(1)**2
        !
    END IF
    !
    ! Third, extract theta from PhiQ_2
    !
    il = iu
    h = 0
    IF (switch_PhiQ1 .EQ. 1) THEN           ! PhiQ_2 lower triangular
        DO ir = 1, num_XC(2)
            DO ic = 1, ir
                h = h+1
                IF (ir .EQ. ic) THEN
                    theta(il+h) = TAN(pisudue*(2.d0*PhiQ_2(ir,ir)-1.d0))
                ELSE IF (ir .GT. ic) THEN
                    theta(il+h) = PhiQ_2(ir,ic)
                END IF
            END DO
        END DO
        iu = iu+num_XC(2)*(num_XC(2)+1)/2
    END IF
    IF (switch_PhiQ2 .EQ. 1) THEN           ! PhiQ_2 full
        !
        ! S1 = lower triangular; N(N+1)/2 parameters
        !
        h = 0
        DO ir = 1, num_XC(2)
            DO ic = 1, ir
                h = h+1
                IF (ir .EQ. ic) THEN
                    theta(il+h) = SQRT(-LOG(PhiQ_2(ir,ic)))
                ELSE IF (ir .GT. ic) THEN
                    theta(il+h) = 0.d0
                END IF
            END DO
        END DO
        !
        ! A = antisymmetric matrix (A = -A'); N(N-1)/2 parameters
        !
        il = il+num_XC(2)*(num_XC(2)+1)/2
        h = 0
        DO ir = 2, num_XC(2)
            DO ic = 1, ir-1
                h = h+1
                theta(il+h) = 0.d0
            END DO
        END DO
        iu = iu+num_XC(2)**2
        !
    END IF
    !
    ! muQ
    !
    il = iu
    IF (switch_muQ1 .EQ. 1) THEN
        theta(il+1:il+num_XC(1)+num_XC(2)) = 0.d0
    END IF
    iu = iu+num_muQ
    !
    ! mu
    !
    il = iu
    IF (switch_mu1 .EQ. 1) THEN
        theta(il+1:il+num_X) = 0.d0
    END IF
    iu = iu+num_mu
    !
    ! Phi
    !
    ! First, create Phi as the PhiQs
    !
    Phi = 0.d0
    DO ir = 1, num_X
        Phi(ir,ir) = 0.5d0+0.5d0*r8_uniform_01()
    END DO
    !
    ! Second, extract theta
    !
    il = iu
    IF (switch_Phi1 .EQ. 1) THEN        ! Phi full
        !
        ! S1 = lower triangular; N(N+1)/2 parameters
        !
        h = 0
        DO ir = 1, num_X
            DO ic = 1, ir
                h = h+1
                IF (ir .EQ. ic) THEN
                    theta(il+h) = SQRT(-LOG(Phi(ir,ic)))
                ELSE IF (ir .GT. ic) THEN
                    theta(il+h) = 0.d0
                END IF
            END DO
        END DO
        !
        ! A = antisymmetric matrix (A = -A'); N(N-1)/2 parameters
        !
        il = il+num_X*(num_X+1)/2
        h = 0
        DO ir = 2, num_X
            DO ic = 1, ir-1
                h = h+1
                theta(il+h) = 0.d0
            END DO
        END DO
        !
    END IF
    IF (switch_Phi2 .EQ. 1) THEN        ! Phi lower triangular
        h = 0
        DO ir = 1, num_X
            DO ic = 1, ir
                ! Cycle if (ir,ic) corresponds to local factors
                ! Only works for models with (XG,1,1) factors
                IF ((ANY(ir .EQ. sel_X(num_XG+1:,num_C))) .AND. (ANY(ic .EQ. sel_X(num_XG+1:,num_C-1)))) CYCLE
                h = h+1
                IF (ir .EQ. ic) THEN
                    theta(il+h) = TAN(pisudue*(2.d0*Phi(ir,ir)-1.d0))
                ELSE IF (ir .GT. ic) THEN
                    theta(il+h) = Phi(ir,ic)
                END IF
            END DO
        END DO
    END IF
    iu = iu+num_Phi
    !
    ! Omega
    !
    il = iu
    theta(il+1:il+num_C) = LOG(1.d-2)
    iu = iu+num_Omega
    !
    ! omega_dex
    !
    il = iu
    theta(il+1) = LOG(1.d-2)
    iu = iu+num_omega_dex
    !
    ! r_inf
    !
    il = iu
    IF (switch_r_inf2 .EQ. 1) THEN
        theta(il+1) = 0.d0
    END IF
    iu = iu+num_r_inf_true
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE HWinival
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE starting_points