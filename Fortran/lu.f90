!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!          LINEAR SYSTEM OF EQUATONS          !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!              LU FACTORIZATION             !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM LU_FACT
	IMPLICIT NONE
	INTEGER, PARAMETER :: N = 4		! DIMENSIONE DEL SISTEMA LINEARE
	REAL, DIMENSION(N,N) :: A=0, L=0,U=0

	WRITE(*,*) ""
	WRITE(*,*) "***********************************************"
	WRITE(*,*) "********    LU MATRIX DECOMPOSITION    ********"
	WRITE(*,*) "***********************************************"
	WRITE(*,*) ""
	WRITE(*,'(A,I1,A,I1)') "   Matrix dimension:",N,"x",N

	! Definition of matrix in a separated subroutine
	CALL readmatrix(A,N)

	
	! A coefficients
	WRITE(*,*) ""
	WRITE(*,*) "  A matrix coefficients:"
	WRITE(*,*) ""
	CALL printmatrix(A,N)
	WRITE(*,*) "-------------------------------------------------"

	! Performing LU factorization without pivoting
	U=A
	CALL lu(A,L,U,N)
	
	! L coefficients
	WRITE(*,*) ""
	WRITE(*,*) "  L matrix coefficients:"
	WRITE(*,*) ""
	CALL printmatrix(L,N)
	WRITE(*,*) "-------------------------------------------------"
	
	! U coefficients
	WRITE(*,*) ""
	WRITE(*,*) "  U matrix coefficients:"
	WRITE(*,*) ""
	CALL printmatrix(U,N)
	WRITE(*,*) "-------------------------------------------------"
	


CONTAINS

! **********************************************************************************
	! Definition of the matrix coefficients
	SUBROUTINE readmatrix(A,N)
		INTEGER :: N		
		REAL, DIMENSION(N,N) :: A
		
		A(1,:) = (/ 1,-2,-4,-3/)
		A(2,:) = (/ 2,-7,-7,-6/)
		A(3,:) = (/ -1,2,6,4/)
		A(4,:) = (/ -4,-1,9,8/)

	END SUBROUTINE readmatrix




! **********************************************************************************
	! Print matrix A
	SUBROUTINE printmatrix(A,N)
		INTEGER :: N, i,j
		REAL, DIMENSION(N,N) :: A

		DO i=1,N
			WRITE(*,'(X*(F10.3))') (A(i,j), j=1,N)
		END DO
		WRITE(*,*) " "


	END SUBROUTINE printmatrix




! **********************************************************************************
	! Implementing LU factorization
	SUBROUTINE lu(A,L,U,N)
		INTEGER :: i,j
		INTEGER :: N
		REAL, DIMENSION(N,N) :: A, L, U
		DO i=1,N
			L(i,i)=1
		END DO


		DO i=1,N
			DO j=i+1,N
				L(j,i)=U(j,i)/U(i,i)
				U(j,:)=U(j,:)-L(j,i)*U(i,:)
			END DO
			
		END DO

	END SUBROUTINE lu	


END PROGRAM LU_FACT
