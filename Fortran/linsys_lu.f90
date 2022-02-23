!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!          LINEAR SYSTEM OF EQUATONS          !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!              LU FACTORIZATION             !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM LIN_SYS_LU

	IMPLICIT NONE
	INTEGER, PARAMETER :: N = 3		! Dimension of linear system
	REAL, DIMENSION(N,N) :: A=0, L=0,U=0
	REAL, DIMENSION(N) :: b, x=0, y=0

	WRITE(*,*) ""
	WRITE(*,*) "*********************************************************"
	WRITE(*,*) "****   LINEAR SYSTEM OF EQUATION: LU FACTORIZATION   ****"
	WRITE(*,*) "*********************************************************"
	WRITE(*,*) ""
	WRITE(*,'(A,I1,A,I1)') "   Input system dimension: ",N,"x",N

	! Definition of the linear system in a separated subroutine
	CALL read_system(A,b,N)

	! Performing LU factorization without pivoting
	U=A
	CALL lu(A,L,U,N)

	! Solution of the linear system Ax=b
	CALL solution(L,U,b,N,x,y)
	

	! Print results and problem constants:
	CALL print_result(A,L,U,b,x,N)

	WRITE(*,*) "-------------------------------------------------"
	


CONTAINS

! *******************************************************************************
	! Definition of the matrix coefficients
	SUBROUTINE read_system(A,b,N)
		INTEGER :: N		
		REAL, DIMENSION(N,N) :: A
		REAL, DIMENSION(N) :: b
		
		A(1,:) = (/ 1,1,1/)
		A(2,:) = (/ 4,3,-1/)
		A(3,:) = (/ 3,5,3/)

		b = (/ 1,6,4/)

	END SUBROUTINE read_system




! *******************************************************************************
	! Print matrix A
	SUBROUTINE printmatrix(A,N)
		INTEGER :: N,i,j
		REAL, DIMENSION(N,N) :: A

		DO i=1,N
			WRITE(*,'(X*(F10.3))') (A(i,j), j=1,N)
		END DO
		
	END SUBROUTINE printmatrix




! *******************************************************************************
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




! *******************************************************************************
	! Implementing solution of the linear system Ax = LUx = b
	! In other words Ux=y -> Ly=b
	SUBROUTINE solution(L,U,b,N,x,y)
		INTEGER :: i,j
		INTEGER :: N
		REAL, DIMENSION(N,N) :: L, U
		REAL, DIMENSION(N) :: b,x,y
		
		
		! Solve for lower-triangular by forward-substitution
		CALL tril(L,y,b,N)

		! Solve for upper-triangular by back-substitution
		CALL triu(U,x,y,N)


	END SUBROUTINE solution




! *******************************************************************************
	! Lower-triangular by forward-substitution algorithm
	SUBROUTINE tril(L,y,b,N)
		INTEGER :: i,j
		INTEGER :: N
		REAL, DIMENSION(N,N) :: L
		REAL, DIMENSION(N) :: b,y
		REAL :: sum1=0
		
		y(1)=b(1)/L(1,1)

		DO i=2,N
			sum1=0
			DO j=1,(i-1)
				sum1=sum1-L(i,j)*y(j)
			END Do
			y(i)=1/L(i,i)*(b(i)+sum1)
		END DO				
	END SUBROUTINE tril




! *******************************************************************************
	! Upper-triangular by backward-substitution algorithm
	SUBROUTINE triu(U,x,y,N)
		INTEGER :: i,j
		INTEGER :: N
		REAL, DIMENSION(N,N) :: U
		REAL, DIMENSION(N) :: x,y
		REAL :: sum2=0
		
		x(N)=y(N)/U(N,N)
		
		i=N-1
		DO WHILE(i.GE.1)
			sum2=0
			j=i+1
			DO WHILE (j.LE.N)
				sum2=sum2-U(i,j)*x(j)
				j=j+1
			END DO
			x(i)=1/U(i,i)*( y(i)+sum2 )
			i=i-1
		END DO

	END SUBROUTINE triu




! *******************************************************************************
	! Upper-triangular by backward-substitution algorithm
	SUBROUTINE print_result(A,L,U,b,x,N)
		INTEGER :: N,j
		REAL, DIMENSION(N,N) :: A,L,U
		REAL, DIMENSION(N) :: x,b
		
		! System defintion
		WRITE(*,*) "-------------------------------------------------"
		WRITE(*,*) "  A matrix coefficients:"
		CALL printmatrix(A,N)
		WRITE(*,*) "-------------------------------------------------"
		WRITE(*,*) "  Constant terms b coefficients:"
		WRITE(*,'(X*(F10.3))') (b(j), j=1,N)
		WRITE(*,*) "-------------------------------------------------"
		WRITE(*,*) "  Solution of the system is:"
		WRITE(*,'(X*(F10.3))') (x(j), j=1,N)
		WRITE(*,*) "-------------------------------------------------"


		! L coefficients
		WRITE(*,*) " "
		WRITE(*,*) " "
		WRITE(*,*) " "
		WRITE(*,*) "  L matrix coefficients:"
		CALL printmatrix(L,N)
		WRITE(*,*) "-------------------------------------------------"
		
		! U coefficients
		WRITE(*,*) ""
		WRITE(*,*) "  U matrix coefficients:"
		WRITE(*,*) ""
		CALL printmatrix(U,N)

	END SUBROUTINE print_result

END PROGRAM LIN_SYS_LU
