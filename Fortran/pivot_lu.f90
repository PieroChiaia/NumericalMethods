! =================================================
! =========   LINEAR SYSTEM OF EQUATION   =========
! =========     LU FACT WITH PIVOTING     =========
! =================================================

PROGRAM LIN_SYS_LU

	IMPLICIT NONE
	! DEFINE THE DIMENSION OF THE SYSTEM
	INTEGER, PARAMETER :: N = 5		
	REAL, DIMENSION(N,N) :: A=0, L=0,U=0
	INTEGER, DIMENSION(N,N) :: P=0
	REAL, DIMENSION(N) :: b, x=0, y=0
	
	CALL SYSTEM('clear')
	WRITE(*,*) ""
	WRITE(*,*) "======================================================="
	WRITE(*,*) "===========    LINEAR SYSTEM OF EQUATION    ==========="
	WRITE(*,*) "===========      LU FACT WITH PIVOTING      ==========="
	WRITE(*,*) "======================================================="
	WRITE(*,*) ""
	WRITE(*,'(A,I1,A,I1)') "      Input system dimension: ",N,"x",N
	
	
	
	! DEFINITION OF LINER SYSTEM EQUATIONS BY MATRIX AND B COEFFICIENTS
	CALL read_system(A,b,N)

	! PERFORM LU FACTORIZATION WITH ROW PIVOTING
	U=A
	CALL lu(L,U,P,N,b)

	! SOLUTION OF PAx=Pb
	CALL solution(L,U,b,N,x,y,P)
	
	! SCREEN PRINT
	CALL print_result(A,L,U,b,x,N,P)

	WRITE(*,*) " "
	WRITE(*,*) "======================================================="
	


CONTAINS



! ======================================================================================
	! DEFINITION OF THE SYSTEM Ax=B
	SUBROUTINE read_system(A,b,N)
		INTEGER :: N		
		REAL, DIMENSION(N,N) :: A
		REAL, DIMENSION(N) :: b
		
		A(1,:) = (/ 14,14,9,3,5/)
		A(2,:) = (/ 14,52,-15,2,-32/)
		A(3,:) = (/ 9,15,36,5,16/)
		A(4,:) = (/ 3,2,5,47,49/)
		A(5,:) = (/ 5,32,16,49,79/)

		b = (/ -15,-100,106,329,463/)

	END SUBROUTINE read_system


! ======================================================================================
	! CONSOLE PRINT MATRIX
	SUBROUTINE printmatrix(A,N)
		INTEGER :: N,i,j
		REAL, DIMENSION(N,N) :: A
		WRITE(*,*) " "
		DO i=1,N
			WRITE(*,'(X*(f10.4))') (A(i,j), j=1,N)
		END DO
		
	END SUBROUTINE printmatrix


! ======================================================================================
	! CONSOLE PRINT MATRIX
	SUBROUTINE printpermatrix(A,N)
		INTEGER :: N,i,j
		INTEGER, DIMENSION(N,N) :: A
		WRITE(*,*) " "
		DO i=1,N
			WRITE(*,'(X*(I5))') (A(i,j), j=1,N)
		END DO
		
	END SUBROUTINE printpermatrix


! ======================================================================================
	! LU FACTORIZATION WITH PARTIAL PIVOTING BY RAWS
	SUBROUTINE lu(L,U,P,N,b)
		INTEGER :: i,j,k,r
		INTEGER :: N
		REAL, DIMENSION(N,N) :: L, U
		INTEGER, DIMENSION(N,N) :: P
		REAL, DIMENSION(N) :: b, temp
		REAL :: m=0
		
		! L AND P DIAGONALS
		DO i=1,N
			L(i,i)=1
			P(i,i)=1
		END DO
		
		! IMPLEMENTING LU
		DO k=1,(N-1)

			! FIND MAXIMUM ON COLOUMN I
			r=0
			m=U(k,k)
			DO j=k,N
				IF (abs(U(j,k)).GE.m) THEN
					m=abs(U(j,k))
					r=j
				END IF
			END DO

			! EXCHANGES RAWS OF P, U, L, b
			temp=U(k,:)
			U(k,:)=U(r,:)
			U(r,:)=temp
			
			temp=0

			temp=P(k,:)
			P(k,:)=P(r,:)
			P(r,:)=temp
			
			temp=0

			temp=L(k,:)
			L(k,:)=L(r,:)
			L(r,:)=temp
			L(k,k)=1
			
			
			! COMPUTES L AND U COEFFICIENTS AT NEW STEP
			DO i=(k+1),N
				L(i,k)=U(i,k)/U(k,k)
				U(i,:)=U(i,:)-L(i,k)*U(k,:)
			END DO
			
		END DO

	END SUBROUTINE lu	


! ======================================================================================
	! SOLUTION OF PAx = LUx = Pb
	SUBROUTINE solution(L,U,b,N,x,y,P)
		INTEGER :: i,j
		INTEGER :: N
		REAL, DIMENSION(N,N) :: L, U
		INTEGER, DIMENSION(N,N) :: P
		REAL, DIMENSION(N) :: b,x,y, Pb
		REAL :: temp
		
		DO i=1,N
			temp=0
			Pb(i)=0
			DO j=1,N
				Pb(i)=Pb(i)+P(i,j)*b(j)
			END DO
		END DO
		b=Pb
		
		! LOWER TRIANG BY FORWARD SUBSTITUTION
		CALL tril(L,y,b,N)

		! UPPER TRIANG BY BACK SUBSTITUTION
		CALL triu(U,x,y,N)


	END SUBROUTINE solution


! ======================================================================================
	! LOWER-TRIANG BY FORWARD SUBSTITUTIONS
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
			y(i)=b(i)+sum1
		END DO				
	END SUBROUTINE tril


! ======================================================================================
	! UPPER-TRIANG BY BACKWARD SUBSTITUTIONS
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


! ======================================================================================
	! PRINT TO SCREEN
	SUBROUTINE print_result(A,L,U,b,x,N,P)
		INTEGER :: N,j
		REAL, DIMENSION(N,N) :: A,L,U
		INTEGER, DIMENSION(N,N) :: P
		REAL, DIMENSION(N) :: x,b
		
		! System defintion
		WRITE(*,*) "-------------------------------------------------------"
		WRITE(*,*) " "
		WRITE(*,*) "  A matrix coefficients:"
		CALL printmatrix(A,N)
		WRITE(*,*) "-------------------------------------------------------"
		WRITE(*,*) "  Constant terms b coefficients (after pivoting):"
		WRITE(*,'(X*(F10.3))') (b(j), j=1,N)
		WRITE(*,*) "-------------------------------------------------------"
		WRITE(*,*) "  Solution of the system is:"
		WRITE(*,'(X*(F10.3))') (x(j), j=1,N)
		WRITE(*,*) "-------------------------------------------------------"


		! L coefficients
		WRITE(*,*) " "
		WRITE(*,*) "  L matrix coefficients:"
		CALL printmatrix(L,N)
		WRITE(*,*) "-------------------------------------------------------"
		! U coefficients
		WRITE(*,*) "  U matrix coefficients:"
		CALL printmatrix(U,N)
		WRITE(*,*) "-------------------------------------------------------"
		! P coefficients
		WRITE(*,*) "  P matrix coefficients:"
		CALL printpermatrix(P,N)

	END SUBROUTINE print_result

END PROGRAM LIN_SYS_LU
