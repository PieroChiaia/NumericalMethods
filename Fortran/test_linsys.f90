! TESTED LINEAR SYSTEM FOR PROGRAM "pivot_lu.f90"

! CASE 1
! =======================================================================
	! DEFINITION OF THE SYSTEM Ax=B
	SUBROUTINE read_system(A,b,N)
		INTEGER :: N		
		REAL, DIMENSION(N,N) :: A
		REAL, DIMENSION(N) :: b
		
		A(1,:) = (/ 9,3,2,0,7/)
		A(2,:) = (/ 7,6,9,6,4/)
		A(3,:) = (/ 2,7,7,8,2/)
		A(4,:) = (/ 0,9,7,2,2/)
		A(5,:) = (/ 7,3,6,4,3/)

		b = (/ 35,58,53,37,39/)

	END SUBROUTINE read_system
	


! CASE 2
! =======================================================================
	! DEFINITION OF THE SYSTEM Ax=B
	SUBROUTINE read_system(A,b,N)
		INTEGER :: N		
		REAL, DIMENSION(N,N) :: A
		REAL, DIMENSION(N) :: b
		
		A(1,:) = (/ 1,0,2,3/)
		A(2,:) = (/ -1,2,2,-3/)
		A(3,:) = (/ 0,1,1,4/)
		A(4,:) = (/ 6,2,2,4/)

		b = (/ 1,-1,2,1/)

	END SUBROUTINE read_system
	


! CASE 3
! =======================================================================
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
