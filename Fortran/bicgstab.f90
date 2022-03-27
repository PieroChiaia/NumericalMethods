        !=======================================================================================!
        !                          NUMERICAL METHODS FOR LINEAR SYSTEMS                         !
        !=======================================================================================!


        !-----------------------------------------------------!
        !          BiCGStab implemented in Fortran90          !
        !-----------------------------------------------------!


MODULE BiCGStab_mod
    implicit none

    INTEGER(kind=2 ), PARAMETER :: sp = kind(1.000)
    INTEGER(kind=2 ), PARAMETER :: dp = kind(1.0d0)
    INTEGER(kind=2 ), PARAMETER :: qp = kind(1.0q0)
    INTEGER(kind=2 ), PARAMETER :: wp = dp

CONTAINS

    SUBROUTINE BiCGStab(A,b,x)
        implicit none
        !--------------------------PARAMETER AND VARIABLE-------------------------------!
        REAL(kind=wp), INTENT(in)                    :: A (:,:)
        REAL(kind=wp), INTENT(in)                    :: b ( : )
        REAL(kind=wp), DIMENSION(1:size(b, dim=1))   :: x, x0

        REAL(kind=wp), DIMENSION(1:size(b, dim=1))   :: r, rp1, rstart
        REAL(kind=wp), DIMENSION(1:size(b, dim=1))   :: p,s,res,v,h,t,rf

        REAL(kind=wp), parameter                     :: eps = 1.e-8
        REAL(kind=wp)                                :: alpha, omega, beta, norm_res
        INTEGER                                      :: k = 0


        alpha=1.0_wp
        omega=1.0_wp
        beta=1.0_wp

        x0=0.0_wp
        

       !First guess: r0* can be chosen arbitrarly
        rstart = b-matmul(A,x0)
        rf=rstart

        r = rstart
        p = r

        s = 0.0_wp

        norm_res=sqrt(dot_product(rstart,rstart))


        DO WHILE (norm_res .GE. eps)
           !BiCGStab standard algorrithm without preconditioner

            alpha = dot_product(r,rf)/(dot_product(matmul(A,p),rf))
                s = r - alpha*matmul(A,p)
            omega = dot_product(matmul(A,s),s)/dot_product(matmul(A,s),matmul(A,s))
                x = x+alpha*p+omega*s
              rp1 = s-omega*matmul(A,s)
             beta = alpha/omega*dot_product(rp1,rf)/dot_product(r,rf)
                p = rp1+beta*(p-omega*matmul(A,p))
                r = rp1

           !Residual of the iteration
            norm_res=sqrt(dot_product(r,r))
            k=k+1

            WRITE(*,'(A,I2)') "Iteration: ", k
            WRITE(*,'(A,E10.4)') "   Res: ", norm_res


        END DO

    END SUBROUTINE BiCGStab   

END MODULE BiCGStab_mod



!==============================================================================================



MODULE Matrices_mod
    USE BiCGStab_mod
    implicit none

    CONTAINS
    
    ! DEFINITION OF THE SYSTEM Ax=B
    SUBROUTINE read_system(A,b)   
        REAL(kind=wp), INTENT(out)   ::   A (:,:)
        REAL(kind=wp), INTENT(out)   ::   b ( : )
        
        A(1,:) = (/ 14,14,9,3,5/)
        A(2,:) = (/ 14,52,-15,2,-32/)
        A(3,:) = (/ 9,15,36,5,16/)
        A(4,:) = (/ 3,2,5,47,49/)
        A(5,:) = (/ 5,32,16,49,79/)

        b = (/ -15,-100,106,329,463/)

    END SUBROUTINE read_system


    ! CONSOLE PRINT MATRIX
    SUBROUTINE printmatrix(A,N)
        INTEGER :: N,i,j
        REAL (kind=wp), DIMENSION(N,N) :: A
        WRITE(*,*) " "
        DO i=1,N
            WRITE(*,'(X*(f10.4))') (A(i,j), j=1,N)
        END DO
        
    END SUBROUTINE printmatrix

END MODULE Matrices_mod



!==============================================================================================



PROGRAM main
use BiCGStab_mod
use Matrices_mod

    implicit none

    INTEGER, parameter             :: m=5, n=5
    REAL (kind=wp), dimension(1:m,1:n)    :: A
    REAL (kind=wp), dimension(1:m    )    :: x_calculated,b
    INTEGER :: j

!-----------------------------A,b DEFINITION--------------------------------!
    CALL read_system(A,b)


!--------------------------------SOLUTION-----------------------------------!
    WRITE(*,*) "BiCGStab without preconditioner"

    CALL BiCGStab(A,b,x_calculated)


    WRITE(*,*) "-----------------------------------------"
    WRITE(*,*) "BiCGStab solution   ="
    WRITE(*,'(X*(F10.3))') (x_calculated(j), j=1,N)


END PROGRAM main
