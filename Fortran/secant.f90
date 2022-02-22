! SOLUZIONE ALLA EQUAZIONE NON LINEARE
! f=@(x) x-log(2-x+x.^2) IN [0,1], zero=0.561
! f=@(x) x.*exp(3*x)-1-x IN [-2,-1], zero=-1.045
! f=@(x) log(2+sin(x))-x IN [1,2], zero=1.054
! f=@(x) 2*sin(x)+x-1 IN [0,1], zero=0.337
! f=@(x) (sin(x)).^2-x-1 IN [-1,0], zero=-0.641

PROGRAM secant
	IMPLICIT NONE
	INTEGER, PARAMETER :: imax = 100		! Parameter = costante
	INTEGER, PARAMETER :: tol = 1.e-8		! Parameter = costante
	INTEGER :: i					! Indice
	REAL :: a, b, c, fa, fb, fc, res		! Dichiaro le variabili
	
	i=0
	res=0
	a=0
	b=1
	c=0

	DO WHILE (res >= tol)
		fa=f_eval(a)
		fb=f_eval(b)
		c=b-(b-a)/(fb-fa)*fb;
		fc=f_eval(c)
		IF (fa*fb<0) THEN
			b=c
		ELSE
			a=c
		END IF
		res=f_eval(c)

		IF (res >= tol) THEN
			i=i+1
			WRITE(*,*) i, c
		ELSE
			EXIT
		END IF
		
	END DO
	
	
CONTAINS


	FUNCTION f_eval( x )
		! Input function
		REAL :: f_eval
		REAL :: x

		f_eval = 2*sin(x)+x-1
	END FUNCTION f_eval


END PROGRAM secant
