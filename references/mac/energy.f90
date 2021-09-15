SUBROUTINE bw98_stretch(val,X,invUV,area) 
IMPLICIT NONE 
REAL(KIND=8) val(1,1) 
val(1,1) = 1.0E+0*area(1,1)*(sqrt((invUV(2,1)*(X(3,3)+tt3)+invUV(&
&1,1)*(X(3,2)+tt3))**2+(invUV(2,1)*(X(2,3)+tt2)+invUV(1,1)*(X(2,2)&
&+tt2))**2+((X(1,3)+tt1)*invUV(2,1)+invUV(1,1)*(X(1,2)+tt1))**2)-1&
&.0E+0)**2
END 
