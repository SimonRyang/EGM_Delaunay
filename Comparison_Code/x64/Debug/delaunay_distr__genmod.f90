        !COMPILER-GENERATED INTERFACE MODULE: Thu Jul 04 13:55:57 2013
        MODULE DELAUNAY_DISTR__genmod
          INTERFACE 
            SUBROUTINE DELAUNAY_DISTR(ADASH,HDASH,WEIGHT,J,K)
              REAL(KIND=8), INTENT(IN) :: ADASH
              REAL(KIND=8), INTENT(IN) :: HDASH
              REAL(KIND=8), INTENT(OUT) :: WEIGHT(3)
              INTEGER(KIND=4), INTENT(OUT) :: J(3)
              INTEGER(KIND=4), INTENT(OUT) :: K(3)
            END SUBROUTINE DELAUNAY_DISTR
          END INTERFACE 
        END MODULE DELAUNAY_DISTR__genmod