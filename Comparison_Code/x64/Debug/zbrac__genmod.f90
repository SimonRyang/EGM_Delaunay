        !COMPILER-GENERATED INTERFACE MODULE: Tue Sep 04 09:13:47 2012
        MODULE ZBRAC__genmod
          INTERFACE 
            SUBROUTINE ZBRAC(FUNC,X1,X2,SUCCES)
              INTERFACE 
                FUNCTION FUNC(X)
                  REAL(KIND=4), INTENT(IN) :: X
                  REAL(KIND=4) :: FUNC
                END FUNCTION FUNC
              END INTERFACE 
              REAL(KIND=4), INTENT(INOUT) :: X1
              REAL(KIND=4), INTENT(INOUT) :: X2
              LOGICAL(KIND=4), INTENT(OUT) :: SUCCES
            END SUBROUTINE ZBRAC
          END INTERFACE 
        END MODULE ZBRAC__genmod
