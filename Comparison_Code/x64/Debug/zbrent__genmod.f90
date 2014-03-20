        !COMPILER-GENERATED INTERFACE MODULE: Tue Sep 04 09:13:45 2012
        MODULE ZBRENT__genmod
          INTERFACE 
            FUNCTION ZBRENT(FUNC,X1,X2,TOL)
              INTERFACE 
                FUNCTION FUNC(X)
                  REAL(KIND=4), INTENT(IN) :: X
                  REAL(KIND=4) :: FUNC
                END FUNCTION FUNC
              END INTERFACE 
              REAL(KIND=4), INTENT(IN) :: X1
              REAL(KIND=4), INTENT(IN) :: X2
              REAL(KIND=4), INTENT(IN) :: TOL
              REAL(KIND=4) :: ZBRENT
            END FUNCTION ZBRENT
          END INTERFACE 
        END MODULE ZBRENT__genmod
