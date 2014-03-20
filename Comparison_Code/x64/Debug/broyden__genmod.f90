        !COMPILER-GENERATED INTERFACE MODULE: Mon Nov 12 14:10:22 2012
        MODULE BROYDEN__genmod
          INTERFACE 
            SUBROUTINE BROYDEN(F_EXO,FVEC,X,N,DF,INTLJ,REEVALJ,CHECK,   &
     &MAXSTP,TOLF)
              INTEGER(KIND=4) :: N
              EXTERNAL F_EXO
              REAL(KIND=8) :: FVEC(N)
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: DF(N,N)
              LOGICAL(KIND=4) :: INTLJ
              LOGICAL(KIND=4) :: REEVALJ
              LOGICAL(KIND=4) :: CHECK
              REAL(KIND=8) :: MAXSTP
              REAL(KIND=8) :: TOLF
            END SUBROUTINE BROYDEN
          END INTERFACE 
        END MODULE BROYDEN__genmod
