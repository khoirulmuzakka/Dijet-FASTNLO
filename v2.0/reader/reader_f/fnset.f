      SUBROUTINE FNSET(VAR,IVAL)
***********************************************************************
*     
*     fastNLO routines to select the contributions to be computed
*     in subsequent calls
*     
*     Input:
*     ------
*     VAR   Contribution name --> translated into integer of sim. name
*     IVAL  Selection flag (0: not selected, 1 selected)
*     
***********************************************************************
      IMPLICIT NONE
      CHARACTER*(*) VAR
      INTEGER IVAL
      INCLUDE "fnx9999.inc"
      DATA IFNFIRST/0/,
     >     PREFTAB/0/,
     >     PORDPTHY/0/, PTHRESHCOR/0/, 
     >     PQUARKCOMPOSITENESS/0/, PADDLED/0/, PTEVED/0/
     >     PNPCOR/0/, PUEVENT/0/, PDATA/0/
      
*---  Reset
      IF (VAR.EQ."P_RESET") THEN
         PREFTAB             = 0
         PORDPTHY            = 0
         PTHRESHCOR          = 0
         PNPCOR              = 0
         PQUARKCOMPOSITENESS = 0
         PADDLED             = 0
         PTEVED              = 0
         PDATA               = 0
      ENDIF

*---  Evaluate standard tables or reference tables
      IF (VAR.EQ."P_REFTAB") PREFTAB = IVAL

*---  Perturbative Contributions - Fixed Orders
      IF (VAR.EQ."P_ORDPTHY") PORDPTHY = IVAL
      
*---  Perturbative Contributions - Corrections
      IF (VAR.EQ."P_THRESHCOR") PTHRESHCOR = IVAL

*---  Multiplicative Corrections
      IF (VAR.EQ."P_NPCOR") PNPCOR = IVAL

*---  New Physics Contributions
      IF (VAR.EQ."P_QUARKCOMPOSITENESS") PQUARKCOMPOSITENESS = IVAL
      IF (VAR.EQ."P_ADDLED") PADDLED = IVAL
      IF (VAR.EQ."P_TEVED") PTEVED = IVAL

*---  Data
      IF (VAR.EQ."P_DATA") PDATA = IVAL

      RETURN
      END
