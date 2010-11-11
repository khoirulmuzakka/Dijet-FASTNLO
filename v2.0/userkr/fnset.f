
**********************************************************************
* fastNLO routines to select the contributions to be computed
* in subsequent calls
*
* notice: any changes in the commonblock need to be copied into
*         the include file fnx9999.inc 
*
*
* question: how can we later determine if the usercode was called 
*           earlier with exactly the same settings?
* possible solution: increment a variable with each call (up to 37k)
*                    if identical in a subsequent call, one can assume 
*                    that FNSET was not called in between
*                                  ? ? ? 
* MW 06/11/2007
**********************************************************************
      Subroutine fnset(var,ival)
      Implicit None
      Character*(*) var
      Integer ival
      Integer IFNfirst,
     +     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
     +     PHADRCOR, PUEVENT
      Common /cfastnlo/ IFNfirst,
     +     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
     +     PHADRCOR, PUEVENT
      Data IFNfirst/0/,
     +     PORDPTHY/0/, PTHRESHCOR/0/, 
     +     PQUARKCOMPOSITENESS/0/, PADDLED/0/, PTEVED/0/
      
c --- Perturbative Contributions - Fixed Orders
      If (var.eq."P_ORDPTHY") PORDPTHY = ival

c --- Perturbative Contributions - Corrections
      If (var.eq."P_THRESHCOR") PTHRESHCOR = ival

c --- New Physics Contributions
      If (var.eq."P_QUARKCOMPOSITENESS") PQUARKCOMPOSITENESS = ival
      If (var.eq."P_ADDLED") PADDLED = ival
      If (var.eq."P_TEVED") PTEVED = ival

      Return
      End
**********************************************************************
