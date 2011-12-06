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
Comment:       Integer IFNfirst, Preftab,
Comment:      +     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
Comment:      +     PNPCOR, PUEVENT, PDATA
Comment:       Common /cfastnlo/ IFNfirst, Preftab,
Comment:      +     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
Comment:      +     PNPCOR, PUEVENT, PDATA
      Include "fnx9999.inc"
      Data IFNfirst/0/,
     +     Preftab/0/,
     +     PORDPTHY/0/, PTHRESHCOR/0/, 
     +     PQUARKCOMPOSITENESS/0/, PADDLED/0/, PTEVED/0/
     +     PNPCOR/0/, PUEVENT/0/, PDATA/0/
      
c --- evaluate 'standard' tables or reference tables
      If (var.eq."P_REFTAB") PREFTAB = ival

c --- Perturbative Contributions - Fixed Orders
      If (var.eq."P_ORDPTHY") PORDPTHY = ival

c --- Perturbative Contributions - Corrections
      If (var.eq."P_THRESHCOR") PTHRESHCOR = ival

c --- Multiplicative Corrections
      If (var.eq."P_NPCOR") PNPCOR = ival

c --- New Physics Contributions
      If (var.eq."P_QUARKCOMPOSITENESS") PQUARKCOMPOSITENESS = ival
      If (var.eq."P_ADDLED") PADDLED = ival
      If (var.eq."P_TEVED") PTEVED = ival

c --- Data
      If (var.eq."P_DATA") PDATA = ival

      Return
      End
**********************************************************************
