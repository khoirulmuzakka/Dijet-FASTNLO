***********************************************************************
*
*     I/O routines to read or write variables from/to fastNLO tables
*
*     Initial version: M. Wobisch, 2007
*     Updated for v2:  K. Rabbertz, 2011
*
*     The following are implemented:
*     ------------------------------
*     FNIODBL   double precision
*     FNIOINT   integer
*     FNIOLINT  long integer
*     FNIOCHAR  character variables (strings)
*
*     In addition:
*     ------------
*     FNIOISEP  reads or writes the "separator" (=1234567890)
*
*     All routines have the input arguments:
*     --------------------------------------
*     CRW    (character)  Can be 'read' or 'write'
*     NUIT   (integer)    Fortran unit to be used to read/write
*     DVAR,IVAR,CVAR      Variable to be read or written,
*                         type depends on the routine
*                         (not used in fniosep)
*     LPRINT (logical)    Print value read from table
*     DSCRPTN (character) Description of table value
*
***********************************************************************
      SUBROUTINE FNIODBL(CRW,NUNIT,DVAR,LPRINT,DSCRPTN)
      IMPLICIT NONE
      CHARACTER*(*) CRW, DSCRPTN
      CHARACTER*40 CHTMP
      INTEGER NUNIT
      DOUBLE PRECISION DVAR
      LOGICAL LPRINT

      IF (CRW.EQ.'read') THEN
         READ(NUNIT,*) DVAR
         IF (LPRINT) THEN
            READ(DSCRPTN,'(A)'),CHTMP
            WRITE(*,'(A,G10.4)'),CHTMP,DVAR
         ENDIF
      ELSE
         IF (DVAR.EQ.0D0) THEN
            WRITE(NUNIT,1100) '0'
         ELSE
            WRITE(NUNIT,1000) DVAR
         ENDIF
      ENDIF

 1000 FORMAT (D13.6)
 1100 FORMAT (A)

      RETURN
      END

***********************************************************************

      SUBROUTINE FNIOINT(CRW,NUNIT,IVAR,LPRINT,DSCRPTN)
      IMPLICIT NONE
      CHARACTER*(*) CRW, DSCRPTN
      CHARACTER*40 CHTMP
      INTEGER NUNIT, IVAR, L, N
      LOGICAL LPRINT
      CHARACTER*5 F(20)
      DATA F/"(I1)","(I2)","(I3)","(I4)","(I5)","(I6)","(I7)","(I8)",
     +     "(I9)","(I10)","(I11)","(I12)","(I13)","(I14)","(I15)",
     +     "(I16)","(I17)","(I18)","(I19)","(I20)"/

      IF (CRW.EQ.'read') THEN
         READ(NUNIT,*) IVAR
         IF (LPRINT) THEN
            READ(DSCRPTN,'(A)'),CHTMP
            WRITE(*,'(A,I10)'),CHTMP,IVAR
         ENDIF
      ELSE
         IF (IVAR.EQ.0) THEN
            L = 1
         ELSE
            L = INT(LOG10(DBLE(ABS(IVAR))))+1
         ENDIF
         IF (IVAR.LT.0) L=L+1
         WRITE(NUNIT,F(L)) IVAR
      ENDIF

      RETURN
      END

***********************************************************************

      SUBROUTINE FNIOLINT(CRW,NUNIT,IVAR,LPRINT,DSCRPTN)
      IMPLICIT NONE
      CHARACTER*(*) CRW, DSCRPTN
      CHARACTER*40 CHTMP
      INTEGER NUNIT, L, N
      INTEGER*8 IVAR
      LOGICAL LPRINT
      CHARACTER*5 F(20)
      DATA F/"(I1)","(I2)","(I3)","(I4)","(I5)","(I6)","(I7)","(I8)",
     +     "(I9)","(I10)","(I11)","(I12)","(I13)","(I14)","(I15)",
     +     "(I16)","(I17)","(I18)","(I19)","(I20)"/

      IF (CRW.EQ.'read') THEN
         READ(NUNIT,*) IVAR
         IF (LPRINT) THEN
            READ(DSCRPTN,'(A)'),CHTMP
            WRITE(*,'(A,I16)'),CHTMP,IVAR
         ENDIF
      ELSE
         IF (IVAR.EQ.0) THEN
            L = 1
         ELSE
            L = INT(LOG10(DBLE(ABS(IVAR))))+1
         ENDIF
         IF (IVAR.LT.0) L=L+1
         WRITE(NUNIT,F(L)) IVAR
      ENDIF

      RETURN
      END

***********************************************************************

      SUBROUTINE FNIOCHAR(CRW,NUNIT,CVAR,LPRINT,DSCRPTN)
      IMPLICIT NONE
      CHARACTER*(*) CRW,CVAR,DSCRPTN
      CHARACTER*40 CHTMP
      INTEGER NUNIT
      LOGICAL LPRINT

      IF (CRW.EQ.'read') THEN
         READ(NUNIT,'(A)') CVAR
         IF (LPRINT) THEN
            READ(DSCRPTN,'(A)'),CHTMP
            WRITE(*,'(A,A)'),CHTMP,CVAR(1:LEN_TRIM(CVAR))
         ENDIF
      ELSE
         WRITE(NUNIT,1000) CVAR(1:LEN_TRIM(CVAR))
      ENDIF

 1000 FORMAT (A)
      RETURN
      END

***********************************************************************

      SUBROUTINE FNIOISEP(CRW,NUNIT,LPRINT,DSCRPTN)
      IMPLICIT NONE
      CHARACTER*(*) CRW, DSCRPTN
      CHARACTER*40 CHTMP
      INTEGER NUNIT, ISEP
      LOGICAL LPRINT

      IF (CRW.EQ.'read') THEN
         READ(NUNIT,*) ISEP
         IF (LPRINT) THEN
            READ(DSCRPTN,'(A)'),CHTMP
            WRITE(*,'(A,I10)'),CHTMP,ISEP
         ENDIF
         IF (ISEP.NE.1234567890) THEN
            WRITE(*,*)'fastNLO: ERROR! '//
     +           'Expected separator line 1234567890 in table and not',
     +           ISEP
            STOP
         ENDIF
      ELSE
         WRITE(NUNIT,"(I10)") 1234567890
      ENDIF

      RETURN
      END
