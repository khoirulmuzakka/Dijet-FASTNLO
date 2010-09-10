C --------------------------------------------------------------
C   M. Wobisch   07/22/2007
C   i/o routines to read or write variables from/to fastNLO tables 
C
C   the following are implemented
C     fniodbl   double precision
C     fnioiny   integer
C     fniolint  long integer
C     fniochar  character variables (strings)
C
C   in addition:
C     fnioisep  reads or writes the "separator" (=1234567890)
C
C   all routines have the input arguments
C     crw   (character)   can be 'read' or 'write'
C     nuit  (integer)     Fortran unit to be used to read/write
C     dvar,ivar,cvar      variable to be read or written, 
C                         type depends on the routine
C                         (not used in fniosep)
C
C        1         2         3         4         5         6         7 
C 3456789012345678901234567890123456789012345678901234567890123456789012
C -------------------------------------------------------------
      Subroutine fniodbl(crw,nunit,dvar)
      Implicit None
      Character*(*) crw
      Integer nunit
      Double Precision dvar

      If (crw.eq.'read') Then
         Read(nunit,*) dvar
      Else 
         If (dvar.eq.0d0) Then
            Write(nunit,1100) '0'
         Else
            Write(nunit,1000) dvar
         Endif
      Endif

c     1000 Format (D24.17)           ! for coefficients
 1000 Format (D13.6)            ! for coefficients
 1100 Format (A) 

      Return
      End
C -------------------------------------------------------------
      Subroutine fnioint(crw,nunit,ivar)
      Implicit None
      Character*(*) crw
      Integer nunit, ivar, l,n
      Character*5 f(20)
      Data f/"(I1)","(I2)","(I3)","(I4)","(I5)","(I6)","(I7)","(I8)",
     +     "(I9)","(I10)","(I11)","(I12)","(I13)","(I14)","(I15)",
     +     "(I16)","(I17)","(I18)","(I19)","(I20)"/

      If (crw.eq.'read') Then
         Read(nunit,*) ivar
      Else 
         If (ivar.eq.0) Then
            l = 1
         Else
            l = int(log10(dble(abs(ivar))))+1
         Endif
         If (ivar.lt.0) l=l+1
         Write(nunit,f(l)) ivar
      Endif

      Return
      End
C -------------------------------------------------------------
      Subroutine fniolint(crw,nunit,ivar)
      Implicit None
      Character*(*) crw
      Integer nunit, l,n
      Integer*8 ivar
      Character*5 f(20)
      Data f/"(I1)","(I2)","(I3)","(I4)","(I5)","(I6)","(I7)","(I8)",
     +     "(I9)","(I10)","(I11)","(I12)","(I13)","(I14)","(I15)",
     +     "(I16)","(I17)","(I18)","(I19)","(I20)"/

      If (crw.eq.'read') Then
         Read(nunit,*) ivar
      Else 
         If (ivar.eq.0) Then
            l = 1
         Else
            l = int(log10(dble(abs(ivar))))+1
         Endif
         If (ivar.lt.0) l=l+1
         Write(nunit,f(l)) ivar
      Endif

      Return
      End
C -------------------------------------------------------------
      Subroutine fniochar(crw,nunit,cvar)
      Implicit None
      Character*(*) crw,cvar
      Integer nunit
      
      If (crw.eq.'read') Then
         Read(nunit,*) cvar
      Else 
         Write(nunit,1000) cvar
      Endif

 1000 Format (A) 
      Return
      End
C -------------------------------------------------------------
      Subroutine fnioisep(crw,nunit)
      Implicit None
      Character*(*) crw
      Integer nunit, isep
      

      If (crw.eq.'read') Then
         Read(nunit,*) isep
         If (isep.ne.1234567890) Then
            Write(*,*) 'fastNLO: ERROR! '//
     +           'Expected separator line 1234567890 in table and not',
     +           isep
            Stop
         Endif
      Else 
         Write(nunit,"(I10)") 1234567890
      Endif

      Return
      End
C -------------------------------------------------------------
