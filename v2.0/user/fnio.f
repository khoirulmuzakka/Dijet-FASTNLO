
C -------------------------------------------------------------
C   M. Wobisch   07/22/2007
C   i/o routines to read/write 
C   double precision, integer, long integer, character variables
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

c 1000 FORMAT (D24.17)           ! for coefficients
 1000 FORMAT (D13.6)           ! for coefficients
 1100 FORMAT (A) 

      Return
      End
C -------------------------------------------------------------
      Subroutine fnioint(crw,nunit,ivar)
      Implicit None
      Character*(*) crw
      Integer nunit, ivar, l,n
      Character*5 f(20)
      Data f/"(I1)","(I2)","(I3)","(I4)","(I5)","(I6)","(I7)","(I8)",
     +     "(I9)","(I10)","(I11)","(I12)","(I13)","(I14)","(I15)","(I16)",
     +     "(I17)","(I18)","(I19)","(I20)"/

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
     +     "(I9)","(I10)","(I11)","(I12)","(I13)","(I14)","(I15)","(I16)",
     +     "(I17)","(I18)","(I19)","(I20)"/

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

 1000 FORMAT (A) 
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
            Write(*,*) ' fastNLO: Iseparator in table not consistent'
            Stop
         Endif
      Else 
         Write(nunit,"(I10)") 1234567890
      Endif

      Return
      End
C -------------------------------------------------------------

