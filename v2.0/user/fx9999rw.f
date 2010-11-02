      Subroutine fx9999rw(crw,filename,ictrb)
* ------------------------------------------------------------
*  fastNLO usercode v2.0                       MW 07/23/2007
* 
*  fx9999rw reads or writes fastNLO v2 tables
*  input:
*     character   crw:        'read' or 'write'
*     character   filename:    name of table
*     integer     ictrb(30):   pointers for writing selected contributions
*                              only contributions n with ictrb(n)<>0 are 
*                              written - reading is not affected
C
C        1         2         3         4         5         6         7 
C 3456789012345678901234567890123456789012345678901234567890123456789012
* ------------------------------------------------------------
      Implicit None
      Character*(*) crw,filename
      Integer Nunit, Ifile,ic,i,j,k,l,n,m,nxmax, nselctrb
      Integer ictrb(30)
      Include 'fnx9999.inc'

      Nunit=2

      If (crw.eq.'write') Then
         Open(unit=Nunit,file=filename,status='unknown')
      Else
         OPEN(Nunit,STATUS='OLD',FILE=FILENAME,IOSTAT=IFILE)
         IF (Ifile .ne. 0) THEN
            WRITE(*,*)"FX9999RW: ERROR! Table file not found, "//
     +           "stopped! IOSTAT = ",Ifile
            STOP
         Endif
      Endif

c --- when writing: count No. of selected contributions to be written
c                   (those with ictrb(ic)=1)
      nselctrb = 0
      If (crw .eq. 'write') Then         
         Do ic=1,NContrib
            If (ictrb(ic).ne.0) Then
               nselctrb = nselctrb + 1
               Write(*,*) '         ... writing contribution ',ic
            Endif
         Enddo
         Write(*,*) ' writing a total of ',nselctrb,' contributions',
     +        ' into: ',filename
      Endif


c --- fastNLO table
c - block A1
      Call fnioisep(crw,nunit)
      Call fnioint(crw,nunit, Itabversion)
      Call fniochar(crw,nunit, ScenName)
      If (crw .eq. 'write') Then
         Call fnioint(crw,nunit, nselctrb)
      Else
         Call fnioint(crw,nunit, Ncontrib)
      Endif
      Call fnioint(crw,nunit, Nmult)
      Call fnioint(crw,nunit, Ndata)
      Call fnioint(crw,nunit, NuserString)
      If (NuserString .gt. MxUser) Then
         Write(*,*) ' NuserString too large ',NuserString,'<=',MxUser
      Endif
      Do i=1,NuserString
        Call fniochar(crw,nunit, UserString(i))
      Enddo
      Call fnioint(crw,nunit, NuserInt)
      If (NuserInt .gt. MxUser) Then
         Write(*,*) ' NuserInt too large ',NuserInt,'<=',MxUser
      Endif
      Do i=1,NuserInt
        Call fnioint(crw,nunit, UserInt(i))
      Enddo
      Call fnioint(crw,nunit, NuserFloat)
      If (NuserFloat .gt. MxUser) Then
         Write(*,*) ' NuserFloat too large ',NuserFloat,'<=',MxUser
      Endif
      Do i=1,NuserFloat
        Call fniodbl(crw,nunit, UserFloat(i))
      Enddo
      Call fnioint(crw,nunit, Imachine) 
      If (Imachine.ne.0) Then
         Write(*,*) ' Imachine different from zero - ',
     +        'not yet implemented: ',Imachine
         Stop
      Endif
c - block A2
      Call fnioisep(crw,nunit)
      Call fnioint(crw,nunit, IpublUnits) 
      Call fnioint(crw,nunit, NscDescript)
      If (NscDescript .gt. MxScDescript) then
         write(*,*) ' NscDescript too large ',NscDescript,'<=',MXScDescript
      Endif
      Do i=1,NscDescript
        Call fniochar(crw,nunit, ScDescript(i))
      Enddo
      Call fniodbl(crw,nunit, Ecms)
      Call fnioint(crw,nunit, ILOord)
      Call fnioint(crw,nunit, NobsBin)
      Call fnioint(crw,nunit, NDim)
      Do i=1,NDim
         Call fniochar(crw,nunit, DimLabel(i))
      Enddo
      Do i=1,NDim
         Call fnioint(crw,nunit, IDiffBin(i))
      Enddo
      Do i=1,NObsBin
         Do j=1,NDim
            Call fniodbl(crw,nunit, LoBin(i,j))
            If (IDiffBin(j).eq.2) Call fniodbl(crw,nunit, UpBin(i,j))
         Enddo
      Enddo
      Do i=1,NObsBin
         Call fniodbl(crw,nunit, BinSize(i))
      Enddo
      Call fnioint(crw,nunit, INormFlag)
      If (INormFlag.gt.1) Call fniochar(crw,nunit, DenomTable)
      If (INormFlag.gt.0) Then
         Do i=1,NObsBin
            Call fnioint(crw,nunit, IDivLoPointer(i))
            Call fnioint(crw,nunit, IDivUpPointer(i))
         Enddo
      Endif

c - block B
      Do ic=1,NContrib
c -      write only selected contributions - those with ictrb(ic)=1         
         If (crw .eq. 'write' .and. ictrb(ic).eq.0) Goto 100

         Call fnioisep(crw,nunit)
         Call fnioint(crw,nunit, IXsectUnits(ic))
         Call fnioint(crw,nunit, IDataFlag(ic))
         Call fnioint(crw,nunit, IAddMultFlag(ic))
         Call fnioint(crw,nunit, IContrFlag1(ic))
         Call fnioint(crw,nunit, IContrFlag2(ic))
         Call fnioint(crw,nunit, IContrFlag3(ic))
         Call fnioint(crw,nunit, NContrDescr(ic))
         Do i=1,NContrDescr(ic)
            Call fniochar(crw,nunit, CtrbDescript(ic,i))
         Enddo
         Call fnioint(crw,nunit, NCodeDescr(ic))
         Do i=1,NCodeDescr(ic)
            Call fniochar(crw,nunit, CodeDescript(ic,i))
         Enddo

c --------------------------- Idata ?????
         If (IDataFlag(ic).eq.1) Then
            Write(*,*) "   Data Blocks can not yet be read"
            STOP
            Goto 100
         Endif

c --------------------------- IAddMult ?????
         If (IAddMultFlag(ic).eq.1) Then
            Write(*,*) "   Multiplicative Blocks can not yet be read"
            STOP
            Goto 100
         Endif

c --- coefficient block
         Call fnioint(crw,nunit, IRef(ic))
         Call fnioint(crw,nunit, IScaleDep(ic))
         Call fniolint(crw,nunit, Nevt(ic))
         Call fnioint(crw,nunit, Npow(ic))
         Call fnioint(crw,nunit, NPDF(ic))
         Do i=1,NPDF(ic)
            Call fnioint(crw,nunit, NPDFPDG(ic,i))
         Enddo
         Call fnioint(crw,nunit, NPDFDim(ic))
         Call fnioint(crw,nunit, NFragFunc(ic))
         Do i=1,NFragFunc(ic)
            Call fnioint(crw,nunit, NFFPDG(ic,i))
         Enddo
         Call fnioint(crw,nunit, NFFDim(ic))
         Call fnioint(crw,nunit, NSubproc(ic))
         Call fnioint(crw,nunit, IPDFdef(ic,1))   
         Call fnioint(crw,nunit, IPDFdef(ic,2))  
         Call fnioint(crw,nunit, IPDFdef(ic,3))  

c --- check consistency between Nsubproc and IPDFdef1,2,
         If (IPDFdef(ic,1).eq.0) Then ! - no predefined set of PDF coefficients
            Continue
         Elseif (IPDFdef(ic,1).eq.1) Then ! --- e+e-
            Continue
         Elseif (IPDFdef(ic,1).eq.2) Then ! --- DIS
            Continue
         Elseif (IPDFdef(ic,1).eq.3) Then ! --- hh/hhbar
            If (IPDFdef(ic,2).eq.1) Then 
               If (IPDFdef(ic,3).eq.1) Then
                  If (Nsubproc(ic).ne.6) Then
                     Write(*,*) " IPDFdef(3)=1 requires ",
     +                    "6 subprocesses for hh"
                     STOP
                  Endif
               Elseif (IPDFdef(ic,3).eq.2) Then 
                  If (Nsubproc(ic).ne.7) Then
                     Write(*,*) " IPDFdef(3)=2 requires ",
     +                    "7 subprocesses for hh"
                     STOP
                  Endif
               Endif
            Else
               Write(*,*) " case IPDFdef(2)<>1 n/a for hh/hhbar"
               STOP
            Endif

         Endif

         If (IPDFdef(ic,1).eq.0) Then ! - no predefined set of PDF coefficients
            Write(*,*) " case IPDFdef(1)=0 not yet implemented"
            STOP
         Endif
         If (NPDF(ic).gt.0) Then
            Do i=1,NObsBin
               Call fnioint(crw,nunit, Nxtot(ic,1,i))
               Do j=1,Nxtot(ic,1,i)
                  Call fniodbl(crw,nunit, XNode1(ic,i,j))
               Enddo
               If (crw.eq.'read') Then
                  Hxlim1(ic,i) = -sqrt(-log10(XNode1(ic,i,1))) ! for HH
               Endif

            Enddo 
            If (NPDFDim(ic).eq.2) Then
               Do i=1,NObsBin
                  Call fnioint(crw,nunit, Nxtot(ic,2,i))
                  Do j=1,Nxtot(ic,2,i)
                     Call fniodbl(crw,nunit, XNode2(ic,i,j))
                  Enddo
               Enddo 
            Endif
         Endif
         IF (NFragFunc(ic).gt.0) then ! - no FFs so far
            write(*,*) " fastNLO: no FragFuncs so far"
            STOP
         Endif

         Call fnioint(crw,nunit, NScales(ic))
         Call fnioint(crw,nunit, NScaleDim(ic))
         Do i=1,NScales(ic)
            Call fnioint(crw,nunit, IScale(ic,i))
         Enddo
         Do i=1,NScaleDim(ic)
            Call fnioint(crw,nunit, NScaleDescript(ic,i))
            Do j=1,NScaleDescript(ic,i)
               Call fniochar(crw,nunit, ScaleDescript(ic,i,j))
            Enddo
         Enddo
         Do i=1,NScaleDim(ic)
            Call fnioint(crw,nunit, NScaleVar(ic,i))
            Call fnioint(crw,nunit, NScaleNode(ic,i))
         Enddo
c         Write(*,*) ic,Nscales(ic),Nscaledim(ic)
c         Write(*,*) ' ',Nscalevar(ic,1),NScaleNode(ic,1)
         Do i=1,NScaleDim(ic)
            Do j=1,NScaleVar(ic,i)
               Call fniodbl(crw,nunit, ScaleFac(ic,i,j))
            Enddo
         Enddo

         Do i=1,NObsBin
            Do j=1,NScaleDim(ic)
               Do k=1,NScaleVar(ic,j)
                  Do l=1,NScaleNode(ic,j)
                     Call fniodbl(crw,nunit, ScaleNode(ic,i,j,k,l))
                     If (crw.eq.'read') Then
                        HScaleNode(ic,i,j,k,l) = 
     +                       log(log(ScaleNode(ic,i,j,k,l)/0.25d0))
                     Endif
                  Enddo
               Enddo
            Enddo
         Enddo
         
         Do i=1,NObsBin
            Do k=1,NScaleVar(ic,1)
               Do l=1,NScaleNode(ic,1)
c --- here we assume NFragFunc=0
                  If (NFragFunc(ic).gt.0) then
                     write(*,*) " NFragFunc>0 not yet implemented"
                     STOP
                  Endif
                  If (NPDFdim(ic).eq.0) Then
                     nxmax = Nxtot(ic,1,i)
                  Elseif (NPDFdim(ic).eq.1) Then
                     nxmax = (Nxtot(ic,1,i)**2+Nxtot(ic,1,i))/2
                  Elseif (NPDFdim(ic).eq.2) Then 
                     nxmax = 
     +                    Nxtot(ic,1,i)*Nxtot(ic,2,i)
                     write(*,*) '  NPDFdim = 2 not yet enabled'
                     Stop
                  Else
                     write(*,*) '  NPDFdim > 2 not enabled'
                     Stop
                  Endif
                  Do m=1,nxmax
                     Do n=1,NSubProc(ic)
                        Call fniodbl(crw,nunit,
     +                       SigmaTilde(ic,i,1,k,l,m,n))
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo
         
 100     Continue               ! - end of block
      Enddo

c - end of table
      Call fnioisep(crw,nunit)
      Call fnioisep(crw,nunit)
      Close(2)

      Return
      End
C -------------------------------------------------------------
