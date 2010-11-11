      Subroutine fx9999rw(crw,filename)
* ------------------------------------------------------------
*  fastNLO usercode v2.0 reads or writes a v2 table
*
*  MW 07/23/2007
*
*  character   crw:        'read' or 'write'
*  character   filename:    name of table
* ------------------------------------------------------------
      Implicit None
      Character*(*) crw,filename
      Integer Nunit, Ifile,ic,i,j,k,l,n,m,nxmax
      Include 'fnx9999.inc'

ckr      write(*,*)"RW: Filename: ",FILENAME

      Nunit=2

      If (crw.eq.'write') Then
         Open(unit=Nunit,file=filename,status='unknown')
      Else
         OPEN(Nunit,STATUS='OLD',FILE=FILENAME,IOSTAT=IFILE)
         IF (Ifile .ne. 0) THEN
            WRITE(*,*)"FX9999RW: ERROR! Table file not found, "//
     >           "stopped! IOSTAT = ",Ifile
            STOP
         Endif
      Endif

c --- fastNLO table
c - block A1
      Call fnioisep(crw,nunit)
      Call fnioint(crw,nunit, Itabversion)
      Call fniochar(crw,nunit, ScenName)
      Call fnioint(crw,nunit, Ncontrib)
      Call fnioint(crw,nunit, Nmult)
      Call fnioint(crw,nunit, Ndata)
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
ckr Define counters and pointers to changes of borders in i'th dimension
ckr (RapIndex, NRapidity, NPt).
ckr Useful in case of npt pt bins in nrap rapidity bins etc.
ckr With LoBin(i,j) works only in two dimensions! 
      Do i=1,NDim
         IDimPointer(1,i) = 1
         NDimCounter(i)   = 1
      Enddo
ckr
      Do i=1,NObsBin
         Do j=1,NDim
            Call fniodbl(crw,nunit, LoBin(i,j))
            If (IDiffBin(j).eq.2) Call fniodbl(crw,nunit, UpBin(i,j))
         Enddo
ckr
         If (i.gt.1) Then
            If (LoBin(i-1,1).ne.LoBin(i,1)) Then
               NDimCounter(1) = NDimCounter(1) + 1  
               IDimPointer(NDimCounter(1),1) = i
            Endif
            If (LoBin(i-1,2).ne.LoBin(i,2)) Then
               NDimCounter(2) = NDimCounter(2) + 1  
               IDimPointer(NDimCounter(2),2) = i
            Endif
         Endif
Comment:          write(*,*)"iobs, ndimcounter1, ndimcounter2",
Comment:      >        i,NDimCounter(1),NDimCounter(2)
Comment:          write(*,*)"iobs, idimpointer1, idimpointer2",i,
Comment:      >        IDimPointer(NDimCounter(1),1),
Comment:      >        IDimPointer(NDimCounter(2),2)
ckr
      Enddo
ckr DEBUG
      Do i=1,NDim
         Write(*,*)"FX9999RW: INFO: Counted ",NDimCounter(i),
     >        " border changes in dimension ",i
         Write(*,*)"          The associated pointers are: "
         Do j=1,NDimCounter(i)
            Write(*,*)"          counter = ",j,
     >           ", pointer = ",IDimPointer(j,i)
         Enddo
      Enddo
ckr End DEBUG
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

         IF (IPDFdef(ic,1).eq.0) then ! - no predefined set of PDF coefficients
            write(*,*) " case IPDFdef(1)=0 not yet implemented"
            STOP
         Endif
         If (NPDF(ic).gt.0) Then
            Do i=1,NObsBin
               Call fnioint(crw,nunit, Nxtot(ic,1,i))
               Do j=1,Nxtot(ic,1,i)
                  Call fniodbl(crw,nunit, XNode1(ic,i,j))
               Enddo
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
         Do i=1,NScaleDim(ic)
            Do j=1,NScaleVar(ic,i)
               Call fniodbl(crw,nunit, ScaleFac(ic,i,j))
cdebug
               Write(*,*)"FX9999RW: IC,IScaleDim,IScaleVar,ScaleFac",
     >              ic,i,j,ScaleFac(ic,i,j)
cdebug
            Enddo
         Enddo

         Do i=1,NObsBin
            Do j=1,NScaleDim(ic)
               Do k=1,NScaleVar(ic,j)
                  Do l=1,NScaleNode(ic,j)
                     Call fniodbl(crw,nunit, ScaleNode(ic,i,j,k,l))
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
     >                       SigmaTilde(ic,i,1,k,l,m,n))
Comment:                         write(*,*)"ic,iobs,iscv,iscn,ix,isub,st",
Comment:      >                       ic,i,k,l,m,n,
Comment:      >                       SigmaTilde(ic,i,1,k,l,m,n)
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo
         


 100     Continue
      Enddo

c      Call fnioint(crw,nunit, )
c      Call fniodbl(crw,nunit, )
c      Call fniochar(crw,nunit, )



c - end of table
      Call fnioisep(crw,nunit)
      Call fnioisep(crw,nunit)
      Close(2)

      Return
      End
