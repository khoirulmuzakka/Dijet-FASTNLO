*******************************************************************
*******************************************************************
      SUBROUTINE FX9999RD(FILENAME,IFLG)
*-----------------------------------------------------------------
* M. Wobisch   read ASCII table of perturbative coefficients
*
* input: FILENAME  name of table
*        IFLG =0 read header (A) 
*             =1 read header plus LO contribution
*             =2 read header plus LO plus NLO contribution
*             =99 read whole table
*
* To Do:
* - problem: Nevt too long?
* - read data blocks
* - read multiplicative factor-blocks
*
* MW 31/05/2007 completely rewritten for v1.5
*-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) FILENAME
      CHARACTER*255 BUFFER
      INTEGER IFLG, IFIRST, IFILE, IC,I,J,K,L,M,N,   nxmax
      INCLUDE 'fnx9999.inc'

      DATA IFIRST/0/
      SAVE IFIRST

      OPEN(2,STATUS='OLD',FILE=FILENAME,IOSTAT=IFILE)
      IF (IFILE .ne. 0) THEN
         WRITE(*,*) '          fastNLO:  table file not found ',
     +        '  -  IOSTAT = ',IFILE
         STOP
      ENDIF
      
c ---------------------------- block A1
      Read(2,*) i
      If (i.ne.Iseparator) Goto 999
      Read(2,*) Itabversion
      Read(2,*) ScenName
      Read(2,*) Ncontrib
      Read(2,*) Nmult
      Read(2,*) Ndata
c ---------------------------- block A2
      Read(2,*) i
      If (i.ne.Iseparator) Goto 999
      Read(2,*) IpublUnits
      Read(2,*) NscDescript
      If (NscDescript .gt. MxScDescript) then
         write(*,*) ' NscDescript too large ',NscDescript,'<=',MXScDescript
      Endif
      Do i=1,NscDescript
         Read(2,*) ScDescript(i)
      Enddo
      Read(2,*) Ecms
      Read(2,*) ILOord
      Read(2,*) NobsBin
      Read(2,*) NDim
      Do i=1,NDim
         Read(2,*) DimLabel(i)
      Enddo
      Do i=1,NDim
         Read(2,*) IDiffBin(i) 
      Enddo
      Do i=1,NObsBin
         Do j=1,NDim
            Read(2,*) LoBin(i,j)
            If (IDiffBin(j).eq.2) Read(2,*) Upbin(i,j)
         Enddo
      Enddo
      Read(2,*) INormFlag
      If (INormFlag.gt.1) then
         Read(2,*) DenomTable
         Do i=1,NObsBin
            Read(2,*) IDivPointer(i)
         Enddo
      Endif
      
c ---------------------------- block B  <<<< need to read this multiple times
      Do ic=1,NContrib
         Read(2,*) i
         If (i.ne.Iseparator) Goto 999
         Read(2,*) IXsectUnits(ic)
         Read(2,*) IDataFlag(ic)
         Read(2,*) IAddMultFlag(ic)
         Read(2,*) IContrFlag1(ic)
         Read(2,*) IContrFlag2(ic)
         Read(2,*) IContrFlag3(ic)
         Read(2,*) NContrDescr(ic)
         Do i=1,NContrDescr(ic)
            Read(2,*) CtrbDescript(ic,i)
         Enddo
         Read(2,*) NCodeDescr(ic)
         Do i=1,NCodeDescr(ic)
            Read(2,*) CodeDescript(ic,i)
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
         Read(2,*) IRef(ic)
         Read(2,*) IScaleDep(ic)
         Read(2,*) Nevt(ic)
         Read(2,*) Npow(ic)
         Read(2,*) NPDF(ic)
         Do i=1,NPDF(ic)
            Read(2,*) NPDFPDG(ic,i)  
         Enddo
         Read(2,*) NPDFDim(ic)
         Read(2,*) NFragFunc(ic)
         Do i=1,NFragFunc(ic)
            Read(2,*) NFFPDG(ic,i)
         Enddo
         Read(2,*) NFFDim(ic) 
         Read(2,*) NSubproc(ic)
         Read(2,*) IPDFdef(ic,1)   
         Read(2,*) IPDFdef(ic,2)   
         Read(2,*) IPDFdef(ic,3)   
         
         IF (IPDFdef(ic,1).eq.0) then ! - no predefined set of PDF coefficients
            write(*,*) " case IPDFdef(1)=0 not yet implemented"
            STOP
         Endif
         If (NPDF(ic).gt.0) Then
            Do i=1,NObsBin
               Read(2,*) Nxtot(ic,1,i)
               Do j=1,Nxtot(ic,1,i)
                  Read(2,*) XNode1(ic,i,j)
               Enddo
            Enddo 
            If (NPDFDim(ic).eq.2) Then
               Do i=1,NObsBin
                  Read(2,*) Nxtot(ic,2,i)
                  Do j=1,Nxtot(ic,2,i)
                     Read(2,*) XNode2(ic,i,j)
                  Enddo
               Enddo 
            Endif
         Endif
         IF (NFragFunc(ic).gt.0) then ! - no FFs so far
            write(*,*) " fastNLO: no FragFuncs so far"
            STOP
         Endif
         Read(2,*) NScales(ic)
         Read(2,*) NScaleDim(ic)
         Do i=1,NScales(ic)
            Read(2,*) IScale(ic,i)
         Enddo
         Do i=1,NScaleDim(ic)
            Read(2,*) NScaleDescript(ic,i)
            Do j=1,NScaleDescript(ic,i)
               Read(2,*) ScaleDescript(ic,i,j)
            Enddo
         Enddo
         Do i=1,NScaleDim(ic)
            Read(2,*) NScaleVar(ic,i)
            Read(2,*) NScaleNode(ic,i)
         Enddo
         Do i=1,NScaleDim(ic)
            Do j=1,NScaleVar(ic,i)
               Read(2,*) ScaleFac(ic,i,j)
            Enddo
         Enddo

         Do i=1,NObsBin
            Do j=1,NScaleDim(ic)
               Do k=1,NScaleVar(ic,j)
                  Do l=1,NScaleNode(ic,j)
                     Read(2,*)ScaleNode(ic,i,j,k,l)
                  Enddo
               Enddo
            Enddo
         Enddo

         Do i=1,NObsBin
c            Do j=1,NScaleDim(ic)     ! <<<< old, wrong concept
            Do k=1,NScaleVar(ic,1)
               Do l=1,NScaleNode(ic,1)
c            Do k=1,NScaleVar(ic,j)
c               Do l=1,NScaleNode(ic,j)
c               Do k1=1,NScaleVar(ic,j)
c                  Do l1=1,NScaleNode(ic,j)
c               Do k2=1,NScaleVar(ic,j)
c                  Do l2=1,NScaleNode(ic,j)
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
                        Read(2,*) SigmaTilde(ic,i,1,k,l,m,n)
c                           Read(2,*) SigmaTilde(ic,i,k1,l1,k2,l2,m,n)
                     Enddo
                  Enddo
               Enddo
            Enddo
         Enddo

 100     Continue
      Enddo

      Read(2,*) i
      If (i.ne.Iseparator) Goto 999
      Read(2,*) i
      If (i.ne.Iseparator) Goto 999
      CLOSE(2)


      RETURN
 999  continue

      close (2) 
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      write(*,*) " >>>>>   fastNLO error in table format "
      write(*,*) " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
      stop
      RETURN

 5000 FORMAT (A,I12,A,A64) 
      END

