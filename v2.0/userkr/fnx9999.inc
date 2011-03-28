* --- force emacs to fortran mode -*-fortran-*-
c ------------------------------------------------------------
c
c  Commonblock definition for fastNLO  v2.0
c
c   >> this code should not be edited
c   >> if you have questions, please contact 
c      the authors at:  fastnlo@cedar.ac.uk
c 
c MW 06/01/07 new design for v2.0
c KR 19.09.09 add second dimension of xsect(,MxCtrb) as before
c                 MW -> remove this 2nd dimension again
c KR 19.09.09 add dimension counters and pointers
c             (can serve as NRapidity, RapIndex, NPt ...)
c ------------------------------------------------------------

c - global fastNLO commonblock
      Integer IFNfirst,
     +     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
     +     PHADRCOR, PUEVENT
      Common /cfastnlo/ IFNfirst,
     +     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
     +     PHADRCOR, PUEVENT

c - fundamental parameters - scenario dependent
      Integer MxObsBin, MxCtrb, MxNxtot, MxScales, MxScaleNode, 
     +     MxSubProc, MxNxMax
      Parameter (MxObsBin=200, MxCtrb=4, MxNxtot=20, MxScales=2, 
     +     MxScaleNode=15, MxSubProc=7, MxNxMax=210)
c - technical parameters
      Integer MxUser, MxScDescript, MxDim, MxContrDescr,
     +     MxCodeDescr, MxPDF, MxFF, MxScaleDim, MxScaleDescript,
     +     MxScaleVar
      Parameter (MxUser=10,MxScDescript=5,MxDim=3,MxContrDescr=10,
     +     MxCodeDescr=10,MxPDF=2,MxFF=1,MxScaleDim=1,MxScaleDescript=5,
     +     MxScaleVar=4)

c --- Block A1 (header part 1)
      Integer Iseparator, Itabversion, Ncontrib,Nmult,Ndata,
     +     NuserString,NuserInt,NuserFloat, UserInt(MxUser), Imachine
      Character*24 ScenName
      Character*64 UserString(MxUser)
      Double Precision UserFloat(MxUser)
      Parameter (Iseparator=1234567890)

c --- Block A2 (header part 2)
      Integer Ipublunits, NscDescript, ILOord, Nobsbin, NDim,
     +     IDiffBin(MxDim), INormFlag, 
     +     IDivLoPointer(MxObsBin), IDivUpPointer(MxObsBin)
      Character*64 ScDescript(MxScDescript), DimLabel(MxDim),
     +     DenomTable
      Double Precision Ecms,LoBin(MxObsBin,MxDim),UpBin(MxObsBin,MxDim),
     +     BinSize(MxObsBin)

c --- Block B
      Integer IXsectUnits(MxCtrb), IDataFlag(MxCtrb), 
     +     IAddMultFlag(MxCtrb), IContrFlag1(MxCtrb),
     +     IContrFlag2(MxCtrb), IContrFlag3(MxCtrb), 
     +     NContrDescr(MxCtrb), NCodeDescr(MxCtrb)
      Character*64 CtrbDescript(MxCtrb,MxContrDescr), 
     +     CodeDescript(MxCtrb,MxCodeDescr)
      Integer*8 Nevt(MxCtrb)
      Integer IRef(MxCtrb), IScaleDep(MxCtrb), Npow(MxCtrb), 
     +     NPDF(MxCtrb), NPDFPDG(MxCtrb,MxPDF),NPDFDim(MxCtrb),
     +     NFragFunc(MxCtrb), NFFPDG(MxCtrb,MxFF), NFFDim(MxCtrb),
     +     NSubproc(MxCtrb), IPDFdef(MxCtrb,3), Nxtot(MxCtrb,2,MxObsBin)
      Double Precision XNode1(MxCtrb,MxObsBin,MxNxtot),
     +     XNode2(MxCtrb,MxObsBin,MxNxtot),
     +     Hxlim1(MxCtrb,MxObsBin)
      Integer NScales(MxCtrb),NScaleDim(MxCtrb), IScale(MxCtrb,MxScales), 
     +     NScaleDescript(MxCtrb,MxScaleDim)
      Character*64 ScaleDescript(MxCtrb,MxScaleDim,MxScaleDescript)
      Integer NScaleVar(MxCtrb,MxScaleDim),NScaleNode(MxCtrb,MxScaleDim)
      Double Precision ScaleFac(MxCtrb,MxScaleDim,MxScaleVar),
     +     ScaleNode(MxCtrb,MxObsBin,MxScaleDim,MxScaleVar,MxScaleNode),
     +     HScaleNode(MxCtrb,MxObsBin,MxScaleDim,MxScaleVar,MxScaleNode),
     +     SigmaTilde(MxCtrb,MxObsBin,1,MxScaleVar,MxScaleNode,
     +     MxNxMax,MxSubProc)

c - output array
ckr      Double Precision xsect(MxObsBin,0:MxCtrb),
ckr     +     xsect2(MxObsBin,0:MxCtrb)
      Double Precision xsect(MxObsBin),
     +     xsect2(MxObsBin)

c - internal variables
ckr Add DimPointers/Counters
      Integer NDimCounter(MxDim),IDimPointer(MxObsBin,MxDim)
      Integer IContr,IContrPointer(MxCtrb),IScalePointer(MxCtrb)
cMW      Integer IContr,IContrPointer(0:MxCtrb),IScalePointer(MxCtrb)
ckr      Double Precision Result(MxObsBin,MxSubproc,MxCtrb),
      Double Precision Result(MxObsBin,0:MxSubproc,0:MxCtrb),
     +     pdf(MxObsBin,MxScaleNode,MxNxMax,MxSubproc)
      

c - put everything into common block
      Common /CFN9999/
c - dbl 
     +     UserFloat,
     +     Ecms,LoBin,UpBin,BinSize, XNode1, XNode2,Hxlim1,
     +     ScaleFac, ScaleNode, HScaleNode, SigmaTilde,
     +     result,pdf,
c - int*8
     +     Nevt,
c - int
     +     Itabversion, Ncontrib,Nmult,Ndata,
     +     NuserString,NuserInt,NuserFloat,UserInt,Imachine,
     +     Ipublunits, NscDescript, ILOord, Nobsbin, NDim,
     +     IDiffBin, INormFlag, IDivLoPointer, IDivUpPointer,
     +     IXsectUnits, IDataFlag, IAddMultFlag, 
     +     IContrFlag1, IContrFlag2, IContrFlag3, 
     +     NContrDescr, NCodeDescr, IRef, IScaleDep, Npow, 
     +     NPDF, NPDFPDG, NPDFDim,
     +     NFragFunc, NFFPDG, NFFDim, NSubproc, IPDFdef, Nxtot,
     +     NScales, NScaleDim, IScale, NScaleDescript,
     +     NScaleVar, NScaleNode,
ckr
     +     NDimCounter,IDimPointer,
     +     IContr,IContrPointer,IScalePointer,
c - char
     +     ScenName,UserString,ScDescript, DimLabel, DenomTable,
     +     CtrbDescript, CodeDescript, ScaleDescript

ckr Add temporarily ?
ckr
      INTEGER IRAPDIM,NRAPIDITY,IPTDIM,NPT(MXOBSBIN)
      COMMON/TFN9999/IRAPDIM,NRAPIDITY,IPTDIM,NPT
ckr
