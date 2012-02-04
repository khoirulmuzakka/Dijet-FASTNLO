* --- force emacs to fortran mode -*-fortran-*-
c ------------------------------------------------------------
c
c  Commonblock definition for fastNLO  v2.0
c
c   >> this code should not be edited
c   >> if you have questions, please contact 
c      the authors at:  fastnlo@cedar.ac.uk
c 
c ------------------------------------------------------------

c - global fastNLO commonblock
      Integer IFNfirst, Preftab,
     +     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
     +     PNPCOR, PUEVENT, PDATA
      Common /cfastnlo/ IFNfirst, Preftab,
     +     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
     +     PNPCOR, PUEVENT, PDATA

c - define numbering scheme of contribution types
      Integer ILO, INLO, ITHC1L, ITHC2L, INPC1, IDATA 
      Parameter(ILO=1,INLO=2,ITHC1L=3,ITHC2L=4,INPC1=5,IDATA=6)
      
c - fundamental parameters - scenario dependent
      Integer MxObsBin, MxCtrb, MxNxtot, MxScales, MxScaleNode, 
     +     MxSubProc, MxNxMax
      Parameter (MxObsBin=200, MxCtrb=6, MxNxtot=20, MxScales=2, 
     +     MxScaleNode=15, MxSubProc=7, MxNxMax=210)
c - technical parameters
      Integer MxUser, MxScDescript, MxDim, MxContrDescr,
     +     MxCodeDescr, MxPDF, MxFF, MxScaleDim, MxScaleDescript,
     +     MxScaleVar, MxNMult, MxNMUnc, MxNMCor, MxNDUnc, MxNDCor
      Parameter (MxUser=10,MxScDescript=5,MxDim=3,MxContrDescr=10,
     +     MxCodeDescr=10,MxPDF=2,MxFF=1,MxScaleDim=1,MxScaleDescript=5,
     +     MxScaleVar=4,MxNMult=2,MxNMUnc=3,MxNMCor=3,
     +     MxNDUnc=10,MxNDCor=10)

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
     +     IContrFlag2(MxCtrb), NScaleDep(MxCtrb), 
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
      Integer NMUncorrel(MxNMult),NMCorrel(MxNMult)
      Integer NDUncorrel,NDCorrel,NDErrMatrix
      Character*64 MUncDescript(MxNMult,MxNMUnc)
      Character*64 MCorDescript(MxNMult,MxNMCor)
      Character*64 DUncDescript(MxNDUnc)
      Character*64 DCorDescript(MxNDCor)
      Double Precision MFact(MxNMult,MxObsBin)
      Double Precision DxVal(MxObsBin),DyVal(MxObsBin)
      Double Precision MUnCorLo(MxNMult,MxObsBin,MxNMUnc)
      Double Precision MUnCorUp(MxNMult,MxObsBin,MxNMUnc)
      Double Precision MCorLo(MxNMult,MxObsBin,MxNMCor)
      Double Precision MCorUp(MxNMult,MxObsBin,MxNMCor)
      Double Precision DUnCorLo(MxObsBin,MxNDUnc)
      Double Precision DUnCorUp(MxObsBin,MxNDUnc)
      Double Precision DCorLo(MxObsBin,MxNDCor)
      Double Precision DCorUp(MxObsBin,MxNDCor)

c - output arrays
      Double Precision Xsect(MxObsBin),Xsnorm(MxObsBin)
      Double Precision XsUnCor(MxObsBin,2),XsCor(MxObsBin,2)

c - internal variables
ckr Add DimPointers/Counters
      Integer NDimCounter(MxDim),IDimPointer(MxObsBin,MxDim)
ckr      Integer IContr
      Integer IContrPointer(MxCtrb),NContrCounter(MxCtrb)
      Integer IContrSelector(MxCtrb),IScalePointer(MxCtrb)
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
     +     MFact, DxVal, DyVal,
     +     MUnCorLo, MUnCorUp, MCorLo, MCorUp,
     +     DUnCorLo, DUnCorUp, DCorLo, DCorUp,
c - int*8
     +     Nevt,
c - int
     +     Itabversion, Ncontrib,Nmult,Ndata,
     +     NuserString,NuserInt,NuserFloat,UserInt,Imachine,
     +     Ipublunits, NscDescript, ILOord, Nobsbin, NDim,
     +     IDiffBin, INormFlag, IDivLoPointer, IDivUpPointer,
     +     IXsectUnits, IDataFlag, IAddMultFlag, 
     +     IContrFlag1, IContrFlag2, NScaleDep, 
     +     NContrDescr, NCodeDescr, IRef, IScaleDep, Npow, 
     +     NPDF, NPDFPDG, NPDFDim,
     +     NFragFunc, NFFPDG, NFFDim, NSubproc, IPDFdef, Nxtot,
     +     NScales, NScaleDim, IScale, NScaleDescript,
     +     NScaleVar, NScaleNode,
     +     NMUncorrel, NMCorrel,
     +     NDUncorrel, NDCorrel, NDErrMatrix,
ckr
     +     NDimCounter,IDimPointer,
ckr     +     IContr,
     +     IContrPointer,NContrCounter,IContrSelector,
     +     IScalePointer,
c - char
     +     ScenName,UserString,ScDescript, DimLabel, DenomTable,
     +     CtrbDescript, CodeDescript, ScaleDescript,
     +     MUncDescript, MCorDescript,
     +     DUncDescript, DCorDescript

ckr Add temporarily ?
ckr
      INTEGER IRAPDIM,NRAPIDITY,IPTDIM,NPT(MXOBSBIN)
      COMMON/TFN9999/IRAPDIM,NRAPIDITY,IPTDIM,NPT
ckr
ckr
ckr Moved from strings.inc
      CHARACTER*2 Cunits(0:15)
      Data Cunits/' b','--','--','mb','--','--','ub','--','--','nb',
     +     '--','--','pb','--','--','fb'/
