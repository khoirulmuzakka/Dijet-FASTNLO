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
c ------------------------------------------------------------

c - global fastNLO commonblock
      Integer IFNfirst,
     +     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
     +     PHADRCOR, PUEVENT
      Common /cfastnlo/ IFNfirst,
     +     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
     +     PHADRCOR, PUEVENT

c - fundamental parameters - scenario dependent
      Integer MxObsBin, MxCtrb, MxNxtot, MXScales, MxScaleNode, 
     +     MxSubProc, MxNxMax
      Parameter (MxObsBin=200,MxCtrb=3, MxNxtot=20, MxScales=2, 
     +     MxScaleNode=5, MxSubProc=7, MxNxMax=210)
c - technical parameters
      Integer MxScDescript, MxDim, MxContrDescr,
     +     MxCodeDescr, MxPDF, MxFF, MxScaleDim, MxScaleDescript,
     +     MxScaleVar
      Parameter (MxScDescript=5,MxDim=3,MxContrDescr=10,
     +     MxCodeDescr=10,MxPDF=2,MxFF=0,MxScaleDim=1,MxScaleDescript=5,
     +     MxScaleVar=4)

c --- Block A1 (header part 1)
      Integer Iseparator, Itabversion, Ncontrib,Nmult,Ndata
      Character*24 ScenName
      Parameter (Iseparator=1234567890)

c --- Block A2 (header part 2)
      Integer Ipublunits, NscDescript, ILOord, Nobsbin, NDim,
     +     IDiffBin(MxDim), INormFlag, IDivPointer(MxObsBin)
      Character*64 ScDescript(MxScDescript), DimLabel(MxDim),
     +     DenomTable
      Double Precision Ecms,LoBin(MxObsBin,MxDim),UpBin(MxObsBin,MxDim)

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
     +     NSubproc(MxCtrb), IPDFCoeff(MxCtrb),Nxtot(MxCtrb,2)
      Double Precision XNode1(MxCtrb,MxObsBin,MxNxtot),
     +     XNode2(MxCtrb,MxObsBin,MxNxtot)
      Integer NScales(MxCtrb),NScaleDim(MxCtrb), IScale(MxCtrb,MxScales), 
     +     NScaleDescript(MxCtrb,MxScaleDim)
      Character*64 ScaleDescript(MxCtrb,MxScaleDim,MxScaleDescript)
      Integer NScaleVar(MxCtrb,MxScaleDim),NScaleNode(MxCtrb,MxScaleDim)
      Double Precision ScaleFac(MxCtrb,MxScaleDim,MxScaleVar),
     +     ScaleNode(MxCtrb,MxObsBin,MxScaleDim,MxScaleVar,MxScaleNode),
     +     SigmaTilde(MxCtrb,MxObsBin,MxScaleDim,MxScaleVar,MxScaleNode,
     +     MxNxMax,MxSubProc)

c - output array
      Double Precision xsect(MxObsBin), xsect2(MxObsBin)

c - internal variables
      Integer IContrib, IContrPoint(MxCtrb),IScalePoint(MxCtrb)
      Double Precision result(MxObsBin,MxSubproc,MxCtrb)


c - put everything into common block
      Common /CFN9999/
c - dbl 
     +     Ecms,LoBin,UpBin, XNode1, XNode2,
     +     ScaleFac, ScaleNode, SigmaTilde,
     +     result,
c - int*8
     +     Nevt,
c - int
     +     Itabversion, Ncontrib,Nmult,Ndata,
     +     Ipublunits, NscDescript, ILOord, Nobsbin, NDim,
     +     IDiffBin, INormFlag, IDivPointer,
     +     IXsectUnits, IDataFlag, IAddMultFlag, 
     +     IContrFlag1, IContrFlag2, IContrFlag3, 
     +     NContrDescr, NCodeDescr, IRef, IScaleDep, Npow, 
     +     NPDF, NPDFPDG, NPDFDim,
     +     NFragFunc, NFFPDG, NFFDim, NSubproc, IPDFCoeff, Nxtot,
     +     NScales, NScaleDim, IScale, NScaleDescript,
     +     NScaleVar, NScaleNode,
     +     IContrib, IContrPoint, IScalePoint,
c - char
     +     ScenName, ScDescript, DimLabel, DenomTable,
     +     CtrbDescript, CodeDescript, ScaleDescript


c ------------------ old stuff
c
      CHARACTER*17 CIREACTION(4)
      CHARACTER*17 CIPROC(4)
      CHARACTER*28 CIALGO(4)
      CHARACTER*5 CJETRES1(4)
      CHARACTER*12 CJETRES2(4)

      DATA CIREACTION/'electron-proton  ','proton-proton    ',
     +     'proton-antiproton','photon-proton'/
      DATA CIPROC/'inclusive jets   ','dijets','three-jets','four-jets'/
      DATA CIALGO/'kT algorithm','midpoint cone algorithm',
     +     'rsep cone algorithm','searchcone algorithm'/
      DATA CJETRES1/'D','Rcone','Rcone','Rcone'/
      DATA CJETRES2/'n/a','f_overlap','Rsep','R_searchcone'/

      CHARACTER*64 CHEADER(22)
      DATA CHEADER/
     +     '################################################################',
     +     ' ',
     +     '  fastNLO - version 2.0         June 3, 2007',
     +     ' ',
     +     '  Thomas Kluge, Klaus Rabbertz, Markus Wobisch',
     +     ' ',
     +     '  if you use this code, please cite as reference:',
     +     '      T. Kluge, K. Rabbertz and M. Wobisch,',
     +     '         publication in preparation,',
     +     '         (hep-ex No. to be added)',
     +     ' ',
     +     '     .... now reading the coefficient table:',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' ',
     +     ' '/

      CHARACTER*64 CNLOJET(4)
      DATA CNLOJET/
     +     '  for NLOJET++ please cite:',
     +     '      Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),',
     +     '      Z. Nagy, Phys. Rev. D68, 094002 (2003).',
     +     ' '/

      CHARACTER*64 CTHRCOR(5)
      DATA CTHRCOR/
     +     '  the 2-loop threshold corrections for the inclusive jet',
     +     '  cross section in pp and ppbar have been computed by:',
     +     '  N. Kidonakis and J.F. Owens - please cite',
     +     '     N. Kidonakis, J.F. Owens, Phys. Rev. D63, 054019 (2001).',
     +     ' '/
