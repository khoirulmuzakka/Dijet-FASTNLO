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
c KR 19.09.09 add some derived properties of table
c             NRapidity, RapIndex, NPt ...
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
      Parameter (MxObsBin=200, MxCtrb=3, MxNxtot=20, MxScales=2, 
     +     MxScaleNode=5, MxSubProc=7, MxNxMax=210)
c - technical parameters
      Integer MxScDescript, MxDim, MxContrDescr,
     +     MxCodeDescr, MxPDF, MxFF, MxScaleDim, MxScaleDescript,
     +     MxScaleVar
      Parameter (MxScDescript=5,MxDim=3,MxContrDescr=10,
     +     MxCodeDescr=10,MxPDF=2,MxFF=1,MxScaleDim=1,MxScaleDescript=5,
     +     MxScaleVar=4)

c --- Block A1 (header part 1)
      Integer Iseparator, Itabversion, Ncontrib,Nmult,Ndata
      Character*24 ScenName
      Parameter (Iseparator=1234567890)

c --- Block A2 (header part 2)
      Integer Ipublunits, NscDescript, ILOord, Nobsbin, NDim,
     +     IDiffBin(MxDim), INormFlag, 
     +     IDivLoPointer(MxObsBin), IDivUpPointer(MxObsBin)
      Character*64 ScDescript(MxScDescript), DimLabel(MxDim),
     +     DenomTable
      Double Precision Ecms,LoBin(MxObsBin,MxDim),UpBin(MxObsBin,MxDim)
ckr
      Integer NRapidity,RapIndex(MxObsBin),NPt(MxObsBin)
ckr

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
     +     XNode2(MxCtrb,MxObsBin,MxNxtot)
      Integer NScales(MxCtrb),NScaleDim(MxCtrb), IScale(MxCtrb,MxScales), 
     +     NScaleDescript(MxCtrb,MxScaleDim)
      Character*64 ScaleDescript(MxCtrb,MxScaleDim,MxScaleDescript)
      Integer NScaleVar(MxCtrb,MxScaleDim),NScaleNode(MxCtrb,MxScaleDim)
      Double Precision ScaleFac(MxCtrb,MxScaleDim,MxScaleVar),
     +     ScaleNode(MxCtrb,MxObsBin,MxScaleDim,MxScaleVar,MxScaleNode),
     +     SigmaTilde(MxCtrb,MxObsBin,1,MxScaleVar,MxScaleNode,
     +     MxNxMax,MxSubProc)

c - output array
      Double Precision xsect(MxObsBin,0:MxCtrb),
     +     xsect2(MxObsBin,0:MxCtrb)

c - internal variables
      Integer IContrib, IContrPoint(MxCtrb),IScalePoint(MxCtrb)
      Double Precision result(MxObsBin,MxSubproc,MxCtrb),
     +        pdf(MxObsBin,MxScaleNode,MxNxMax,MxSubproc)


c - put everything into common block
      Common /CFN9999/
c - dbl 
     +     Ecms,LoBin,UpBin, XNode1, XNode2,
     +     ScaleFac, ScaleNode, SigmaTilde,
     +     result,pdf,
c - int*8
     +     Nevt,
c - int
     +     Itabversion, Ncontrib,Nmult,Ndata,
     +     Ipublunits, NscDescript, ILOord, Nobsbin, NDim,
     +     IDiffBin, INormFlag, IDivLoPointer, IDivUpPointer,
     +     IXsectUnits, IDataFlag, IAddMultFlag, 
     +     IContrFlag1, IContrFlag2, IContrFlag3, 
     +     NContrDescr, NCodeDescr, IRef, IScaleDep, Npow, 
     +     NPDF, NPDFPDG, NPDFDim,
     +     NFragFunc, NFFPDG, NFFDim, NSubproc, IPDFdef, Nxtot,
     +     NScales, NScaleDim, IScale, NScaleDescript,
     +     NScaleVar, NScaleNode,
     +     IContrib, IContrPoint, IScalePoint,
c - char
     +     ScenName, ScDescript, DimLabel, DenomTable,
     +     CtrbDescript, CodeDescript, ScaleDescript


c ------------------ old stuff
c
Comment:       CHARACTER*17 CIREACTION(4)
Comment:       CHARACTER*17 CIPROC(4)
Comment:       CHARACTER*28 CIALGO(4)
Comment:       CHARACTER*5 CJETRES1(4)
Comment:       CHARACTER*12 CJETRES2(4)
Comment: 
Comment:       DATA CIREACTION/'electron-proton  ','proton-proton    ',
Comment:      +     'proton-antiproton','photon-proton'/
Comment:       DATA CIPROC/'inclusive jets   ','dijets','three-jets','four-jets'/
Comment:       DATA CIALGO/'kT algorithm','midpoint cone algorithm',
Comment:      +     'rsep cone algorithm','searchcone algorithm'/
Comment:       DATA CJETRES1/'D','Rcone','Rcone','Rcone'/
Comment:       DATA CJETRES2/'n/a','f_overlap','Rsep','R_searchcone'/
Comment: 
Comment:       CHARACTER*64 CHEADER(22)
Comment:       DATA CHEADER/
Comment:      +     '################################################################',
Comment:      +     ' ',
Comment:      +     '  fastNLO - version 2.0         June 3, 2007',
Comment:      +     ' ',
Comment:      +     '  Thomas Kluge, Klaus Rabbertz, Markus Wobisch',
Comment:      +     ' ',
Comment:      +     '  if you use this code, please cite as reference:',
Comment:      +     '      T. Kluge, K. Rabbertz and M. Wobisch,',
Comment:      +     '         publication in preparation,',
Comment:      +     '         (hep-ex No. to be added)',
Comment:      +     ' ',
Comment:      +     '     .... now reading the coefficient table:',
Comment:      +     ' ',
Comment:      +     ' ',
Comment:      +     ' ',
Comment:      +     ' ',
Comment:      +     ' ',
Comment:      +     ' ',
Comment:      +     ' ',
Comment:      +     ' ',
Comment:      +     ' ',
Comment:      +     ' '/
Comment: 
Comment:       CHARACTER*64 CNLOJET(4)
Comment:       DATA CNLOJET/
Comment:      +     '  for NLOJET++ please cite:',
Comment:      +     '      Z. Nagy, Phys. Rev. Lett. 88, 122003 (2002),',
Comment:      +     '      Z. Nagy, Phys. Rev. D68, 094002 (2003).',
Comment:      +     ' '/
Comment: 
Comment:       CHARACTER*64 CTHRCOR(5)
Comment:       DATA CTHRCOR/
Comment:      +     '  the 2-loop threshold corrections for the inclusive jet',
Comment:      +     '  cross section in pp and ppbar have been computed by:',
Comment:      +     '  N. Kidonakis and J.F. Owens - please cite',
Comment:      +     '     N. Kidonakis, J.F. Owens, Phys. Rev. D63, 054019 (2001).',
Comment:      +     ' '/
