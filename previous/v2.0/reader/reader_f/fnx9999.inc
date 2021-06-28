*     --- force emacs to fortran mode -*-fortran-*-
***********************************************************************
*
*     Common block definition for fastNLO v2
*
*     This code should not be edited!
*
***********************************************************************
*---  Global fastNLO commonblock
      INTEGER IFNFIRST, PREFTAB,
     >     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
     >     PNPCOR, PUEVENT, PDATA
      COMMON /CFASTNLO/ IFNFIRST, PREFTAB,
     >     PORDPTHY, PTHRESHCOR, PQUARKCOMPOSITENESS, PADDLED, PTEVED,
     >     PNPCOR, PUEVENT, PDATA

*---  Define numbering and naming scheme of contribution types
      INTEGER NTYPE, ILO, INLO, ITHC1L, ITHC2L, INPC1, IDATA
      PARAMETER(NTYPE=6,ILO=1,INLO=2,ITHC1L=3,ITHC2L=4,INPC1=5,IDATA=6)
      CHARACTER*6  CCONNAME(6)
      CHARACTER*40 CTYPNAME(20)
      CHARACTER*10 CORDNAME(4,4)
      INTEGER ICONTRTYPE(NTYPE),ICONTRORD(NTYPE)
      DATA CCONNAME/"LO","NLO","THC1L","THC2L","NPC1","DATA"/
      DATA CTYPNAME/"Fixed order calculation","Threshold corrections",
     >     "Electroweak corrections","Non-perturbative corrections",
     >     "Undefined", "Data", "Undefined", "Undefined",
     >     "Undefined", "Undefined", "Undefined",
     >     "Quark compositeness", "ADD-LED", "TeV 1-ED", "Unknown",
     >     "Unknown", "Unknown", "Unknown", "Unknown", "Unknown"/
      DATA CORDNAME/"LO","NLO","NNLO","N3LO",
     >     "1-loop","2-loop","3-loop","4-loop",
     >     "Undef","Undef","Undef","Undef",
     >     "LO MC","NLO MC","NNLO MC","N3LO MC"/
      DATA ICONTRTYPE/1,1,2,2,4,6/
      DATA ICONTRORD/1,2,1,2,1,0/

*---  Fundamental parameters - scenario dependent
      INTEGER MXOBSBIN, MXCTRB, MXNXTOT, MXSCALES, MXSCALENODE,
     >     MXSUBPROC, MXNXMAX
      PARAMETER (MXOBSBIN=216, MXCTRB=6, MXNXTOT=20, MXSCALES=2,
     >     MXSCALENODE=25, MXSUBPROC=7, MXNXMAX=325)
ckr     >     MXSCALENODE=15, MXSUBPROC=7, MXNXMAX=210)

*---  Technical parameters
      INTEGER MXUSER, MXSCDESCRIPT, MXDIM, MXCONTRDESCR,
     >     MXCODEDESCR, MXPDF, MXFF, MXSCALEDIM, MXSCALEDESCRIPT,
     >     MXSCALEVAR, MXNMULT, MXNMUNC, MXNMCOR, MXNDUNC, MXNDCOR
      PARAMETER (MXUSER=10,MXSCDESCRIPT=20,MXDIM=3,MXCONTRDESCR=10,
     >     MXCODEDESCR=10,MXPDF=2,MXFF=1,MXSCALEDIM=1,MXSCALEDESCRIPT=5,
     >     MXSCALEVAR=4,MXNMULT=2,MXNMUNC=3,MXNMCOR=3,
     >     MXNDUNC=10,MXNDCOR=10)

*---  Block A1 (header part 1)
      INTEGER ISEPARATOR, ITABVERSION, NCONTRIB,NMULT,NDATA,
     >     NUSERSTRING,NUSERINT,NUSERFLOAT, USERINT(MXUSER), IMACHINE
      CHARACTER*24 SCENNAME
      CHARACTER*80 USERSTRING(MXUSER)
      DOUBLE PRECISION USERFLOAT(MXUSER)
      PARAMETER (ISEPARATOR=1234567890)

*---  Block A2 (header part 2)
      INTEGER IPUBLUNITS, NSCDESCRIPT, ILOORD, NOBSBIN, NDIM,
     >     IDIFFBIN(MXDIM), INORMFLAG,
     >     IDIVLOPOINTER(MXOBSBIN), IDIVUPPOINTER(MXOBSBIN)
      CHARACTER*80 SCDESCRIPT(MXSCDESCRIPT), DIMLABEL(MXDIM),
     >     DENOMTABLE
      DOUBLE PRECISION ECMS,LOBIN(MXOBSBIN,MXDIM),UPBIN(MXOBSBIN,MXDIM),
     >     BINSIZE(MXOBSBIN)

*---  Block B
      INTEGER IXSECTUNITS(MXCTRB), IDATAFLAG(MXCTRB),
     >     IADDMULTFLAG(MXCTRB), ICONTRFLAG1(MXCTRB),
     >     ICONTRFLAG2(MXCTRB), NSCALEDEP(MXCTRB),
     >     NCONTRDESCR(MXCTRB), NCODEDESCR(MXCTRB)
      CHARACTER*80 CTRBDESCRIPT(MXCTRB,MXCONTRDESCR),
     >     CODEDESCRIPT(MXCTRB,MXCODEDESCR)
      INTEGER*8 NEVT(MXCTRB)
      DOUBLE PRECISION DEVT(MXCTRB)
      INTEGER IREF(MXCTRB), ISCALEDEP(MXCTRB), NPOW(MXCTRB),
     >     NPDF(MXCTRB), NPDFPDG(MXCTRB,MXPDF),NPDFDIM(MXCTRB),
     >     NFRAGFUNC(MXCTRB), NFFPDG(MXCTRB,MXFF), NFFDIM(MXCTRB),
     >     NSUBPROC(MXCTRB), IPDFDEF(MXCTRB,3), NXTOT(MXCTRB,2,MXOBSBIN)
      DOUBLE PRECISION XNODE1(MXCTRB,MXOBSBIN,MXNXTOT),
     >     XNODE2(MXCTRB,MXOBSBIN,MXNXTOT),
     >     HXLIM1(MXCTRB,MXOBSBIN)
      INTEGER NSCALES(MXCTRB),NSCALEDIM(MXCTRB), ISCALE(MXCTRB,MXSCALES)
     >     ,NSCALEDESCRIPT(MXCTRB,MXSCALEDIM)
      CHARACTER*80 SCALEDESCRIPT(MXCTRB,MXSCALEDIM,MXSCALEDESCRIPT)
      INTEGER NSCALEVAR(MXCTRB,MXSCALEDIM),NSCALENODE(MXCTRB,MXSCALEDIM)
      DOUBLE PRECISION SCALEFAC(MXCTRB,MXSCALEDIM,MXSCALEVAR),
     >     SCALENODE(MXCTRB,MXOBSBIN,MXSCALEDIM,MXSCALEVAR,MXSCALENODE),
     >     HSCALENODE(MXCTRB,MXOBSBIN,MXSCALEDIM,MXSCALEVAR,MXSCALENODE)
     >     ,SIGMATILDE(MXCTRB,MXOBSBIN,1,MXSCALEVAR,MXSCALENODE,MXNXMAX
     >     ,MXSUBPROC)
      INTEGER NMUNCORREL(MXNMULT),NMCORREL(MXNMULT)
      INTEGER NDUNCORREL,NDCORREL,NDERRMATRIX
      CHARACTER*80 MUNCDESCRIPT(MXNMULT,MXNMUNC)
      CHARACTER*80 MCORDESCRIPT(MXNMULT,MXNMCOR)
      CHARACTER*80 DUNCDESCRIPT(MXNDUNC)
      CHARACTER*80 DCORDESCRIPT(MXNDCOR)
      DOUBLE PRECISION MFACT(MXNMULT,MXOBSBIN)
      DOUBLE PRECISION DXVAL(MXOBSBIN),DYVAL(MXOBSBIN)
      DOUBLE PRECISION MUNCORLO(MXNMULT,MXOBSBIN,MXNMUNC)
      DOUBLE PRECISION MUNCORUP(MXNMULT,MXOBSBIN,MXNMUNC)
      DOUBLE PRECISION MCORLO(MXNMULT,MXOBSBIN,MXNMCOR)
      DOUBLE PRECISION MCORUP(MXNMULT,MXOBSBIN,MXNMCOR)
      DOUBLE PRECISION DUNCORLO(MXOBSBIN,MXNDUNC)
      DOUBLE PRECISION DUNCORUP(MXOBSBIN,MXNDUNC)
      DOUBLE PRECISION DCORLO(MXOBSBIN,MXNDCOR)
      DOUBLE PRECISION DCORUP(MXOBSBIN,MXNDCOR)

*---  Output arrays
      DOUBLE PRECISION XSECT(MXOBSBIN),XSNORM(MXOBSBIN)
      DOUBLE PRECISION XSCALE(MXOBSBIN)
      DOUBLE PRECISION XSUNCOR(MXOBSBIN,2),XSCOR(MXOBSBIN,2)

*---  Internal variables
*---  Add DimPointers/Counters
      INTEGER NBINCOUNTER(MXDIM)
      INTEGER NDIVCOUNTER(MXOBSBIN)
      INTEGER IOBSPOINTER(MXOBSBIN,MXOBSBIN)

      INTEGER ICONTRPOINTER(MXCTRB),NCONTRCOUNTER(MXCTRB)
      INTEGER ICONTRSELECTOR(MXCTRB),ISCALEPOINTER(MXCTRB)
      LOGICAL LCONTR(MXCTRB)
      DOUBLE PRECISION RESULT(MXOBSBIN,0:MXSUBPROC,0:MXCTRB),
     >     RESMUR(MXOBSBIN,0:MXSUBPROC,0:MXCTRB),
     >     PDF(MXOBSBIN,MXSCALENODE,MXNXMAX,MXSUBPROC)

*---  Put everything into common block
      COMMON /CFN9999/

*---  dbl
     >     USERFLOAT,
     >     ECMS,LOBIN,UPBIN,BINSIZE, XNODE1, XNODE2,HXLIM1,
     >     SCALEFAC, SCALENODE, HSCALENODE, SIGMATILDE,
     >     RESULT, RESMUR, PDF,
     >     MFACT, DXVAL, DYVAL,
     >     MUNCORLO, MUNCORUP, MCORLO, MCORUP,
     >     DUNCORLO, DUNCORUP, DCORLO, DCORUP, DEVT,
*---  int*8
     >     NEVT,
*---  int
     >     ITABVERSION, NCONTRIB,NMULT,NDATA,
     >     NUSERSTRING,NUSERINT,NUSERFLOAT,USERINT,IMACHINE,
     >     IPUBLUNITS, NSCDESCRIPT, ILOORD, NOBSBIN, NDIM,
     >     IDIFFBIN, INORMFLAG, IDIVLOPOINTER, IDIVUPPOINTER,
     >     IXSECTUNITS, IDATAFLAG, IADDMULTFLAG,
     >     ICONTRFLAG1, ICONTRFLAG2, NSCALEDEP,
     >     NCONTRDESCR, NCODEDESCR, IREF, ISCALEDEP, NPOW,
     >     NPDF, NPDFPDG, NPDFDIM,
     >     NFRAGFUNC, NFFPDG, NFFDIM, NSUBPROC, IPDFDEF, NXTOT,
     >     NSCALES, NSCALEDIM, ISCALE, NSCALEDESCRIPT,
     >     NSCALEVAR, NSCALENODE,
     >     NMUNCORREL, NMCORREL,
     >     NDUNCORREL, NDCORREL, NDERRMATRIX,
     >     NBINCOUNTER,NDIVCOUNTER,IOBSPOINTER,
     >     ICONTRPOINTER,NCONTRCOUNTER,ICONTRSELECTOR,
     >     ISCALEPOINTER,
*---  log
     >     LCONTR,
*---  char
     >     SCENNAME,USERSTRING,SCDESCRIPT, DIMLABEL, DENOMTABLE,
     >     CTRBDESCRIPT, CODEDESCRIPT, SCALEDESCRIPT,
     >     MUNCDESCRIPT, MCORDESCRIPT,
     >     DUNCDESCRIPT, DCORDESCRIPT

C---  Tests
C---  INTEGER IRAPDIM,NRAPIDITY,IPTDIM,NPT(MXOBSBIN)
C---  COMMON/TFN9999/IRAPDIM,NRAPIDITY,IPTDIM,NPT

*---  Moved from strings.inc
      CHARACTER*2 CUNITS(0:15)
      DATA CUNITS/' b','--','--','mb','--','--','ub','--','--','nb',
     >     '--','--','pb','--','--','fb'/