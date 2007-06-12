#ifndef __physconst__
#define __physconst__

#define MZ       91.1882 // Z mass, review of particle physics 2000
#define PI        3.14159265358979323846264
#define PIHALF    PI/2.0
#define TWOPI     6.28318530717958647692528
#define TWOPISQR 39.47841760435743447533796
#define NF    5 // number of active flavours

#define BETA0  (11. - 2./3.*NF) // The beta coefficients of the QCD beta function
#define BETA1  (51. - 19./3.*NF)

#define QCDB         BETA0/(4.*PI)           // The b coefficients 
#define QCDBPRIME 2.*BETA1/(4.*TWOPISQR*QCDB) 
#define LAMBDA12  1.147686435 //Factor between definitions (1) and (2) of LambdaQCD (see p35 Ellis,Sterling,Webber)
#define RADTODEG  57.2957795131

#endif
