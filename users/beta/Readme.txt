
************************************************************************
**********    fastNLO     - beta-release    Setember 13, 2005
************************************************************************


The fastNLO beta release includes everything needed to compute
the inclusive jet cross section, as measured at the Tevatron
in Run I by the CDF and the D0 Collaborations, in NLO perturbative QCD 

  B. Abbott et al. (D0 Collaboration), Phys. Rev. Lett. {86} 1707 (2001).
  T. Affolder et al. (CDF Collaboration), Phys. Rev. D64, 032001 (2001).



* -------------------------------------------------------------
*     Files
* -------------------------------------------------------------
Two files are needed: 
   fbt1001.tab.gz         the table with the precomputed 
                          perturbative coefficients
   fastNLO-beta.tar.gz    the user code to read the table, extract 
                          the PDF information from LHAPDF, 
                          and compute the cross sections 

The following files are stored in fastNLO-beta.tar.gz:
   Readme.txt          this file
   Makefile            (you need to set the paths to LHAPDF/CERNLIB)
   example01.f         an example how to call the fastNLO code
                       (you need to set the path to the LHAPDF grid file)

   fbt1001.f           code, specific for observable
   fbt1001.inc         includes variable definitions
   fn-pdf.f            compute linear combinations of PDFs
   fn-alphas-demo.f    example routine to compute alpha_s

   gz.F                three files, needed to read gzipped tables
   ftn_gzio.c                    - " -
   gzio.inc                      - " -
   
   fn-interface.f      user-interface to PDFs and alpha_s
                       (to be adopted for individual needs)


* -------------------------------------------------------------
*     Installation 
* -------------------------------------------------------------

(1) LHAPDF:
   The fastNLO example uses LHAPDF to access the PDFs.
   We recommend to use LHAPDF and install the latest version from:
      http://hepforge.cedar.ac.uk/lhapdf/
   The installation procedure is very easy, following the guide at:
      http://hepforge.cedar.ac.uk/lhapdf/install


(2) fastNLO
   unpack the code:
      tar -xvzf fastNLO-beta.tar.gz
   This extracts the files into a new directory: "fastNLO-beta"
   Now go into this directory and edit the "Makefile":
     - set the path to your LHAPDF library file "libLHAPDF.a"
     - set the path to your CERNLIB
   Edit the file "example1.f"
     The grid files for the PDFs are stored in the
     subdirectory "PDFsets" of LHAPDF.
     In the call of "InitPDFset" (in example1.f) you need
     to set the path to the CTEQ6.1 grid "cteq61.LHgrid"
     Example:
       call InitPDFset('/path/lhapdf-4.1/PDFsets/cteq61.LHgrid')

   Copy the table file "fbt1001.tab.gz" from the webpage
   into your fastNLO-beta directory (this is 36M large!)
   You have two options:
     - Either you keep the table file gzipped.
       This needs less storage, but takes much longer to read
     - You gunzip the table file (which increases its size to 84MB)
       This makes it much faster to read
   In example1.f you set "FILENAME" correspondingly:
      FILENAME = 'fbt1001.tab.gz'
   or
      FILENAME = 'fbt1001.tab'

   Now you type
      make example1
   This creates your executable  "example1"
   Now, run the example:
      ./example1


* -------------------------------------------------------------
*     About the example
* -------------------------------------------------------------
 
The example "example1.f" is very easy to read:
 - First, LHAPDF is initialized for the CTEQ6.1M PDFs (the central
   result, iset=0).

 - The cross section is computed for all 123 D0 and CDF analysis bins.
   The D0 measurement was made in five regions of pseudorapidity, and 
   the CDF measurement was made in one region.
   The output is ordered as follows:
    - first come all D0 cross sections in the central pseudorapidity 
      bin, in increasing pT; then the other D0 results in increasing
      pseudorapidity.
    - the last numbers are the CDF results at central pseudorapidities

   The first column is the LO result - the second column is the 
   NLO correction - and the last column is the sum of the two, i.e.
   the cross section up to NLO.

   The ouput is made for five different settings of the renormalization 
   and factorization scales
    1: mu/ET = 1/2 (=central choice)
    2: mu/ET = 1/4
    3: mu/ET = 1
    4: mu/ET = 2

 - In the last step we loop over all 40 CTEQ6.1M PDF error sets.
   Here we only compute the cross section for the central scale 
   choice: mu=ET/2 
   After each step the cross section results are printed.

 - The cross sections can also be accessed through the  
   array XST1001(n,iord,iscale) where:
        n: continuous bin number for all D0/CDF bins (from 1 to 123,
           in the same order as they are printed)
     iord: order  1 LO, 2 NLO correction (add both to get the NLO x-section)
   iscale: scale setting for mu_r,mu_f (1-4; see output for values)

note: all cross sections are divided by the bin width in ET and eta
      and the results are presented in pico barn.


* -------------------------------------------------------------
*     Including fastNLO in a PDF fit
* -------------------------------------------------------------

To include the fastNLO code in a PDF fit, the user needs to modify the
interface in  "fn-interface.f" to read his/her own PDFs/alpha_s 
parametrizations.
The code is easy to understand, but be careful that:
 - alpha_s needs to be divided by 2Pi
 - the PDFs are multiplied by the momentum fraction x (i.e.: xg(x))

The cross sections are then computed, as in the example, by calling
   call FB0001CC(FILENAME, 1, 1, XST0001)
 - the table is only read during the first call
 - if only the results for the central scale (=ET/2) are needed, the 
   second argument should be set to Zero (this saves time!)
 - if no ASCII output is desired, the third argument should 
   be set to Zero.
 - the cross section results can be accessed through the array
   XST0001, as described above


* -------------------------------------------------------------
*     Remarks:
* -------------------------------------------------------------

(1)
The CDF measurement is not published in bins. Only differential values
at the bin centers are published. The integration technique in NLOJET++,
however, requires to integrate over bins.
Therefore we have used the bin centers, together with the smooth ansatz 
function, both published by CDF, to determine bins for which the 
bin-averaged cross section agrees better than 0.01% with the differential 
value at the bin center. (this is the reason for the strange values of 
the CDF ET bins)


(2) 
The cone jet algorithms used in Run I were not infrared safe. Therefore
it is not possible to make pQCD predictions for the original algorithms.
A usual choice is to compute the cross sections for the infrared safe
midpoint cone algorithm, but reduce the phase space for the clustering
of partons by adding the Rsep parameter (usually with Rsep=1.3).
This procedure is very questionable, but we have adopted it here
because this was the choice in the original publications and also
in the CTEQ and MRST PDF fits.
A table that was computed without the Rsep parameter is available
on request from the authors.


(3)
The renormalization and factorization scales within each bin 
are proportional to a fixed ET value inside the bin. This value 
is given by:     ET = ETmin  + 0.3* (ETmax-ETmin)
This is close to the bin-center.
Therefore the LO results may differ by less than 1% from the
results obtained using a variable scale, proportional to ET. 
The differences for the NLO results are much smaller.


(4) 
The table is computed to a very high statistical precision
using more than 2000 CPU hours. Please note that a much larger 
time would be needed to reach this precision using the stand-alone 
NLOJET++ code if the internal histogramming package is used.
All NLO cross sections have a statistical precision much better 
than 1 percent. Exact numbers will be given later in the official 
release.


(5)
On the author's desktop PC (not the most recent model) the example
job takes 3:50 minutes.
It takes 25 second to read the table during the first call and
3:25 minutes to do the all computations of the cross sections for
the 41 CTEQ6.1 PDF sets.
Please note that the time for the latter is totally dominated by 
the time to access the PDFs in LHAPDF


(6)
For further information on updates and new observables,
please join the fastNLO users mailing list, as described on 
the fastNLO webpage:
      http://hepforge.cedar.ac.uk/fastnlo/

(7) 
If you have questions, don't hesitate to contact the
authors at   fastnlo@cedar.ac.uk


(8) 
Please remember that this is only a test release and the
results should not yet be used in publications. If you want
to publish results, obtained using the fastNLO code, please
contact the authors before.
We would appreciate if you give us all kinds of feedback:
what you liked, what you would like to be improved, etc.
The first official release will come on a short timescale 
(1-2 months from today).


(9)
The basic concepts for fastNLO are described in a talk given 
at the TeV4LHC meeting in April at CERN:
http://hepforge.cedar.ac.uk/fastnlo/docs/talk-050429tev4lhc.pdf
A publication is in preparation.


September 13, 2005

Thomas Kluge, Klaus Rabbertz, Markus Wobisch