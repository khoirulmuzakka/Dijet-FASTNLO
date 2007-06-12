**********************************************************************
**********************************************************************
**  
**  09/19/2001 M. Wobisch
**
**  code for the calculation of NLO predictions 
**  for H1 jet cross sections from convolution
**  of arbitrary PDFs and alpha_s with precalculated 
**  coefficient functions (from DISENT), stored
**  in external grids.
**
* routines:  
*
* TSH0001  user interface - returns two arrays of NLO cross sections
* MW_XS_DIJET   computes the dijet cross section
* MW_XS_INCJET  computes the inclusive jet cross section
* MW_GRID_DIJ   reads the pre-computed perturbative coefficients
*               for the dijet cross section
* MW_GRID_INC   reads the pre-computed perturbative coefficients
*               for the inclusive jet cross section
**********************************************************************
**********************************************************************

      SUBROUTINE TSH0001(XINCL,XDIJET)
**********************************************************************
* 09/19/2001  M. Wobisch
*   user interface
*
* input:
*
* output (double precision):
*             xincl(4,4)  inclusive jet cross section at NLO
*            xdijet(4,5)  dijet cross section at NLO
*                  for both: first index is Q2 bin number
*                            second index is ET or Xi bin number
**********************************************************************
      IMPLICIT NONE
      double precision xincl(4,4) , xdijet(4,5)
      INTEGER IFIRST
      SAVE IFIRST
      DATA IFIRST/0/ 

      IF (IFIRST .eq. 0) THEN
         IFIRST = 1
         WRITE(*,*) ' +-------------------------------------------------'  
         WRITE(*,*) ' | fastNLO        - scenario FNH0001'
         WRITE(*,*) ' |   compute the NLO pQCD predictions for:' 
         WRITE(*,*) ' |      H1 Collaboration, C. Adloff et al.,'
         WRITE(*,*) ' |               Eur. Phys. J. C19 (2001) 289.'
         WRITE(*,*) ' |                            (hep-ex/0010054)'
         WRITE(*,*) ' |   implemented are those bins for the inclusive'
         WRITE(*,*) ' |      jet and dijet cross sections that have been'
         WRITE(*,*) ' |      used in the gluon fit in the H1 publication'
         WRITE(*,*) ' |      (section 5.4, Fig 18)'
         WRITE(*,*) ' |   the computation is made for'
         WRITE(*,*) ' |      mu_r=ET'
         WRITE(*,*) ' |      mu_f=sqrt(200)GeV <- the average ET of the jets'
         WRITE(*,*) ' |   the scale depencence and the hadronization'
         WRITE(*,*) ' |      corrections are published in the H1 paper'
         WRITE(*,*) ' |   the pQCD coefficients are computed using DISENT'
         CALL MW_GRID_DIJ
         CALL MW_GRID_INC 
      ENDIF

      CALL MW_XS_DIJET(xdijet)
      CALL MW_XS_INCJET(xincl)

      RETURN 
      END

C -------------------------------------------------------------------

      SUBROUTINE MW_XS_DIJET(xdijet)
      IMPLICIT NONE
      DOUBLE PRECISION  xdijet(4,5)
      DOUBLE PRECISION grid_dij(4,5,3,4,4,55)

      DOUBLE PRECISION OBS(6,41)
      DOUBLE PRECISION ALPS, G_HI_BINCENT(4)
      INTEGER  IDX,IMUR,IOAS,IFLAV,IBIN,IOBS,  I,J 
      DOUBLE PRECISION DX ,MUF,Q0(4),PDF(3),AS(4), XCONT
      DOUBLE PRECISION pi, FNALPHAS
      DOUBLE PRECISION DXPDF(-6:6)
      INTEGER IFIRST
      SAVE IFIRST, Q0
      DATA IFIRST/0/
      COMMON /CH0001b/ grid_dij

      IF (IFIRST .eq. 0) THEN
         IFIRST = 1
         pi = 3.14159265358979323846d0
         G_HI_BINCENT(1) = 100D0 
         G_HI_BINCENT(2) = 170D0
         G_HI_BINCENT(3) = 300D0
         G_HI_BINCENT(4) = 700D0
c --- get mean values of Mu_r regions
         DO I=1,4
            Q0(i) = SQRT(G_HI_BINCENT(i))
c            WRITE(*,*) '  *  HDJ ren. scale bin-center ',i,q0(i)**2
         ENDDO
      ENDIF

c --- reset observable-array
      DO I=1,5
         DO J=1,4
            xdijet(j,i) = 0d0
         ENDDO
      ENDDO

c --- start calculation of cross section
      DO IDX = 2,55             ! loop over x values - read pdf
         DX = 10D0** ( DBLE(1-IDX) / 25D0 )
         call FNEPPDF(dx,sqrt(200d0), PDF)   ! 1g  2sigma  3delta

         DO IMUR = 1,4 ! loop over Mu_r ranges - get alpha_s
            AS(1) = FNALPHAS(Q0(IMUR))*2d0*Pi
            AS(2) = AS(1) ** 2 
            AS(3) = AS(1) ** 3 
            AS(4) = AS(1) ** 4 

            DO IOAS = 1,4
               DO IFLAV = 1,3
                  DO IOBS = 1,4 ! Q2 bin
                     DO IBIN = 1,5 ! Xi bin
                        XCONT =  AS(IOAS) * PDF(IFLAV) * grid_dij(
     +                       IOBS,IBIN,IFLAV,IOAS,IMUR,IDX)
                        xdijet(IOBS,IBIN) = XDIJET(IOBS,IBIN) + XCONT
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END

C --------------------------------------------------------------------

      SUBROUTINE MW_XS_INCJET(xincl)
      IMPLICIT NONE
      DOUBLE PRECISION  xincl(4,4), grid_inc(4,4,3,4,75)
      INTEGER IDX,IMUR,IOAS,IFLAV,IBIN,IOBS,  I,J
      DOUBLE PRECISION DX,MUF,Q0(4),SQ0(4),PDF(3), AS(4), XCONT
      DOUBLE PRECISION pi, FNALPHAS
      INTEGER IFIRST
      COMMON /CH0001a/ grid_inc
      SAVE IFIRST, Q0, SQ0
      DATA Q0/  70D0, 200D0, 500D0 , 1500D0 /
      DATA IFIRST/0/

      IF (IFIRST .eq. 0) THEN
         IFIRST = 1
         pi = 3.14159265358979323846d0
         DO I=1,4
            SQ0(i) = SQRT(Q0(i))
         ENDDO
      ENDIF

c --- reset observable-array
      DO I=1,4                  ! Q2 bin
         DO J=1,4               ! ET bin
            xincl(I,J) = 0D0
         ENDDO
      ENDDO

c --- start calculation of cross section
      DO IDX = 2,75             ! loop over x values - read pdf
         DX = 10D0** ( DBLE(1-IDX) / 25D0 )
         call FNEPPDF(dx,sqrt(200d0), PDF)   ! 1g  2sigma  3delta

         DO IMUR = 1,4          ! loop over Mu_r ranges - get alpha_s
                                ! for inclusive jets:  imur = ibin (=ETbin)
            IBIN = IMUR
            AS(1) = FNALPHAS(SQ0(IMUR))*2d0*Pi
            AS(2) = AS(1) ** 2 
            AS(3) = AS(1) ** 3 
            AS(4) = AS(1) ** 4 
            DO IOAS = 1,4
               DO IFLAV = 1,3
                  DO IOBS = 1,4 ! Q2 bin
                     XCONT =  AS(IOAS) * PDF(IFLAV) * grid_inc(
     +                    IOBS,IBIN,IFLAV,IOAS,IDX)
                     xincl(IOBS,IBIN) = xincl(IOBS,IBIN) + XCONT
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
C --------------------------------------------------------------------

      SUBROUTINE MW_GRID_DIJ
*******************************************************************
* read coefficient grid for high Q2 dijet data from file in /newgrid
*
*******************************************************************
      IMPLICIT NONE
      INTEGER IFIRST, IFILE, I        ,J,K,L,M,N
      DOUBLE PRECISION grid_dij(4,5,3,4,4,55)

      DATA IFIRST/0/
      SAVE IFIRST
      COMMON /CH0001b/ grid_dij

      IF (IFIRST.EQ.1) RETURN
      IFIRST = 1

      OPEN(2,STATUS='OLD',FILE='path/fnh0001dij.dat' ,IOSTAT=IFILE)

      IF (IFILE .ne. 0) THEN
         WRITE(*,*) ' NLOH1JET:  file fnh0001dij.dat not found',
     +        '  -  IOSTAT = ',IFILE
         STOP
      ENDIF

      READ(2,'(4(1X,D17.10))') GRID_DIJ
      CLOSE(2)

      RETURN
      END
C------------------------------------------------------------------

      SUBROUTINE MW_GRID_INC
*******************************************************************
* read coefficient grid for high Q2 incl. jet data from file in /newgrid
*
*******************************************************************
      IMPLICIT NONE
      INTEGER IFIRST, IFILE, I
      DOUBLE PRECISION grid_inc(4,4,3,4,75)
      DATA IFIRST/0/
      SAVE IFIRST
      COMMON /CH0001a/ grid_inc

      IF (IFIRST.EQ.1) RETURN
      IFIRST = 1

      OPEN(2,STATUS='OLD',FILE='path/fnh0001inc.dat' ,IOSTAT=IFILE)
      IF (IFILE .ne. 0) THEN
         WRITE(*,*) ' NLOH1JET:  file fnh0001inc.dat not found',
     +        '  -  IOSTAT = ',IFILE
         STOP
      ENDIF

      DO I=1,13
         READ(2,*)
      ENDDO

      READ(2,'(4(1X,D17.10))') grid_inc
      CLOSE(2)

      RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE TSH0001out(XINCL,XDIJET)
**********************************************************************
* 07/26/2005  M. Wobisch
*
*  write output
**********************************************************************
      IMPLICIT NONE
      double precision xincl(4,4) , xdijet(4,5)
      INTEGER I,J

      write (*,*) ' '
      write (*,*) ' # fastNLO - results for: H1 DIS jets'
      write (*,*) '   FNH0001: H1 inclusive jet cross section (Q2,ET)'
      do i=1,4                  ! Q2 Bins
         do j=1,4               ! ET Bins
            write(*,*) '    Q2 range ',i,'  ET bin ',j,' - ',xincl(i,j)
         enddo
      enddo
      write (*,*) '   FNH0001:  H1 dijet cross section (Q2,Xi)'
      do i=1,4                  ! Q2 Bins
         do j=1,5               ! Xi Bins
            write(*,*) '    Q2 range ',i,'  Xi bin ',j,' - ',xdijet(i,j)
         enddo
      enddo
      write (*,*) '             -------------------'

      RETURN 
      END
