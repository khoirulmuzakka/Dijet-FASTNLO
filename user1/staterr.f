      PROGRAM STATERR
* -------------------------------------------------------------------
* M. Wobisch 12/27/2005
*
* fastNLO - program to compute statistical errors from a large set
*           of LO/NLO tables
*
* -------------------------------------------------------------------
      implicit none
c - Attention!!! - this must be declared consistent with its 
c                  definition in the commonblock!!!!!
      double precision xst1001(286,3,4) 

      CHARACTER*255 FILENAME, LOFILE, NLOFILE
      CHARACTER*255 HISTOFILE
      integer i,j,k,nborn,nnlo
      double precision sum(900,4),sum2(900,4),nevt1(900,4),nevt2(900,4),
     +     mean1(900,4),sigma1(900,4),mean2(900,4),sigma2(900,4),
     +     lmax1(900,4),lmin1(900,4), lmax2(900,4),lmin2(900,4),
     +     nv, NEVTS, val
      CHARACTER*2 NO(0:99)
      DATA NO/'00','01','02','03','04','05','06','07','08','09',
     +     '10','11','12','13','14','15','16','17','18','19',
     +     '20','21','22','23','24','25','26','27','28','29',
     +     '30','31','32','33','34','35','36','37','38','39',
     +     '40','41','42','43','44','45','46','47','48','49',
     +     '50','51','52','53','54','55','56','57','58','59',
     +     '60','61','62','63','64','65','66','67','68','69',
     +     '70','71','72','73','74','75','76','77','78','79',
     +     '80','81','82','83','84','85','86','87','88','89',
     +     '90','91','92','93','94','95','96','97','98','99'/

c --- parse command line - only histo-filename
      IF ( IARGC().LT.1)  HISTOFILE= 'fastnlo-stat.hbk'
      IF ( IARGC().GT.0)  CALL GETARG(1,HISTOFILE)
      IF ( IARGC().GT.1)  THEN
         write(*,*) 'fastNLO: too many arguments given. Stopping'
         RETURN
      ENDIF
      
c - Initialize LHAPDF    - for CTEQ6.1M   
      call InitPDFset('/disk2/work/wobisch/lhapdf-4.1/PDFsets/cteq61.LHgrid')
c
c      call InitPDFset('/h1/h1gen/lhapdf/LHAPDFv4/PDFsets/cteq61.LHgrid')
c - initialize one member, 0=best fit member
      call InitPDF(0)




c - Number of files   <<<<<<<< to be edited by user!!!!!!!
      nborn = 4
      nnlo  = 18

c - >>>> note: the user als needs to edit the filenames 
c              in the LO and in the NLO loops below



      do j=1,900                ! bins
         do k=1,4               ! scales
            nevt1(j,k) = 0d0
            mean1(j,k) = 0d0
            sigma1(j,k)= 0d0
            nevt2(j,k) = 0d0
            mean2(j,k) = 0d0
            sigma2(j,k)= 0d0
            lmax1(j,k) = 0d0
            lmax2(j,k) = 0d0
            lmin1(j,k) = 1d99
            lmin2(j,k) = 1d99
         enddo
      enddo


c ==================================================================
c === for LO ======================================================
c ==================================================================
      write(*,*) '     fastNLO: compute stat errors for LO'

      do j=1,286
         do k =1,4
            sum(j,k)  = 0d0
            sum2(j,k) = 0d0
         enddo
      enddo

c -   loop over files
      do i=0,nborn
         write(*,*) ' ###############################################'
         write(*,*) ' ############  next table  #####################'
         write(*,*) ' ###############################################'

         filename = 'path/fnt1001/fn01ml-born-'//NO(i)//'.stc' ! <<<<< user!!!
         write(*,*) filename,'<<<<<'
         call FX9999CC(FILENAME, -1.0d0 , 1, 0 , XST1001)

         do j=1,286
            do k =1,4
               val = xst1001(j,1,k)
               if (val.lt.lmin1(j,k)) lmin1(j,k)=val
               if (val.gt.lmax1(j,k)) lmax1(j,k)=val
               sum(j,k)  = sum(j,k)  + val 
               sum2(j,k) = sum2(j,k) + val*val
               nv = NEVTS(1)
               nevt1(j,k)  = nevt1(j,k) + nv
            enddo
         enddo
         write (*,*) " >>>>> N=",nv
      enddo

c -   extract mean values and standard deviations
      do k=1,4                  ! scales
         write(*,*) ' next scale: ',k
         do j=1,123             ! 286  ! 900
            mean1(j,k) = sum(j,k) / (nborn+1)
            sigma1(j,k)= sqrt((sum2(j,k)-(sum(j,k)*sum(j,k))/(nborn+1))/nborn)
            write(*,900) j, mean1(j,k),
     +           100d0*sigma1(j,k)/mean1(j,k),
     +           100d0*(lmin1(j,k)-mean1(j,k))/mean1(j,k),
     +           100d0*(lmax1(j,k)-mean1(j,k))/mean1(j,k)
         enddo
      enddo



c ==================================================================
c === for NLO ======================================================
c ==================================================================
      write(*,*)
      write(*,*) ' ########################################################'
      write(*,*) ' ########################################################'
      write(*,*) ' ########################################################'
      write(*,*) ' ########################################################'
      write(*,*) ' ########################################################'
      write(*,*) ' ########################################################'
      write(*,*) '     fastNLO: compute stat errors for NLO'

      do j=1,286
         do k =1,4
            sum(j,k)  = 0d0
            sum2(j,k) = 0d0
         enddo
      enddo

c -   loop over files
      do i=0,nnlo
         write(*,*) ' ###############################################'
         write(*,*) ' ############  next table  #####################'
         write(*,*) ' ###############################################'

         filename = 'path/fnt1001/fn01mn-nlo-'//NO(i)//'.stc' ! <<<<< user!!!
         write(*,*) filename,'<<<<<'
         call FX9999CC(FILENAME, -1.0d0 , 1, 0 , XST1001)

         do j=1,286
            do k =1,4
               val = xst1001(j,1,k)+xst1001(j,2,k)
               if (val.lt.lmin2(j,k)) lmin2(j,k)=val
               if (val.gt.lmax2(j,k)) lmax2(j,k)=val
               sum(j,k)  = sum(j,k)  + val 
               sum2(j,k) = sum2(j,k) + val*val
               nv = NEVTS(1)
               nevt2(j,k)  = nevt2(j,k) + nv
            enddo
         enddo
         write (*,*) " >>>>> N=",nv
      enddo

c -   extract mean values and standard deviations
      do k=1,4                  ! scales
         write(*,*) ' next scale: ',k
         do j=1,123             ! 286  ! 900
            mean2(j,k) = sum(j,k) / (nnlo+1)
            sigma2(j,k)= sqrt((sum2(j,k)-(sum(j,k)*sum(j,k))/(nnlo+1))/nnlo)
            write(*,900) j, mean2(j,k),
     +           100d0*sigma2(j,k)/mean2(j,k),
     +           100d0*(lmin2(j,k)-mean2(j,k))/mean2(j,k),
     +           100d0*(lmax2(j,k)-mean2(j,k))/mean2(j,k)
         enddo
      enddo


 900   FORMAT (I4,F13.5,"  in %:",3F8.3)
       write(*,*) 'bin       mean          sigma     min     max'


c - book, fill and store histograms - works now! (Sept 15, 2005 MW)
c - >> works not yet for threshold corrections
c      call FNHIST(HISTOFILE,1)  ! to plot x-sect histos

      END


C-------------------------------------------------------------------
      Double Precision Function NEVTS(i)
      implicit none
      integer i
      INCLUDE 'fnx9999.inc'

      NEVTS=dble(NEVT(i))

      RETURN
      END
