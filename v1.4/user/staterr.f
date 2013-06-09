      PROGRAM STATERR
* -------------------------------------------------------------------
* M. Wobisch 12/27/2005
*
* fastNLO - program to compute statistical errors from a large set
*           of LO/NLO tables
*
* the user needs to edit four lines
*       - No of available tables in LO
*       - No of available tables in NNLO
*       - name of LO table files
*       - name of NLO table files
*
* MW Jan 24, 2006  compute weighted mean,rms according to No. events
* MW Feb 1, 2006  rewrite for new usercode v14a with single scale
*                 reorganzize to allow simple loops over many scenarios
* -------------------------------------------------------------------
      implicit none

c - Initialize LHAPDF    - for CTEQ6.1M
      call InitPDFset('/work/shootingstar-clued0/wobisch/lhapdf500/share/lhapdf/PDFsets/cteq61.LHgrid')
c      call InitPDFset('/disk2/work/wobisch/lhapdf-4.1/PDFsets/cteq61.LHgrid')
c      call InitPDFset('/h1/h1gen/lhapdf/LHAPDFv4/PDFsets/cteq61.LHgrid')
c - initialize one member, 0=best fit member
      call InitPDF(0)


c --- call statistical error-code for each scenario
c
c     arguments:
c      - No of Born tables minus 1 (start tables from 00)
c      - name of Born tables  (+ 2digit number + '.stc')
c      - No of NLO tables minus 1 (start tables from 00)
c      - name of NLO tables  (+ 2digit number + '.stc')
c      - name of histogram file for output

      goto 544

      CALL STATCODE(   4, 'path/fnt1001/ft11ml-born-',
     +     12, 'path/fnt1001/ft11mn-nlo-',     'fnt1001-stat.hbk')


      CALL STATCODE(4, 'path/fnt1003/ft13ml-born-',
     +     19, 'path/fnt1003/ft13mn-nlo-',     'fnt1003-stat.hbk')



 545  continue
      CALL STATCODE(7, 'path/fnt2003plus/ft23pl-born-',
     +     25, 'path/fnt2003plus/ft23pn-nlo-',     'fnt2003plus-stat.hbk')

 546  continue
      CALL STATCODE(7, 'path/fnt2003plus/ft23pl-born-',
     +     25, 'path/fnt2003plus/ft23p5n-nlo-',     'fnt2003plusD05-stat.hbk')

 547  continue
      CALL STATCODE(7, 'path/fnt2003plus/ft23pl-born-',
     +     25, 'path/fnt2003plus/ft23p1n-nlo-',     'fnt2003plusD10-stat.hbk')

      CALL STATCODE(9, 'path/fnt1007/ft17rl-born-',
     +     77, 'path/fnt1007/ft17rn-nlo-',     'fnt1007rsep-stat.hbk')
c          78/79 are copied to 11,11


      CALL STATCODE(4, 'path/fnt1002prod/ft12ml-born-',
     +     59, 'path/fnt1002prod/ft12mn-nlo-',     'fnt1002-stat.hbk')

 543  continue
      CALL STATCODE(19, 'path/fnt200a/ft2aml-born-',
     +     43, 'path/fnt200a/ft2amn-nlo-',     'fnt200a-stat.hbk')



      CALL STATCODE(4, 'path/fnt1004/ft14ml-born-',
     +     19, 'path/fnt1004/ft14mn-nlo-',     'fnt1004midp-stat.hbk')
      CALL STATCODE(4, 'path/fnt1004/ft14rl-born-',
     +     19, 'path/fnt1004/ft14rn-nlo-',     'fnt1004rsep-stat.hbk')

      CALL STATCODE(4, 'path/fnt1004/ft15ml-born-',
     +     19, 'path/fnt1004/ft15mn-nlo-',     'fnt1005midp-stat.hbk')
      CALL STATCODE(4, 'path/fnt1004/ft15rl-born-',
     +     19, 'path/fnt1004/ft15rn-nlo-',     'fnt1005rsep-stat.hbk')

      CALL STATCODE(9, 'path/fnt1007/ft17rl-born-',
     +     58, 'path/fnt1007/ft17mn-nlo-',     'fnt1007midp-stat.hbk')

      CALL STATCODE(5, 'path/fnt1008/ft18ml-born-',
     +     59, 'path/fnt1008/ft18rn-nlo-',     'fnt1008rsep-stat.hbk')

      CALL STATCODE(5, 'path/fnt1008/ft18ml-born-',
     +     59, 'path/fnt1008/ft18mn-nlo-',     'fnt1008midp-stat.hbk')

      CALL STATCODE(11, 'path/fnl0001/fl01kl-born-',
     +     136, 'path/fnl0001/fl01kn-nlo-',     'fnl0001-stat.hbk')

      CALL STATCODE(9, 'path/fnr0001/fr01l-born-',
     +     39, 'path/fnr0001/fr01n-nlo-',     'fnr0001-stat.hbk')

 544  continue


c ->
      CALL STATCODE(29, 'path/fnt2001/ft21dl-born-',
     +     282, 'path/fnt2001/ft21dn-nlo-',     'fnt2001-stat.hbk')




c---------------------------------------

c      CALL STATCODE(9, 'path/fnl000a/fl0aml-born-',
c     +     21, 'path/fnl000a/fl0amn-nlo-',     'fnl000a-stat.hbk')

      write(*,*)
      write(*,*)"now move histogram files:"
      write(*,*)"mv fn*stat*k /work/joker-clued0/wobisch/fastNLO/hbook-v1c/stat"
      RETURN
      END

C -----------------------------------------------------------------

      SUBROUTINE STATCODE(nborn,cborntab,nnlo,cnlotab,chist)
      implicit none
c - Attention!!! - this must be declared consistent with its
c                  definition in the commonblock!!!!!
      double precision xsect(900,3)

      CHARACTER*(*) cborntab
      CHARACTER*(*) cnlotab
      CHARACTER*(*) chist

      CHARACTER*255 FILENAME, LOFILE, NLOFILE
      CHARACTER*255 HISTOFILE
      integer i,j,k,nborn,nnlo,  ntotal, nmax, ICYCLE, istat2
      double precision sum(900,4),sum2(900,4),nevt1(900,4),nevt2(900,4),
     +     mean1(900,4),sigma1(900,4),mean2(900,4),sigma2(900,4),
     +     lmax1(900,4),lmin1(900,4), lmax2(900,4),lmin2(900,4),
     +     NEVTS, val,  nv, nvtot
      Integer NJmin(900),NJmax(900)
      CHARACTER*2 NO(0:289)
      Double Precision mu(4)
      DATA mu/0.25d0, 0.5d0, 1.0d0, 2.0d0/

c - HBOOK common
      INTEGER NWPAWC
      PARAMETER (NWPAWC=2500000)
      REAL HMEMOR(NWPAWC)
      COMMON /PAWC/ HMEMOR

      DATA NJMIN/900*-1/, NJMAX/900*-1/

      DATA NO/'00','01','02','03','04','05','06','07','08','09',
     +     '10','11','12','13','14','15','16','17','18','19',
     +     '20','21','22','23','24','25','26','27','28','29',
     +     '30','31','32','33','34','35','36','37','38','39',
     +     '40','41','42','43','44','45','46','47','48','49',
     +     '50','51','52','53','54','55','56','57','58','59',
     +     '60','61','62','63','64','65','66','67','68','69',
     +     '70','71','72','73','74','75','76','77','78','79',
     +     '80','81','82','83','84','85','86','87','88','89',
     +     '90','91','92','93','94','95','96','97','98','99',
     +     'a0','a1','a2','a3','a4','a5','a6','a7','a8','a9',
     +     'b0','b1','b2','b3','b4','b5','b6','b7','b8','b9',
     +     'c0','c1','c2','c3','c4','c5','c6','c7','c8','c9',
     +     'd0','d1','d2','d3','d4','d5','d6','d7','d8','d9',
     +     'e0','e1','e2','e3','e4','e5','e6','e7','e8','e9',
     +     'f0','f1','f2','f3','f4','f5','f6','f7','f8','f9',
     +     'g0','g1','g2','g3','g4','g5','g6','g7','g8','g9',
     +     'h0','h1','h2','h3','h4','h5','h6','h7','h8','h9',
     +     'i0','i1','i2','i3','i4','i5','i6','i7','i8','i9',
     +     'j0','j1','j2','j3','j4','j5','j6','j7','j8','j9',   ! -> 199
     +     'k0','k1','k2','k3','k4','k5','k6','k7','k8','k9',
     +     'l0','l1','l2','l3','l4','l5','l6','l7','l8','l9',
     +     'm0','m1','m2','m3','m4','m5','m6','m7','m8','m9',
     +     'n0','n1','n2','n3','n4','n5','n6','n7','n8','n9',
     +     'o0','o1','o2','o3','o4','o5','o6','o7','o8','o9',
     +     'p0','p1','p2','p3','p4','p5','p6','p7','p8','p9',   ! -> 259
     +     'q0','q1','q2','q3','q4','q5','q6','q7','q8','q9',
     +     'r0','r1','r2','r3','r4','r5','r6','r7','r8','r9',   ! -> 279
     +     's0','s1','s2','s3','s4','s5','s6','s7','s8','s9',   ! -> 289
     +     't0','t1','t2','t3','t4','t5','t6','t7','t8','t9'    ! -> 299
     +     /

      HISTOFILE= chist
      do j=1,900                ! bins
         do k=1,4               ! scales
            nevt1(j,k) = 0d0
            mean1(j,k) = 0d0
            sigma1(j,k)= 0d0
            nevt2(j,k) = 0d0
            mean2(j,k) = 0d0
            sigma2(j,k)= 0d0
            lmax1(j,k) = -1d99
            lmax2(j,k) = -1d99
            lmin1(j,k) = 1d99
            lmin2(j,k) = 1d99
         enddo
      enddo

c ==================================================================
c === for LO ======================================================
c ==================================================================
      write(*,*) '     fastNLO: compute stat errors for LO'

      do j=1,900
         do k =1,4
            sum(j,k)  = 0d0
            sum2(j,k) = 0d0
         enddo
      enddo

c -   loop over files
      nvtot = 0d0
      do i=0,nborn
         write(*,*) ' ###############################################'
         write(*,*) ' ############  next table  #####################'
         write(*,*) ' ###############################################'

         filename = cborntab//NO(i)//'.stc'
         write(*,*) filename,'<<<<<'

c - loop over scales
         do k=1,4
            call FX9999CC(FILENAME, mu(k) , mu(k), 0 , XSECT)
            if (i.eq.0 .and. k.eq.1) then
               ntotal = nmax(1)
               write(*,*) ' this file has ',ntotal,' bins'
            endif

            nv = NEVTS(1)/1000d0
            if (k.eq.1) nvtot = nvtot + nv
            do j=1,ntotal
               nevt1(j,k)  = nevt1(j,k) + nv ! ?? don't need ??
               val = xsect(j,1)
               if (val.lt.lmin1(j,k)) lmin1(j,k)=val
               if (val.gt.lmax1(j,k)) lmax1(j,k)=val
               sum(j,k)  = sum(j,k)  + val * nv
               sum2(j,k) = sum2(j,k) + val*val * nv
            enddo
         enddo
         write (*,*) " >>>>> N=",nv
      enddo

c -   extract mean values and standard deviations
      do k=1,4
         write(*,*) ' next scale: ',k,"    .... N=",nv
         do j=1,ntotal
            mean1(j,k) = sum(j,k) / nvtot
            sigma1(j,k)= sqrt((sum2(j,k)-(sum(j,k)*sum(j,k)/nvtot))/nvtot)
            sigma1(j,k)= sigma1(j,k)* sqrt(1d0/nborn)
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

      do j=1,ntotal
         do k =1,4
            sum(j,k)  = 0d0
            sum2(j,k) = 0d0
         enddo
      enddo

c -   loop over files
      nvtot = 0d0
      do i=0,nnlo
         write(*,*) ' ###############################################'
         write(*,*) ' ############  next table  #####################'
         write(*,*) ' ###############################################'

c         filename = 'path/fnt1001/ft11mn-nlo-'//NO(i)//'.stc' !<<<< user!!
         filename = cnlotab//NO(i)//'.stc' ! <<<user!!!
         write(*,*) filename,'<<<<<'

c - loop over scales
         do k=1,4
            call FX9999CC(FILENAME, mu(k), mu(k), 0 , XSECT)

            nv = NEVTS(2)/1000d0
            if (k.eq.1) nvtot = nvtot + nv
            do j=1,ntotal
               nevt2(j,k)  = nevt2(j,k) + nv ! ?? dont need this ??
               val = xsect(j,1)+xsect(j,2)
               if (val.lt.lmin2(j,k)) then
                  lmin2(j,k)=val
                  NJmin(j) = i
               endif
               if (val.gt.lmax2(j,k)) then
                  lmax2(j,k)=val
                  NJmax(j) = i
               endif
               sum(j,k)  = sum(j,k)  + val * nv
               sum2(j,k) = sum2(j,k) + val*val * nv
            enddo
         enddo
         write (*,*) " >>>>> N=",nv
      enddo

c -   extract mean values and standard deviations
      do k=1,4                  ! scales
         write(*,*) ' next scale: ',k
         do j=1,ntotal
            mean2(j,k) = sum(j,k) / nvtot
            sigma2(j,k)= sqrt((sum2(j,k)-(sum(j,k)*sum(j,k)/nvtot))/nvtot)
            sigma2(j,k)= sigma2(j,k)* sqrt(1d0/nnlo)
            write(*,901) j, mean2(j,k),
     +           100d0*sigma2(j,k)/mean2(j,k),
     +           100d0*(lmin2(j,k)-mean2(j,k))/mean2(j,k),
     +           100d0*(lmax2(j,k)-mean2(j,k))/mean2(j,k),
     +           NJmin(j),NJmax(j)
         enddo
      enddo


 900   FORMAT (I4,F16.5,"  in %:",3F10.3)
 901   FORMAT (I4,F16.5,"  in %:",3F9.2,2X,2I5)
       write(*,*) 'bin       mean          sigma     min     max'


c ====== book, fill and store histograms for stat errors ================

       WRITE(*,*) ' ------------------ histo filling -------------'

c - open & book
       CALL HLIMIT(NWPAWC)
       CALL HROPEN(11,'fastNLO',histofile,'N',1024,ISTAT2)
       if (ISTAT2.NE.0) then
          WRITE(*,*) ' FNHBOOK: could not open histofile ',istat2
       endif

c  -  LO only at single scale -  NLO at four scales
       call hbook1(101,'bin number',ntotal,0.5,real(ntotal+0.5),0)
       do i=1,4  ! scales
          call hbook1(200+i,'bin number',ntotal,0.5,real(ntotal+0.5),0)
       enddo

       do j=1,ntotal
          call hfill(101,real(j),0.0,real(100d0*sigma1(j,1)/mean1(j,1)))
          do k=1,4
             call hfill(200+k, real(j) ,0.0,
     +            real(100d0*sigma2(j,k)/mean2(j,k)))
          enddo
       enddo
       write(*,*) ' store in histos 101,201-204 in: ',histofile
       CALL HROUT (0,ICYCLE,' ')
       CALL HREND ('fastNLO')

       END


C-------------------------------------------------------------------
      Double Precision Function NEVTS(i)
      implicit none
      integer i
      INCLUDE 'fnx9999.inc'

      NEVTS=dble(NEVT(i))

      RETURN
      END
C-------------------------------------------------------------------
      INTEGER Function NMAX(i)
      implicit none
      integer i,j1,j2, ntotal
      INCLUDE 'fnx9999.inc'

      ntotal = 0
      do j1=1,nrapidity          ! Rapidity Bins
         do j2=1,NPT(j1)          ! pT Bins
            ntotal = ntotal + 1
         enddo
      enddo
      nmax= ntotal

      RETURN
      END
