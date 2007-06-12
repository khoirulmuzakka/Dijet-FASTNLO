      PROGRAM EXAMPLE02
* -------------------------------------------------------------------
* M. Wobisch 08/30/2005
*
* fastNLO - example program to compute statitical errors
*           from a large number of result tables
*
* the program requires to have 
*   - the total result stored in table   path/fnt1001.tab
*   - 5 LO tables to be stored as        path/fnt1001-[ij].tab
*                  where [ij] is 'lo01',...,'lo05'
*   - 20 NLO tables to be stored as      path/fnt1001-[ij].tab
*                  where [ij] is 'nlo01',...,'nlo20'
*
* -------------------------------------------------------------------
      implicit none
      CHARACTER*4 LOTAB(5)
      CHARACTER*5 NLOTAB(20)
      CHARACTER*40 filename
      double precision xst1001(123,2,5),reference(123,2,5), val
      double precision lmean(123,5),lrms(123,5),lmax(123,5),lmin(123,5)
      double precision nmean(123,5),nrms(123,5),nmax(123,5),nmin(123,5)
      integer imufflag,  numlo,numnlo, i,j,k
      DATA LOTAB/'lo01','lo02','lo03','lo04','lo05'/
      DATA NLOTAB/'nlo01','nlo02','nlo03','nlo04','nlo05','nlo06','nlo07',
     +     'nlo08','nlo09','nlo10','nlo11','nlo12','nlo13','nlo14','nlo15',
     +     'nlo16','nlo17','nlo18','nlo19','nlo20'/

      numlo =5
      numnlo=20
      imufflag = 0              ! 1: vary muf

c - Initialize LHAPDF    - for CTEQ6.1M   
      call InitPDFset('/disk2/work/wobisch/LHAPDFv4/PDFsets/cteq61.LHgrid')
c - initialize one member, 0=best fit member
      call InitPDF(0)


c - compute the central cross sections (=reference values)
      write(*,*) '     fastNLO: compute the reference values'
      call FNCALC('path/fnt1001.tab', imufflag, reference)




c ------- determine stat. errors
      write(*,*) '     fastNLO: compute the statistical fluctuations'
c - LO
      write(*,*)
      do k=1,5 
         do j=1,123
            lmean(j,k)= 0d0
            lrms(j,k) = 0d0
            lmax(j,k) = -9d40
            lmin(j,k) = 9d40
         enddo
      enddo
      do i=1,numlo
         filename='path/fnt1001-'//LOTAB(i)//'.tab'
         write(*,*) ' >>',filename,'<<  '
         call FNCALC(filename, imufflag, XST1001)
         do j=1,123
            do k=1,5
               val = XST1001(j,1,k)
               lmean(j,k)=lmean(j,k) + val/numlo
               lrms(j,k)=lrms(j,k) + (val**2)/numlo
               if (val.gt.lmax(j,k)) lmax(j,k) = val
               if (val.lt.lmin(j,k)) lmin(j,k) = val
            enddo
         enddo
      enddo

      write(*,*) '  --------------------------------------------'

      do k=1,1                  !,5
         write(*,*)
         write(*,*) ' ----------- LO:  scale variation ',k
         write(*,*) '      everything *1000 -> in per mil !!!!!'
         write(*,*) '   rms/tot (rms/tot)/sqrt(n) (highest-tot)/tot',
     +        ' (smallest-tot)/tot'
         do j=1,123
            val = reference(j,1,k)
            write(*,*) j,'  ',
     +           1000d0*sqrt(lrms(j,k)-lmean(j,k)**2)/val ,'  ',
     +           1000d0*sqrt(lrms(j,k)-lmean(j,k)**2)/val/sqrt(dble(numlo))
     +           ,'  ',1000d0*(lmax(j,k)-val)/val
     +           ,'  ',1000d0*(lmin(j,k)-val)/val
         enddo
      enddo



c - NLO
      write(*,*)
      do k=1,5 
         do j=1,123
            lmean(j,k)= 0d0
            lrms(j,k) = 0d0
            lmax(j,k) = -9d40
            lmin(j,k) = 9d40
         enddo
      enddo
      do i=1,numnlo
         filename='path/fnt1001-'//NLOTAB(i)//'.tab'
         write(*,*) ' >>',filename,'<<  '
         call FNCALC(filename, imufflag, XST1001)
         do j=1,123
            do k=1,5
               val = XST1001(j,1,k)+XST1001(j,2,k)
               lmean(j,k)=lmean(j,k) + val/numnlo
               lrms(j,k)=lrms(j,k) + (val**2)/numnlo
               if (val.gt.lmax(j,k)) lmax(j,k) = val
               if (val.lt.lmin(j,k)) lmin(j,k) = val
            enddo
         enddo
      enddo

      do k=1,1                  !,5
         write(*,*)
         write(*,*) ' ----------- NLO:  scale variation ',k
         write(*,*) '      everything *1000 -> in per mil !!!!!'
         write(*,*) '   rms/tot (rms/tot)/sqrt(n) (highest-tot)/tot',
     +        ' (smallest-tot)/tot'
         do j=1,123
            val = reference(j,1,k)+reference(j,2,k)
            write(*,*) j,'  ',
     +           1000d0*sqrt(lrms(j,k)-lmean(j,k)**2)/val ,'  ',
     +           1000d0*sqrt(lrms(j,k)-lmean(j,k)**2)/val/sqrt(dble(numnlo))
     +           ,'  ',1000d0*(lmax(j,k)-val)/val
     +           ,'  ',1000d0*(lmin(j,k)-val)/val
         enddo
      enddo





      RETURN
      END
