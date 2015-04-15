      INTEGER FUNCTION IDim0Bin(iObsBin)
*     Returns bin number in first dimension
*     Valid for up to triple differential binnings
*     There always must be at least one bin!
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INTEGER iObsBin,i0bin,iobs
      DOUBLE PRECISION lo0bin

      if ( NOBSBIN.EQ.0 ) then
         WRITE(*,*)'ALLUNC:GetIDim0Bin: ERROR! '//
     >        "No observable bins defined, aborted!"
         STOP
      endif
      if ( iObsBin.GT.NObsBin ) then
         WRITE(*,*)'ALLUNC:GetIDim0Bin: ERROR! '//
     >        "Observable bin out of range, aborted!"
         STOP
      endif

      i0bin = 1
      lo0bin = LOBIN(1,1)
      IDim0Bin = 0
      DO iobs=1,NOBSBIN
         if ( lo0bin.LT.LOBIN(iobs,1) ) then
            lo0bin = LOBIN(iobs,1)
            i0bin = i0bin + 1;
         endif
         if ( iobs.EQ.iObsBin ) then
            IDim0Bin = i0bin
            return
         endif
      enddo
      WRITE(*,*)'ALLUNC:GetIDim0Bin: ERROR! '//
     >     'Observable bin not found. This should never happen, '//
     >     'aborted!'
      STOP

      END FUNCTION



      INTEGER FUNCTION IDim1Bin(iObsBin)
* Returns bin number in second dimension
* Valid for up to triple differential binnings
* 1d binning --> error exit
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INTEGER iObsBin,i0bin,i1bin,iobs
      DOUBLE PRECISION lo0bin,lo1bin
      if ( NDIM < 2 ) then
         WRITE(*,*)'ALLUNC:GetIDim1Bin: ERROR! '//
     >        "No second dimension available, aborted!"
         STOP
      endif
      if ( NOBSBIN.EQ.0 ) then
         WRITE(*,*)'ALLUNC:GetIDim1Bin: ERROR! '//
     >        "No observable bins defined, aborted!"
         STOP
      endif
      if ( iObsBin.GT.NObsBin ) then
         WRITE(*,*)'ALLUNC:GetIDim1Bin: ERROR! '//
     >        "Observable bin out of range, aborted!"
         STOP
      endif

      i0bin = 1
      i1bin = 1
      lo0bin = LOBIN(1,1)
      lo1bin = LOBIN(1,2)
      DO iobs=1,NOBSBIN
         if ( lo0bin.LT.LOBIN(iobs,1) ) then
            lo0bin = LOBIN(iobs,1)
            lo1bin = LOBIN(iobs,2)
            i0bin = i0bin + 1
            i1bin = 1
         else if ( lo1bin < LOBIN(iobs,2) ) then
            lo1bin = LOBIN(iobs,2)
            i1bin = i1bin + 1
         endif
         if ( iobs.EQ.iObsBin ) then
            IDim1Bin = i1bin
            return
         endif
      enddo

      WRITE(*,*)'ALLUNC:GetIDim1Bin: ERROR! '//
     >     'Observable bin not found. This should never happen, '//
     >     'aborted!'
      STOP

      END FUNCTION



      INTEGER FUNCTION IDim2Bin(iObsBin)
* Returns bin number in third dimension
* Valid for up to triple differential binnings
* 1d, 2d binning --> error exit
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INTEGER iObsBin,i0bin,i1bin,i2bin,iobs
      DOUBLE PRECISION lo0bin,lo1bin,lo2bin
      if ( NDIM < 3 ) then
         WRITE(*,*)'ALLUNC:GetIDim2Bin: ERROR! '//
     >        "No third dimension available, aborted!"
         STOP
      endif
      if ( NOBSBIN.EQ.0 ) then
         WRITE(*,*)'ALLUNC:GetIDim2Bin: ERROR! '//
     >        "No observable bins defined, aborted!"
         STOP
      endif
      if ( iObsBin.GT.NObsBin ) then
         WRITE(*,*)'ALLUNC:GetIDim2Bin: ERROR! '//
     >        "Observable bin out of range, aborted!"
         STOP
      endif

      i0bin = 1
      i1bin = 1
      i2bin = 1
      lo0bin = LOBIN(1,1)
      lo1bin = LOBIN(1,2)
      lo2bin = LOBIN(1,3)
      DO iobs=1,NOBSBIN
         if ( lo0bin.LT.LOBIN(iobs,1) ) then
            lo0bin = LOBIN(iobs,1)
            lo1bin = LOBIN(iobs,2)
            lo2bin = LOBIN(iobs,3)
            i0bin = i0bin + 1
            i1bin = 1
            i2bin = 1
         else if ( lo1bin.LT.LOBIN(iobs,2) ) then
            lo1bin = LOBIN(iobs,2)
            lo2bin = LOBIN(iobs,3)
            i1bin = i1bin + 1
            i2bin = 1
         else if ( lo2bin.LT.LOBIN(iobs,3) ) then
            lo2bin = LOBIN(iobs,3)
            i2bin = i2bin + 1
         endif

         if ( iobs.EQ.iObsBin ) then
            IDim2Bin = i2bin
            return
         endif
      enddo

      WRITE(*,*)'ALLUNC:GetIDim2Bin: ERROR! '//
     >     'Observable bin not found. This should never happen, '//
     >     'aborted!'
      STOP

      END FUNCTION



      INTEGER FUNCTION NDim0Bins()
*     Returns number of bins in first dimension
*     Valid for up to triple differential binnings
*     There always must be at least one bin!
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INTEGER IDim0Bin
      NDim0Bins = IDim0Bin(NOBSBIN)
      return

      END FUNCTION



      INTEGER FUNCTION NDim1Bins(i0Bin)
*     Returns number of bins in second dimension for i0Bin in first
*     Valid for up to triple differential binnings
*     1d binning --> error exit
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INTEGER IDim0Bin,IDim1Bin,i0Bin,iobs
      if ( NDIM < 2 ) then
         WRITE(*,*)'ALLUNC:NDim1Bins: ERROR! '//
     >        "No second dimension available, aborted!"
         STOP
      endif

      NDim1Bins = 0
      DO iobs=1,NOBSBIN
         if ( IDim0Bin(iobs).EQ.i0Bin +1) then
            NDim1Bins = IDim1Bin(iobs-1)
            exit
         else if ( iobs == NOBSBIN ) then
            NDim1Bins = IDim1Bin(iobs)
            exit
         endif
      enddo
      return

      END FUNCTION



      INTEGER FUNCTION NDim2Bins(i0Bin,i1Bin)
*     Returns number of bins in third dimension for i0Bin in first
*     and i1Bin in second dimension
*     Valid for up to triple differential binnings
*     1d, 2d binning --> error exit
      IMPLICIT NONE
      INCLUDE "fnx9999.inc"
      INTEGER IDim0Bin,IDim1Bin,IDim2Bin,i0Bin,i1bin,iobs
      if ( NDIM < 3 ) then
         WRITE(*,*)'ALLUNC:NDim2Bins: ERROR! '//
     >        "No third dimension available, aborted!"
         STOP
      endif

      NDim2Bins = 0
      DO iobs=1,NOBSBIN
         if ( IDim0Bin(iobs).EQ.i0Bin.AND.
     >        IDim1Bin(iobs).EQ.i1Bin+1) then
            NDim2Bins = IDim2Bin(iobs-1)
            exit
         else if ( IDim0Bin(iobs).EQ.i0Bin+1.AND.
     >           IDim1Bin(iobs-1).EQ.i1Bin) then
            NDim2Bins = IDim2Bin(iobs-1)
            exit
         else if ( iobs == NOBSBIN ) then
            NDim2Bins = IDim2Bin(iobs)
            exit
         endif
      enddo
      return

      END FUNCTION
