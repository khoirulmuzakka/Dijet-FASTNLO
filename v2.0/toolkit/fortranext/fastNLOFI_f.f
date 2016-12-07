      SUBROUTINE fastnlo_initpdf
         REAL*8 A1, A2, A3
         COMMON A1, A2, A3
         A1 = 0.5
         A2 = 1
         A3 = 2
      END

      subroutine fastnlo_getxfx(RESULT, X, MUF, F)
         REAL*8 RESULT, X, MUF
         INTEGER F

         REAL*8 A1, A2, A3
         COMMON A1, A2, A3

         RESULT = A1 * X**A2 * (1 - X)**A3
      end

      subroutine fastnlo_evolve_as(RESULT, Q)
         REAL*8 RESULT, Q
         RESULT = 0.118 / ( 1 - 0.118 / 9.423 * dlog( Q**2 / 91.18**2))
      end

      PROGRAM main
         INTEGER CTX
         INTEGER IDX, BINS, DIMS
         REAL*8 XS(1024)
         REAL*8 BININFO(4)
         CHARACTER*65 TABLE
         REAL*8 MUR, MUF
         TABLE = "../data/check/InclusiveNJets_fnr0001midpHT_I723509"//
     _        "_v23_fix.tab.gz"//CHAR(0)

         call fastnlo_create(CTX, TABLE)
         write(*,*) CTX
         call fastnlo_getcrosssection(CTX, XS, BINS)
         write(*,*) CTX, BINS
         do IDX = 1,BINS
           call fastnlo_getobsbindimbounds(CTX, IDX - 1, BININFO, DIMS)
           write(*,*) IDX,BININFO,"->",XS(IDX)
         enddo

Comment:          MUR = 2E0
Comment:          MUF = 2E0
Comment:          call fastnlo_setscalefactorsmurmuf(CTX, MUR, MUF)
Comment:          call fastnlo_getcrosssection(CTX, XS, BINS)
Comment:          do IDX = 1,BINS
Comment:             write(*,*) IDX,BININFO,"->",XS(IDX)
Comment:          enddo

         call fastnlo_destroy(CTX)

      END PROGRAM
