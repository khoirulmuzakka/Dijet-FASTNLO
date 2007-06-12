c       a simple fortran wrapper around LHAPDFv2,3 for easy interfacing
c       adapted for LHAPDFv4 by Mike Whalley
c       to g++       
c       Stefan Gieseke 2004
	subroutine finitpdfset(name)
	character name*200
	call InitPDFset(name)
	end
c
	subroutine finitpdf(mem)
	integer mem
	call InitPDF(mem)
	end
c
	subroutine fnumberpdf(set)
	integer set
	call numberPDF(set)
	end
c
	subroutine fevolvepdf(x,Q,f)
	real*8 x,Q
	real*8 f(-6:6)
	call evolvePDF(x,Q,f)
	end
c
	subroutine fevolvepdfp(x,Q2,P2,ip,f)
	real*8 x,Q2,P2
	real*8 f(-6:6)
	integer ip
	call evolvePDFp(x,Q2,P2,ip,f)
	end
c
	subroutine falphaspdf(Q,ans)
	real*8 ans,Q
	ans=alphasPDF(Q)
	end
c
	subroutine fgetorderpdf(order)
	integer order
	call getorderpdf(order)
	end
c
	subroutine fgetorderas(order)
	integer order
	call getorderas(order)
	end
c
	subroutine fgetdesc()
	call getdesc()
	end
c
	subroutine fgetqmass(nf,mass)
	integer nf
	real*8 mass
	call getqmass(nf,mass)
	end
c
	subroutine fgetthreshold(nf,Q)
	integer nf
	real*8 Q
	call getthreshold(nf,Q)
	end
c
	subroutine fgetnf(nfmax)
	integer nfmax
	call getnf(nfmax)
	end
c
       subroutine fgetlam4(mem,xlam4)
       integer mem
       real*8 xlam4
       call getlam4(mem,xlam4)
       end
c
       subroutine fgetlam5(mem,xlam5)
       integer mem
       real*8 xlam5
       call getlam5(mem,xlam5)
       end
c
