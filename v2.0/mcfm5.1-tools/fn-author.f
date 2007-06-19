*********************************************************************
*                                         M. Wobisch  06/12/2007
*
*  fastNLO interface to MCFM code and authoring tools for 
*  fastNLO tables
*
*********************************************************************
      Subroutine FnInterface(wt)
      Implicit None
      Integer i,j,k
      Double Precision wt

      include 'fastnlo.f'
      include 'qcdcouple.f'
      character*4 part
      common/part/part

c --- two cross checks - in 0th and 1st order ... 
c     can we reproduce the original MCFM results?
      fnsum(1) = fnsum(1) + wt  ! simply add the MCFM user weights
      fntmpsum = 0d0            ! multiply the extracted coefficients & PDFs 
      Do i=-6,6
         Do j=-6,6
            fntmpsum = fntmpsum + fnsigma*fncoeff(i,j)*fnpdf(i,j)
         Enddo
      Enddo
      If (part .eq. 'virt') then
         Do i=-6,6
            Do j=-6,6
               fntmpsum = fntmpsum+fnsigma*fncfv2(i,j)*fnpdf2(i,j)
               fntmpsum = fntmpsum+fnsigma*fncfv3(i,j)*fnpdf3(i,j)
            Enddo
         Enddo
      Endif
      fnsum(2) = fnsum(2) + fntmpsum

c      write(*,*) 'fastNLO interface: wt,fnsum=',wt,fntmpsum,fnsigma
c      write(*,*) '        ',fnsum(1),fnsum(2)
c      write(*,*) 'x12 ',sqrt(fnx(1)*fnx(2))*1960d0,fnx(1),fnx(2)
c      write(*,*) 'QCD coupling gsq,as ',gsq,as

      Return
      End
*********************************************************************
      Subroutine FnInit
      Implicit None
      Integer i

      do i=1,40
         write(*,*) '   fastNLO interface: initialization'
      Enddo

      Return
      End
*********************************************************************
      Subroutine FnTerm
      Implicit None
      Include 'fastnlo.f'

      write(*,*) '   fastNLO interface: termination:'

      write(*,*) ' '
      write(*,*) ' fastNLO test output'
      write(*,*) ' result (sum wgt) = ',fnsum(1)
      write(*,*) ' result (sum x*P) = ',fnsum(2)

      Return
      End
