c
c This program computes the price of a European option (call and put) 
c using the simple jump diffusion model of Kou (2002) discussed in Chapter 6 
c of the book Tsay (2002).
c
c The program is intended for use in my course ``Analysis of Financial Time 
c Series'' at GSB. I tested the program, but errors may still exist. Please 
c bring any error to my attention. 
c
c You may use the program at your own risk!  Ruey S. Tsay, December 2000.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit double precision(a-h,o-z)
c      
      write(6,1)
 1    format(1x,'Current and strike prices: ',$)
      read(5,*) pr, sk
      write(6,2)
 2    format(1x,'time to expiration (in years): ',$)
      read(5,*) dt
      write(6,3)
 3    format(1x,'Interest rate and volatility (std): ',$)
      read(5,*) free, sig
      if(sig .le. 0.0d0)then
       print*,'volatility must be positive!'
       stop
      endif
      write(6,4)
 4    format(1x,'lambda, kappa, eta: ',$)
      read(5,*) poi, rka, eta
      if((eta.le.0.0d0).or.(eta.gt.1.0d0))then
       print*,'eta is out of range!'
       stop
      endif
c
      nterm = 10
      write(6,5)
 5    format(1x,'terms of truncation? (default = 10): ',$)
      read(5,*) nterm
      if(nterm .le. 10) nterm = 10
      if(nterm .gt. 16) nterm = 16
c
      psi = dexp(rka)/(1.0d0-eta*eta)-1.0d0
      prsk = dlog(pr/sk)
      sigsrt = sig*dsqrt(dt)
      sig2h = sig*sig/2.0d0
      tlamb = poi*dt
      rtexp = dexp(-free*dt)
      ome1 = -prsk+tlamb*psi-(free-sig2h)*dt
      w1 = poi*psi*dt
      sigt = dexp(sig*sig*dt/(2.0d0*eta*eta))
      st2pi = dsqrt(6.28318530718d0)
      steta = sigsrt/eta
c
c----- compute the portion of geometric Brownian motion
      hp = (prsk+(free+sig2h-poi*psi)*dt)/sigsrt
      hm = (prsk+(free-sig2h-poi*psi)*dt)/sigsrt
      call gauss(hp,cdf1)
      call gauss(hm,cdf2)
      cgeo = dexp(-tlamb)*(pr*dexp(-w1)*cdf1-sk*rtexp*cdf2)
      print*,'Price due to the geometric part: ', cgeo
c---- compute the price for the jump process portion
c
      cj = 0.0d0
c
      do 100 n = 1, nterm
       tmp = dfloat(n)*rka/sigsrt
       bp = hp+tmp
       bm = hm+tmp
       ome = ome1-dfloat(n)*rka
       tmp1 = ome/sigsrt
       cp = steta+tmp1
       cm = steta-tmp1
       omexp = dexp(ome/eta)
       omexpi = dexp(-ome/eta)
       a11 = dexp(-w1+dfloat(n)*rka)*pr/2.0d0
       a22 = rtexp*omexpi*sigt*sk/2.0d0
       a33 = rtexp*omexp*sigt*sk/2.0d0
c
       call gauss(bp,cdf1)
       call gauss(bm,cdf2)
c
       fac = factor(n)
       npow = 2**(2*n-1)
       coef = (tlamb**n)*dexp(-tlamb)/fac
c
       a = 0.0d0
       do 50 j=1, n
        m = 2*n-j-1
        call comb(m,n-1,icnt)
        icnt = icnt*(2**j)
        coef1 = coef*dfloat(icnt)/dfloat(npow)
        a1 = (1.0d0/(1.0d0-eta)**j+1.0d0/(1.0d0+eta)**j)*a11*cdf1 
     &   -rtexp*sk*cdf2
        a2 = 0.0d0
        a3 = 0.0d0
        do 30 ii=1, j
         i = ii-1
          w2 = 1.0d0/(1.0d0-eta)**(j-i)-1.0d0
          w2 = w2*(steta**i)*hh(cm,i)/st2pi
          w3 = 1.0d0-1.0d0/(1.0d0+eta)**(j-i)
          w3 = w3*(steta**i)*hh(cp,i)/st2pi
          a2 = a2 + w2
          a3 = a3 + w3
 30      continue
         a2 = a2*a22
         a3 = a3*a33
         a = coef1*(a1+a2+a3)+a
 50       continue
         cj = cj + a
 100      continue
c
      print*,'price due to jump component: ', cj
      c = cj + cgeo
      print*,'Price of a call: ', c
c
c---- put option
c
      p = c+sk*dexp(-free*dt)-pr
c 
      print*,'Price of a put: ', p
c
      stop
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine comb(m,n,icnt)
c----- combites the combinations      
      implicit double precision(a-h,o-z)
      if((m .lt. n).or.(m .lt. 0))then
       icnt = 0
       return
      endif
      if(m.eq.n)then
       icnt = 1
       return
      endif
c
      n1 = 1
      do 5 j=n, 1, -1
 5       n1 = n1*j
c
      icnt = 1
      do 10 i=m, (m-n+1), -1
 10      icnt =icnt*i
      icnt = icnt/n1
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function factor(n)
      implicit double precision(a-h,o-z)
c
      if(n .lt. 0)then
        factor = 0.0d0
        return
       endif
       if(n.le. 1)then
        factor = 1.0d0
        return
       endif
       a = 1.0d0
       do 10 i=n, 1, -1
 10       a = a*dfloat(i)
       factor = a
c
       return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double precision function hh(x,n)
c
      implicit double precision(a-h,o-z)
      dimension h(100)
      if(n .gt. 100)then
       print*,'dimension in Hn function is too low!'
       stop
      endif
c
      pisqt = dsqrt(6.28318530718d0)
      h(1) = dexp(-x*x/2.0d0)
      if(n .eq. -1)then
       hh = h(1)
       return
      endif
      tmp = -x
      call gauss(tmp,cdf)
      h(2) = pisqt*cdf
      if(n .eq. 0)then
        hh = h(2)
        return
      endif
      do 40 j=3, 2+n
 40        h(j) = (h(j-2)-x*h(j-1))/dfloat(j-2)
 50        continue
c
       hh = h(n+2)
c
       return
       end
ccccccccccccccccccccccccccccccccccccccc
       subroutine gauss(x,cdf)
c
       implicit double precision(a-h,o-z)
       dimension a(5)
c
       a(1) = 0.319381530d0
       a(2) = -0.356563782d0
       a(3) = 1.781477937d0
       a(4) = -1.821255978d0
       a(5) = 1.330274429d0
c
c---------------- n: is used to indicate negative value.
       n = 0
c
       if(x .lt. 0.0d0) then
        x = -x 
        n = 1
       endif
       xk = 1.0d0/(1.0d0+0.2316419d0*x)
       xf = dexp(-x*x/2.0d0)/dsqrt(6.28318530718d0)
       tmp = 0.0d0
       do 10 i=1, 5
 10       tmp = tmp + a(i)*xk**i
       cdf = 1.0d0-xf*tmp
       if(n.eq.1)cdf = 1.0d0-cdf
c
       return
       end
cccccccccccccccccccccccccccccccccccccc

     
