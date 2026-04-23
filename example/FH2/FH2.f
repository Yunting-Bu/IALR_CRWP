c     -----------------------------------------------------------------
c     POT0 wrapper for the Stark-Werner F+H2 potential energy surface
c
c     Interface:  call pot0(nPES, diaV, bond)
c
c     Input:
c       nPES     - integer, PES index (1 = ground state SW surface)
c       bond(3)  - double precision, internuclear distances in bohr
c                  bond(1) = A-B distance (AB, F-H)
c                  bond(2) = B-C distance (BC, H-H)
c                  bond(3) = A-C distance (AC, F-H)
c
c     Output:
c       diaV     - double precision, potential energy in Hartree
c                  (relative to the F+H2 asymptote)
c
c     Data files required in working directory:
c       SW.3p  - three-body fit parameters
c       SW.2p  - two-body fit parameters
c
c     Reference:
c       K. Stark and H-J. Werner, J. Chem. Phys. 104, 6515 (1996)
c     -----------------------------------------------------------------
c
      subroutine pot0(nPES, diaV, bond)
      implicit double precision (a-h,o-z)
      integer nPES
      integer b, c, d
      dimension bond(3)
      parameter (nmx=200)
      dimension a(nmx), b(nmx), c(nmx), d(nmx), p(nmx)
      save icall, a, b, c, d, p, ma, np
      data icall /0/
c
      if (icall .eq. 0) then
         !write (6,61)
         open (unit=1, file='SW.3p', status='old')
         i = 1
   2     read (1,*,end=3) nparm, b(i), c(i), d(i), a(i)
            i = i + 1
            goto 2
   3     ma = i - 1
         close (unit=1)
         open (unit=1, file='SW.2p', status='old')
         i = 1
   4     read (1,*,end=5) p(i)
            i = i + 1
            goto 4
   5     np = i - 1
         close (unit=1)
         icall = 1
      endif
c
c     Aguado-Paniagua type fit
c
      x = min(bond(1), bond(3))
      y = bond(2)
      z = max(bond(1), bond(3))
      b1 = a(ma-5)
      b2 = a(ma-4)
      b3 = a(ma-3)
      x0 = a(ma-2)
      y0 = a(ma-1)
      z0 = a(ma)
      fit = 0.0d0
      do i = 1, ma-6
         expon = b(i)*b1*(x-x0) + c(i)*b2*(y-y0) + d(i)*b3*(z-z0)
         fex = exp(-expon)
         fxy = ((x**b(i))*(y**c(i))*(z**d(i)))*fex
         fit = fit + a(i)*fxy
      enddo
      xr = x - p(3)
      yr = y - p(9)
      zr = z - p(3)
      fx = dexp(-p(2)*xr)
      fy = dexp(-p(8)*yr)
      fz = dexp(-p(2)*zr)
      xval = -p(1)*(1.0d0+p(2)*xr+p(4)*xr**2+p(5)*xr**3)*fx + p(6)
      yval = -p(7)*(1.0d0+p(8)*yr+p(10)*yr**2+p(11)*yr**3)*fy + p(12)
      zval = -p(1)*(1.0d0+p(2)*zr+p(4)*zr**2+p(5)*zr**3)*fz + p(6)
c
c     Return potential in Hartree (no eV conversion)
c
      diaV = fit + xval + yval + zval
      return
  61  format(/1x,
     +'This calculation is using the SW F+H2 PES'/1x,
     +'Please cite: ',
     +'K.Stark and H-J.Werner, J. Chem. Phys. 104, 6515 (1996)')
      end
