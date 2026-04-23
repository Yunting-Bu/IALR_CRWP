        subroutine POT0(nPES,Vtot,bond)
          implicit none
  
          integer, intent(in) :: nPES 
          real(8), intent(out) :: Vtot(1,1) 
          real(8), intent(in) :: bond(3)
          real(8) :: Z, r, th
          real(8) :: vint, vHF
          real(8), external :: vr

          call bond_to_jacobi(bond, r, Z, th)
          call interaction_PES(r,Z,th,vint)
          vHF = vr(r)
          Vtot(1,1) = vint+vHF


        end subroutine

        subroutine bond_to_jacobi(bond,r,RR,th)
        implicit none
        real(8), intent(in)  :: bond(3)     ! bond(1)=Ar-H, bond(2)=H-F, bond(3)=Ar-F
        real(8), intent(out) :: r           ! Jacobi r (H–F)
        real(8), intent(out) :: RR          ! Jacobi R
        real(8), intent(out) :: th          ! Jacobi angle (degree)
        real(8) :: rArH, rHF, rArF, a1, a2, r1, r2
        real(8) :: pi, mH, mF, mAr

        mH  = 1.00784d0
        mF  = 18.998403d0
        mAr = 39.948d0

        pi = dacos(-1.d0)

        rArH = bond(1)   ! Ar–H
        rHF = bond(2)   ! H–F
        rArF = bond(3)   ! Ar–F

        r = rHF
        r1 = mF * r / (mH+mF)
        r2 = mH * r / (mH+mF)
        
        a1 = (rArH**2+r**2-rArF**2) / (2.d0*r*rArH)
        RR = dsqrt(rArH**2+r1**2-2.d0*r1*rArH*a1)
        a2 = (RR**2+r1**2-rArH**2) / (2.d0*r1*RR)
        if (a2 > 1.d0) a2 = 1.d0
        if (a2 < -1.d0) a2 = -1.d0     
        th = dacos(a2)*180.d0/pi
        if (RR == 0.d0) th = 90.d0


        end subroutine bond_to_jacobi

