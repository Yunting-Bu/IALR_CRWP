module m_InitWP
    use m_MachinaBasic, only : f8, c8 
    use m_gPara
    use m_Basis
    use m_Potent
    implicit none

    public
    private :: GaussianWavePacket, Z_initGaussWP
    
contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine Z_initGaussWP()
        implicit none
        real(f8) :: an
        integer :: iZ 

!> Gaussian wave-packet in DVR representation
        do iZ = 1, IALR%nZ_IALR
            initGaussWP(iZ) = GaussianWavePacket(initWP%Ec,initWP%delta,initWP%Zc,massTot,Z_IALR(iZ)) &
                              * dsqrt(Z_IALR(2)-Z_IALR(1))
        end do 

        an = 0.0_f8
        do iZ = 1, IALR%nZ_IALR
            an = an + initGaussWP(iZ)*initGaussWP(iZ)
        end do 

        !if (an > 0.0_f8) initGaussWP(:) = initGaussWP(:) / dsqrt(an)

    end subroutine Z_initGaussWP
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setInitTotWP()
        implicit none
        real(f8) :: normWPTot, normWPZ
        integer :: iZ, iPES, j, ir, ir0
        integer :: Kmax, K, ijKl, ijKa, ijKi
        real(f8) :: delta, CGtmp, fact
        real(f8), external :: CG 
        real(f8), allocatable :: initWP_BLK(:)

!> WF at Z in DVR
        allocate(initGaussWP(IALR%nZ_IALR))
        call Z_initGaussWP()

        Kmax = min(initWP%j0, initWP%Jtot)
        allocate(initWP_BLK(initWP%Kmin:Kmax))
        allocate(int_TDWP(IALR%nZ_I,IALR%vint,int_njK,nPES), &
                int_TDWPm(IALR%nZ_I,IALR%vint,int_njK,nPES))
        allocate(asy_TDWP(IALR%nZ_IA,IALR%vasy,asy_njK,nPES), &
                asy_TDWPm(IALR%nZ_IA,IALR%vasy,asy_njK,nPES))
        allocate(lr_TDWP(IALR%nZ_IALR,IALR%vlr,lr_njK,1), &
                lr_TDWPm(IALR%nZ_IALR,IALR%vlr,lr_njK,1))

!> SF to BF transMat
!> BLK = sqrt(2-delta_K0)*sqrt((2L+1)/(2J+1))*<j,K,L,0|J,K>
        do K = initWP%Kmin, Kmax
            delta = merge(1.0_f8, dsqrt(2.0_f8), K == 0)
            fact = dsqrt((2.0_f8*initWP%l0+1.0_f8)/(2.0_f8*initWP%Jtot+1.0_f8))
            CGtmp = CG(initWP%j0,K,initWP%l0,0,initWP%Jtot)
            initWP_BLK(K) = delta*fact*CGtmp
        end do

!> Construct initial WP in diabatic representation for lr
        lr_TDWP = 0.0_f8; lr_TDWPm = 0.0_f8
        normWPTot = 0.0_f8
        do K = initWP%Kmin, Kmax
            if (K <= initWP%j0) then 
                ijKl = lr_seqjK(initWP%j0,K)
                if (ijKl == 0) then 
                    write(outFileUnit,*) 'Error in find ijK for initial WP construction in lr'
                    write(outFileUnit,*) 'POSITION: initWP.f90, subroutine getInitTotWP()'
                end if
                ir0 = initWP%v0 + 1
                if (ir0 < 1 .or. ir0 > IALR%vlr) then
                    write(outFileUnit,*) 'Error in initial v0 index for WP construction in lr'
                    write(outFileUnit,*) 'POSITION: initWP.f90, subroutine getInitTotWP()'
                else
                    ir = ir0
                    do iZ = 1, IALR%nZ_IALR
                        lr_TDWP(iZ,ir,ijKl,1) = initWP_BLK(K) * initGaussWP(iZ)
                        normWPTot = normWPTot + dabs(lr_TDWP(iZ,ir,ijKl,1))**2
                    end do
                end if
            end if
        end do

!> Rebuild lr_TDWP to asy_TDWP
        asy_TDWP = 0.0_f8; asy_TDWPm = 0.0_f8
        !do iPES = 1, nPES
            do K = initWP%Kmin, Kmax
                if (K <= initWP%j0) then 
                    ijKl = lr_seqjK(initWP%j0,K)
                    ijKa = asy_seqjK(initWP%j0,K)
                    if (ijKa == 0 .or. ijKl == 0) then 
                        write(outFileUnit,*) 'Error in find ijK for initial WP construction in asy'
                        write(outFileUnit,*) 'POSITION: initWP.f90, subroutine getInitTotWP()'
                    end if
                    asy_TDWP(:,1:IALR%vlr,ijKa,initWP%PES0) = lr_TDWP(1:IALR%nZ_IA,1:IALR%vlr,ijKl,1)
                    !asy_TDWP(:,1:IALR%vlr,ijKa,iPES) = lr_TDWP(1:IALR%nZ_IA,1:IALR%vlr,ijKl,1)
                end if
            end do
        !end do 

!> Rebuild asy_TDWP to int_TDWP
        int_TDWP = 0.0_f8; int_TDWPm = 0.0_f8
        do iPES = 1, nPES
            do K = initWP%Kmin, Kmax
                do j = initWP%jmin, IALR%jasy, initWP%jinc 
                    if (K <= j) then 
                        ijKa = asy_seqjK(j,K)
                        ijKi = int_seqjK(j,K)
                        if (ijKa == 0 .or. ijKi == 0) then 
                            write(outFileUnit,*) 'Error in find ijK for initial WP construction in int'
                            write(outFileUnit,*) 'POSITION: initWP.f90, subroutine getInitTotWP()'
                        end if
                        int_TDWP(:,1:IALR%vasy,ijKi,iPES) = asy_TDWP(1:IALR%nZ_I,1:IALR%vasy,ijKa,iPES)
                    end if
                end do
            end do
        end do

        write(outFileUnit,'(1x,a)') "=====> Initial WP infromation <====="
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a,f15.9)') 'Z0 = ', initWP%Zc, ' , E0 (eV) = ', initWP%Ec*au2ev, ' , delta = ', initWP%delta
        normWPZ = sum( abs(initGaussWP)**2 )
        write(outFileUnit,'(1x,a,f15.9)') 'Initial Gaussian wave-packet normalization check: ', normWPZ
        write(outFileUnit,'(1x,a,f15.9)') 'Initial adiabatic total wave-packet normalization check: ', normWPTot
        write(outFileUnit,'(1x,a)') ''

    end subroutine setInitTotWP
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setEnergyAM()
        implicit none
        !> For Ricatti-Bessel function 
        real(f8) :: rbZ, rbZU, rbZP, rbZUP 
        real(f8) :: fact, wZ
        !!> Phase for trans the Bessel to Hankel
        !complex(c8) :: phase
        integer :: iZ, iEtot, chkUnit

        fact = dsqrt(massTot/(2.0_f8*pi))
        !phase = exp(-img*pi*0.5_f8*initWP%l0)
        wZ = dsqrt(Z_IALR(2)-Z_IALR(1))

        !> rbZ - Ricatti-Bessel j_l
        !> rbZU - Ricatti-Bessel y_l
        !> Second kind Ricatti-Hankel: h^(2)_l = j_l - i*y_l 
        !> a_i(E) = <i\sqrt{mu_R/(2\pi*k_{l0})}h^(2)_{l0}(k_{l0}R)|G(R)>
        !> Ref: PHYSICAL REVIEW A 74, 022703 (2006)

        if (IF_modHankel) then
            call set1DVeff()
            call setModifiedHankel()
        else
            do iEtot = 1, nEtot
                initAM(iEtot) = imgZero
                do iZ = 1, IALR%nZ_IALR
                    call rbesjy(initWP%l0*1.0_f8, kReact(iEtot)*Z_IALR(iZ), rbZ, rbZP, rbZU, rbZUP)
                    initAM(iEtot) = initAM(iEtot) + fact/dsqrt(kReact(iEtot)) &
                                    * (-img) * cmplx(rbZ, rbZU, kind=c8) * initGaussWP(iZ) * wZ
                end do 
            end do
        end if

        open(newunit=chkUnit, file='Energy_Info.chk', status='replace', action='write')
        write(chkUnit,'(1x,a)') "# iEtot, Ecol(eV), |a_i(E)|^2"
        do iEtot = 1, nEtot
            write(chkUnit,'(1x,i5,3x,f20.10,3x,f20.10)') &
                    iEtot, Ecol(iEtot)*au2ev, &
                    cdabs(initAM(iEtot))**2*2.0_f8*pi
        end do

        close(chkUnit)

        write(outFileUnit,'(1x,a)') '=====> Energy information <====='
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a,i5)') 'Total energy points: ', nEtot
        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a)') 'Collsion energy interval [ ', Ecol(1)*au2ev, ' eV , ', Ecol(nEtot)*au2ev, ' eV ]'
        write(outFileUnit,'(1x,a)') 'Please check Energy_Info.chk !'
        write(outFileUnit,*) ''


    end subroutine setEnergyAM
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine initDFFT()
        implicit none
        integer :: i, j 
        real(f8) :: Zfact

        lr_nsZ = 3*IALR%nZ_IALR+15
        allocate(lr_WSaveZ(lr_nsZ))
        call dsinti(IALR%nZ_IALR,lr_WSaveZ)

        asy_nsZ = 3*IALR%nZ_IA+15
        allocate(asy_WSaveZ(asy_nsZ,asy_njK))
        do i = 1, asy_njK
            call dsinti(IALR%nZ_IA,asy_WSaveZ(1,i))
        end do

        int_nsZ = 3*IALR%nZ_I+15
        allocate(int_WSaveZ(int_nsZ,int_njK))
        do i = 1, int_njK
            call dsinti(IALR%nZ_I,int_WSaveZ(1,i))
        end do

!> Combine kinMat_r and kinMat_Z 
        allocate(int_kinEigen(IALR%nZ_I,IALR%vint))
        allocate(asy_kinEigen(IALR%nZ_IA,IALR%vasy))
        allocate(lr_kinEigen(IALR%nZ_IALR,IALR%vlr))

        !> interaction region
        Zfact = 0.5_f8/(IALR%nZ_I+1.0_f8)
        do j = 1, IALR%vint
            do i = 1, IALR%nZ_I
                int_kinEigen(i,j) = (int_rKinMat(j) + int_ZkinMat(i))*Zfact
            end do
        end do

        !> asymptotic region
        Zfact = 0.5_f8/(IALR%nZ_IA+1.0_f8)
        do j = 1, IALR%vasy
            do i = 1, IALR%nZ_IA
                asy_kinEigen(i,j) = (asy_rKinMat(j) + asy_ZkinMat(i))*Zfact
            end do
        end do

        !> long-range region
        Zfact = 0.5_f8/(IALR%nZ_IALR+1.0_f8)
        do j = 1, IALR%vlr
            do i = 1, IALR%nZ_IALR
                lr_kinEigen(i,j) = (lr_rKinMat(j) + lr_ZkinMat(i))*Zfact
            end do
        end do
    end subroutine initDFFT
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine set1DVeff()
!> Compute 1D effective potential along Z in SF representation:
!> Veff(Z) = <v0 j0 l0 | V(Z,r,theta) - V_BC(r) | v0 j0 l0>
        implicit none
        integer :: iZ, ir, ith, ijK, iv0, K, Kmax
        real(f8) :: Z, r, theta, VBC
        real(f8) :: Vdia(1,1), bond(3)
        real(f8) :: Vth, rWF2, VeffK
        real(f8) :: delta, fact, CGtmp, BLK2
        real(f8), external :: CG
        real(f8), allocatable :: Vpot(:,:,:)

        allocate(Vpot(IALR%nZ_IALR, IALR%vlr, IALR%lr_nA))
        allocate(VeffMat(IALR%nZ_IALR))

!> Compute V(Z,r,theta) - V_BC(r) on the grid
        do iZ = 1, IALR%nZ_IALR
            Z = Z_IALR(iZ)
            do ir = 1, IALR%vlr
                r = r_LR(ir)
                do ith = 1, IALR%lr_nA
                    theta = acos(lr_ANode(ith))
                    call Jacobi2Bond(Z, r, theta, atomMass(2), atomMass(3), bond)
                    call POT0(1, Vdia, bond)
                    call setPot_vj0(r, VBC)
                    Vpot(iZ, ir, ith) = min((Vdia(1,1) - VBC), VMaxCut)
                end do
            end do
        end do

        iv0 = initWP%v0 + 1
        Kmax = min(initWP%j0, initWP%Jtot)

        VeffMat = 0.0_f8
        do K = initWP%Kmin, Kmax
            ijK = lr_seqjK(initWP%j0, K)
!> BF to SF coefficient squared: |BLK(K)|^2
            delta = merge(1.0_f8, dsqrt(2.0_f8), K == 0)
            fact = dsqrt((2.0_f8*initWP%l0 + 1.0_f8) / (2.0_f8*initWP%Jtot + 1.0_f8))
            CGtmp = CG(initWP%j0, K, initWP%l0, 0, initWP%Jtot)
            BLK2 = (delta * fact * CGtmp)**2

            do iZ = 1, IALR%nZ_IALR
                VeffK = 0.0_f8
                do ir = 1, IALR%vlr
                    rWF2 = lr_POWF(ir,iv0)**2
                    Vth = 0.0_f8
                    do ith = 1, IALR%lr_nA
                        Vth = Vth+lr_YMat(ith, ijK)**2*Vpot(iZ, ir, ith)*lr_AWeight(ith)**2
                    end do
                    VeffK = VeffK + rWF2 * Vth
                end do
                VeffMat(iZ) = VeffMat(iZ) + BLK2 * VeffK
            end do
        end do

        deallocate(Vpot)

    end subroutine set1DVeff
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setModifiedHankel()
        implicit none
        !> For Ricatti-Bessel function 
        real(f8) :: rbZ, rbZU, rbZP, rbZUP 
        real(f8) :: fact, wZ
        complex(c8), allocatable :: modHankel(:)
        integer :: iZ, iEtot, chkUnit

        fact = dsqrt(massTot/(2.0_f8*pi))
        wZ = dsqrt(Z_IALR(2)-Z_IALR(1))
        allocate(ModHankel(IALR%nZ_IALR))

        do iEtot = 1, nEtot 
            initAM(iEtot) = imgZero
            modHankel = imgZero
            call rbesjy(initWP%l0*1.0_f8, kReact(iEtot)*Z_IALR(IALR%nZ_IALR), rbZ, rbZP, rbZU, rbZUP)
            modHankel(IALR%nZ_IALR) = fact/dsqrt(kReact(iEtot)) * (-img) * cmplx(rbZ, rbZU, kind=c8)
            call rbesjy(initWP%l0*1.0_f8, kReact(iEtot)*Z_IALR(IALR%nZ_IALR-1), rbZ, rbZP, rbZU, rbZUP)
            modHankel(IALR%nZ_IALR-1) = fact/dsqrt(kReact(iEtot)) * (-img) * cmplx(rbZ, rbZU, kind=c8)
            call Numerov(IALR%nZ_IALR, Z_IALR(2)-Z_IALR(1), modHankel, VeffMat, Ecol(iEtot), massTot)
            do iZ = 1, IALR%nZ_IALR
                initAM(iEtot) = initAM(iEtot) + modHankel(iZ)*initGaussWP(iZ)*wZ
            end do
        end do  
        
    end subroutine setModifiedHankel
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine Numerov(nZ, dZ, cmplxY, Vpot, E0, tms)
!> Numerov method for solving the 1D Schrodinger equation
!> Y(i-2) = (2*(1-5/12*dZ^2*k(i-1))*Y(i-1) - (1+dZ^2/12*k(i))*Y(i)) / (1+dZ^2/12*k(i-2))
!> Ref: COMPUTERS IN PHYSICS, VOL. 11, NO. 5, SEP/OCT 1997
        implicit none
        integer, intent(in) :: nZ
        real(f8), intent(in) :: dZ, E0, tms 
        complex(c8), intent(inout) :: cmplxY(nZ)
        real(f8), intent(in) :: Vpot(nZ)
        integer :: iZ 
        real(f8) :: k1, k2, k3

        do iZ = nZ, 3, -1
            k1 = 2.0_f8*(E0-Vpot(iZ-2))*tms
            k2 = 2.0_f8*(E0-Vpot(iZ-1))*tms
            k3 = 2.0_f8*(E0-Vpot(iZ))*tms
            cmplxY(iZ-2) = ((2.0_f8-10.0_f8*dZ**2/12.0_f8*k2)*cmplxY(iZ-1) &
                            - (1.0_f8+dZ**2/12.0_f8*k3)*cmplxY(iZ)) &
                            / (1.0_f8+dZ**2/12.0_f8*k1)
        end do
    
    end subroutine Numerov
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    real(f8) function GaussianWavePacket(E0, delta, Z0, mass, Z) result(WP)
        implicit none
        real(f8), intent(in) :: E0, delta, Z0, mass, Z
        real(f8) :: fact, expPart, cosPart 

        fact = (1.0_f8/(pi*delta*delta))**(0.25)
        expPart = -(Z-Z0)**2 / delta**2
        cosPart = Z * dsqrt(2.0_f8*E0*mass)
        WP = fact * dexp(expPart/2.0_f8)*dcos(cosPart)

        return
    end function GaussianWavePacket
!> ------------------------------------------------------------------------------------------------------------------ <!
    
end module m_InitWP
