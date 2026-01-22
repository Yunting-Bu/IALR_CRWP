module m_InitWP
    use m_MachinaBasic, only : f8, c8 
    use m_gPara
    use m_Basis
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

    end subroutine Z_initGaussWP
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setInitTotWP()
        implicit none
        real(f8) :: normWPTot, normWPZ
        integer :: iZ, iPES, j, ir
        integer :: Kmax, K, ijKl, ijKa, ijKi
        real(f8) :: delta, CGtmp, fact
        real(f8), external :: CG 
        real(f8), allocatable :: initWP_BLK(:)

!> WF at Z in DVR
        allocate(initGaussWP(IALR%nZ_IALR))
        call Z_initGaussWP()

        allocate(initWP_BLK(initWP%Kmin:Kmax))
        allocate(int_TDWP(IALR%nZ_I,IALR%vint,int_njK,nPES), &
                int_TDWPm(IALR%nZ_I,IALR%vint,int_njK,nPES))
        allocate(asy_TDWP(IALR%nZ_IA,IALR%vasy,asy_njK,nPES), &
                asy_TDWPm(IALR%nZ_IA,IALR%vasy,asy_njK,nPES))
        allocate(lr_TDWP(IALR%nZ_IALR,IALR%vlr,lr_njK,1), &
                lr_TDWPm(IALR%nZ_IALR,IALR%vlr,lr_njK,1))

!> SF to BF transMat
!> BLK = sqrt(2-delta_K0)*sqrt((2L+1)/(2J+1))*<j,K,L,0|J,K>
        Kmax = min(initWP%j0, initWP%Jtot)
        do K = initWP%Kmin, Kmax
            delta = merge(1.0_f8, dsqrt(2.0_f8), K == 0)
            fact = dsqrt((2.0_f8*initWP%l0+1)/(2.0_f8*initWP%Jtot+1))
            CGtmp = CG(initWP%j0,K,initWP%l0,0,initWP%Jtot)
            initWP_BLK(K) = delta*fact*CGtmp
        end do

!> Construct initial WP in adiabatic representation for lr
        lr_TDWP = 0.0_f8; lr_TDWPm = 0.0_f8
        normWPTot = 0.0_f8
        do K = initWP%Kmin, Kmax
            if (K <= initWP%j0) then 
                ijKl = lr_seqjK(initWP%j0,K)
                if (ijKl == 0) then 
                    write(outFileUnit,*) 'Error in find ijK for initial WP construction in lr'
                    write(outFileUnit,*) 'POSITION: initWP.f90, subroutine getInitTotWP()'
                end if
                do iZ = 1, IALR%nZ_IALR
                    lr_TDWP(iZ,:,ijKl,1) = initWP_BLK(K) * initGaussWP(iZ)
                    do ir = 1, IALR%vlr
                        normWPTot = normWPTot + dabs(lr_TDWP(iZ,ir,ijKl,1))**2
                    end do
                end do
            end if
        end do

!> Rebuild lr_TDWP to asy_TDWP
        asy_TDWP = 0.0_f8; asy_TDWPm = 0.0_f8
        do iPES = 1, nPES
            do K = initWP%Kmin, Kmax
                if (K <= initWP%j0) then 
                    ijKl = lr_seqjK(initWP%j0,K)
                    ijKa = asy_seqjK(initWP%j0,K)
                    if (ijKa == 0 .or. ijKl == 0) then 
                        write(outFileUnit,*) 'Error in find ijK for initial WP construction in asy'
                        write(outFileUnit,*) 'POSITION: initWP.f90, subroutine getInitTotWP()'
                    end if
                    asy_TDWP(:,1:IALR%vlr,ijKa,iPES) = lr_TDWP(1:IALR%nZ_IA,1:IALR%vlr,ijKl,1)
                end if
            end do
        end do 

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
        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a,f15.9)') 'Z0 = ', initWP%Zc, ' , E0 = ', initWP%Ec, ' , delta = ', initWP%delta
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
        complex(c8) :: phase
        integer :: iZ, iEtot

        fact = dsqrt(massTot/(2.0_f8*pi))
        !phase = exp(-img*pi*0.5_f8*initWP%l0)
        wZ = dsqrt(Z_IALR(2)-Z_IALR(1))

        !> rbZ - Ricatti-Bessel j_l
        !> rbZU - Ricatti-Bessel y_l
        !> Second kind Ricatti-Hankel: h^(2)_l = j_l - i*y_l 
        !> a_i(E) = <i\sqrt{mu_R/(2\pi*k_{l0})}h^(2)_{l0}(k_{l0}R)|G(R)>
        !> Ref: PHYSICAL REVIEW A 74, 022703 (2006)

        do iEtot = 1, nEtot
            initAM(iEtot) = imgZero
            do iZ = 1, IALR%nZ_IALR
                call rbesjy(initWP%l0*1.0_f8, kReact(iEtot)*Z_IALR(iZ), rbZ, rbZP, rbZU, rbZUP)
                initAM(iEtot) = initAM(iEtot) + fact/dsqrt(kReact(iEtot)) &
                                * -img * cmplx(rbZ, rbZU, kind=c8) * initGaussWP(iZ) * wZ
            end do 
        end do

        write(outFileUnit,'(1x,a)') '=====> Energy information <====='
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a)') 'All unit in eV!'
        write(outFileUnit,'(1x,a5,3x,a20,3x,a20,3x,a20,3x,a20)') 'nEtot','Ecol', 'Wave Number', 'Amplitude (real)', 'Amplitude (img)'
        write(outFileUnit,'(1x,5("-"),3x,20("-"),3x,20("-"),3x,20("-"),3x,20("-"))')
        do iEtot = 1, nEtot
            write(outFileUnit,'(1x,i5,3x,f20.10,3x,f20.10,3x,f20.10,3x,f20.10)') iEtot, Ecol(iEtot)*au2ev, kReact(iEtot)*au2ev, &
                                                                                initAM(iEtot)%re*au2ev, initAM(iEtot)%im*au2ev
        end do 
        write(outFileUnit,*) ''

        do iEtot = 1, nEtot
            write(888,*) Ecol(iEtot)*au2ev, cdabs(initAM(iEtot))**2
        end do

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