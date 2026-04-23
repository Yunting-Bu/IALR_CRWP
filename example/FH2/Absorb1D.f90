program Absorb1D
    use m_MachinaBasic, only : f8, c8
    use m_Basis, only : sinDVRGrid
    implicit none
    
    real(f8), parameter :: au2ev = 27.211386245988_f8
    real(f8), parameter :: ev2au = 1.0_f8 / au2ev
    real(f8), parameter :: au2cm = 219474.6313705_f8
    real(f8), parameter :: cm2au = 1.0_f8 / au2cm 
    real(f8), parameter :: au2K = 315777.663_f8
    real(f8), parameter :: K2au = 1.0_f8 / au2K
    real(f8), parameter :: amu2au = 1822.888486209_f8
    real(f8), parameter :: pi = 4.0_f8*atan(1.0_f8)
    complex(c8), parameter :: img = (0.0_f8, 1.0_f8)
    complex(c8), parameter :: imgZero = (0.0_f8, 0.0_f8)
    

    !> For namelist NIP
    real(f8) :: Cabs, AbsRange, covThreshold, fluxPos 
    !> For namelist WP_Grids
    integer :: nZ, energyUnit 
    real(f8) :: Z_range(2), Z0, Ec, delta 
    character(len=2) :: Atoms(3)
    !> For namelist energy_Prop
    integer :: timeStep, timePrint
    real(f8) :: E_range(2), dE 

    namelist /NIP/ Cabs, AbsRange, covThreshold, fluxPos
    namelist /WP_Grids/ nZ, energyUnit, Z_range, Z0, Ec, delta, Atoms
    namelist /energy_Prop/ timeStep, timePrint, E_range, dE

    integer :: inpFileUnit
    !> For energy
    integer :: nEtot
    real(f8), allocatable :: Etot(:)
    real(f8) :: energyUnitTrans(4) = [cm2au, ev2au, K2au, 1.0_f8]
    !> For getMass
    real(f8) :: atomMass(3), massTot
    !> For sin DVR grid
    real(f8), allocatable :: ZGrid(:)
    !> For initial wave packet
    real(f8), allocatable :: initGaussWP(:)

    open(newunit=inpFileUnit, file='Vabs.inf', status='old')
    read(inpFileUnit, nml=NIP)
    read(inpFileUnit, nml=WP_Grids)
    read(inpFileUnit, nml=energy_Prop)
    close(inpFileUnit)


    call energySet()
    call getMass()

    allocate(ZGrid(nZ))
    call sinDVRGrid(Z_range(1), Z_range(2), nZ, ZGrid)
    allocate(initGaussWP(nZ))
    call Z_initGaussWP()



contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine energySet()
        implicit none
        integer :: iEtot

        nEtot = int((E_range(2) - E_range(1))/dE) + 1
        allocate(Etot(nEtot))
        do iEtot = 1, nEtot 
            Etot(iEtot) = E_range(1) + (iEtot-1)*dE 
            if (Etot(iEtot) > E_range(2)) then
                Etot(iEtot) = E_range(2)
            end if 
            !> convert to au
            Etot(iEtot) = Etot(iEtot) * energyUnitTrans(energyUnit)
        end do
        Ec = Ec * energyUnitTrans(energyUnit)

    end subroutine energySet
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine getMass()
        implicit none
        integer, parameter :: nSupportedElements = 10
        character(len=2), parameter :: elemSymbol(nSupportedElements) = &
            ['H ', 'D ', 'He', 'Li', 'N ', 'O ', 'F ', 'S ', 'Cl', 'Ar']
        real(f8), parameter :: elemMass(nSupportedElements) = &
            [1.00782503223_f8, 2.01410177812_f8, 4.002602_f8, 6.938_f8, 14.00307400443_f8, &
            15.99491461957_f8, 18.99840316273_f8, 32.06_f8, 34.968852682_f8, 39.9623831237_f8]
        
        integer :: iatm, ipos

        do iatm = 1, 3 
            ipos = findloc(elemSymbol, Atoms(iatm), dim=1)
            if (ipos == 0) then
                write(outFileUnit,*) "Error: Unsupported element ", Atoms(iatm)
                write(outFileUnit,*) "POSITION: Absorb1D.f90, subroutine getMass()"
                stop
            else
                atomMass(iatm) = elemMass(ipos)*amu2au
            end if
        end do

        massTot = atomMass(1)*(atomMass(2)+atomMass(3))/(atomMass(1)+atomMass(2)+atomMass(3))

    end subroutine getMass
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine Z_initGaussWP()
        implicit none
        real(f8) :: an
        integer :: iZ 

!> Gaussian wave-packet in DVR representation
        do iZ = 1, nZ
            initGaussWP(iZ) = GaussianWavePacket(Ec,delta,Z0,massTot,ZGrid(iZ)) &
                              * dsqrt(ZGrid(2)-ZGrid(1))
        end do 

        an = 0.0_f8
        do iZ = 1, nZ
            an = an + initGaussWP(iZ)*initGaussWP(iZ)
        end do 

    end subroutine Z_initGaussWP
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    real(f8) function GaussianWavePacket(E0, del, Zc, mass, Z) result(WP)
        implicit none
        real(f8), intent(in) :: E0, del, Zc, mass, Z
        real(f8) :: fact, expPart, cosPart 

        fact = (1.0_f8/(pi*del*del))**(0.25)
        expPart = -(Z-Zc)**2 / del**2
        cosPart = Z * dsqrt(2.0_f8*E0*mass)
        WP = fact * dexp(expPart/2.0_f8)*dcos(cosPart)

        return
    end function GaussianWavePacket
!> ------------------------------------------------------------------------------------------------------------------ <!


end program Absorb1D