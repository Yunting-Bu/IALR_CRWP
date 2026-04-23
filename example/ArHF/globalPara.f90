module m_gPara
    use m_MachinaBasic , only : f8, c8
    implicit none
    public

    private :: diatomParity, energySet

!> ========== Constants ==========

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

!> ========== Global Parameters ==========

    type :: initWP_class
        integer :: j0, v0, l0, Jtot, tpar, jpar
        real(f8) :: Zc, delta, Ec
        integer :: PES0
        integer :: Kmin, jmin, jinc 
    end type initWP_class

    type :: IALR_class
        integer :: nZ_IALR, nZ_IA, nZ_I, nr_DVR
        integer :: vint, jint, vasy, jasy, vlr 
        integer :: int_nA, asy_nA, lr_nA
        real(f8) :: Z_range(2), r_range(2)
    end type IALR_class

    type :: Vabs_class
        character(len=3) :: type
        integer :: nZasy, nZlr, nrint
        real(f8) :: Casy, Clr, Cr
        real(f8) :: Zasy_range, Zlr_range, rabs_range
    end type Vabs_class

    type :: channel1_class
        integer :: nr, vmax, jmax, jpar, ichoice 
        real(f8) :: r_range(2), Zinf, massAB, massTot
        integer :: jmin, jinc, Kmin, njK, nIC
    end type channel1_class

    type :: channel2_class
        integer :: nr, vmax, jmax, jpar, ichoice
        real(f8) :: r_range(2), Zinf, massAC, massTot
        integer :: jmin, jinc, Kmin, njK, nIC
    end type channel2_class

!> ========== Parameters from input file ==========

!> Task types and propagation settings
    integer :: nChannel, nPES, energyUnit
    integer :: timeTot, timePrint
    real(f8) :: timeStep
!> Energy input
    integer :: nEtot, outFileUnit
    real(f8) :: E_range(2), dE, TMaxCut, VMaxCut
    real(f8), allocatable :: Etot(:)
    real(f8), allocatable :: Ecol(:)
    real(f8) :: energyUnitTrans(4) = [cm2au, ev2au, K2au, 1.0_f8]
!> Wave number in reactant coordinate, k_i and a_i(E)
    real(f8), allocatable :: kReact(:)
    complex(c8), allocatable :: initAM(:)
!> Atoms and masses
    real(f8) :: atomMass(3), massBC, massTot
    character(len=2) :: Atoms(3)
!> Output and other flags
    character(len=4) :: potentialType
    character(len=50) :: outfile
!> Inelastic
    logical :: IF_inelastic
    real(f8) :: inePos
    integer :: iInePos
    real(f8), allocatable :: Evj_BC(:,:,:)
    real(f8), allocatable :: WFvj_BC(:,:,:,:,:,:)
    complex(c8), allocatable :: ine_TIDWF(:,:,:,:,:)
!> Flux
    character(len=1) :: IDflux
    integer :: iZFluxPos
    real(f8) :: fluxPos
!    real(f8), allocatable :: diffZFlux(:)
!> For new Hankel method
    logical :: IF_modHankel
    real(f8), allocatable :: VeffMat(:)
!> Basis sets
!> theta
    integer :: int_njK
    integer :: asy_njK
    integer :: lr_njK
    integer, allocatable :: int_jKPair(:,:)
    integer, allocatable :: asy_jKPair(:,:)
    integer, allocatable :: lr_jKPair(:,:)
    integer, allocatable :: int_seqjK(:,:)
    integer, allocatable :: asy_seqjK(:,:)
    integer, allocatable :: lr_seqjK(:,:)
    real(f8), allocatable :: int_lo(:)
!> Grids and weights for K independent Gauss-Legendre quadrature
    real(f8), allocatable :: asy_ANode(:), asy_AWeight(:)
    real(f8), allocatable :: int_ANode(:), int_AWeight(:)
    real(f8), allocatable :: lr_ANode(:), lr_AWeight(:)
!> FBR to DVR transformation matrix
    real(f8), allocatable :: int_YMat(:,:), asy_YMat(:,:), lr_YMat(:,:)
!> Z
!> I - interaction region
!> A - asymptotic region
!> LR - long range region
!> Interaction-Asympotic-Lonng-Range (IALR) range
!> r3   |----|
!> r2   |----|----|
!>      |----|----|----|
!> r1   |----|----|----|
!>      Z1   Z2   Z3   Z4
!>  nPODVR in range [r1, r3] = vint, in range [r1, r2] = vasy 
!>  VAbs: Z3 = Zabs_asy_end, Z4 = Zabs_lr_end
    real(f8), allocatable :: Z_IALR(:)
    real(f8), allocatable :: initGaussWP(:)
!> r
    real(f8), allocatable :: r_DVR(:), r_Int(:), r_Asy(:), r_LR(:)
    real(f8), allocatable :: int_POWF(:,:), asy_POWF(:,:), lr_POWF(:,:)
    real(f8), allocatable :: int_PO2FBR(:,:), asy_PO2FBR(:,:), lr_PO2FBR(:,:)
    real(f8), allocatable :: int_POEig(:), asy_POEig(:), lr_POEig(:)
    real(f8), allocatable :: int_PO2DVR(:,:), asy_PO2DVR(:,:)
!> Vabs
    real(f8), allocatable :: asy_ZFabs(:), lr_ZFabs(:), int_rFabs(:)
!> Hamiltonian matrix
    real(f8), allocatable :: int_ZkinMat(:), asy_ZkinMat(:), lr_ZkinMat(:)
    real(f8), allocatable :: int_rKinMat(:), asy_rKinMat(:), lr_rKinMat(:)
    real(f8), allocatable :: int_kinEigen(:,:), asy_kinEigen(:,:), lr_kinEigen(:,:)
    real(f8), allocatable :: int_rotMat(:,:), asy_rotMat(:,:), lr_rotMat(:,:)
    real(f8), allocatable :: IA_CPMat(:,:,:), lr_CPMat(:,:,:)
!> In case nK == 1
    real(f8), allocatable :: int_CPK1(:,:), asy_CPK1(:,:), lr_CPK1(:,:)
!> Interaction potential matrix and diabatic-to-adiabatic transformation matrix
    real(f8), allocatable :: int_Vdia(:,:,:,:,:), asy_Vdia(:,:,:,:,:), lr_Vdia(:,:,:)
!> DFFT 
    integer :: int_nsZ, asy_nsZ, lr_nsZ
    real(f8), allocatable :: int_WSaveZ(:,:), asy_WSaveZ(:,:), lr_WSaveZ(:)
!> Wave packet during propagation
    !> m for k-1, ALR_TD for k
    real(f8), allocatable :: int_TDWP(:,:,:,:)
    real(f8), allocatable :: int_TDWPm(:,:,:,:)
    real(f8), allocatable :: asy_TDWP(:,:,:,:)
    real(f8), allocatable :: asy_TDWPm(:,:,:,:)
    real(f8), allocatable :: lr_TDWP(:,:,:,:)
    real(f8), allocatable :: lr_TDWPm(:,:,:,:)
    complex(c8), allocatable :: intAB_TIDWF(:,:,:,:)
    complex(c8), allocatable :: intAC_TIDWF(:,:,:,:)
!> Chebyshev 
    real(f8) :: Hplus, Hminus
    real(f8), allocatable :: ChebyAngle(:)
!> RCB intermediate coordinates
!> 1 - Z_r/r_r, 2 - rij/Zij, 2 - costheta_r
!> 4 - r_p, 5 - costheta_p, 6 - beta, 7 - gij
    real(f8), allocatable :: InterCoor_AB(:,:), InterCoor_AC(:,:)
!> 1 - iZ/ir, 2 - ith
    integer, allocatable :: idXY_AB(:,:), idXY_AC(:,:)
!> RCB Uij matrices 
    real(f8), allocatable :: Uij_AB(:,:), Uij_AC(:,:)
!> AB and AC jK pair
    integer, allocatable :: AB_jKPair(:,:), AC_jKPair(:,:)
    integer, allocatable :: AB_seqjK(:,:), AC_seqjK(:,:)
!> Product vib-rotational WF in intermediate coordinates, IC for interCoor
    real(f8), allocatable :: ICWFvj_AB(:,:,:,:), ICWFvj_AC(:,:,:,:)
!> AtDMat for product channels
    real(f8), allocatable :: AtDMat_AB(:,:,:), AtDMat_AC(:,:,:)
!> Trans matrix from Ka to Kb
    real(f8), allocatable :: TKMat_AB(:,:,:), TKMat_AC(:,:,:)
!> Orbital angular momentum in product channels
    integer, allocatable :: qn_l_AB(:), qn_l_AC(:)
!> SF to BF transformation matrix in product channels
    real(f8), allocatable :: BLK_AB(:,:), BLK_AC(:,:)
!> S matrix in SF
    complex(c8), allocatable :: Smat_AB(:,:,:,:), Smat_AC(:,:,:,:)
!> Inelastic S matrix in SF
    complex(c8), allocatable :: Smat_BC(:,:,:,:)
!> Orbital angular momentum and BF to SF for inelastic channel
    integer, allocatable :: qn_l_BC(:)
    real(f8), allocatable :: BLK_BC(:,:)
!> Tpyes declarations
    type(initWP_class) :: initWP
    type(IALR_class) :: IALR 
    type(Vabs_class) :: Vabs 
    type(channel1_class) :: channel1
    type(channel2_class) :: channel2 

!> ========== Namelists ==========

    namelist /task/ nChannel, IF_inelastic, IF_modHankel, inePos, IDflux, fluxPos, &
                    Atoms, nPES, energyUnit, potentialType, outfile
    namelist /energy/ E_range, dE, TMaxCut, VMaxCut
    namelist /initWavePacket/ initWP
    namelist /IALRset/ IALR
    namelist /VabsAndDump/ Vabs
    namelist /propagation/ timeTot, timePrint
    namelist /productChannel1/ channel1
    namelist /productChannel2/ channel2


contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine initPara()
        implicit none
        integer :: inpFileUnit

!> ========== Read parameters from input file ==========

        open(newunit=inpFileUnit, file="ABC.inf", status='old')
        read(inpFileUnit, nml=task) 
        read(inpFileUnit, nml=initWavePacket)
        read(inpFileUnit, nml=IALRset)
        read(inpFileUnit, nml=propagation)
        read(inpFileUnit, nml=energy) 
        read(inpFileUnit, nml=VabsAndDump)
        read(inpFileUnit, nml=productChannel1)

        rewind(inpFileUnit)

        if (nChannel == 2) read(inpFileUnit, nml=productChannel2)
        rewind(inpFileUnit)

        close(inpFileUnit)

        open(newunit=outFileUnit, file=trim(outfile)//".out", status='replace')
        !> Should add more output files later

        if (initWP%tpar == 1) then
            initWP%Kmin = 0
        else if (initWP%tpar == -1) then
            initWP%Kmin = 1
        else
            write(outFileUnit,*) "Error: total parity tpar must be 1 or -1."
            write(outFileUnit,*) "POSITION: globalPara.f90, subroutine initPara()"
            stop
        end if

        !> nChannel=0: inelastic only, force IF_inelastic
        if (nChannel == 0) then
            IF_inelastic = .true.
        end if

        if (initWP%tpar /= (-1)**(initWP%Jtot+initWP%j0+initWP%l0)) then 
            write(outFileUnit,*) "Error: total pairty incoorent with l0."
            write(outFileUnit,*) "POSITION: globalPara.f90, subroutine initPara()"
            stop
        end if

        if (initWP%l0 > initWP%j0+initWP%Jtot) then 
            write(outFileUnit,*) "Error: l0 > Jtot + j0."
            write(outFileUnit,*) "POSITION: globalPara.f90, subroutine initPara()"
            stop
        end if

        if (initWP%l0 < abs(initWP%j0-initWP%Jtot)) then 
            write(outFileUnit,*) "Error: l0 < |Jtot - j0|."
            write(outFileUnit,*) "POSITION: globalPara.f90, subroutine initPara()"
            stop
        end if

        call getMass()
        call diatomParity()
        call energySet()

!> ========== Output ==========

        write(outFileUnit,'(1x,a)') " =====> Input parameters <======"
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a,a,a)') "Reaction Channel: ", &
            merge("Inelastic only (A+BC)       ", &
            merge("A + BC -> AB + C            ", "A + BC -> AB + C and AC + B ", nChannel==1), nChannel==0)
        write(outFileUnit,'(1x,6a)') "Atoms A, B, C: ", Atoms(1),', ', Atoms(2),', ', Atoms(3)
        write(outFileUnit,'(1x,a,f15.9)') "Reduced Masses of BC (a.u.):", massBC
        write(outFileUnit,'(1x,a,f15.9)') "Reduced Masses of ABC (a.u.):", massTot
        write(outFileUnit,'(1x,a,3i4)') "v0, j0, l0 : ", initWP%v0, initWP%j0, initWP%l0
        write(outFileUnit,'(1x,a,i4)') "Total angular momentum Jtot: ", initWP%Jtot
        write(outFileUnit,'(1x,a,i4)') "Total parity: ", initWP%tpar
        write(outFileUnit,'(1x,a,i4)') "Number of PES: ", nPES
        write(outFileUnit,'(1x,a,i4)') "Initial PES: ", initWP%PES0
        write(outFileUnit,'(1x,a)') ''

    end subroutine initPara
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine energySet()
        implicit none
        integer :: iEtot
!> ========== Construct collision energy ==========

        nEtot = int((E_range(2) - E_range(1))/dE) + 1
        allocate(Ecol(nEtot),Etot(nEtot))
        allocate(initAM(nEtot))
        allocate(kReact(nEtot))
        allocate(ChebyAngle(nEtot))
        do iEtot = 1, nEtot 
            Ecol(iEtot) = E_range(1) + (iEtot-1)*dE 
            if (Ecol(iEtot) > E_range(2)) then
                Ecol(iEtot) = E_range(2)
            end if 
            !> convert to au
            Ecol(iEtot) = Ecol(iEtot) * energyUnitTrans(energyUnit)
            kReact(iEtot) = dsqrt(2.0_f8*massTot*Ecol(iEtot))
        end do
        initWP%Ec = initWP%Ec * energyUnitTrans(energyUnit)
        TMaxCut = TMaxCut * energyUnitTrans(energyUnit)
        VMaxCut = VMaxCut * energyUnitTrans(energyUnit)

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
                write(outFileUnit,*) "POSITION: globalPara.f90, subroutine getMass()"
                stop
            else
                atomMass(iatm) = elemMass(ipos)*amu2au
            end if
        end do

        !> Reduce mass of BC : mB*mC/(mB+mC)
        !> Reduce mass of ABC : mA*(mB+mC)/(mA+mB+mC)
        massBC = atomMass(2)*atomMass(3)/(atomMass(2)+atomMass(3))
        massTot = atomMass(1)*(atomMass(2)+atomMass(3))/(atomMass(1)+atomMass(2)+atomMass(3))

    end subroutine getMass
!> ------------------------------------------------------------------------------------------------------------------ <!
    
!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine diatomParity()
        implicit none
        !> jpar = -1, 0, 1, jmin = 1, 0, 0, jinc = 2, 1, 2
        integer :: jmin(-1:1) = [1, 0, 0]
        integer :: jinc(-1:1) = [2, 1, 2]
        real(f8) :: eps 
        
        eps = 100.0_f8 * epsilon(atomMass(1))

        !>if (abs(atomMass(2)-atomMass(3)) >= eps) then
            !> BC is heteronuclear
        !>    if (initWP%jpar /= 0) then
        !>        write(outFileUnit,*) "Error: For heteronuclear, BC parity jpar must be 0."
        !>        write(outFileUnit,*) "POSITION: globalPara.f90, subroutine diatomParity()"
        !>        stop
        !>    end if
        !>end if
        initWP%jmin = jmin(initWP%jpar)
        initWP%jinc = jinc(initWP%jpar)

        if (abs(atomMass(1)-atomMass(2)) >= eps) then
            !> A and B are heteronuclear
            if (channel1%jpar /= 0) then
                write(outFileUnit,*) "Error: For heteronuclear, AB parity jpar must be 0."
                write(outFileUnit,*) "POSITION: globalPara.f90, subroutine diatomParity()"
                stop
            end if
        end if
        channel1%jmin = jmin(channel1%jpar)
        channel1%jinc = jinc(channel1%jpar)

        if (nChannel == 2) then
            if (abs(atomMass(1)-atomMass(3)) >= eps) then
                !> A and C are heteronuclear
                if (channel2%jpar /= 0) then
                    write(outFileUnit,*) "Error: For heteronuclear, AC parity jpar must be 0."
                    write(outFileUnit,*) "POSITION: globalPara.f90, subroutine diatomParity()"
                    stop
                end if
            end if
            channel2%jmin = jmin(channel2%jpar)
            channel2%jinc = jinc(channel2%jpar)
        end if

    end subroutine diatomParity
!> ------------------------------------------------------------------------------------------------------------------ <!

end module m_gPara
