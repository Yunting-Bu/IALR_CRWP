module m_Potent
    use m_MachinaBasic, only : f8, BinReadWrite 
    use m_gPara
    implicit none

    public
    private :: setVabs, setPot
    
contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine Jacobi2Bond(Z, r, theta, massB, massC, bond)
        implicit none
        real(f8), intent(in) :: Z, r, theta 
        real(f8), intent(in) :: massB, massC 
        real(f8), intent(inout) :: bond(3)
        real(f8) :: rOB, rOC

!> bond(1) = rAB, bond(2) = rBC, bond(3) = rAC
!> theta in radian
        bond(2) = r
        rOB = massC / (massB+massC) * r 
        rOC = massB / (massB+massC) * r 
        bond(1) = dsqrt(Z*Z+rOB*rOB-2.0_f8*Z*rOB*dcos(theta))
        bond(3) = dsqrt(Z*Z+rOC*rOC+2.0_f8*Z*rOC*dcos(theta))

    end subroutine Jacobi2Bond
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine bond2Jacobi(bond, Z, r, theta, massB, massC)
        implicit none
        real(f8), intent(in) :: bond(3)
        real(f8), intent(in) :: massB, massC 
        real(f8), intent(out) :: Z, r, theta 
        real(f8) :: r1, r2, a1, a2

!> theta in radian
        r = bond(2)
        r1 = massB * r / (massB+massC)
        r2 = massC * r / (massB+massC)

        a1 = (bond(1)**2+r**2-bond(3)**2) / (2.0_f8*r*bond(1))
        Z = dsqrt(bond(1)**2+r1**2-2.0_f8*r1*bond(1)*a1)
        a2 = (Z**2+r2**2-bond(3)**2) / (2.0_f8*r2*Z)

        if (a2 > 1.0_f8) a2 = 1.0_f8
        if (a2 < -1.0_f8) a2 = -1.0_f8
        theta = acos(a2)
        if (Z == 0.0_f8) theta = pi / 2.0_f8
    
    end subroutine bond2Jacobi
!> ------------------------------------------------------------------------------------------------------------------ <!


!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine diagDiaVmat(bond, AtDMat, adiaV)
        implicit none
        real(f8), intent(in) :: bond(3)
        real(f8), intent(inout) :: AtDMat(nPES,nPES), adiaV(nPES)
        real(f8) :: diaV(nPES,nPES)
        real(f8), allocatable :: work(:)
        integer :: lwork, info

!> bond(3) is the three bonds coordinate
!> POT0 is the common PES interface
        call POT0(nPES,diaV,bond)
        if (nPES > 1) then 
            allocate(work(1))
            call dsyev('V', 'U', nPES, diaV, nPES, adiaV, work, -1, info)
            lwork = int(work(1))
            deallocate(work)

            allocate(work(lwork))
            call dsyev('V', 'U', nPES, diaV, nPES, adiaV, work, lwork, info)
            if (info /= 0) then
                write(outFileUnit,*) "Error in DSYEV of diaV, info =", info
                write(outFileUnit,*) "POSITION: potent.f90, subroutine diagDiaVmat()"
                stop
            end if
            deallocate(work)
        else 
            AtDMat = 1.0_f8
            adiaV(1) = diaV(1,1)
        end if 
    
    end subroutine diagDiaVmat
!> ------------------------------------------------------------------------------------------------------------------ <!
    
!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setPot_vj0(r,pot_vj0)
        implicit none
        real(f8), intent(in) :: r 
        real(f8), intent(out) :: pot_vj0
        real(f8) :: bond(3), AtD(nPES,nPES), Vadia(nPES)
        real(f8), parameter :: longDistance = 30.0_f8
        
!> In this situation, bond(2) is the length of BC
!> You should check the PES interface!
        bond(1) = longDistance
        bond(2) = r
        bond(3) = longDistance

        call diagDiaVmat(bond, AtD, Vadia)
        pot_vj0 = Vadia(initWP%PES0) 
        pot_vj0 = pot_vj0 + initWP%j0*(initWP%j0+1.0_f8)/(2.0_f8*massBC*r*r)

    end subroutine setPot_vj0
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setPot_Product(type, rp, Vpot)
        implicit none
        character(len=*), intent(in) :: type
        real(f8), intent(in) :: rp 
        !> Adiabatic potential
        real(f8), intent(out) :: Vpot(nPES)
        real(f8) :: AtD(nPES,nPES)
        real(f8) :: bond(3)
        real(f8), parameter :: longDistance = 30.0_f8
        integer :: IDreac

        bond = longDistance
!> C + AB, bond(1) = rAB
!> B + AC, bond(3) = rAC
        if (type == 'A+BC->C+AB') then 
            IDreac = 1
            bond(IDreac) = rp
        else if (type == 'A+BC->B+AC') then
            IDreac = 3
            bond(IDreac) = rp
        else
            write(outFileUnit,'(1x,a)') 'Error: unknown reaction type in setPotProduct!'
            write(outFileUnit,'(1x,a)') 'POSITION: potent.f90, subroutine setPotProduct()'
            stop
        end if

        call diagDiaVmat(bond, AtD, Vpot)

    end subroutine setPot_Product
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setVabs(range, Cabs, nGird, grid, nabs, VabsMat)
        implicit none
        real(f8), intent(in) :: range, Cabs 
        integer, intent(in) :: nGird 
        real(f8), intent(in) :: grid(nGird)
        integer, intent(out) :: nabs
        real(f8), allocatable, intent(inout) :: VabsMat(:)
        real(f8) :: rangeAll, rangeNoAbs
        integer :: i

!> Since DVR grid don't have the boundary point
        rangeAll = grid(nGird) + grid(2) - grid(1)
        rangeNoAbs = rangeAll - range
        nabs = 0
        do i = 1, nGird
            if (grid(i) >= rangeNoAbs) then
                !VabsMat(i) = dexp(-Cabs * ((grid(i)-rangeNoAbs)/(rangeAll - rangeNoAbs))**2)
                VabsMat(i) = dexp(-Cabs * ((grid(i)-rangeNoAbs)**2))
            else
                VabsMat(i) = 1.0_f8
                nabs = i
            end if
        end do

        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a)') 'Range of Vabs: [', rangeNoAbs, ', ', rangeAll, ' ]' 
        write(outFileUnit,'(1x,a,i5,a)') 'When nDVR >=', nabs, ', the absorbing potential starts to work.'
        write(outFileUnit,'(1x,a)') 'Please ensure that the intermediate coordinate lies outside the absorbing region !!'
        write(outFileUnit,'(1x,a)') ''

    end subroutine setVabs
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine initAllVabs()
        implicit none

        write(outFileUnit,'(1x,a)') '=====> Vasb information <====='
        write(outFileUnit,'(1x,a)') ''
!> Vabs in r
        allocate(int_rFabs(IALR%vint))
        write(outFileUnit,'(1x,a)') 'Absorbing potential in r:'
        call setVabs(Vabs%rabs_range,Vabs%Cr,IALR%vint,r_Int,Vabs%nrint,int_rFabs) 
!> Vabs in Z_int
!        allocate(int_ZFabs(IALR%nZ_I))
!        write(outFileUnit,'(1x,a)') 'Absorbing potential in Z_int: '
!        call setVabs(Vabs%Zint_range,Vabs%Cint,IALR%nZ_I,Z_IALR(1:IALR%nZ_I),Vabs%nZint,int_ZFabs)
!        write(101,*) int_ZFabs
!> Vabs in long-range
        allocate(lr_ZFabs(IALR%nZ_IALR))
        write(outFileUnit,'(1x,a)') 'Absorbing potential in Z_lr:'
        call setVabs(Vabs%Zlr_range,Vabs%Clr,IALR%nZ_IALR,Z_IALR,Vabs%nZlr,lr_ZFabs)
        write(102,*) lr_ZFabs
!> Vabs in asymptotic
        allocate(asy_ZFabs(IALR%nZ_IA))
        write(outFileUnit,'(1x,a)') 'Absorbing potential in Z_asy:'
        write(outFileUnit,'(1x,a)') 'Note that the Vabs only works for the channel with (v, j, iPES) /= (v0, j0, initPES)!'
        call setVabs(Vabs%Zasy_range,Vabs%Casy,IALR%nZ_IA,Z_IALR(1:IALR%nZ_IA),Vabs%nZasy,asy_ZFabs)
        write(outFileUnit,'(1x,a)') ''
    end subroutine initAllVabs
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setPot(region, nZ, ZGrid, nr, rGrid, nth, thGrid, Vpot)
        implicit none
        character(len=*), intent(in) :: region
        integer, intent(in) :: nZ, nr, nth
        real(f8), intent(in) :: ZGrid(nZ), rGrid(nr), thGrid(nth)
        real(f8), intent(out) :: Vpot(nPES,nPES,nZ,nr,nth)
        real(f8):: Vdia(nPES,nPES), Vadia(nPES), AtDMat(nPES,nPES)
        real(f8) :: bond(3), VBC(nr)
        real(f8) :: Z, r, theta
        integer :: iPES, jPES, kPES
        integer :: iZ, ir, ith
        character(len=100) :: fileName

        do ir = 1, nr
            r = rGrid(ir)
            call setPot_vj0(r,VBC(ir))
        end do

        do iZ = 1, nZ
            Z = ZGrid(iZ)
            do ir = 1, nr
                r = rGrid(ir)
                do ith = 1, nth
                    theta = thGrid(ith)
                    call Jacobi2Bond(Z, r, theta, atomMass(2), atomMass(3), bond)
                    call POT0(nPES,Vdia,bond)

                    if (nPES == 1) then 
                        Vpot(1,1,iZ,ir,ith) = min((Vdia(1,1)-VBC(ir)), VMaxCut)
                    else 
                        call diagDiaVmat(bond, AtDMat, Vadia)
                        do iPES = 1, nPES
                            Vadia(iPES) = min((Vadia(iPES)-VBC(ir)), VMaxCut)
                        end do 
                        do iPES = 1, nPES
                            do jPES = 1, nPES
                                Vpot(iPES,jPES,iZ,ir,ith) = 0.0_f8
                                do kPES = 1, nPES
                                    Vpot(iPES,jPES,iZ,ir,ith) = Vpot(iPES,jPES,iZ,ir,ith) + &
                                                                AtDMat(kPES,iPES) * Vadia(kPES) * AtDMat(kPES,jPES)
                                end do
                            end do
                        end do
                    end if

                end do
            end do
        end do

        fileName = trim(region)//'_Vdia_'//trim(outfile)//'.bin'
        call BinReadWrite(fileName, Vpot, 'write')

    end subroutine setPot
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setPot_IALR()
        implicit none
        integer :: asy_nZ, lr_nZ 
        real(f8), allocatable :: Z_asy(:), Z_lr(:)
        character(len=3) :: region
        character(len=100) :: fileName
        integer :: iZ, ir
        real(f8) :: Z, r, theta, bond(3)
        real(f8) :: VBC, Vdia(1,1)
        
        asy_nZ = IALR%nZ_IA - IALR%nZ_I 
        lr_nZ = IALR%nZ_IALR - IALR%nZ_IA 
        allocate(Z_asy(asy_nZ), Z_lr(lr_nZ))
!> Z grid in asymptotic region
        Z_asy = Z_IALR(IALR%nZ_I+1:IALR%nZ_IA)
!> Z grid in long-range region
        Z_lr = Z_IALR(IALR%nZ_IA+1:IALR%nZ_IALR)

        allocate(int_Vdia(nPES,nPES,IALR%nZ_I,IALR%vint,IALR%int_nA))
        allocate(asy_Vdia(nPES,nPES,asy_nZ,IALR%vasy,IALR%asy_nA))
        allocate(lr_Vdia(lr_nZ,IALR%vlr))

!> interaction region
        region = 'int'
        if (trim(potentialType) == 'New') then
            write(outFileUnit,'(1x,a)') 'Calculating interaction potential in the interaction region...'
            call setPot(region, IALR%nZ_I, Z_IALR(1:IALR%nZ_I), IALR%vint, r_Int, IALR%int_nA, int_ANode, int_Vdia)
        else if (trim(potentialType) == 'Read') then
            fileName = trim(region)//'_Vdia_'//trim(outfile)//'.bin'
            call BinReadWrite(fileName, int_Vdia, 'read')
            write(outFileUnit,'(1x,a,a)') 'Read int_Vdia in file: ', trim(fileName)
        else
            write(outFileUnit,'(1x,a)') 'Error: unknown potentialType, must Read or New !'
            write(outFileUnit,'(1x,a)') 'POSITION: potent.f90, subroutine setPot_IALR()'
            stop
        end if
!> asymptotic region
        region = 'asy'
        if (trim(potentialType) == 'New') then
            write(outFileUnit,'(1x,a)') 'Calculating interaction potential in the asymptotic region...'
            call setPot(region, asy_nZ, Z_asy, IALR%vasy, r_Asy, IALR%asy_nA, asy_ANode, asy_Vdia)
        else if (trim(potentialType) == 'Read') then
            fileName = trim(region)//'_Vdia_'//trim(outfile)//'.bin'
            call BinReadWrite(fileName, asy_Vdia, 'read')
            write(outFileUnit,'(1x,a,a)') 'Read asy_Vdia in file: ', trim(fileName)
        else
            write(outFileUnit,'(1x,a)') 'Error: unknown potentialType, must Read or New !'
            write(outFileUnit,'(1x,a)') 'POSITION: potent.f90, subroutine setPot_IALR()'
            stop
        end if
!> long-range region
        region = 'lr'
        fileName = trim(region)//'_Vdia_'//trim(outfile)//'.bin'
        if (trim(potentialType) == 'New') then
            write(outFileUnit,'(1x,a)') 'Calculating interaction potential in the long-range region...'
            do  iZ = 1, lr_nZ
                Z = Z_lr(iZ)
                do ir = 1, IALR%vlr
                    r = r_LR(ir)
                    theta = 0.0_f8
                    call Jacobi2Bond(Z, r, theta, atomMass(2), atomMass(3), bond)
                    call POT0(1, Vdia, bond)
                    call setPot_vj0(r, VBC)
                    lr_Vdia(iZ, ir) = min((Vdia(1,1)-VBC), VMaxCut)
                end do 
            end do
            call BinReadWrite(fileName, lr_Vdia, 'write')
        else if (trim(potentialType) == 'Read') then
            call BinReadWrite(fileName, lr_Vdia, 'read')
            write(outFileUnit,'(1x,a,a)') 'Read lr_Vdia in file: ', trim(fileName)
        else
            write(outFileUnit,'(1x,a)') 'Error: unknown potentialType, must Read or New !'
            write(outFileUnit,'(1x,a)') 'POSITION: potent.f90, subroutine setPot_IALR()'
            stop
        end if
        write(outFileUnit,'(1x,a)') ''

        deallocate(Z_asy, Z_lr)

    end subroutine setPot_IALR
!> ------------------------------------------------------------------------------------------------------------------ <!
end module m_Potent