program CrossSection
    use m_MachinaBasic , only : f8, c8
    use iso_fortran_env, only : iostat_end
    implicit none

    !> Constants
    real(f8), parameter :: au2ev = 27.211386245988_f8
    real(f8), parameter :: ev2au = 1.0_f8 / au2ev
    real(f8), parameter :: au2SI = 6.12612d-9
    real(f8), parameter :: pi = 4.0_f8*atan(1.0_f8)
    real(f8), parameter:: kB = 3.1668d-6
    complex(c8), parameter :: img = (0.0_f8, 1.0_f8)
    complex(c8), parameter :: imgZero = (0.0_f8, 0.0_f8)

    integer, parameter :: maxFinal = 200

    integer :: nEtot, nJtot, JtotMax
    integer :: v0, j0, l0
    integer :: ichnl, nTemp, nFinal
    integer :: nLr, nLp, nLp_max, iPES
    integer :: tpar, jfMax
    integer :: vRead, jRead
    integer :: inpFileUnit, outFileUnit
    integer :: iJtot, iEtot, imj, imjp
    integer :: iLr, iLp, ith, ipar, iml
    integer :: io_stat, iline
    integer :: iTemp, iFinal
    integer :: vf_cur, jf_cur

    integer :: vf(maxFinal), jf(maxFinal)

    real(f8) :: tms, theta, Lr
    real(f8) :: tmp1, tmp2
    real(f8) :: Sreal, Simag
    real(f8) :: temp, logStep
    real(f8) :: tempRange(2)
    real(f8) :: s2, tempFact, BoltzWeight
    real(f8), allocatable :: y2a(:)
    real(f8), allocatable :: Ecol(:), Ecol_au(:), xint(:), yint(:), csint(:)
    real(f8), allocatable :: kvj(:)
    real(f8), allocatable :: ICS(:), prop(:,:)
    real(f8), allocatable :: DCS_SF(:,:)
    real(f8), allocatable :: temper(:), rate(:)
    real(f8), allocatable :: CG1(:,:,:), CG2(:,:,:,:)
    real(f8), allocatable :: Smat(:,:,:,:,:,:)
    real(f8), external :: CG

    complex(c8) :: Tmat
    complex(c8), allocatable :: Ylm(:,:,:)
    complex(c8), allocatable :: amplitude_SF(:,:,:,:)
    complex(c8), external :: spharmonic

    character(len=2) :: reactChannel(0:2) = ['BC','AB','AC']
    character(len=100) :: fileName, dirName, SmatDir, outFile
    character(len=200) :: chartmp, tmpStr

    logical :: f_exists, hasData
    logical :: doAmplitude, doDCS, doProb, doICS, doRate

    namelist /stateInfo/ ichnl, v0, j0, JtotMax, iPES, nEtot, nTemp, tempRange, &
                        nFinal, vf, jf, SmatDir, &
                        doAmplitude, doDCS, doProb, doICS, doRate

    ! Defaults
    nFinal = 1
    vf = 0
    jf = 0
    nTemp = 0
    tempRange = [0.001_f8, 100.0_f8]
    doAmplitude = .true.
    doDCS = .true.
    doProb = .true.
    doICS = .true.
    doRate = .true.

    open(newunit=inpFileUnit, file='CrossSection.inf', status='old')
    read(inpFileUnit, nml=stateInfo)
    close(inpFileUnit)

    ! Handle dependencies: auto-enable prerequisites
    if (doDCS) doAmplitude = .true.
    if (doRate) doICS = .true.
    if (doICS) doProb = .true.

    nJtot = JtotMax
    nLr = nJtot + j0
    jfMax = maxval(jf(1:nFinal))
    nLp_max = nJtot + jfMax

    ! Allocate arrays
    allocate(kvj(nEtot))
    allocate(Ecol(nEtot))
    allocate(Smat(nEtot,0:nJtot,0:nLr,0:nLp_max,2,2))

    if (doAmplitude) then
        allocate(CG1(0:nJtot,-j0:j0,0:nLr))
        allocate(CG2(0:nJtot,-jfMax:jfMax,0:nLp_max,-nLp_max:nLp_max))
        allocate(Ylm(0:180,0:nLp_max,-nLp_max:nLp_max))
        allocate(amplitude_SF(nEtot,0:180,-j0:j0,-jfMax:jfMax))
    end if

    if (doDCS) allocate(DCS_SF(nEtot,0:180))

    if (doProb) then
        allocate(ICS(nEtot))
        allocate(prop(nEtot,0:nJtot))
    end if

    ! ------- First pass: read tms and Ecol from any S-matrix file -------
    tms = 0.0_f8
    Ecol = 0.0_f8
    call readFirstSmat(trim(SmatDir), reactChannel(ichnl), j0, nJtot, iPES, nEtot, tms, Ecol)

    ! Compute kvj (independent of final state)
    do iEtot = 1, nEtot
        kvj(iEtot) = sqrt(2.0_f8 * Ecol(iEtot) * ev2au * tms)
    end do

    write(*,'(a,f12.4,a)') ' Reduced mass (A+BC) = ', tms, ' au'
    write(*,'(a,f8.3,a,f8.3,a)') ' Energy range: ', Ecol(1), ' to ', Ecol(nEtot), ' eV'
    write(*,'(a,i0,a)') ' Computing ', nFinal, ' final state(s)'

    ! ------- Pre-compute CG1 and Ylm (only needed for amplitude/DCS) -------
    if (doAmplitude) then
        CG1 = 0.0_f8
        do iJtot = 0, nJtot
            do imj = -j0, j0
                do iLr = 0, nLr
                    if (iJtot > j0+iLr .or. iJtot < abs(j0-iLr)) cycle
                    CG1(iJtot,imj,iLr) = CG(j0, imj, iLr, 0, iJtot)
                end do
            end do
        end do

        Ylm = imgZero
        do ith = 0, 180
            do iLp = 0, nLp_max
                do iml = -nLp_max, nLp_max
                    if (abs(iml) > iLp) cycle
                    theta = ith * pi / 180.0_f8
                    Ylm(ith,iLp,iml) = spharmonic(theta, 0.0_f8, iLp, iml)
                end do
            end do
        end do
    end if

    ! ------- Pre-compute temperature grid (only needed for rate) -------
    if (doRate .and. nTemp > 0) then
        allocate(temper(nTemp), rate(nTemp))
        allocate(xint(nEtot), yint(nEtot), csint(nEtot))
        allocate(Ecol_au(nEtot))
        allocate(y2a(nEtot))
        Ecol_au = Ecol * ev2au

        if (nTemp == 1) then
            temper(1) = tempRange(1)
        else
            if (tempRange(1) > 0.0_f8 .and. tempRange(2) > tempRange(1)) then
                logStep = log(tempRange(2)/tempRange(1)) / real(nTemp-1, f8)
                do iTemp = 1, nTemp
                    temper(iTemp) = tempRange(1) * exp(logStep * real(iTemp-1, f8))
                end do
            else
                do iTemp = 1, nTemp
                    temper(iTemp) = tempRange(1) + (tempRange(2)-tempRange(1)) * real(iTemp-1, f8) / real(nTemp-1, f8)
                end do
            end if
        end if
    end if

    ! ===================== Loop over final states =====================
    do iFinal = 1, nFinal
        vf_cur = vf(iFinal)
        jf_cur = jf(iFinal)
        nLp = nJtot + jf_cur

        write(*,'(a,i0,a,i0,a,i0,a)') ' --- Final state ', iFinal, ': vf=', vf_cur, ', jf=', jf_cur, ' ---'

        ! ------- Read S-matrix for this (vf, jf) -------
        call read_Smatrix()

        if (.not. hasData) then
            write(*,*) '  No S-matrix data found, skipping.'
            cycle
        end if

        if (doAmplitude) call calc_scatt_amplitude()
        if (doDCS)       call calc_DCS()
        if (doProb)      call calc_probability()
        if (doICS)       call calc_ICS()
        if (doRate .and. nTemp > 0) call calc_rate_constant()

    end do ! iFinal loop

    write(*,*) 'Done.'

contains

    !> Read S-matrix for current (vf_cur, jf_cur) final state
    subroutine read_Smatrix()
        implicit none

        Smat = 0.0_f8
        hasData = .false.

        do tpar = -1, 1, 2
            if (tpar == -1) ipar = 1
            if (tpar == 1) ipar = 2
            write(dirName, '(a,a,i0,a)') trim(SmatDir), '/p', tpar, '/'
            do iJtot = 0, nJtot
                do l0 = abs(iJtot-j0), iJtot+j0
                    write(fileName,'(a,a,i0,a,i0,a,i0,a)') reactChannel(ichnl), '_J', iJtot, 'L', l0, 'iPES', iPES, '.Smat'
                    inquire(file=trim(dirName)//trim(fileName), exist=f_exists)
                    if (.not. f_exists) cycle
                    write(*,'(a,a,a,a)') '  Reading dir: ', trim(dirName), '  file: ', trim(fileName)
                    open(newunit=inpFileUnit, file=trim(dirName)//trim(fileName), status='old')
                    do iline = 1, 2
                        read(inpFileUnit,'(a)') chartmp
                    end do
                    read(inpFileUnit,'(a)') chartmp
                    read(inpFileUnit,'(a)') chartmp
                    do iline = 1, 5
                        read(inpFileUnit,'(a)') chartmp
                    end do
                    do
                        read(inpFileUnit,'(a)', iostat=io_stat) chartmp
                        if (io_stat /= 0) exit
                        if (len_trim(chartmp) == 0) cycle
                        if (chartmp(1:1) == '#') cycle
                        backspace(inpFileUnit)
                        read(inpFileUnit,*, iostat=io_stat) iEtot, tmp1, vRead, jRead, Lr, tmp2, Sreal, Simag, tmp2
                        if (io_stat /= 0) cycle
                        if (vRead == vf_cur .and. jRead == jf_cur) then
                            Smat(iEtot,iJtot,l0,nint(Lr),ipar,1) = Sreal
                            Smat(iEtot,iJtot,l0,nint(Lr),ipar,2) = Simag
                            hasData = .true.
                        end if
                    end do
                    close(inpFileUnit)
                end do
            end do
        end do
    end subroutine read_Smatrix

    !> Compute scattering amplitude in space-fixed frame
    subroutine calc_scatt_amplitude()
        implicit none

        ! Compute CG2 Clebsch-Gordan coefficients for current jf_cur
        CG2 = 0.0_f8
        do iJtot = 0, nJtot
            do imj = -j0, j0
                do imjp = -jf_cur, jf_cur
                    do iLr = 0, nLr
                        do iLp = 0, nLp
                            iml = imj - imjp
                            if (abs(iml) > iLp) cycle
                            if (iJtot > jf_cur+iLp .or. iJtot < abs(jf_cur-iLp)) cycle
                            CG2(iJtot,imjp,iLp,iml) = CG(jf_cur, imjp, iLp, iml, iJtot)
                        end do
                    end do
                end do
            end do
        end do

        ! Compute scattering amplitude
        amplitude_SF = imgZero
        do iEtot = 1, nEtot
            do ith = 0, 180
                do imj = -j0, j0
                    do imjp = -jf_cur, jf_cur
                        iml = imj - imjp
                        amplitude_SF(iEtot,ith,imj,imjp) = imgZero
                        do iJtot = 0, nJtot
                            do ipar = 1, 2
                                do iLr = 0, nLr
                                    do iLp = 0, nLp
                                        if (abs(iml) > iLp) cycle
                                        Tmat%re = -Smat(iEtot,iJtot,iLr,iLp,ipar,1)
                                        Tmat%im = -Smat(iEtot,iJtot,iLr,iLp,ipar,2)
                                        ! Elastic: T = delta_{l,l'} - S; Inelastic: T = -S
                                        if (v0 == vf_cur .and. j0 == jf_cur .and. iLr == iLp) then
                                            Tmat%re = Tmat%re + 1.0_f8
                                        end if
                                        if (Tmat == imgZero) cycle
                                        if (CG1(iJtot,imj,iLr) == 0.0_f8 .or. CG2(iJtot,imjp,iLp,iml) == 0.0_f8) cycle
                                        amplitude_SF(iEtot,ith,imj,imjp) = amplitude_SF(iEtot,ith,imj,imjp) &
                                                                            + dsqrt(2.0_f8*iLr+1.0_f8) * img**(iLr-iLp+1) &
                                                                            * Ylm(ith,iLp,iml) * Tmat * dsqrt(pi) &
                                                                            * CG1(iJtot,imj,iLr) * CG2(iJtot,imjp,iLp,iml)
                                    end do
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine calc_scatt_amplitude

    !> Compute differential cross section and write output
    subroutine calc_DCS()
        implicit none

        do iEtot = 1, nEtot
            do ith = 0, 180
                DCS_SF(iEtot,ith) = 0.0_f8
                do imj = -j0, j0
                    do imjp = -jf_cur, jf_cur
                        DCS_SF(iEtot,ith) = DCS_SF(iEtot,ith) &
                                            + amplitude_SF(iEtot,ith,imj,imjp)%re**2 &
                                            + amplitude_SF(iEtot,ith,imj,imjp)%im**2
                    end do
                end do
                DCS_SF(iEtot,ith) = 2.0_f8*pi*DCS_SF(iEtot,ith) / (kvj(iEtot)**2*(2.0_f8*j0+1.0_f8))
            end do
        end do

        write(outFile,'(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a)') &
            'DCS_',reactChannel(ichnl), 'v', v0, 'j', j0, 'vf', vf_cur, 'jf', jf_cur, 'iPES', iPES, '.dat'
        open(newunit=outFileUnit, file=trim(outFile), status='replace')
        write(outFileUnit,'(a)') '#theta(deg)    DCS_SF(bohr^2/sr)'
        do ith = 0, 180
            write(outFileUnit,'(f15.9,1x,*(f15.9))') ith*1.0_f8, (DCS_SF(iEtot,ith), iEtot=1,nEtot)
        end do
        close(outFileUnit)
        write(*,*) '  DCS -> ', trim(outFile)
    end subroutine calc_DCS

    !> Compute reaction probability from S-matrix
    subroutine calc_probability()
        implicit none

        do iEtot = 1, nEtot
            do iJtot = 0, nJtot
                prop(iEtot,iJtot) = 0.0_f8
                do ipar = 1, 2
                    do iLr = 0, nLr
                        do iLp = 0, nLp
                            if (Smat(iEtot,iJtot,iLr,iLp,ipar,1) == 0.0_f8 .and. Smat(iEtot,iJtot,iLr,iLp,ipar,2) == 0.0_f8) cycle
                            prop(iEtot,iJtot) = prop(iEtot,iJtot) + Smat(iEtot,iJtot,iLr,iLp,ipar,1)**2 + Smat(iEtot,iJtot,iLr,iLp,ipar,2)**2
                        end do
                    end do
                end do
            end do
        end do
    end subroutine calc_probability

    !> Compute integral cross section from probability and write output
    subroutine calc_ICS()
        implicit none

        ICS = 0.0_f8
        do iEtot = 1, nEtot
            do iJtot = 0, nJtot
                ICS(iEtot) = ICS(iEtot) + prop(iEtot,iJtot) * (2.0_f8*iJtot+1.0_f8)
            end do
            ICS(iEtot) = pi*ICS(iEtot) / (kvj(iEtot)**2*(2.0_f8*j0+1.0_f8))
        end do

        write(outFile,'(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a)') &
            'ICS_',reactChannel(ichnl), 'v', v0, 'j', j0, 'vf', vf_cur, 'jf', jf_cur, 'iPES', iPES, '.dat'
        open(newunit=outFileUnit, file=trim(outFile), status='replace')
        write(outFileUnit,'(a)') '#Ecol(eV)    prop(Jmax)    ICS(bohr^2)'
        do iEtot = 1, nEtot
            write(outFileUnit,'(f15.9,1x,f15.9,1x,f15.9)') Ecol(iEtot), prop(iEtot,nJtot), ICS(iEtot)
        end do
        close(outFileUnit)
        write(*,*) '  ICS -> ', trim(outFile)
    end subroutine calc_ICS

    !> Compute thermal rate constant and write output
    subroutine calc_rate_constant()
        implicit none

        xint = Ecol_au
        yint = ICS
        call spline(xint,yint,nEtot,1.0e30_f8,1.0e30_f8,y2a)
        do iEtot = 1, nEtot
            call splint(xint,yint,y2a,nEtot,Ecol_au(iEtot),s2)
            if (s2 < 0.0_f8) s2 = 0.0_f8
            csint(iEtot) = s2
        end do
        do iTemp = 1, nTemp
            temp = temper(iTemp)
            tempFact = 1.0_f8 / (kB*temp) * sqrt(8.0_f8/(pi*tms*kB*temp))
            tmp2 = 0.0_f8
            do iEtot = 1, nEtot
                BoltzWeight = exp(-Ecol_au(iEtot)/kB/temp)
                tmp2 = tmp2 + csint(iEtot) * BoltzWeight
            end do
            rate(iTemp) = tempFact * tmp2 * (Ecol_au(2)-Ecol_au(1)) * au2SI
        end do

        write(outFile,'(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a)') &
            'Rate_',reactChannel(ichnl), 'v', v0, 'j', j0, 'vf', vf_cur, 'jf', jf_cur, 'iPES', iPES, '.dat'
        open(newunit=outFileUnit, file=trim(outFile), status='replace')
        write(outFileUnit,'(a)') '#Temperature(K)    Rate(cm3*s-1*mole-1)'
        do iTemp = 1, nTemp
            write(outFileUnit,'(f8.3,3x,e15.9)') temper(iTemp), rate(iTemp)
        end do
        close(outFileUnit)
        write(*,*) '  Rate -> ', trim(outFile)
    end subroutine calc_rate_constant

end program CrossSection


subroutine readFirstSmat(SmatDir, chnlName, j0, nJtot, iPES, nEtot, tms, Ecol)
    use m_MachinaBasic, only : f8
    use iso_fortran_env, only : iostat_end
    implicit none
    character(len=*), intent(in) :: SmatDir, chnlName
    integer, intent(in) :: j0, nJtot, iPES, nEtot
    real(f8), intent(out) :: tms, Ecol(nEtot)

    integer :: tpar, iJtot, l0, iEtot, io_stat, iline, inpFileUnit
    integer :: vRead, jRead
    real(f8) :: Lr, tmp1, tmp2, Sreal, Simag
    character(len=100) :: fileName, dirName
    character(len=200) :: chartmp
    logical :: f_exists

    tms = 0.0_f8
    Ecol = 0.0_f8

    do tpar = -1, 1, 2
        write(dirName, '(a,a,i0,a)') trim(SmatDir), '/p', tpar, '/'
        do iJtot = 0, nJtot
            do l0 = abs(iJtot-j0), iJtot+j0
                write(fileName,'(a,a,i0,a,i0,a,i0,a)') chnlName, '_J', iJtot, 'L', l0, 'iPES', iPES, '.Smat'
                inquire(file=trim(dirName)//trim(fileName), exist=f_exists)
                if (.not. f_exists) cycle
                write(*,'(a,a,a,a)') '  Reading dir: ', trim(dirName), '  file: ', trim(fileName)
                open(newunit=inpFileUnit, file=trim(dirName)//trim(fileName), status='old')
                do iline = 1, 2
                    read(inpFileUnit,'(a)') chartmp
                end do
                read(inpFileUnit,'(a)') chartmp
                read(chartmp(index(chartmp,':')+1:),*) tms
                read(inpFileUnit,'(a)') chartmp
                do iline = 1, 5
                    read(inpFileUnit,'(a)') chartmp
                end do
                do
                    read(inpFileUnit,'(a)', iostat=io_stat) chartmp
                    if (io_stat /= 0) exit
                    if (len_trim(chartmp) == 0) cycle
                    if (chartmp(1:1) == '#') cycle
                    backspace(inpFileUnit)
                    read(inpFileUnit,*, iostat=io_stat) iEtot, Ecol(iEtot), vRead, jRead, Lr, tmp1, Sreal, Simag, tmp2
                    if (io_stat /= 0) cycle
                end do
                close(inpFileUnit)
                return
            end do
        end do
    end do
end subroutine readFirstSmat
