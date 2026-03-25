program CrossSection
    use m_MachinaBasic , only : f8, c8
    use iso_fortran_env, only : iostat_end
    implicit none
    
    !> Constants
    real(f8), parameter :: au2ev = 27.211386245988_f8
    real(f8), parameter :: ev2au = 1.0_f8 / au2ev
    real(f8), parameter :: pi = 4.0_f8*atan(1.0_f8)
    complex(c8), parameter :: img = (0.0_f8, 1.0_f8)
    complex(c8), parameter :: imgZero = (0.0_f8, 0.0_f8)

    integer :: nEtot, nJtot 
    integer :: v0, j0, vf, jf
    integer :: ichnl
    integer :: nLr, nLp, l0, iPES
    integer :: tpar
    integer :: vRead, jRead
    integer :: inpFileUnit, outFileUnit
    integer :: iJtot, iEtot, imj, imjp
    integer :: iLr, iLp, ith, ipar, iml
    integer :: io_stat, iline

    real(f8) :: tms, theta, Lr
    real(f8) :: tmp1, tmp2
    real(f8) :: Sreal, Simag
    real(f8), allocatable :: Ecol(:)
    real(f8), allocatable :: kvj(:)
    real(f8), allocatable :: ICS(:), prop(:,:)
    real(f8), allocatable :: DCS_SF(:,:)
    real(f8), allocatable :: CG1(:,:), CG2(:,:,:,:)
    real(f8), allocatable :: Smat(:,:,:,:,:)
    real(f8), external :: CG

    complex(c8) :: Tmat
    complex(c8), allocatable :: Ylm(:,:,:)
    complex(c8), allocatable :: amplitude_SF(:,:,:,:)
    complex(c8), external :: spharmonic

    character(len=2) :: reactChannel(2) = ['AB','AC']
    character(len=100) :: fileName, dirName, SmatDir, outFile
    character(len=200) :: chartmp, tmpStr

    logical :: f_exists

    namelist /stateInfo/ ichnl, v0, j0, l0, nJtot, iPES, nEtot, vf, jf, SmatDir

    open(newunit=inpFileUnit, file='CrossSection.inf', status='old')
    read(inpFileUnit, nml=stateInfo)
    close(inpFileUnit)

    nLp = nJtot + jf

    allocate(kvj(nEtot))
    allocate(Ecol(nEtot))
    allocate(ICS(nEtot))
    allocate(prop(nEtot,0:nJtot))
    allocate(DCS_SF(nEtot,0:180))
    allocate(Smat(nEtot,0:nJtot,0:nLp,2,2))
    allocate(Ylm(0:180,0:nLp,-nLp:nLp))
    allocate(CG1(0:nJtot,-j0:j0))
    allocate(CG2(0:nJtot,-jf:jf,0:nLp,-nLp:nLp))
    allocate(amplitude_SF(nEtot,0:180,-j0:j0,-jf:jf))

    do tpar = -1, 1, 2
        if (tpar == -1) ipar = 1
        if (tpar == 1) ipar = 2
        write(dirName, '(a,a,i0,a)') trim(SmatDir), '/p', tpar, '/'
        write(*,*) 'Reading S-matrix from directory: ', trim(dirName)
        do iJtot = 0, nJtot
            write(fileName,'(a,a,i0,a,i0,a,i0,a)') reactChannel(ichnl), '_J', iJtot, 'L', l0, 'iPES', iPES, '.Smat'
            inquire(file=trim(dirName)//trim(fileName), exist=f_exists)
            if (.not. f_exists) cycle 
            open(newunit=inpFileUnit, file=trim(dirName)//trim(fileName), status='old')
            write(*,*) 'Reading file: ', trim(dirName)//trim(fileName)
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
                if (io_stat == iostat_end) exit
                backspace(inpFileUnit)
                read(inpFileUnit,*) iEtot, Ecol(iEtot), vRead, jRead, Lr, tmp1, Sreal, Simag, tmp2
                if (vRead == vf .and. jRead == jf) then
                    kvj(iEtot) = sqrt(2.0_f8*Ecol(iEtot)*ev2au*tms)
                    Smat(iEtot,iJtot,nint(Lr),ipar,1) = Sreal
                    Smat(iEtot,iJtot,nint(Lr),ipar,2) = Simag
                end if
            end do
            close(inpFileUnit)
        end do
    end do
    
    Ylm = imgZero
    CG1 = 0.0_f8
    CG2 = 0.0_f8

    do iJtot = 0, nJtot
        do imj = -j0, j0 
            do imjp = -jf, jf 
                do iLp = 0, nLp 
                    iml = imj - imjp
                    if (abs(iml) > iLp) cycle
                    if ((iJtot > j0+l0 .or. iJtot < abs(j0-l0)) .or. (iJtot > jf+iLp .or. iJtot < abs(jf-iLp))) cycle
                    CG1(iJtot,imj) = CG(j0, imj, l0, 0, iJtot)
                    CG2(iJtot,imjp,iLp,iml) = CG(jf, imjp, iLp, iml, iJtot)
                end do
            end do
        end do
    end do

    do ith = 0, 180
        do iLp = 0, nLp
            do iml = -nLp, nLp
                if (abs(iml) > iLp) cycle
                theta = ith * pi / 180.0_f8
                Ylm(ith,iLp,iml) = spharmonic(theta,0.0_f8, iLp, iml)
            end do
        end do
    end do
    
    write(*,*) 'Calculating scattering amplitudes...'
    do iEtot = 1, nEtot 
        do ith = 0, 180
            do imj = -j0, j0 
                do imjp = -jf, jf 
                    theta = ith * pi / 180.0_f8
                    iml = imj - imjp
                    amplitude_SF(iEtot,ith,imj,imjp) = imgZero
                    do iJtot = 0, nJtot
                        do ipar = 1, 2 
                            do iLp = 0, nLp 
                                if (abs(iml) > iLp) cycle 
                                Tmat%re = -Smat(iEtot,iJtot,iLp,ipar,1)
                                Tmat%im = -Smat(iEtot,iJtot,iLp,ipar,2)
                                if (Tmat == imgZero) cycle
                                if (CG1(iJtot,imj) == 0.0_f8 .or. CG2(iJtot,imjp,iLp,iml) == 0.0_f8) cycle
                                amplitude_SF(iEtot,ith,imj,imjp) = amplitude_SF(iEtot,ith,imj,imjp) &
                                                                    + dsqrt(2.0_f8*l0+1.0_f8) * img**(l0-iLp+1) &
                                                                    * Ylm(ith,iLp,iml) * Tmat * dsqrt(pi) &
                                                                    * CG1(iJtot,imj) * CG2(iJtot,imjp,iLp,iml)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do

    write(*,*) 'Calculating DCS in SF...'
    do iEtot = 1, nEtot
        do ith = 0, 180
            theta = ith * pi / 180.0_f8
            DCS_SF(iEtot,ith) = 0.0_f8
            do imj = -j0, j0
                do imjp = -jf, jf
                    DCS_SF(iEtot,ith) = DCS_SF(iEtot,ith) &
                                        + amplitude_SF(iEtot,ith,imj,imjp)%re**2 &
                                        + amplitude_SF(iEtot,ith,imj,imjp)%im**2
                end do
            end do
            DCS_SF(iEtot,ith) = 2.0_f8*pi*DCS_SF(iEtot,ith) / (kvj(iEtot)**2*(2.0_f8*j0+1.0_f8))
        end do
    end do

    write(outFile,'(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'DCS_',reactChannel(ichnl), 'v', v0, 'j', j0, 'vf', vf, 'jf', jf, 'iPES', iPES, '.dat'
    open(newunit=outFileUnit, file=trim(outFile), status='replace')
    write(outFileUnit,'(a)') '#theta(deg)    DCS_SF(bohr^2/sr)'
    do ith = 0, 180
        write(outFileUnit,'(f15.9,1x,*(f12.9))') ith*1.0_f8, (DCS_SF(iEtot,ith), iEtot=1,nEtot)
    end do
    close(outFileUnit)
    write(*,*) 'DCS written to file: ', trim(outFile)

    do iEtot = 1, nEtot
        do iJtot = 0, nJtot
            prop(iEtot,iJtot) = 0.0_f8
            do ipar = 1, 2
                do iLp = 0, nLp
                    if (Smat(iEtot,iJtot,iLp,ipar,1) == 0.0_f8 .and. Smat(iEtot,iJtot,iLp,ipar,2) == 0.0_f8) cycle
                    prop(iEtot,iJtot) = prop(iEtot,iJtot) + Smat(iEtot,iJtot,iLp,ipar,1)**2 + Smat(iEtot,iJtot,iLp,ipar,2)**2
                end do
            end do
        end do 
    end do

    write(*,*) 'Writing ICS...'
    ICS = 0.0_f8
    do iEtot = 1, nEtot
        ICS(iEtot) = 0.0_f8
        do iJtot = 0, nJtot
            ICS(iEtot) = ICS(iEtot) + prop(iEtot,iJtot) * (2.0_f8*iJtot+1.0_f8)
        end do
        ICS(iEtot) = pi*ICS(iEtot) / (kvj(iEtot)**2*(2.0_f8*j0+1.0_f8))
    end do

    write(outFile,'(a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a)') 'ICS_',reactChannel(ichnl), 'v', v0, 'j', j0, 'vf', vf, 'jf', jf, 'iPES', iPES, '.dat'
    open(newunit=outFileUnit, file=trim(outFile), status='replace')
    write(outFileUnit,'(a)') '#Ecol(eV)    prop(Jmax)    ICS(bohr^2)'
    do iEtot = 1, nEtot
        write(outFileUnit,'(f8.3,1x,f12.5,1x,f12.5)') Ecol(iEtot), prop(iEtot,nJtot), ICS(iEtot)
    end do
    close(outFileUnit)
    write(*,*) 'ICS written to file: ', trim(outFile)

    
end program CrossSection