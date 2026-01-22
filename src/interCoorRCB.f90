module m_InterCoorRCB
    use m_MachinaBasic, only : f8
    use m_gPara
    use m_Basis
    use m_Potent
    implicit none

!> Mass constants in RCB, beta = C+AB and gamma = B+AC
    real(f8), private :: A_beta, B_beta
    real(f8), private :: A_gamma, B_gamma
!> DVR grids
    real(f8), allocatable, private :: rp1(:), rp2(:)
!> DVR to FBR transformation matrices
    real(f8), allocatable, private :: DVR2FBR1(:,:), DVR2FBR2(:,:)
!> Product diatomic potential 
    real(f8), allocatable, private :: VAB(:,:), VAC(:,:)
!> Product diatomic vib-rot energy and wavefunction in product coordinates
    real(f8), allocatable, private :: Evj_AB(:,:,:), Evj_AC(:,:,:)
    real(f8), allocatable, private :: WFvj_AB(:,:,:,:), WFvj_AC(:,:,:,:)
!> Boundary of RCB coordinates
    real(f8), private :: RCB_ABZ_range(2), RCB_ABr_range(2) 
    real(f8), private :: RCB_ACZ_range(2), RCB_ACr_range(2) 

    public
    private :: findjMax, DVR_VibRotp
    private :: setUij, interCoorWFvj, setAtDMat
    private :: setNumberInterCoor, setTransKMat
    private :: RCB_rmid, RCB_Zmid
    private :: rJacobiFactor, ZJacobiFactor

contains


!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine initRCB()
        implicit none
        real(f8) :: Vadia(nPES)
        integer :: ir, vmax, jmax, Kmax, ijK, j, K
        character(len=10) :: type 

!> Channel 1
        type = 'A+BC->C+AB'
!> Set number of intermediate coordinates
        call setNumberInterCoor(type, channel1%ichoice, channel1%Zinf, channel1%r_range, atomMass, channel1%nIC, interCoor_AB, idXY_AB)
!> Set Uij
        call setUij(channel1%ichoice, channel1%nIC, interCoor_AB, idXY_AB, Uij_AB)
        channel1%massAB = atomMass(1)*atomMass(2)/(atomMass(1)+atomMass(2))
        channel1%massTot = atomMass(3)*(atomMass(2)+atomMass(1))/(atomMass(1)+atomMass(2)+atomMass(3))
!> Set DVR grids
        allocate(rp1(channel1%nr))
        call sinDVRGrid(channel1%nr, channel1%r_range(1), channel1%r_range(2), rp1)
!> Set transformation matrix from DVR to FBR
        allocate(DVR2FBR1(channel1%nr,channel1%nr))
        call setDVR2FBR(channel1%nr, DVR2FBR1)
!> Set diatomic potential
        allocate(VAB(channel1%nr,nPES))
        do ir = 1, channel1%nr
            call setPot_Product(type, rp1(ir), Vadia)
            VAB(ir,:) = Vadia
        end do
!> Find vmax for channel 1
        call findvMax(channel1%massAB, channel1%nr, rp1, VAB, vmax)
        if (channel1%vmax < vmax) then
            write(outFileUnit,'(1x,a,i4,i4)') 'Warning: channel 1 vmax increased from ', &
                                            channel1%vmax, ' to ', vmax
            channel1%vmax = vmax
        end if
!> Find jmax for channel 1
        call findjMax(channel1%massAB, channel1%nr, rp1, VAB, channel1%jmin, channel1%jinc, jmax)
        if (channel1%jmax < jmax) then
            write(outFileUnit,'(1x,a,i4,i4)') 'Warning: channel 1 jmax increased from ', &
                                            channel1%jmax, ' to ', jmax
            channel1%jmax = jmax
        end if
!> AB vib-rot energy and wavefunction
        allocate(Evj_AB(0:channel1%vmax,channel1%jmin:channel1%jmax,nPES))
        allocate(WFvj_AB(channel1%nr,0:channel1%vmax,channel1%jmin:channel1%jmax,nPES))
        call DVR_VibRotp(type)
!> Set jK pair
        channel1%Kmin = merge(0, 1, initWP%tpar == 1)
        call setjKPair(channel1%jmin, channel1%jmax, channel1%jinc, channel1%Kmin, channel1%njK, AB_jKPair)
        Kmax = min(channel1%jmax, initWP%Jtot)
        allocate(AB_seqjK(channel1%jmin:channel1%jmax,channel1%Kmin:Kmax))
        AB_seqjK = -1
        do ijK = 1, channel1%njK 
            j = AB_jKPair(ijK,1)
            K = AB_jKPair(ijK,2)
            AB_seqjK(j,K) = ijK
        end do
!> WF_AB in intermediate coordinates
        call interCoorWFvj(type, channel1%ichoice)
!> Transformation matrix from Ka to Kb
        call setTransKMat(type)
!> AtD matrix for channel 1
        call setAtDMat(channel1%nIC,channel1%ichoice,atomMass(1),atomMass(2),InterCoor_AB,AtDMat_AB)

        if (nChannel == 2) then
!> Channel 2
            type = 'A+BC->B+AC'
!> Set number of intermediate coordinates
            call setNumberInterCoor(type,channel2%ichoice,channel2%Zinf,channel2%r_range,atomMass,channel2%nIC, interCoor_AC, idXY_AC)
!> Set Uij
            call setUij(channel2%ichoice, channel2%nIC, interCoor_AC, idXY_AC, Uij_AC)
            channel2%massAC = atomMass(1)*atomMass(3)/(atomMass(1)+atomMass(3))
            channel2%massTot = atomMass(2)*(atomMass(3)+atomMass(1))/(atomMass(1)+atomMass(2)+atomMass(3))
!> Set DVR grids
            allocate(rp2(channel2%nr))
            call sinDVRGrid(channel2%nr, channel2%r_range(1), channel2%r_range(2), rp2)
!> Set transformation matrix from DVR to FBR
            allocate(DVR2FBR2(channel2%nr,channel2%nr))
            call setDVR2FBR(channel2%nr, DVR2FBR2)
!> Set diatomic potential
            allocate(VAC(channel2%nr,nPES))
            do ir = 1, channel2%nr
                call setPot_Product(type, rp2(ir), Vadia)
                VAC(ir,:) = Vadia
            end do
!> Find vmax for channel 2
            call findvMax(channel2%massAC, channel2%nr, rp2, VAC, vmax)
            if (channel2%vmax < vmax) then
                write(outFileUnit,'(1x,a,i4,i4)') 'Warning: channel 2 vmax increased from ', &
                                                channel2%vmax, ' to ', vmax
                channel2%vmax = vmax
            end if
!> Find jmax for channel 2
            call findjMax(channel2%massAC, channel2%nr, rp2, VAC, channel2%jmin, channel2%jinc, jmax)
            if (channel2%jmax < jmax) then
                write(outFileUnit,'(1x,a,i4,i4)') 'Warning: channel 2 jmax increased from ', &
                                                channel2%jmax, ' to ', jmax
                channel2%jmax = jmax
            end if 
!> AC vib-rot energy and wavefunction
            allocate(Evj_AC(0:channel2%vmax,channel2%jmin:channel2%jmax,nPES))
            allocate(WFvj_AC(channel2%nr,0:channel2%vmax,channel2%jmin:channel2%jmax,nPES))
            call DVR_VibRotp(type)
!> Set jK pair
            channel2%Kmin = merge(0, 1, initWP%tpar == 1)
            call setjKPair(channel2%jmin, channel2%jmax, channel2%jinc, channel2%Kmin, channel2%njK, AC_jKPair)
            Kmax = min(channel2%jmax, initWP%Jtot)
            allocate(AC_seqjK(channel2%jmin:channel2%jmax,channel2%Kmin:Kmax))
            AC_seqjK = -1
            do ijK = 1, channel2%njK 
                j = AC_jKPair(ijK,1)
                K = AC_jKPair(ijK,2)
                AC_seqjK(j,K) = ijK
            end do
!> WF_AC in intermediate coordinates
            call interCoorWFvj(type, channel2%ichoice)
!> Transformation matrix from Ka to Kb
            call setTransKMat(type)
!> AtD matrix for channel 2
            call setAtDMat(channel2%nIC,channel2%ichoice,atomMass(1),atomMass(3),InterCoor_AC,AtDMat_AC)
        end if

        write(outFileUnit,'(1x,a)') '=====> Intermediate coordinate RCB information <====='
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a)') 'Intermediate coordinate ichoice = 1, R_r -> Z_p, X, Y = Z_r, theta_r.'
        write(outFileUnit,'(1x,a)') 'Intermediate coordinate ichoice = 2, Z_r -> Z_p, X, Y = r_r, theta_r.'
        write(outFileUnit,'(1x,a,i2)') 'Number of product channels: ', nChannel
        if (nChannel == 0) then 
            write(outFileUnit,'(1x,a)') 'No reactive channels.'
        else 
            write(outFileUnit,'(1x,a)') 'Product channel 1 : A + BC -> C + AB :'
            write(outFileUnit,'(1x,a,f15.9)') "Reduced Masses of AB (a.u.):", channel1%massAB
            write(outFileUnit,'(1x,a,f15.9)') "Total Masses of C + AB (a.u.):", channel1%massTot
            write(outFileUnit,'(1x,a,i3,a,i3,a,i2)') 'vmax = ', channel1%vmax, ' jmax = ', channel1%jmax, ' jpar = ', channel1%jpar
            write(outFileUnit,'(1x,a,f15.9,a)') 'Evj of AB (v0,jmin) = ', Evj_AB(0,channel1%jmin,initWP%PES0)*au2ev, ' eV.'
            write(outFileUnit,'(1x,a,f15.9,a)') 'Evj of AB (v0,jmin) = ', Evj_AB(0,channel1%jmin,initWP%PES0)*au2cm, ' cm-1.'
            write(outFileUnit,'(1x,a)') 'Please check the energy!'
            write(outFileUnit,'(1x,a,i4)') 'Number of jK pairs in channel 1: ', channel1%njK
            write(outFileUnit,'(1x,a,i2)') 'ichoice = ', channel1%ichoice
            write(outFileUnit,'(1x,a,i4)') 'Number of intermediate coordinates in channel 1: ', channel1%nIC
            write(outFileUnit,'(1x,a)') 'Please check InterCoorWF_AB.chk for the wavefunction normalization!'
            write(outFileUnit,'(1x,a,f15.9,a,f15.9)') 'Zmin: ', RCB_ABZ_range(1), ' Zmax: ', RCB_ABZ_range(2)
            write(outFileUnit,'(1x,a,f15.9,a,f15.9)') 'Rmin: ', RCB_ABr_range(1), ' Rmax: ', RCB_ABr_range(2)
            write(outFileUnit,'(1x,a)') ''
            if (nChannel == 2) then
                write(outFileUnit,'(1x,a)') 'Product channel 2 : A + BC -> B + AC :'
                write(outFileUnit,'(1x,a,f15.9)') "Reduced Masses of AC (a.u.):", channel2%massAC
                write(outFileUnit,'(1x,a,f15.9)') "Total Masses of B + AC (a.u.):", channel2%massTot
                write(outFileUnit,'(1x,a,i3,a,i3,a,i2)') 'vmax = ', channel2%vmax, ' jmax = ', channel2%jmax, ' jpar = ', channel2%jpar
                write(outFileUnit,'(1x,a,f15.9,a)') 'Evj of AC (v0,jmin) = ', Evj_AC(0,channel2%jmin,initWP%PES0)*au2ev, ' eV.'
                write(outFileUnit,'(1x,a,f15.9,a)') 'Evj of AC (v0,jmin) = ', Evj_AC(0,channel2%jmin,initWP%PES0)*au2cm, ' cm-1.'
                write(outFileUnit,'(1x,a)') 'Please check the energy!'
                write(outFileUnit,'(1x,a,i4)') 'Number of jK pairs in channel 2: ', channel2%njK
                write(outFileUnit,'(1x,a,i2)') 'ichoice = ', channel2%ichoice
                write(outFileUnit,'(1x,a,i4)') 'Number of intermediate coordinates in channel 2: ', channel2%nIC
                write(outFileUnit,'(1x,a)') 'Please check InterCoorWF_AC.chk for the wavefunction normalization!'
                write(outFileUnit,'(1x,a,f15.9,a,f15.9)') 'Zmin: ', RCB_ACZ_range(1), ' Zmax: ', RCB_ACZ_range(2)
                write(outFileUnit,'(1x,a,f15.9,a,f15.9)') 'Rmin: ', RCB_ACr_range(1), ' Rmax: ', RCB_ACr_range(2)
                write(outFileUnit,'(1x,a)') ''
            end if
        end if
    end subroutine initRCB
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine findvMax(mass, nr, rp, Vpot, vmax)
        implicit none
        real(f8), intent(in) :: mass
        integer, intent(in) :: nr
        real(f8), intent(in) :: rp(nr)
        real(f8), intent(in) :: Vpot(nr,nPES)
        integer, intent(out) :: vmax
        real(f8) :: DVREig(nr), DVRWF(nr,nr), Vtmp(nr)
        real(f8) :: range, EtotMax
        integer :: v, vmaxtmp, vlimit, ir, iPES

        vlimit = 41
        vmaxtmp = -1
        range = rp(nr) - rp(1)
        EtotMax = maxval(Ecol)+int_POEig(initWP%v0+1)

        do iPES = 1, nPES
            Vtmp(:) = Vpot(:,iPES)
            do v = 0, vlimit
                call DVR_vib(nr,mass,range,Vtmp,DVREig,DVRWF)
                if (EtotMax <= DVREig(v+1)) then 
                    vmaxtmp = v - 1
                    exit
                end if
            end do
        end do
        if (vmaxtmp == -1) then
            write(outFileUnit,*) 'Error in findvMax: vmax not found!'
            write(outFileUnit,*) 'POSITION: interCoorRCB.f90, subroutine findvMax'
            stop
        end if
        vmax = vmaxtmp
    end subroutine findvMax
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine findjMax(mass, nr, rp, Vpot, jmin, jinc, jmax)
        implicit none
        real(f8), intent(in) :: mass
        integer, intent(in) :: nr
        real(f8), intent(in) :: rp(nr)
        real(f8), intent(in) :: Vpot(nr,nPES)
        integer, intent(in) :: jmin, jinc
        integer, intent(out) :: jmax
        real(f8) :: Vtmp(nr)
        real(f8) :: DVREig(nr), DVRWF(nr,nr)
        real(f8) :: range, EtotMax
        integer :: j, jlimit, jmaxtmp, ir, iPES

        jlimit = 233
        jmaxtmp = -1
        range = rp(nr) - rp(1)
        EtotMax = maxval(Ecol)+int_POEig(initWP%v0+1)
        do iPES = 1, nPES 
            do j = jmin, jlimit, jinc
                do ir = 1, nr 
                    Vtmp(ir) = Vpot(ir,iPES) + real(j*(j+1),f8)/(2.0_f8*mass*rp(ir)**2)
                end do 
                call DVR_vib(nr,mass,range,Vtmp,DVREig,DVRWF)
                if (EtotMax <= DVREig(1)) then 
                    jmaxtmp = j - 1
                    exit
                end if
            end do
        end do
        if (jmaxtmp == -1) then
            write(outFileUnit,*) 'Error in findjMax: jmax not found!'
            write(outFileUnit,*) 'POSITION: interCoorRCB.f90, subroutine findjMax'
            stop
        end if
        jmax = jmaxtmp
    end subroutine findjMax
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine DVR_VibRotp(type)
        implicit none
        character(len=*), intent(in) :: type
        !> type: 'A+BC->C+AB' or 'A+BC->B+AC'
        integer :: ir, v, j, i, iPES
        real(f8) :: mass, range
        real(f8), allocatable :: DVREig(:), DVRWF(:,:)
        real(f8), allocatable :: Vpot(:)

        if (type == 'A+BC->C+AB') then 
            mass = channel1%massAB
            range = channel1%r_range(2) - channel1%r_range(1)
            allocate(DVREig(channel1%nr))
            allocate(DVRWF(channel1%nr,channel1%nr))
            allocate(Vpot(channel1%nr))
            WFvj_AB = 0.0_f8
            do iPES = 1, nPES
                do j = channel1%jmin, channel1%jmax, channel1%jinc
                    do ir = 1, channel1%nr 
                        Vpot(ir) = VAB(ir,iPES) + real(j*(j+1),f8)/(2.0_f8*mass*rp1(ir)**2)
                    end do 
                    call DVR_vib(channel1%nr,mass,range,Vpot, DVREig,DVRWF)
                    do v = 0, channel1%vmax
                        Evj_AB(v,j,iPES) = DVREig(v+1)
                        do ir = 1, channel1%nr
                            do i = 1, channel1%nr
                                WFvj_AB(ir,v,j,iPES) = WFvj_AB(ir,v,j,iPES) + DVR2FBR1(ir,i)*DVRWF(i,v+1)
                            end do
                        end do
                    end do
                end do
            end do
            deallocate(DVREig)
            deallocate(DVRWF)
            deallocate(Vpot)
        else if (type == 'A+BC->B+AC') then
            mass = channel2%massAC
            range = channel2%r_range(2) - channel2%r_range(1)
            allocate(DVREig(channel2%nr))
            allocate(DVRWF(channel2%nr,channel2%nr))
            allocate(Vpot(channel2%nr))
            WFvj_AC = 0.0_f8
            do iPES = 1, nPES
                do j = channel2%jmin, channel2%jmax, channel2%jinc
                    do ir = 1, channel2%nr 
                        Vpot(ir) = VAC(ir,iPES) + real(j*(j+1),f8)/(2.0_f8*mass*rp2(ir)**2)
                    end do 
                    call DVR_vib(channel2%nr,mass,range,Vpot,DVREig,DVRWF)
                    do v = 0, channel2%vmax
                        Evj_AC(v,j,iPES) = DVREig(v+1)
                        do ir = 1, channel2%nr
                            do i = 1, channel2%nr
                            WFvj_AC(ir,v,j,iPES) = WFvj_AC(ir,v,j,iPES) + DVR2FBR2(ir,i)*DVRWF(i,v+1)
                            end do
                        end do
                    end do
                end do
            end do
            deallocate(DVREig)
            deallocate(DVRWF)
            deallocate(Vpot)
        end if

    end subroutine DVR_VibRotp
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine interCoorWFvj(type, ichoice)
        implicit none
        character(len=*), intent(in) :: type
        integer, intent(in) :: ichoice
        real(f8) :: gij, djFact, tmp, rp, norm, dpos, tmpW
        real(f8), allocatable :: pos(:)
        real(f8), external :: spgndr
        integer :: iPos, nPos, ith, v, j, iPES, iDVR, ijK, K, iIC
        integer :: normUnit, Kmax
        character(len=20) :: fileName 

        djFact = merge(0.5_f8, 1.0_f8, initWP%jinc==2)
        Kmax = min(IALR%jint, initWP%Jtot)

        if (ichoice == 1) then 
            nPos = IALR%nZ_I
            allocate(pos(nPos))
            pos(:) = Z_IALR(1:nPos)
            dpos = merge(Z_IALR(2)-Z_IALR(1), 1.0_f8, nPos > 1)
        else if (ichoice == 2) then
            nPos = IALR%vint
            allocate(pos(nPos))
            pos(:) = r_Int(1:nPos)
            dpos = merge(r_Int(2)-r_Int(1), 1.0_f8, nPos > 1)
        else
            write(outFileUnit,*) 'Error in interCoorWFvj: unknown ichoice ', ichoice
            write(outFileUnit,*) 'POSITION: interCoorRCB.f90, subroutine interCoorWFvj'
            stop
        end if

        if (type == 'A+BC->C+AB') then 
            RCB_ABZ_range(1) = 100.0_f8; RCB_ABZ_range(2) = 0.0_f8
            RCB_ABr_range(1) = 100.0_f8; RCB_ABr_range(2) = 0.0_f8
            allocate(ICWFvj_AB(channel1%nIC,0:channel1%vmax,channel1%njK,nPES))
            allocate(intAB_TIDWF(channel1%nIC,initWP%Kmin:Kmax,nPES,nEtot))
            intAB_TIDWF = imgZero
            ICWFvj_AB = 0.0_f8
            do iPES = 1, nPES 
                do iIC = 1, channel1%nIC
                    iPos = idXY_AB(iIC,1)
                    ith = idXY_AB(iIC,2)
                    rp = interCoor_AB(iIC,4)
                    gij = interCoor_AB(iIC,7)/channel1%Zinf/rp*pos(iPos)*interCoor_AB(iIC,2)
                    tmpW = gij*dsqrt(int_AWeight(ith)*djFact)*dsqrt(dpos)*dsqrt(2.0_f8/(channel1%r_range(2)-channel1%r_range(1)))
                    do ijK = 1, channel1%njK
                        j = AB_jKPair(ijK,1)
                        K = AB_jKPair(ijK,2)
                        do v = 0, channel1%vmax
                            tmp = 0.0_f8
                            do iDVR = 1, channel1%nr 
                                tmp = tmp + WFvj_AB(iDVR,v,j,iPES) &
                                    * dsin(iDVR*pi*(rp-channel1%r_range(1))/(channel1%r_range(2)-channel1%r_range(1)))
                            end do
                            if (abs(tmp) > 0.01_f8) then 
                                if (ichoice == 1) then 
                                    RCB_ABZ_range(1) = min(RCB_ABZ_range(1), pos(iPos))
                                    RCB_ABZ_range(2) = max(RCB_ABZ_range(2), pos(iPos))
                                    RCB_ABr_range(1) = min(RCB_ABr_range(1), interCoor_AB(iIC,2))
                                    RCB_ABr_range(2) = max(RCB_ABr_range(2), interCoor_AB(iIC,2))
                                else if (ichoice == 2) then
                                    RCB_ABZ_range(1) = min(RCB_ABZ_range(1), interCoor_AB(iIC,2))
                                    RCB_ABZ_range(2) = max(RCB_ABZ_range(2), interCoor_AB(iIC,2))
                                    RCB_ABr_range(1) = min(RCB_ABr_range(1), pos(iPos))
                                    RCB_ABr_range(2) = max(RCB_ABr_range(2), pos(iPos))
                                end if
                            end if
                            ICWFvj_AB(iIC,v,ijK,iPES) = tmpW*tmp*spgndr(j,K,interCoor_AB(iIC,5))
                        end do
                    end do
                end do 
            end do 
            !> Normal check
            fileName = 'InterCoorWF_AB.chk'
            open(unit=normUnit, file=fileName, status='replace', action='write')
            write(normUnit,'(1x,a)') '# iPES  v  j  K  Evj(eV)  Norm'
            do iPES = 1, nPES 
                do ijK = 1, channel1%njK
                    j = AB_jKPair(ijK,1)
                    K = AB_jKPair(ijK,2)
                    do v = 0, channel1%vmax
                        norm = 0.0_f8
                        do iIC = 1, channel1%nIC
                            norm = norm + ICWFvj_AB(iIC,v,ijK,iPES)**2
                        end do
                        write(normUnit,'(1x,i4,1x,i4,1x,i4,1x,i4,1x,f15.9,1x,f15.9)') iPES, v, j, K, Evj_AB(v,j,iPES)*au2ev,norm
                    end do
                end do
            end do
            close(normUnit)
        else if (type == 'A+BC->B+AC') then
            RCB_ACZ_range(1) = 100.0_f8; RCB_ACZ_range(2) = 0.0_f8
            RCB_ACr_range(1) = 100.0_f8; RCB_ACr_range(2) = 0.0_f8
            allocate(ICWFvj_AC(channel2%nIC,0:channel2%vmax,channel2%njK,nPES))
            allocate(intAC_TIDWF(channel2%nIC,initWP%Kmin:Kmax,nPES,nEtot))
            intAC_TIDWF = imgZero
            ICWFvj_AC = 0.0_f8
            do iPES = 1, nPES 
                do iIC = 1, channel2%nIC
                    iPos = idXY_AC(iIC,1)
                    ith = idXY_AC(iIC,2)
                    rp = interCoor_AC(iIC,4)
                    gij = interCoor_AC(iIC,7)/channel2%Zinf/rp*pos(iPos)*interCoor_AC(iIC,2)
                    tmpW = gij*dsqrt(int_AWeight(ith)*djFact)*dsqrt(dpos)*dsqrt(2.0_f8/(channel2%r_range(2)-channel2%r_range(1)))
                    do ijK = 1, channel2%njK
                        j = AC_jKPair(ijK,1)
                        K = AC_jKPair(ijK,2)
                        do v = 0, channel2%vmax
                            tmp = 0.0_f8
                            do iDVR = 1, channel2%nr 
                                tmp = tmp + WFvj_AC(iDVR,v,j,iPES) &
                                    * dsin(iDVR*pi*(rp-channel2%r_range(1))/(channel2%r_range(2)-channel2%r_range(1)))
                            end do
                            if (abs(tmp) > 0.01_f8) then 
                                if (ichoice == 1) then 
                                    RCB_ACZ_range(1) = min(RCB_ACZ_range(1), pos(iPos))
                                    RCB_ACZ_range(2) = max(RCB_ACZ_range(2), pos(iPos))
                                    RCB_ACr_range(1) = min(RCB_ACr_range(1), interCoor_AC(iIC,2))
                                    RCB_ACr_range(2) = max(RCB_ACr_range(2), interCoor_AC(iIC,2))
                                else if (ichoice == 2) then
                                    RCB_ACZ_range(1) = min(RCB_ACZ_range(1), interCoor_AC(iIC,2))
                                    RCB_ACZ_range(2) = max(RCB_ACZ_range(2), interCoor_AC(iIC,2))
                                    RCB_ACr_range(1) = min(RCB_ACr_range(1), pos(iPos))
                                    RCB_ACr_range(2) = max(RCB_ACr_range(2), pos(iPos))
                                end if
                            end if
                            ICWFvj_AC(iIC,v,ijK,iPES) = tmpW*tmp*spgndr(j,K,interCoor_AC(iIC,5))
                        end do
                    end do
                end do
            end do
            !> Normal check
            fileName = 'InterCoorWF_AC.chk'
            open(unit=normUnit, file=fileName, status='replace', action='write')
            write(normUnit,'(1x,a)') '# iPES  v  j  K  Evj(eV)  Norm'
            do iPES = 1, nPES 
                do ijK = 1, channel2%njK    
                    j = AC_jKPair(ijK,1)
                    K = AC_jKPair(ijK,2)
                    do v = 0, channel2%vmax
                        norm = 0.0_f8
                        do iIC = 1, channel2%nIC
                            norm = norm + ICWFvj_AC(iIC,v,ijK,iPES)**2
                        end do
                        write(normUnit,'(1x,i4,1x,i4,1x,i4,1x,i4,1x,f15.9,1x,f15.9)') iPES, v, j, K, Evj_AC(v,j,iPES)*au2ev,norm
                    end do
                end do
            end do
            close(normUnit)
        end if

    end subroutine interCoorWFvj
!> ------------------------------------------------------------------------------------------------------------------ <!
    
!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setTransKMat(type)
        implicit none
        character(len=*), intent(in) :: type
        real(f8), external :: dm 
        real(f8) :: beta, fact
        real(f8), external :: CG
        integer :: iIC, Kr, Kp, Krmax, Kpmax
        integer :: ijK, l, j

        Krmax = min(IALR%jint, initWP%Jtot)

        if (type == 'A+BC->C+AB') then 
            Kpmax = min(channel1%jmax, initWP%Jtot)

            !> Transformation matrix Ka -> Kb 
            allocate(TKMat_AB(channel1%nIC,initWP%Kmin:Krmax,channel1%Kmin:Kpmax))
            do Kp = channel1%Kmin, Kpmax
                do Kr = initWP%Kmin, Krmax
                    do iIC = 1, channel1%nIC
                        beta = interCoor_AB(iIC,6)
                        TKMat_AB(iIC,Kr,Kp) = dm(initWP%Jtot, Kp, Kr, beta) &
                                            + dm(initWP%Jtot, Kp, -Kr, beta) &
                                            * initWP%tpar*(-1.0_f8)**(Kr)
                    end do
                    if (Kr == 0) TKMat_AB(:,Kr,Kp) = TKMat_AB(:,Kr,Kp)/dsqrt(2.0_f8)
                    if (Kp == 0) TKMat_AB(:,Kr,Kp) = dsqrt(2.0_f8)*TKMat_AB(:,Kr,Kp)
                end do
            end do

            !> Qn_l in AB
            allocate(qn_l_AB(channel1%njK))
            do ijK = 1, channel1%njK
                j = AB_jKPair(ijK,1)
                do l = abs(initWP%Jtot - j), initWP%Jtot + j
                    if ((-1)**(l+j+initWP%Jtot) == initWP%tpar) then
                        qn_l_AB(ijK) = l
                    end if
                end do
            end do

            !> BLK in AB
            allocate(BLK_AB(channel1%Kmin:Kpmax,channel1%njK))
            BLK_AB = 0.0_f8
            do ijK = 1, channel1%njK
                j = AB_jKPair(ijK,1)
                l = qn_l_AB(ijK)
                fact = dsqrt((2.0_f8*l+1.0_f8)/(2.0_f8*initWP%Jtot+1.0_f8))
                do Kp = channel1%Kmin, Kpmax
                    if (Kp <= j) then 
                        BLK_AB(Kp,ijK) = fact*CG(j,Kp,l,0,initWP%Jtot)
                        if (Kp /= 0) BLK_AB(Kp,ijK) = BLK_AB(Kp,ijK)*dsqrt(2.0_f8)
                    end if
                end do
            end do

        else if (type == 'A+BC->B+AC') then
            Kpmax = min(channel2%jmax, initWP%Jtot)

            !> Transformation matrix Ka -> Kc
            allocate(TKMat_AC(channel2%nIC,initWP%Kmin:Krmax,channel2%Kmin:Kpmax))
            do Kp = channel2%Kmin, Kpmax
                do Kr = initWP%Kmin, Krmax
                    do iIC = 1, channel2%nIC
                        beta = interCoor_AC(iIC,6)
                        TKMat_AC(iIC,Kr,Kp) = dm(initWP%Jtot, Kp, Kr, beta) &
                                            + dm(initWP%Jtot, Kp, -Kr, beta) &
                                            * initWP%tpar*(-1.0_f8)**(Kr)
                    end do
                    if (Kr == 0) TKMat_AC(:,Kr,Kp) = TKMat_AC(:,Kr,Kp)/dsqrt(2.0_f8)
                    if (Kp == 0) TKMat_AC(:,Kr,Kp) = dsqrt(2.0_f8)*TKMat_AC(:,Kr,Kp)
                end do
            end do

            !> Qn_l in AC
            allocate(qn_l_AC(channel2%njK))
            do ijK = 1, channel2%njK
                j = AC_jKPair(ijK,1)
                do l = abs(initWP%Jtot - j), initWP%Jtot + j
                    if ((-1)**(l+j+initWP%Jtot) == initWP%tpar) then
                        qn_l_AC(ijK) = l
                    end if
                end do
            end do

            !> BLK in AC
            allocate(BLK_AC(channel2%Kmin:Kpmax,channel2%njK))
            BLK_AC = 0.0_f8
            do ijK = 1, channel2%njK
                j = AC_jKPair(ijK,1)
                l = qn_l_AC(ijK)
                fact = dsqrt((2.0_f8*l+1.0_f8)/(2.0_f8*initWP%Jtot+1.0_f8))
                do Kp = channel2%Kmin, Kpmax
                    if (Kp <= j) then 
                        BLK_AC(Kp,ijK) = fact*CG(j,Kp,l,0,initWP%Jtot)
                        if (Kp /= 0) BLK_AC(Kp,ijK) = BLK_AC(Kp,ijK)*dsqrt(2.0_f8)
                    end if
                end do
            end do
        end if

    end subroutine setTransKMat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setAtDMat(nIC, ichoice, massB, massC, interCoor, AtDMat)
        implicit none
        integer, intent(in) :: nIC
        integer, intent(in) :: ichoice
        !> AB or AC
        real(f8), intent(in) :: massB, massC
        real(f8), intent(in) :: interCoor(nIC,7)
        real(f8), allocatable, intent(out) :: AtDMat(:,:,:)
        real(f8) :: bond(3), AtD(nPES,nPES), Vadia(nPES)
        real(f8) :: Z, rij, r, Zij, th  
        integer :: iIC

        allocate(AtDMat(nPES,nPES,nIC))
        do iIC = 1, nIC 
            if (ichoice == 1) then
                Z = interCoor(iIC,1)
                rij = interCoor(iIC,2)
                th = acos(interCoor(iIC,3))
                call Jacobi2Bond(Z,rij,th,massB,massC,bond)
            else if (ichoice == 2) then
                r = interCoor(iIC,1)
                Zij = interCoor(iIC,2)
                th = acos(interCoor(iIC,3))
                call Jacobi2Bond(Zij,r,th,massB,massC,bond)
            end if
            call diagDiaVmat(bond, AtD, Vadia)
            AtDMat(:,:,iIC) = AtD(:,:)
        end do

    end subroutine setAtDMat
!> ------------------------------------------------------------------------------------------------------------------ <!
    
!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setNumberInterCoor(type, ichoice, Z0p, r_range, mass, nIC, interCoor, idXY)
        implicit none
        character(len=*), intent(in) :: type
        !> ichoice: 1 for (Z,theta), 2 for (r,theta)
        integer, intent(in) :: ichoice
        real(f8), intent(in) :: Z0p
        real(f8), intent(in) :: r_range(2)
        !> Mass in A, B, C order
        real(f8), intent(in) :: mass(3)
        integer, intent(out) :: nIC
        real(f8), allocatable, intent(out) :: interCoor(:,:)
        integer, allocatable, intent(out) :: idXY(:,:)
        real(f8) :: R_Jacobi(3), P_Jacobi(3)
        real(f8) ::  costh, Z, r, rp, beta
        real(f8) :: rij, Zij, gij
        integer :: ir, iZ, ith, iIC 

        A_beta = (mass(1)+mass(2))/mass(1)-mass(3)/(mass(2)+mass(3))
        B_beta = (mass(1)+mass(2))/mass(1)
        A_gamma = (mass(1)+mass(3))/mass(1)-mass(2)/(mass(2)+mass(3))
        B_gamma = (mass(1)+mass(3))/mass(1)

        nIC = 0
        if (ichoice == 1) then
            do iZ = 1, IALR%nZ_I
                do ith = 1, IALR%int_nA
                    costh = int_ANode(ith)
                    Z = Z_IALR(iZ)
                    if (type == 'A+BC->C+AB') then
                        rij = RCB_rmid(A_beta, B_beta, Z0p, Z, costh)
                    else if (type == 'A+BC->B+AC') then
                        rij = RCB_rmid(A_gamma, B_gamma, Z0p, Z, -costh)
                    end if
                    if (rij > r_Int(1) .and. rij < r_Int(IALR%vint)) then
                        R_Jacobi(1) = rij 
                        R_Jacobi(2) = Z
                        R_Jacobi(3) = costh
                        call Jacobi_R2P(type, mass, R_Jacobi, P_Jacobi)
                        rp = P_Jacobi(1)
                        if (rp > r_range(1) .and. rp < r_range(2)) then
                            nIC = nIC + 1
                        end if
                    end if
                end do
            end do
        else if (ichoice == 2) then
            do ir = 1, IALR%vint
                do ith = 1, IALR%int_nA
                    costh = int_ANode(ith)
                    r = r_Int(ir)
                    if (type == 'A+BC->C+AB') then
                        Zij = RCB_Zmid(A_beta, B_beta, Z0p, r, costh)
                    else if (type == 'A+BC->B+AC') then
                        Zij = RCB_Zmid(A_gamma, B_gamma, Z0p, r, -costh)
                    end if
                    if (Zij > Z_IALR(1) .and. Zij < Z_IALR(IALR%nZ_I)) then
                        R_Jacobi(1) = r
                        R_Jacobi(2) = Zij
                        R_Jacobi(3) = costh
                        call Jacobi_R2P(type, mass, R_Jacobi, P_Jacobi)
                        rp = P_Jacobi(1)
                        if (rp > r_range(1) .and. rp < r_range(2)) then
                            nIC = nIC + 1
                        end if
                    end if
                end do
            end do
        else
            write(outFileUnit,*) 'Error in setNumberInterCoor: unknown ichoice ', ichoice
            write(outFileUnit,*) 'POSITION: interCoorRCB.f90, subroutine setNumberInterCoor'
            stop
        end if

        if (nIC == 0) then
            write(outFileUnit,*) 'Error in setNumberInterCoor: nIC = 0!'
            write(outFileUnit,*) 'POSITION: interCoorRCB.f90, subroutine setNumberInterCoor'
            stop
        end if

        allocate(interCoor(nIC,7))
        allocate(idXY(nIC,2))

        iIC = 0
        if (ichoice == 1) then
            do iZ = 1, IALR%nZ_I
                do ith = 1, IALR%int_nA
                    costh = int_ANode(ith)
                    Z = Z_IALR(iZ)
                    if (type == 'A+BC->C+AB') then
                        rij = RCB_rmid(A_beta, B_beta, Z0p, Z, costh)
                        gij = rJacobiFactor(A_beta, B_beta, Z0p, Z, costh)
                    else if (type == 'A+BC->B+AC') then
                        rij = RCB_rmid(A_gamma, B_gamma, Z0p, Z, -costh)
                        gij = rJacobiFactor(A_gamma, B_gamma, Z0p, Z, -costh)
                    end if
                    if (rij > r_Int(1) .and. rij < r_Int(IALR%vint)) then
                        R_Jacobi(1) = rij 
                        R_Jacobi(2) = Z
                        R_Jacobi(3) = costh
                        call Jacobi_R2P(type, mass, R_Jacobi, P_Jacobi, beta)
                        rp = P_Jacobi(1)
                        if (rp > r_range(1) .and. rp < r_range(2)) then
                            iIC = iIC + 1
                            interCoor(iIC,1) = Z 
                            interCoor(iIC,2) = rij 
                            interCoor(iIC,3) = costh
                            interCoor(iIC,4) = rp 
                            interCoor(iIC,5) = P_Jacobi(3)
                            interCoor(iIC,6) = beta
                            interCoor(iIC,7) = gij
                            idXY(iIC,1) = iZ
                            idXY(iIC,2) = ith
                        end if
                    end if
                end do
            end do
        else if (ichoice == 2) then
            do ir = 1, IALR%vint
                do ith = 1, IALR%int_nA
                    costh = int_ANode(ith)
                    r = r_Int(ir)
                    if (type == 'A+BC->C+AB') then
                        Zij = RCB_Zmid(A_beta, B_beta, Z0p, r, costh)
                        gij = ZJacobiFactor(A_beta, B_beta, Z0p, r, costh)
                    else if (type == 'A+BC->B+AC') then
                        Zij = RCB_Zmid(A_gamma, B_gamma, Z0p, r, -costh)
                        gij = ZJacobiFactor(A_gamma, B_gamma, Z0p, r, -costh)
                    end if
                    if (Zij > Z_IALR(1) .and. Zij < Z_IALR(IALR%nZ_I)) then
                        R_Jacobi(1) = r
                        R_Jacobi(2) = Zij
                        R_Jacobi(3) = costh
                        call Jacobi_R2P(type, mass, R_Jacobi, P_Jacobi, beta)
                        rp = P_Jacobi(1)
                        if (rp > r_range(1) .and. rp < r_range(2)) then
                            iIC = iIC + 1
                            interCoor(iIC,1) = r 
                            interCoor(iIC,2) = Zij 
                            interCoor(iIC,3) = costh
                            interCoor(iIC,4) = rp
                            interCoor(iIC,5) = P_Jacobi(3)
                            interCoor(iIC,6) = beta
                            interCoor(iIC,7) = gij
                            idXY(iIC,1) = ir
                            idXY(iIC,2) = ith
                        end if
                    end if
                end do
            end do
        end if

    end subroutine setNumberInterCoor
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setUij(ichoice, nIC, interCoor, idXY, Uij)
        implicit none
        !> ichoice: 1 for (Z,theta), 2 for (r,theta)
        integer, intent(in) :: ichoice
        integer, intent(in) :: nIC
        real(f8), intent(in) :: interCoor(nIC,7)
        integer, intent(in) :: idXY(nIC,2)
        real(f8), allocatable, intent(out) :: Uij(:,:)
        real(f8) :: rFact, Zfact
        integer :: iZ, ir, iU, iIC

        rFact = 2.0_f8/dsqrt((IALR%r_range(2)-IALR%r_range(1))*(IALR%vint+1.0_f8))
        Zfact = 2.0_f8/dsqrt((Z_IALR(IALR%nZ_I)-Z_IALR(1))*(IALR%nZ_I+1.0_f8))

!> Set the Uij transformation matrix
        if (ichoice == 1) then
            !> Uij(r_ij,nIC)
            allocate(Uij(IALR%vint,nIC))
        else if (ichoice == 2) then
            !> Uij(Z_ij,nIC)
            allocate(Uij(IALR%nZ_I,nIC))
        else
            write(outFileUnit,*) 'Error in setUij: unknown ichoice ', ichoice
            write(outFileUnit,*) 'POSITION: interCoorRCB.f90, subroutine setUij'
            stop
        end if

        do iIC = 1, nIC
            if (ichoice == 1) then
                do iU = 1, IALR%vint
                    Uij(iU,iIC) = 0.0_f8
                    do ir = 1, IALR%vint 
                        Uij(iU,iIC) = Uij(iU,iIC) &
                                    + dsin(ir*pi*(interCoor(iIC,2)-IALR%r_range(1))/(IALR%r_range(2)-IALR%r_range(1))) &
                                    * dsin(iU*ir*pi/(IALR%vint+1.0_f8)) 
                    end do
                    Uij(iU,iIC) = Uij(iU,iIC) * rFact * interCoor(iIC,7)
                end do
            else if (ichoice == 2) then
                do iU = 1, IALR%nZ_I 
                    Uij(iU,iIC) = 0.0_f8
                    do iZ = 1, IALR%nZ_I 
                        Uij(iU,iIC) = Uij(iU,iIC) &
                                    + dsin(iZ*pi*(interCoor(iIC,1)-Z_IALR(1))/(Z_IALR(IALR%nZ_I)-Z_IALR(1))) &
                                    * dsin(iU*iZ*pi/(IALR%nZ_I+1.0_f8)) 
                    end do
                    Uij(iU,iIC) = Uij(iU,iIC) * Zfact * interCoor(iIC,7)
                end do
            end if
        end do

    end subroutine setUij
!> ------------------------------------------------------------------------------------------------------------------ <!
        
!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine Jacobi_R2P(type, mass, R_Jacobi, P_Jacobi, beta)
        implicit none
        character(len=*), intent(in) :: type
        !> Mass in A, B, C order
        real(f8), intent(in) :: mass(3)
        !> Jacobi in (r,Z,costheta) order
        real(f8), intent(in) :: R_Jacobi(3)
        real(f8), intent(out) :: P_Jacobi(3)
        real(f8), intent(out), optional :: beta
        real(f8) :: a1, a2, a3 

!> Transform Jacobi coordinates from reactant to product
!> Ref: J. Chem. Phys. 106 (5), 1 February 1997
        if (type == 'A+BC->C+AB') then 
            a1 = mass(3)/(mass(2)+mass(3))
            a2 = mass(1)/(mass(1)+mass(2))
            a3 = a1*a2 - 1.0_f8
        else if (type == 'A+BC->B+AC') then 
            a1 = -mass(2)/(mass(2)+mass(3))
            a2 = -mass(1)/(mass(1)+mass(3))
            a3 = a1*a2 - 1.0_f8
        else
            write(outFileUnit,*) 'Error in Jacobi_R2P: unknown type ', type
            write(outFileUnit,*) 'POSITION: interCoorRCB.f90, subroutine Jacobi_R2P'
            stop
        end if

        P_Jacobi(1) = dsqrt(R_Jacobi(2)**2 & 
                            + (a1*R_Jacobi(1))**2 &
                            + 2.0_f8*a1*R_Jacobi(1)*R_Jacobi(2)*R_Jacobi(3))
        P_Jacobi(2) = dsqrt((a2*R_Jacobi(2))**2 & 
                            + (a3*R_Jacobi(1))**2 &
                            + 2.0_f8*a2*a3*R_Jacobi(1)*R_Jacobi(2)*R_Jacobi(3))
        P_Jacobi(3) = (a2*R_Jacobi(2)**2 &
                        + a1*a3*R_Jacobi(1)**2 &
                        + (a1*a2+a3)*R_Jacobi(1)*R_Jacobi(2)*R_Jacobi(3)) &
                        / (P_Jacobi(1)*P_Jacobi(2))

        if (present(beta)) then
            if (type == 'A+BC->C+AB') then 
                beta = (-a2*R_Jacobi(2)-a3*R_Jacobi(1)*R_Jacobi(3)) / P_Jacobi(2)
            else if (type == 'A+BC->B+AC') then
                beta = (a2*R_Jacobi(2)+a3*R_Jacobi(1)*R_Jacobi(3)) / P_Jacobi(2)
            end if

            if (dabs(beta) > 1.0_f8) then 
                write(outFileUnit,*) 'Error in Jacobi_R2P: |cos(beta)| > 1 ', beta
                write(outFileUnit,*) 'POSITION: interCoorRCB.f90, subroutine Jacobi_R2P'
                stop
            end if
            beta = dacos(beta)
        end if

    end subroutine Jacobi_R2P
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    real(f8) function RCB_Zmid(A, B, Z0P, r, costh) result(Zmid)
        implicit none
        real(f8), intent(in) :: A, B
        real(f8), intent(in) :: Z0P, r, costh

!> A+BC -> C+AB, th = theta 
!> A+BC -> B+AC, th = pi-theta, should intent -costh

        Zmid = A*r*costh &
                + dsqrt((A*r)**2*(costh**2-1.0_f8) + (B*Z0P)**2)

    end function RCB_Zmid 
!> ------------------------------------------------------------------------------------------------------------------ <!
    
!> ------------------------------------------------------------------------------------------------------------------ <!
    real(f8) function RCB_rmid(A, B, Z0P, Z, costh) result(rmid)
        implicit none
        real(f8), intent(in) :: A, B
        real(f8), intent(in) :: Z0P, Z, costh

        rmid = (Z*costh  &
                + dsqrt(Z*Z*(costh**2-1.0_f8) + (B*Z0P)**2)) / A   
    
    end function RCB_rmid
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    real(f8) function rJacobiFactor(A, B, Z0P, Z, costh) result(r_gij)
        implicit none
        real(f8), intent(in) :: A, B
        real(f8), intent(in) :: Z0P, Z, costh

        r_gij = B*B*Z0P &
                / A / dsqrt(Z*Z*(costh**2-1.0_f8) + (B*Z0P)**2)
        r_gij = dsqrt(r_gij)

    end function rJacobiFactor
!> ------------------------------------------------------------------------------------------------------------------ <!
    
!> ------------------------------------------------------------------------------------------------------------------ <!
    real(f8) function ZJacobiFactor(A, B, Z0P, r, costh) result(Z_gij)
        implicit none
        real(f8), intent(in) :: A, B
        real(f8), intent(in) :: Z0P, r, costh

        Z_gij = B*B*Z0P &
                / dsqrt((A*r)**2*(costh**2-1.0_f8) + (B*Z0P)**2)
        Z_gij = dsqrt(Z_gij)

    end function ZJacobiFactor



end module m_InterCoorRCB