module m_InterCoorRCB
    use m_MachinaBasic, only : f8
    use m_gPara
    use m_Basis
    use m_Potent
    implicit none

!> Mass constants in RCB, beta = C+AB and gamma = B+AC
    real(f8), private :: A_beta, B_beta
    real(f8), private :: A_gamma, B_gamma
!> Intermidiate coordinates rij/Zij and gij
    real(f8), allocatable, private :: interCoor_AB(:,:,:), interCoor_AC(:,:,:)
!> DVR grids
    real(f8), allocatable, private :: rp1(:), rp2(:)
!> DVR to FBR transformation matrices
    real(f8), allocatable, private :: DVR2FBR1(:,:), DVR2FBR2(:,:)
!> Product diatomic potential 
    real(f8), allocatable, private :: VAB(:,:), VAC(:,:)
!> Product diatomic vib-rot energy and wavefunction in product coordinates
    real(f8), allocatable, private :: Evj_AB(:,:,:), Evj_AC(:,:,:)
    real(f8), allocatable, private :: WFvj_AB(:,:,:,:), WFvj_AC(:,:,:,:)

    public
    private :: findjMax, DVR_VibRotp
    private :: setUij, interCoorWFvj
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
        !> Set Uij
        call setUij(type,channel1%ichoice,atomMass,channel1%Zinf,interCoor_AB,Uij_AB)
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

        if (nChannel == 2) then
            !> Channel 2
            type = 'A+BC->B+AC'
            !> Set Uij
            call setUij(type,channel2%ichoice,atomMass,channel2%Zinf,interCoor_AC,Uij_AC)
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
            write(outFileUnit,'(1x,a)') 'Please check InterCoorWF_AB.chk for the wavefunction normalization!'
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
                write(outFileUnit,'(1x,a)') 'Please check InterCoorWF_AC.chk for the wavefunction normalization!'
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
                ! Ensure accumulators are zeroed before +=
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
                ! Ensure accumulators are zeroed before +=
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
        real(f8) :: Jacobi_AB(3), Jacobi_AC(3), Jacobi_BC(3)
        real(f8), allocatable :: pos(:)
        real(f8), external :: spgndr
        integer :: iPos, nPos, ith, v, j, iPES, iDVR, ijK, K
        integer :: normUnit
        character(len=20) :: fileName 

        djFact = merge(0.5_f8, 1.0_f8, initWP%jinc==2)

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
            allocate(ICWFvj_AB(nPos,IALR%int_nA,0:channel1%vmax,channel1%njK,nPES))
            ICWFvj_AB = 0.0_f8
            do iPES = 1, nPES 
                do ith = 1, IALR%int_nA
                    ! For AB channel, use cos(pi - theta) = -cos(theta)
                    Jacobi_BC(3) = int_ANode(ith) 
                    do iPos = 1, nPos
                        if (ichoice == 1) then 
                            Jacobi_BC(1) = interCoor_AB(iPos,ith,1)
                            Jacobi_BC(2) = pos(iPos)
                        else if (ichoice == 2) then
                            Jacobi_BC(1) = pos(iPos)
                            Jacobi_BC(2) = interCoor_AB(iPos,ith,1)
                        end if
                        if (Jacobi_BC(1) > 0.0_f8 .and. Jacobi_BC(2) > 0.0_f8) then 
                            call Jacobi_R2P(type, atomMass, Jacobi_BC, Jacobi_AB)
                            rp = Jacobi_AB(1)
                            write(9876,*) rp
                            if (rp < channel1%r_range(1) .or. rp > channel1%r_range(2)) cycle
                            gij = interCoor_AB(iPos,ith,2)/channel1%Zinf/rp*pos(iPos)*interCoor_AB(iPos,ith,1)
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
                                    ICWFvj_AB(iPos,ith,v,ijK,iPES) = tmpW*tmp*spgndr(j,K,Jacobi_AB(3))
                                end do
                            end do
                        end if 
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
                        do ith = 1, IALR%int_nA
                            do iPos = 1, nPos
                                norm = norm + ICWFvj_AB(iPos,ith,v,ijK,iPES)**2
                            end do
                        end do
                        write(normUnit,'(1x,i4,1x,i4,1x,i4,1x,i4,1x,f15.9,1x,f15.9)') iPES, v, j, K, Evj_AB(v,j,iPES)*au2ev,norm
                    end do
                end do
            end do
            close(normUnit)
        else if (type == 'A+BC->B+AC') then
            allocate(ICWFvj_AC(nPos,IALR%int_nA,0:channel2%vmax,channel2%njK,nPES))
            ICWFvj_AC = 0.0_f8
            do iPES = 1, nPES 
                do ith = 1, IALR%int_nA
                    Jacobi_BC(3) = int_ANode(ith) 
                    do iPos = 1, nPos
                        if (ichoice == 1) then 
                            Jacobi_BC(1) = interCoor_AC(iPos,ith,1)
                            Jacobi_BC(2) = pos(iPos)
                        else if (ichoice == 2) then
                            Jacobi_BC(1) = pos(iPos)
                            Jacobi_BC(2) = interCoor_AC(iPos,ith,1)
                        end if
                        if (Jacobi_BC(1) > 0.0_f8 .and. Jacobi_BC(2) > 0.0_f8) then
                            call Jacobi_R2P(type, atomMass, Jacobi_BC, Jacobi_AC)
                            rp = Jacobi_AC(1)
                            if (rp < channel2%r_range(1) .or. rp > channel2%r_range(2)) cycle
                            gij = interCoor_AC(iPos,ith,2)/channel2%Zinf/rp*pos(iPos)*interCoor_AC(iPos,ith,1)
                            write(9876,*) interCoor_AC(iPos,ith,1)
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
                                    ICWFvj_AC(iPos,ith,v,ijK,iPES) = tmp*tmpW*spgndr(j,K,Jacobi_AC(3))
                                end do 
                            end do
                        end if
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
                        do ith = 1, IALR%int_nA
                            do iPos = 1, nPos
                                norm = norm + ICWFvj_AC(iPos,ith,v,ijK,iPES)**2
                            end do
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
    subroutine setUij(type, ichoice, mass, Z0P, interCoor, Uij)
        implicit none
        character(len=*), intent(in) :: type
        !> ichoice: 1 for (Z,theta), 2 for (r,theta)
        integer, intent(in) :: ichoice
        !> Mass in A, B, C order
        real(f8), intent(in) :: mass(3)
        real(f8), intent(in) :: Z0P
        !> For rij or Zij and gij
        real(f8), allocatable, intent(out) :: interCoor(:,:,:)
        real(f8), allocatable, intent(out) :: Uij(:,:,:)
        real(f8) :: r_gij, Z_gij
        real(f8) :: costh, r, Z, rmid, Zmid
        real(f8) :: rFact, Zfact, A, B
        integer :: iZ, ith, ir, iU

        A_beta = (mass(1)+mass(2))/mass(1)-mass(3)/(mass(2)+mass(3))
        B_beta = (mass(1)+mass(2))/mass(1)
        A_gamma = (mass(1)+mass(3))/mass(1)-mass(2)/(mass(2)+mass(3))
        B_gamma = (mass(1)+mass(3))/mass(1)
        rFact = 2.0_f8/dsqrt((IALR%r_range(2)-IALR%r_range(1))*(IALR%vint+1.0_f8))
        Zfact = 2.0_f8/dsqrt((Z_IALR(IALR%nZ_I)-Z_IALR(1))*(IALR%nZ_I+1.0_f8))

        if (type == 'A+BC->C+AB') then
            A = A_beta
            B = B_beta
        else if (type == 'A+BC->B+AC') then
            A = A_gamma
            B = B_gamma
        else
            write(outFileUnit,*) 'Error in setUij: unknown type ', type
            write(outFileUnit,*) 'POSITION: interCoorRCB.f90, subroutine setUij'
            stop
        end if

!> Set the Uij transformation matrix
        if (ichoice == 1) then
            !> Uij(r_ij,Z_i,theta_j)
            allocate(Uij(IALR%vint,IALR%nZ_I,IALR%int_nA))
            allocate(interCoor(IALR%nZ_I,IALR%int_nA,2))
        else if (ichoice == 2) then
            !> Uij(Z_ij,r_i,theta_j)
            allocate(Uij(IALR%nZ_I,IALR%vint,IALR%int_nA))
            allocate(interCoor(IALR%vint,IALR%int_nA,2))
        else
            write(outFileUnit,*) 'Error in setUij: unknown ichoice ', ichoice
            write(outFileUnit,*) 'POSITION: interCoorRCB.f90, subroutine setUij'
            stop
        end if
        do ith = 1, IALR%int_nA 
            if (type == 'A+BC->C+AB') then
                ! For AB channel, geometry implies pi - theta
                costh = int_ANode(ith)
            else if (type == 'A+BC->B+AC') then
                costh = -int_ANode(ith)
            end if

            if (ichoice == 1) then
                do iZ = 1, IALR%nZ_I 
                    Z = Z_IALR(iZ)
                    r_gij = rJacobiFactor(A, B, Z0P, Z, costh)
                    rmid = RCB_rmid(A, B, Z0P, Z, costh)
                    if (rmid < r_Int(1) .or. rmid > r_Int(IALR%vint)) rmid = 0.0_f8
                    interCoor(iZ,ith,1) = rmid
                    interCoor(iZ,ith,2) = r_gij

                    do iU = 1, IALR%vint
                        Uij(iU,iZ,ith) = 0.0_f8
                        do ir = 1, IALR%vint 
                            Uij(iU,iZ,ith) = Uij(iU,iZ,ith) &
                                            + dsin(ir*pi*(rmid-IALR%r_range(1))/(IALR%r_range(2)-IALR%r_range(1))) &
                                            * dsin(iU*ir*pi/(IALR%vint+1.0_f8)) 
                        end do
                        Uij(iU,iZ,ith) = Uij(iU,iZ,ith) * rFact * r_gij
                    end do
                end do
            else if (ichoice == 2) then
                do ir = 1, IALR%vint
                    r = r_Int(ir)
                    Z_gij = ZJacobiFactor(A, B, Z0P, r, costh)
                    Zmid = RCB_Zmid(A, B, Z0P, r, costh)
                    if (Zmid < Z_IALR(1) .or. Zmid > Z_IALR(IALR%nZ_I)) Zmid = 0.0_f8
                    interCoor(ir,ith,1) = Zmid
                    interCoor(ir,ith,2) = Z_gij

                    do iU = 1, IALR%nZ_I 
                        Uij(iU,ir,ith) = 0.0_f8
                        do iZ = 1, IALR%nZ_I 
                            Uij(iU,ir,ith) = Uij(iU,ir,ith) &
                                            + dsin(iZ*pi*(Zmid-Z_IALR(1))/(Z_IALR(IALR%nZ_I)-Z_IALR(1))) &
                                            * dsin(iU*iZ*pi/(IALR%nZ_I+1.0_f8)) 
                        end do
                        Uij(iU,ir,ith) = Uij(iU,ir,ith) * Zfact * Z_gij
                    end do
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
                beta = (a2*R_Jacobi(2)+a3*R_Jacobi(1)*R_Jacobi(3)) / P_Jacobi(2)
            else if (type == 'A+BC->B+AC') then
                beta = (-a2*R_Jacobi(2)-a3*R_Jacobi(1)*R_Jacobi(3)) / P_Jacobi(2)
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