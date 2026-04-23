module m_Basis
    use m_MachinaBasic, only : f8, c8
    use m_gPara
    use m_Potent
    implicit none
    public

    private :: setANodeAndWeight

contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine sinDVRGrid(nGrid, rMin, rMax, Grid)
        implicit none
        integer, intent(in) :: nGrid
        real(f8), intent(in) :: rMin, rMax
        real(f8), intent(inout) :: Grid(:)
        integer :: i
    
        do i = 1, nGrid
            Grid(i) = rMin + i*(rMax - rMin)/(nGrid + 1)
        end do
    end subroutine sinDVRGrid
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setDVR2FBR(nZ, BMat)
        implicit none
        integer, intent(in) :: nZ
        real(f8), intent(out) :: BMat(nZ,nZ)
        integer :: i, j
        real(f8) :: fact 

        fact = dsqrt(2.0_f8/(nZ+1))
        do i = 1, nZ
            do j = 1, i 
                BMat(i,j) = fact * dsin(pi * i * j / (nZ + 1))
                BMat(j,i) = BMat(i,j)
            end do
        end do

    end subroutine setDVR2FBR
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setjKPair(jmin, jmax, jinc, Kmin, njK, jKPair)
        implicit none
        integer, intent(in) :: jmin, jmax, jinc, Kmin
        integer, intent(out) :: njK
        integer, allocatable, intent(out) :: jKPair(:,:)
        integer :: j, K, Kmax
        integer :: idx

        njK = 0
        Kmax = min(jmax, initWP%Jtot)
        do K = Kmin, Kmax
            do j = jmin, jmax, jinc
                if (j >= K) then
                    njK = njK + 1
                end if
            end do
        end do
        allocate(jKPair(njK, 2))

        idx = 0
        do K = Kmin, Kmax
            do j = jmin, jmax, jinc
                if (j >= K) then
                    idx = idx + 1
                    jKPair(idx, 1) = j
                    jKPair(idx, 2) = K
                end if
            end do
        end do

    end subroutine setjKPair
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setANodeAndWeight(jpar, nth, nodes, weights)
        implicit none
        integer, intent(in) :: jpar, nth
        real(f8), intent(inout) :: nodes(:), weights(:)
        integer :: i, ntmp 

!> Take advantage of the symmetry of potential the grid used for angle could be used half.
        if (jpar /= 0) then 
            ntmp = 2*nth 
        else
            ntmp = nth 
        end if 

!> First call ntmp = 2*nth for the system has symmetry.
        block 
            real(f8) :: tmpNode(ntmp), tmpWeight(ntmp)
            
            call GAULEG(-1.0_f8, 1.0_f8, tmpNode, tmpWeight, ntmp)
!> Take half of grid
            do i = 1, nth 
                nodes(i) = tmpNode(i)
                if (jpar /= 0) then
                    weights(i) = 2.0_f8*tmpWeight(i)
                else 
                    weights(i) = tmpWeight(i)
                end if 
            end do 
        end block

    end subroutine setANodeAndWeight
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setA_DVR2FBR(njK, nth, jKPair, ANode, AWeight, A_DVR2FBR)
        implicit none
        integer, intent(in) :: njK, nth
        integer, intent(in) :: jKPair(:,:)
        real(f8), intent(in) :: ANode(:), AWeight(:)
        real(f8), intent(out) :: A_DVR2FBR(:,:)
        integer :: ijK, ith, j, K
        real(f8), external :: spgndr 

        do ijK = 1, njK
            do ith = 1, nth
                j = jKPair(ijK,1)
                K = jKPair(ijK,2)
                if (K > j) then
                    write(outFileUnit,*) 'Error in setA_DVR2FBR: K > j'
                    write(outFileUnit,*) 'POSITION: basis.f90, subroutine setA_DVR2FBR()'
                    stop
                end if
                A_DVR2FBR(ith,ijK) = spgndr(j,K,ANode(ith)) * dsqrt(AWeight(ith))
            end do
        end do
    end subroutine setA_DVR2FBR
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine findPositionID(nGrid, Grid, pos, posID)
        implicit none
        integer, intent(in) :: nGrid
        real(f8), intent(in) :: Grid(nGrid)
        real(f8), intent(in) :: pos
        integer, intent(out) :: posID
        integer :: i

        posID = 0
        do i = 1, nGrid
            if (Grid(i) >= pos) then
                posID = i
                exit
            end if
        end do
        if (posID == 0) then
            write(outFileUnit,*) 'Error in findPositionID: position out of grid range'
            write(outFileUnit,*) 'POSITION: basis.f90, subroutine findPositionID()'
            stop
        end if

    end subroutine findPositionID
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine PODVR_BCvj()
        implicit none
        integer :: ir, i, j, v, qn_j, K, ijK, ith, iPES, Kmax
        real(f8) :: mass, r_range, tmp
        real(f8) :: Vadia(nPES), Vpot(IALR%nr_DVR), VBC(IALR%nr_DVR, nPES)
        real(f8) :: DVREig(IALR%nr_DVR), DVRWF(IALR%nr_DVR, IALR%nr_DVR)
        real(f8) :: HRefMat(IALR%vasy,IALR%vasy), EMat(IALR%vasy,IALR%vasy)
        real(f8) :: XMat(IALR%vasy,IALR%vasy), poTransMat(IALR%nr_DVR, IALR%vasy)
        real(f8) :: poEig(IALR%vasy), poGrid(IALR%vasy)
        real(f8), allocatable :: work(:)
        integer :: lwork, info
        character(len=10) :: type

        mass = massBC
        r_range = IALR%r_range(2) - IALR%r_range(1)
        Kmax = min(IALR%jasy, initWP%Jtot)
        type = 'A+BC->A+BC'

        allocate(Evj_BC(0:IALR%vasy-1,initWP%jmin:IALR%jasy,nPES))
        allocate(WFvj_BC(IALR%vasy,IALR%asy_nA,0:IALR%vasy-1, &
                        initWP%jmin:IALR%jasy,initWP%Kmin:Kmax,nPES))
        
        Evj_BC = 0.0_f8
        WFvj_BC = 0.0_f8

        do ir = 1, IALR%nr_DVR
            call setPot_Product(type, r_DVR(ir), Vadia)
            VBC(ir,:) = Vadia
        end do

        do iPES = 1, nPES 
            Vpot(:) = VBC(:,iPES)
            call DVR_vib(IALR%nr_DVR,mass,r_range,Vpot,DVREig,DVRWF)

            do i = 1, IALR%nr_DVR
                do j = 1, IALR%vasy
                    poTransMat(i,j) = DVRWF(i,j)
                end do 
            end do

            do i = 1, IALR%vasy
                call phaseTrans(IALR%nr_DVR, poTransMat(:,i))
            end do

            do i = 1, IALR%vasy
                do j = i, IALR%vasy
                    tmp = 0.0_f8
                    do ir = 1, IALR%nr_DVR
                        tmp = tmp + poTransMat(ir,i) * r_DVR(ir) * poTransMat(ir,j)
                    end do
                    XMat(i,j) = tmp
                    XMat(j,i) = tmp
                end do
            end do

            allocate(work(1))
            call dsyev('V', 'U', IALR%vasy, XMat, IALR%vasy, poGrid, work, -1, info)
            lwork = int(work(1))
            deallocate(work)

            allocate(work(lwork))
            call dsyev('V', 'U', IALR%vasy, XMat, IALR%vasy, poGrid, work, lwork, info)
            if (info /= 0) then
                write(outFileUnit,*) "Error in DSYEV of PODVR XMat, info =", info
                write(outFileUnit,*) "POSITION: basis.f90, subroutine PODVR_BCvj()"
                stop
            end if
            deallocate(work)

            do i = 1, IALR%vasy
                call phaseTrans(IALR%vasy, XMat(:,i))
            end do

            do i = 1, IALR%vasy
                do j = 1, i
                    tmp = 0.0_f8
                    do ir = 1, IALR%vasy
                        tmp = tmp + XMat(ir,i) * DVREig(ir) * XMat(ir,j)
                    end do
                    HRefMat(i,j) = tmp
                    HRefMat(j,i) = tmp
                end do
            end do

            do qn_j = initWP%jmin, IALR%jasy, initWP%jinc
                EMat = HRefMat
                do i = 1, IALR%vasy
                    EMat(i,i) = EMat(i,i) + &
                                (qn_j*(qn_j+1))/(2.0_f8*mass*poGrid(i)**2)
                end do

                allocate(work(1))
                call dsyev('V', 'U', IALR%vasy, EMat, IALR%vasy, poEig, work, -1, info)
                lwork = int(work(1))
                deallocate(work)

                allocate(work(lwork))
                call dsyev('V', 'U', IALR%vasy, EMat, IALR%vasy, poEig, work, lwork, info)
                if (info /= 0) then
                    write(outFileUnit,*) "Error in DSYEV of PODVR HRefMat, info =", info
                    write(outFileUnit,*) "POSITION: basis.f90, subroutine PODVR_BCvj()"
                    stop
                end if
                deallocate(work)

                do i = 1, IALR%vasy 
                    call phaseTrans(IALR%vasy, EMat(:,i))
                end do 

                do v = 0, IALR%vasy-1
                    Evj_BC(v, qn_j, iPES) = poEig(v+1)
                end do 

                do i = 1, IALR%vasy 
                    do v = 0, IALR%vasy-1
                        do K = initWP%Kmin, min(qn_j, initWP%Jtot)
                            ijK = asy_seqjK(qn_j,K)
                            if (ijK > 0) then
                                do ith = 1, IALR%asy_nA
                                    WFvj_BC(i,ith,v,qn_j,K,iPES) = EMat(i,v+1)*asy_YMat(ith,ijK)
                                end do
                            end if
                        end do
                    end do
                end do
            end do
        end do 

    end subroutine PODVR_BCvj
!> ------------------------------------------------------------------------------------------------------------------ <!

    subroutine setBasis()
        implicit none
        real(f8) :: range
        real(f8), allocatable :: DVRGrid(:), DVREig(:), DVRWF(:,:)
        real(f8), allocatable :: VBC(:)
        real(f8) :: r, Ev0j0
        integer :: ir, Kmax, K, j, ijK, jpar_eff
        integer :: EvjUnit
        
!> Z in sin-DVR : |Z_l>
        allocate(Z_IALR(IALR%nZ_IALR))
        call sinDVRGrid(IALR%nZ_IALR, IALR%Z_range(1), IALR%Z_range(2), Z_IALR)
        
!> Find the position ID for inelastic
        if (IF_inelastic) then
            call findPositionID(IALR%nZ_IALR, Z_IALR, inePos, iInePos)
            if (iInePos <= IALR%nZ_I .or. iInePos > IALR%nZ_IA) then
                write(outFileUnit,*) 'Error in setBasis: inelastic position out of asymptotic range'
                write(outFileUnit,*) 'POSITION: basis.f90, subroutine setBasis()'
                stop
            end if
        end if

!> r in PODVR : |r_p>
        IALR%vlr = initWP%v0+1
        allocate(r_Int(IALR%vint), r_Asy(IALR%vasy), r_LR(IALR%vlr))
        allocate(int_POWF(IALR%vint,IALR%vint), asy_POWF(IALR%vasy,IALR%vasy), lr_POWF(IALR%vlr, IALR%vlr))
        allocate(int_PO2FBR(IALR%vint,IALR%vint), asy_PO2FBR(IALR%vasy,IALR%vasy), lr_PO2FBR(IALR%vlr, IALR%vlr))
        allocate(int_POEig(IALR%vint), asy_POEig(IALR%vasy), lr_POEig(IALR%vlr))
        allocate(r_DVR(IALR%nr_DVR))
        allocate(DVRGrid(IALR%nr_DVR), DVREig(IALR%nr_DVR), DVRWF(IALR%nr_DVR, IALR%nr_DVR))
        allocate(VBC(IALR%nr_DVR))
        allocate(int_PO2DVR(IALR%nr_DVR,IALR%vint))
        allocate(asy_PO2DVR(IALR%nr_DVR,IALR%vasy))
        call sinDVRGrid(IALR%nr_DVR, IALR%r_range(1), IALR%r_range(2),DVRGrid)
        r_DVR = DVRGrid
        range = IALR%r_range(2) - IALR%r_range(1)
        do ir = 1, IALR%nr_DVR
            r = DVRGrid(ir)
            call setPot_vj0(r,VBC(ir))
        end do
        call DVR_vib(IALR%nr_DVR,massBC,range,VBC,DVREig,DVRWF)
        do ir = 1, IALR%nr_DVR
            call phaseTrans(IALR%nr_DVR, DVRWF(:,ir))
        end do
!> long range region
        call PODVR(IALR%vlr,IALR%nr_DVR,DVRWF,DVRGrid,DVREig,r_LR,lr_POEig,lr_POWF,lr_PO2FBR)
!> asymptotic region
        call PODVR(IALR%vasy,IALR%nr_DVR,DVRWF,DVRGrid,DVREig,r_Asy,asy_POEig,asy_POWF,asy_PO2FBR)
!> interaction region
        int_POEig = 0.0_f8; int_POWF = 0.0_f8; int_PO2FBR = 0.0_f8
        asy_PO2DVR = DVRWF(:,1:IALR%vasy)
        int_PO2DVR = DVRWF(:,1:IALR%vint)
        call PODVR(IALR%vint,IALR%nr_DVR,DVRWF,DVRGrid,DVREig,r_Int,int_POEig,int_POWF,int_PO2FBR)

        write(outFileUnit,'(1x,a)') '=====> Initail ro-vibrational state information <====='
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a)') 'lr-PODVR grids range: [', r_LR(1), ', ', r_LR(IALR%vlr), '  ] a.u.'
        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a)') 'asy-PODVR grids range: [', r_Asy(1), ', ', r_Asy(IALR%vasy), '  ] a.u.'
        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a)') 'int-PODVR grids range: [', r_Int(1), ', ', r_Int(IALR%vint), '  ] a.u.'
        write(outFileUnit,'(1x,a)') "Please ensure that the PODVR grid covers the relevant region of the BC potential!"
        write(outFileUnit,*) ''

!> theta in FBR : |JMjK>
!> long range region
        call setjKPair(initWP%j0, initWP%j0, 1, initWP%Kmin, lr_njK, lr_jKPair)
        Kmax = min(initWP%j0, initWP%Jtot)
        allocate(lr_seqjK(initWP%j0:initWP%j0,initWP%Kmin:Kmax))
        lr_seqjK = -1
        do ijK = 1, lr_njK
            j = lr_jKPair(ijK,1)
            K = lr_jKPair(ijK,2)
            lr_seqjK(j,K) = ijK
        end do
!> asymptotic region
        call setjKPair(initWP%jmin, IALR%jasy, initWP%jinc, initWP%Kmin, asy_njK, asy_jKPair)
        Kmax = min(IALR%jasy, initWP%Jtot)
        allocate(asy_seqjK(initWP%jmin:IALR%jasy,initWP%Kmin:Kmax))
        asy_seqjK = -1
        do ijK = 1, asy_njK
            j = asy_jKPair(ijK,1)
            K = asy_jKPair(ijK,2)
            asy_seqjK(j,K) = ijK
        end do
!> interaction region
        call setjKPair(initWP%jmin, IALR%jint, initWP%jinc, initWP%Kmin, int_njK, int_jKPair)
        Kmax = min(IALR%jint, initWP%Jtot)
        allocate(int_seqjK(initWP%jmin:IALR%jint,initWP%Kmin:Kmax))
        int_seqjK = -1
        do ijK = 1, int_njK
            j = int_jKPair(ijK,1)
            K = int_jKPair(ijK,2)
            int_seqjK(j,K) = ijK
        end do

!> theta grid and weight
        Kmax = min(IALR%jasy, initWP%Jtot)
!> long range region
        IALR%lr_nA = initWP%j0/initWP%jinc + 1
        allocate(lr_ANode(IALR%lr_nA), lr_AWeight(IALR%lr_nA))
        call setANodeAndWeight(initWP%jpar, IALR%lr_nA, lr_ANode, lr_AWeight)
        allocate(lr_YMat(IALR%lr_nA, lr_njK))
        call setA_DVR2FBR(lr_njK, IALR%lr_nA, lr_jKPair, lr_ANode, lr_AWeight, lr_YMat)
!> asymptotic region
        IALR%asy_nA = IALR%jasy/initWP%jinc + 1
        allocate(asy_ANode(IALR%asy_nA), asy_AWeight(IALR%asy_nA))
        call setANodeAndWeight(initWP%jpar, IALR%asy_nA, asy_ANode, asy_AWeight)
        allocate(asy_YMat(IALR%asy_nA, asy_njK))
        call setA_DVR2FBR(asy_njK, IALR%asy_nA, asy_jKPair, asy_ANode, asy_AWeight, asy_YMat)
!> interaction region
        IALR%int_nA = IALR%jint/initWP%jinc + 1
        allocate(int_ANode(IALR%int_nA), int_AWeight(IALR%int_nA))
        call setANodeAndWeight(initWP%jpar, IALR%int_nA, int_ANode, int_AWeight)
        allocate(int_YMat(IALR%int_nA, int_njK))
        call setA_DVR2FBR(int_njK, IALR%int_nA, int_jKPair, int_ANode, int_AWeight, int_YMat)

        if (IF_inelastic) call PODVR_BCvj()

        open(newunit=EvjUnit, file='BC_roVibEvj0.chk', status='replace', action='write')
        write(EvjUnit,'(a)') '# Ro-vibrational eigenvalues of BC at j0, unit in eV'
        write(EvjUnit,'(a)') '# v, Evj0'
        if (IF_inelastic) then
            do ir = 1, IALR%vasy
                write(EvjUnit,'(i5,2x,f15.9)') ir-1, Evj_BC(ir-1,initWP%j0,initWP%PES0)*au2ev
            end do
        else
            do ir = 1, IALR%vint
                write(EvjUnit,'(i5,2x,f15.9)') ir-1, int_POEig(ir)*au2ev
            end do
        end if
        close(EvjUnit)

        if (IF_inelastic) then
            Ev0j0 = Evj_BC(initWP%v0,initWP%j0,initWP%PES0)
        else
            Ev0j0 = int_POEig(initWP%v0+1)
        end if
        write(outFileUnit,'(1x,a,2i2,a,i2)') 'Initial ro-vibrational energy of (v0, j0) = (', initWP%v0, initWP%j0, ') at iPES = ', initWP%PES0
        write(outFileUnit,'(1x,a,f15.9,a)') 'Evj of BC = ', Ev0j0*au2ev, ' eV.'
        write(outFileUnit,'(1x,a,f15.9,a)') 'Evj of BC = ', Ev0j0*au2cm, ' cm-1.'
        write(outFileUnit,'(1x,a)') 'Please check the ro-vibrational eigenvalues of BC at j0 in BC_roVibEvj0.chk!'
        write(outFileUnit,*) ''

!> Output basis information
        write(outFileUnit,'(1x,a)') '=====> Basis information <====='
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a)') 'IALR basis set, Z in sin-DVR |Z_l>, r in PODVR |r_p>, theta in FBR |JMjK>'
        write(outFileUnit,'(1x,a,i6)') 'Number of Z grid points : ', IALR%nZ_IALR
        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a)') 'Long-range Z grid range: [ ', Z_IALR(IALR%nZ_IA+1), ',', Z_IALR(IALR%nZ_IALR), ' ] a.u.'
        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a)') 'Asymptotic Z grid range: [ ', Z_IALR(IALR%nZ_I+1), ',', Z_IALR(IALR%nZ_IA), ' ] a.u.'
        write(outFileUnit,'(1x,a,f15.9,a,f15.9,a)') 'Interaction Z grid range: [ ', Z_IALR(1), ',', Z_IALR(IALR%nZ_I), ' ] a.u.'
        if (IF_inelastic) then
            write(outFileUnit,'(1x,a,f15.9,a)') 'Inelastic position: ', inePos, ' a.u.'
        end if
        write(outFileUnit,'(1x,a,i6)') 'Number of r PODVR points in interaction region : ', IALR%vint
        write(outFileUnit,'(1x,a,i6)') 'Number of r PODVR points in asymptotic region : ', IALR%vasy
        write(outFileUnit,'(1x,a,i6)') 'Number of r PODVR points in long range region : ', IALR%vlr
        write(outFileUnit,'(1x,a,i6)') 'Number of (j, K) pairs in interaction region : ', int_njK
        write(outFileUnit,'(1x,a,i6)') 'Number of (j, K) pairs in asymptotic region : ', asy_njK
        write(outFileUnit,'(1x,a,i6)') 'Number of (j, K) pairs in long range region : ', lr_njK
        write(outFileUnit,'(1x,a,i6)') 'Number of Gaussian-Legendre quadrature in interaction region : ', IALR%int_nA
        write(outFileUnit,'(1x,a,i6)') 'Number of Gaussian-Legendre quadrature in asymptotic region : ', IALR%asy_nA
        write(outFileUnit,'(1x,a,i6)') 'Number of Gaussian-Legendre quadrature in long range region : ', IALR%lr_nA

        Kmax = min(IALR%jint, initWP%Jtot)
        write(outFileUnit,'(1x,a,3i5)') 'Maximum numbers of (v, j, K) in interaction region = ', IALR%vint, IALR%jint, Kmax
        Kmax = min(IALR%jasy, initWP%Jtot)
        write(outFileUnit,'(1x,a,3i5)') 'Maximum numbers of (v, j, K) in asymptotic region = ', IALR%vasy, IALR%jasy, Kmax
        write(outFileUnit,'(1x,a)') ''

        deallocate(DVRGrid, DVREig, DVRWF, VBC)

    end subroutine setBasis
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine DVR_vib(nDVR, mass, range, Vdiatom, eigVal, eigVec)
!> Calculate the DVR eigenvalues and eigenvectors for a given diatomic potential
!> Ref: J. Chem. Phys. Vol. 96 (3), 1 Feb 1992, pp 1982-1991
        implicit none
        integer, intent(in) :: nDVR
        real(f8), intent(in) :: mass, range
!> Vdiatom: nDVR
        real(f8), intent(in) :: Vdiatom(:)
!> Only calculate vibrational states
        real(f8), intent(inout) :: eigVal(:)
        real(f8), intent(inout) :: eigVec(:,:) 
        integer :: i, j, info, lwork
        real(f8), allocatable :: work(:)
        real(f8) :: fact

        associate(n => nDVR, C => eigVec, V => Vdiatom)
        fact = pi*pi/(4.0_f8*mass*range*range)
!> Since diagonal T matrix is C matrix
            do i = 1, n 
                do j = 1, i-1
                    C(i,j) = fact  * &
                             ((dsin(pi * (i-j) / (2.0_f8 * (n+1))))**(-2) - &
                             (dsin(pi * (i+j) / (2.0_f8 * (n+1))))**(-2)) * (-1.0_f8) **(i-j)
                end do
                C(i,i) = fact * &
                         ((2.0_f8 * (n+1)**2 +1.0_f8) / 3.0_f8 - &
                         (dsin(pi * i / (n+1)))**(-2)) + V(i)
            end do 
            do i = 1, n
                do j = i+1, n
                    C(i,j) = C(j,i)
                end do
            end do
!> Test for the need space of work 
            allocate(work(1))
            call dsyev('V', 'U', n, C, n, eigVal, work, -1, info)
            lwork = int(work(1))
            deallocate(work)

            allocate(work(lwork))
            call dsyev('V', 'U', n, C, n, eigVal, work, lwork, info)
            if (info /= 0) then
                write(outFileUnit,*) "Error in DSYEV of DVR eigen, info =", info
                write(outFileUnit,*) "POSITION: basis.f90, subroutine DVR_vib()"
                stop
            end if
            deallocate(work)


        end associate

    end subroutine DVR_vib
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine PODVR(nPODVR, nDVR, DVRWF, DVRGrid, DVREig, POGrid, POEig, POWF, PO2FBR)
!> Calculate the PODVR basis and eigenvalues/eigenfunctions for given DVR basis
!> See Chemical Physics Letters 1992, 190 (3–4), 225–230.
        implicit none
        integer, intent(in) :: nPODVR, nDVR
        real(f8), intent(in) :: DVRWF(:,:)
        real(f8), intent(in) :: DVRGrid(:)
        real(f8), intent(in) :: DVREig(:)
        real(f8), intent(out) :: POGrid(nPODVR)
        real(f8), intent(out) :: POEig(nPODVR)
        real(f8), intent(out) :: POWF(nPODVR,nPODVR)
        real(f8), intent(out) :: PO2FBR(nPODVR,nPODVR)
        real(f8) :: tmp
        real(f8) :: POTransMat(nDVR,nPODVR), HRefMat(nPODVR,nPODVR), Xmat(nPODVR,nPODVR)
        real(f8), allocatable :: work(:)
        integer :: i, j, l, info, lwork

!> DVR to PODVR transformation matrix
        do i = 1, nDVR
            do j = 1, nPODVR
                POTransMat(i,j) = DVRWF(i,j)
            end do
        end do

        !do i = 1, nPODVR
        !    call phaseTrans(nDVR,POTransMat(:,i))
        !end do

!> X_PO = TransMat * X_DVR * TransMat^T
!> <F|X|F'>
        do i = 1, nPODVR
            do j = i, nPODVR
                tmp = 0.0_f8
                do l = 1, nDVR
                    tmp = tmp + POTransMat(l,i) * DVRGrid(l) * POTransMat(l,j)
                end do
                Xmat(i,j) = tmp
                Xmat(j,i) = tmp
            end do
        end do

        allocate(work(1))
        call dsyev('V', 'U', nPODVR, Xmat, nPODVR, POGrid, work, -1, info)
        lwork = int(work(1))
        deallocate(work)
!> Diagonalize Xmat, Eigenvalues are POGrid, eigenvectors are TransMat between FBR and PODVR, named Xmat
        allocate(work(lwork))
        call dsyev('V', 'U', nPODVR, Xmat, nPODVR, POGrid, work, lwork, info)
        if (info /= 0) then
            write(outFileUnit,*) "Error in DSYEV of PODVR, info =", info
            write(outFileUnit,*) "POSITION: basis.f90, subroutine PODVR()"
            stop
        end if
        deallocate(work)

        do i = 1, nPODVR
            call phaseTrans(nPODVR, Xmat(:,i))
        end do
        !call phase(nPODVR, Xmat)
        PO2FBR = Xmat

!> HRef(i,j) = Sum_l Xmat(l,i) * DVREig(l) * Xmat(l,j)
        HRefMat = 0.0_f8
        do i = 1, nPODVR
            do j = 1, i 
                tmp = 0.0_f8
                do l = 1, nPODVR
                    tmp = tmp + Xmat(l,i) * DVREig(l) * Xmat(l,j)
                end do
                HRefMat(i,j) = tmp
                HRefMat(j,i) = tmp
            end do
        end do

        allocate(work(1))
        call dsyev('V', 'U', nPODVR, HRefMat, nPODVR, POEig, work, -1, info)
        lwork = int(work(1))
        deallocate(work)

        allocate(work(lwork))
        call dsyev('V', 'U', nPODVR, HRefMat, nPODVR, POEig, work, lwork, info)
        if (info /= 0) then
            write(outFileUnit,*) "Error in DSYEV of PODVR HRefMat, info =", info
            write(outFileUnit,*) "POSITION: basis.f90, subroutine PODVR()"
            stop
        end if
        deallocate(work)

        do i = 1, nPODVR
            call phaseTrans(nPODVR, HRefMat(:,i))
        end do
        !call phase(nPODVR, HRefMat) 
        POWF = HRefMat

    end subroutine PODVR
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine phase(n, A)
        implicit none
        integer, intent(in) :: n 
        real(f8), intent(inout) :: A(:,:)
        integer :: i, j, id 
        real(f8) :: tmp

        do i = 1, n 
            tmp = 0.0_f8
            do j = 1, n
                if (abs(A(j,i)) > tmp) then 
                    tmp = abs(A(j,i))
                    id = j
                end if
            end do
            if (A(id,i) < 0.0_f8) then
                do j = 1, n 
                    A(j,i) = -A(j,i)
                end do
            end if
        end do
    end subroutine phase
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine phaseTrans(n, vec)
!> Phase, the sign is choosen such that the first non-zero element is positive 
!> See J. Chem. Phys. 158, 054801 (2023), Sec. II H
        implicit none
        integer, intent(in) :: n
        real(f8), intent(inout) :: vec(n)
        integer :: i
        logical :: found
        found = .false.

        do i = 1, n 
            if (abs(vec(i)) > 1.0e-4_f8) then 
                found = .true.
                exit
            end if
        end do

        if (.not. found) then 
            write(outFileUnit,*) "Error: All vector elements near 0."
            write(outFileUnit,*) "POSITION: basis.f90, subroutine phaseTrans()"
            stop
        end if

        if (vec(i) < 0.0_f8) then
            vec = -vec
        end if

    end subroutine phaseTrans
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine rBasisTrans(type, nZ, nr, njK, BMat, auxWP, WP, WPm)
        implicit none
        character(len=*), intent(in) :: type
        integer, intent(in) :: nZ, nr, njK
        real(f8), intent(in) :: BMat(nr,nr)
        real(f8), intent(inout) :: auxWP(:,:,:,:), WP(:,:,:,:), WPm(:,:,:,:)
        integer :: ijK, iPES

        if (trim(type) == 'DVR2FBR') then 
!$OMP parallel do default(shared) private(ijK,iPES) collapse(2)
            do iPES = 1, nPES
                do ijK = 1, njK
                    call dgemm('N', 'T', nZ, nr, nr, &
                                1.0_f8, WP(1,1,ijK,iPES), nZ, BMat, nr, &
                                0.0_f8, auxWP(1,1,ijK,iPES), nZ)
                    call dcopy(nr*nZ, auxWP(1,1,ijK,iPES), 1, WP(1,1,ijK,iPES), 1)
                    call dgemm('N', 'T', nZ, nr, nr, &
                                1.0_f8, WPm(1,1,ijK,iPES), nZ, BMat, nr, &
                                0.0_f8, auxWP(1,1,ijK,iPES), nZ)
                    call dcopy(nr*nZ, auxWP(1,1,ijK,iPES), 1, WPm(1,1,ijK,iPES), 1)
                end do
            end do
!$OMP end parallel do
        else if (trim(type) == 'FBR2DVR') then
!$OMP parallel do default(shared) private(ijK,iPES) collapse(2)
            do iPES = 1, nPES
                do ijK = 1, njK
                    call dgemm('N', 'N', nZ, nr, nr, &
                                1.0_f8, WP(1,1,ijK,iPES), nZ, BMat, nr, &
                                0.0_f8, auxWP(1,1,ijK,iPES), nZ)
                    call dcopy(nr*nZ, auxWP(1,1,ijK,iPES), 1, WP(1,1,ijK,iPES), 1)
                    call dgemm('N', 'N', nZ, nr, nr, &
                                1.0_f8, WPm(1,1,ijK,iPES), nZ, BMat, nr, &
                                0.0_f8, auxWP(1,1,ijK,iPES), nZ)
                    call dcopy(nr*nZ, auxWP(1,1,ijK,iPES), 1, WPm(1,1,ijK,iPES), 1)
                end do
            end do
!$OMP end parallel do
        else
            write(outFileUnit,*) "Error in rBasisTrans: Unknown type ", trim(type)
            write(outFileUnit,*) "POSITION: propagation.f90, subroutine rBasisTrans()"
            stop
        end if  

    end subroutine rBasisTrans
!> ------------------------------------------------------------------------------------------------------------------ <!
    

end module m_Basis
