module m_Prop
    use m_MachinaBasic, only : f8, c8
    use m_gPara
    use m_InitWP, only : initDFFT
    use m_Potent, only : setPot_IALR
    use m_Basis, only : rBasisTrans
    implicit none

!> Hamiltonian scale
    real(f8), private :: TZMax, TrMax, rotMax, CPMax, VMax 
    real(f8), private :: TZMin, TrMin, rotMin, CPMin, VMin
    integer, private :: nKBlock
!> Auxiliary wave-packet arrays for propagation
    real(f8), allocatable, private :: int_auxWP(:,:,:,:)
    real(f8), allocatable, private :: asy_auxWP(:,:,:,:)
    real(f8), allocatable, private :: lr_auxWP(:,:,:,:)
!> VAction tmp arrays
    real(f8), allocatable, private :: asy_Vtmp(:,:), int_Vtmp(:,:)

    public
    private :: ChebyshevRecursion, setTIDWF, totHamAction
    private :: setKinMat, setRotMat, setCPMat
    private :: kinAction, rotAction, CPAction, VpotAction
    private :: lambdaPlus, lambdaMinus

contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine propProcess()
        implicit none
        real(f8) :: lr_norm, asy_norm, int_norm, tot_norm
        real(f8) :: cpu_start, cpu_now, cpu_dt
        integer :: int_njMax, asy_njMax
        integer :: iStep, iPrint
        integer :: iZ, ir, ith

        int_njMax = max(int_njK, IALR%int_nA)
        asy_njMax = max(asy_njK, IALR%asy_nA)
        allocate(int_auxWP(IALR%nZ_I, IALR%vint, int_njMax, nPES))
        allocate(asy_auxWP(IALR%nZ_IA, IALR%vasy, asy_njMax, nPES))
        allocate(lr_auxWP(IALR%nZ_IALR, IALR%vlr, lr_njK, 1))
        allocate(int_Vtmp(nPES,IALR%int_nA))
        allocate(asy_Vtmp(nPES,IALR%asy_nA))
        if (IF_inelastic) then
            allocate(ine_TIDWF(nEtot,IALR%vasy,IALR%asy_nA,initWP%Kmin:min(initWP%Jtot,IALR%jasy),nPES))
            ine_TIDWF = imgZero
        end if

!> k = 0, |\phi>_0 = |initWP>
        lr_TDWPm = lr_TDWP; lr_TDWP = 0.0_f8
        asy_TDWPm = asy_TDWP; asy_TDWP = 0.0_f8
        int_TDWPm = int_TDWP; int_TDWP = 0.0_f8

!> Propagation loop
        iPrint = 0
        iStep = 0
        write(outFileUnit,'(1x,a)') '=====> Chebyshev propagation <====='
        write(outFileUnit,'(1x,a)') ''
!> Print title header once
        write(outFileUnit,'(1x,a8,4a14,a12)') 'Step', 'intNorm', 'asyNorm', 'lrNorm', 'totNorm', 'CPU(s)'
        write(outFileUnit,'(1x,a)') repeat('-', 76)
        flush(outFileUnit)
!> Print initial norm (Step 0)
        lr_norm = sum(abs(lr_TDWPm(IALR%nZ_IA+1:IALR%nZ_IALR,:,:,1))**2)
        asy_norm = sum(abs(asy_TDWPm(IALR%nZ_I+1:IALR%nZ_IA,:,:,:))**2)
        int_norm = sum(abs(int_TDWPm(1:IALR%nZ_I,:,:,:))**2)
        tot_norm = lr_norm + asy_norm + int_norm
        write(outFileUnit,'(1x,i8,4f14.5,a12)') iStep, int_norm, asy_norm, lr_norm, tot_norm, '--'

        call cpu_time(cpu_start)
        flush(outFileUnit)
        do iStep = 1, timeTot
            call ChebyshevRecursion(iStep)
            if (IF_inelastic) call ineSetTIDWF(iStep)
            if (nChannel >= 1) call setTIDWF(iStep,'A+BC->C+AB',channel1%ichoice)
            if (nChannel == 2) then
                call setTIDWF(iStep,'A+BC->B+AC',channel2%ichoice)
            end if
            iPrint = iPrint + 1
            if (iPrint == timePrint) then
                call cpu_time(cpu_now)
                cpu_dt = cpu_now - cpu_start
                lr_norm = sum(abs(lr_TDWPm(IALR%nZ_IA+1:IALR%nZ_IALR,:,:,1))**2)
                asy_norm = sum(abs(asy_TDWPm(IALR%nZ_I+1:IALR%nZ_IA,:,:,:))**2)
                int_norm = sum(abs(int_TDWPm(1:IALR%nZ_I,:,:,:))**2)
                tot_norm = lr_norm + asy_norm + int_norm
                write(outFileUnit,'(1x,i8,4f14.5,f12.3)') iStep, int_norm, asy_norm, lr_norm, tot_norm, cpu_dt
                flush(outFileUnit)
                cpu_start = cpu_now
                iPrint = 0
            end if
        end do
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a)') 'Propagation done!'
        write(outFileUnit,'(1x,a)') ''
        flush(outFileUnit)

    end subroutine propProcess
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine ChebyshevRecursion(iStep)
        implicit none
        integer, intent(in) :: iStep
        real(f8) :: WFtmp, fact, VabsTmp 
        integer :: iPES, ir, iZ, ijK, j 
        logical :: is_init_channel

!> Hamiltonian action
        call totHamAction()

!> Recursion relation
        if (iStep == 1) then
            fact = 1.0_f8
        else
            fact = 0.5_f8
        end if
!> Long range region
!$OMP parallel do default(shared) private(ijK,ir,iZ,WFtmp)
        do ijK = 1, lr_njK
            do ir = 1, IALR%vlr
                do iZ = IALR%nZ_IA+1, IALR%nZ_IALR
                    WFtmp = fact*lr_TDWPm(iZ,ir,ijK,1)
                    lr_TDWPm(iZ,ir,ijK,1) = 2.0_f8*lr_TDWP(iZ,ir,ijK,1)*lr_ZFabs(iZ)
                    lr_TDWP(iZ,ir,ijK,1) = -WFtmp*lr_ZFabs(iZ)
                end do
            end do
        end do
!$OMP end parallel do

!> Asymptotic region
        do iPES = 1, nPES
!$OMP parallel do default(shared) private(ijK,ir,iZ,j,VabsTmp,WFtmp,is_init_channel)
            do ijK = 1, asy_njK
                j = asy_jKPair(ijK,1)
                is_init_channel = (iPES == initWP%PES0 .and. j == initWP%j0)
                do ir = 1, IALR%vasy
                    do iZ = IALR%nZ_I+1, IALR%nZ_IA
                        if (is_init_channel .and. ir == initWP%v0+1) then
                            VabsTmp = 1.0_f8
                        else
                            VabsTmp = asy_ZFabs(iZ)
                        end if
                        WFtmp = fact*asy_TDWPm(iZ,ir,ijK,iPES)
                        asy_TDWPm(iZ,ir,ijK,iPES) = 2.0_f8*asy_TDWP(iZ,ir,ijK,iPES)*VabsTmp
                        asy_TDWP(iZ,ir,ijK,iPES) = -WFtmp*VabsTmp
                    end do
                end do
            end do
!$OMP end parallel do
        end do
!> Interaction region
        do iPES = 1, nPES
!$OMP parallel do default(shared) private(ijK,ir,iZ,j,VabsTmp,WFtmp,is_init_channel)
            do ijK = 1, int_njK
                do ir = 1, IALR%vint
                    do iZ = 1, IALR%nZ_I
                        WFtmp = fact*int_TDWPm(iZ,ir,ijK,iPES)
                        int_TDWPm(iZ,ir,ijK,iPES) = 2.0_f8*int_TDWP(iZ,ir,ijK,iPES)
                        int_TDWP(iZ,ir,ijK,iPES) = -WFtmp
                    end do
                end do
            end do
!$OMP end parallel do
        end do
    end subroutine ChebyshevRecursion
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setTIDWF(iStep,type,ichoice)
        implicit none
        integer, intent(in) :: iStep
        character(len=*), intent(in) :: type
        integer, intent(in) :: ichoice
        complex(c8) :: FourierFactor(nEtot)
        real(f8) :: an
        integer :: iZ, ir, ith, iPES, iEtot, iIC, nIC
        integer :: K, Kmax
        integer :: ntk, nti, j
        real(f8), external :: ddot
        external :: zaxpy

        FourierFactor(:) = cdexp(cmplx(0.0_f8,-iStep*ChebyAngle(:), kind=c8))

        if (type == 'A+BC->C+AB') then
            nIC = channel1%nIC
        else if (type == 'A+BC->B+AC') then
            nIC = channel2%nIC
        end if

        Kmax = min(IALR%jasy, initWP%Jtot)

        nti = 0
        do K = initWP%Kmin, Kmax
            !> Count number of (j,K) pairs with K=ka
            ntk = 0
            do j = initWP%jmin, IALR%jint, initWP%jinc
                if (j >= K) ntk = ntk + 1
            end do

            do iPES = 1, nPES
                !> Perform transformation: FBR -> DVR for angle
                call dgemm('N','T', IALR%vint*IALR%nZ_I, IALR%int_nA, ntk, &
                            1.0_f8, int_TDWPm(1,1,nti+1,iPES), IALR%vint*IALR%nZ_I, &
                            int_YMat(1,nti+1), IALR%int_nA, &
                            0.0_f8, int_auxWP(1,1,1,iPES), IALR%vint*IALR%nZ_I)
            end do

            if (ichoice == 1) then
!$OMP parallel do default(shared) private(iPES,iIC,iZ,ith,an) collapse(2)
                do iIC = 1, nIC
                    do iPES = 1, nPES
                        if (type == 'A+BC->C+AB') then
                            iZ = idXY_AB(iIC,1)
                            ith = idXY_AB(iIC,2)
                            an = ddot(IALR%vint, Uij_AB(1,iIC), 1, int_auxWP(iZ,1,ith,iPES), IALR%nZ_I)
                            call zaxpy(nEtot, cmplx(an,0.0_f8,kind=c8), FourierFactor, 1, &
                                        intAB_TIDWF(1,iIC,K,iPES), 1)
                        else if (type == 'A+BC->B+AC') then
                            iZ = idXY_AC(iIC,1)
                            ith = idXY_AC(iIC,2)
                            an = ddot(IALR%vint, Uij_AC(1,iIC), 1, int_auxWP(iZ,1,ith,iPES), IALR%nZ_I)
                            call zaxpy(nEtot, cmplx(an,0.0_f8,kind=c8), FourierFactor, 1, &
                                        intAC_TIDWF(1,iIC,K,iPES), 1)
                        end if
                    end do
                end do
!$OMP end parallel do
            else if (ichoice == 2) then
!$OMP parallel do default(shared) private(iPES,iIC,ir,ith,an) collapse(2)
                do iIC = 1, nIC
                    do iPES = 1, nPES
                        if (type == 'A+BC->C+AB') then
                            ir = idXY_AB(iIC,1)
                            ith = idXY_AB(iIC,2)
                            an = ddot(IALR%nZ_I, Uij_AB(1,iIC), 1, int_auxWP(1,ir,ith,iPES), 1)
                            call zaxpy(nEtot, cmplx(an,0.0_f8,kind=c8), FourierFactor, 1, &
                                        intAB_TIDWF(1,iIC,K,iPES), 1)
                        else if (type == 'A+BC->B+AC') then
                            ir = idXY_AC(iIC,1)
                            ith = idXY_AC(iIC,2)
                            an = ddot(IALR%nZ_I, Uij_AC(1,iIC), 1, int_auxWP(1,ir,ith,iPES), 1)
                            call zaxpy(nEtot, cmplx(an,0.0_f8,kind=c8), FourierFactor, 1, &
                                        intAC_TIDWF(1,iIC,K,iPES), 1)
                        end if
                    end do
                end do
!$OMP end parallel do
            end if
            nti = nti + ntk
        end do

    end subroutine setTIDWF
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine ineSetTIDWF(iStep)
!> Accumulate inelastic TIDWF at Z = Z_IALR(iInePos) only
        implicit none
        integer, intent(in) :: iStep
        integer :: nti, K, Kmax, ntk, j
        integer :: ir, ijK, iPES, iEtot, ith
        complex(c8) :: FourierFactor(nEtot)
        real(f8) :: rthTmp(IALR%vasy, IALR%asy_nA)
        real(f8) :: rthSlice(IALR%vasy, IALR%asy_nA)
        real(f8) :: invSqrtDZ

        FourierFactor(:) = cdexp(cmplx(0.0_f8,-iStep*ChebyAngle(:), kind=c8))
        invSqrtDZ = 1.0_f8 / dsqrt(Z_IALR(2) - Z_IALR(1))

        Kmax = min(IALR%jasy, initWP%Jtot)
        nti = 0
        do K = initWP%Kmin, Kmax
            ntk = 0
            do j = initWP%jmin, IALR%jasy, initWP%jinc
                if (j >= K) ntk = ntk + 1
            end do

            do iPES = 1, nPES
                !> FBR -> DVR for angle (r remains in PODVR-FBR)
                call dgemm('N','T', IALR%vasy*IALR%nZ_IA, IALR%asy_nA, ntk, &
                            1.0_f8, asy_TDWPm(1,1,nti+1,iPES), IALR%vasy*IALR%nZ_IA, &
                            asy_YMat(1,nti+1), IALR%asy_nA, &
                            0.0_f8, asy_auxWP(1,1,1,iPES), IALR%vasy*IALR%nZ_IA)
            end do

!> Extract Z-slice at iInePos and transform r from PODVR-FBR to PODVR-DVR
            do iPES = 1, nPES
                do ith = 1, IALR%asy_nA
                    do ir = 1, IALR%vasy
                        rthTmp(ir, ith) = asy_auxWP(iInePos, ir, ith, iPES) * invSqrtDZ
                    end do
                end do
                call dgemm('T', 'N', IALR%vasy, IALR%asy_nA, IALR%vasy, &
                            1.0_f8, asy_PO2FBR, IALR%vasy, &
                            rthTmp, IALR%vasy, &
                            0.0_f8, rthSlice, IALR%vasy)
!$OMP parallel do default(shared) private(ith,ir) collapse(2)
                do ith = 1, IALR%asy_nA
                    do ir = 1, IALR%vasy
                        call zaxpy(nEtot, cmplx(rthSlice(ir,ith),0.0_f8,kind=c8), FourierFactor, 1, &
                                    ine_TIDWF(1,ir,ith,K,iPES), 1)
                    end do
                end do
!$OMP end parallel do
            end do
            nti = nti + ntk
        end do

    end subroutine ineSetTIDWF
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine totHamAction()
        implicit none
        integer :: iPES, ir, ijKi, iZ


!>  r in FBR and Z in DVR now

        call kinAction()
        call CPAction()

!> Trans r to DVR
        call rBasisTrans('FBR2DVR',IALR%nZ_I,IALR%vint,int_njK, &
                        int_PO2FBR,int_auxWP,int_TDWP,int_TDWPm)
        call rBasisTrans('FBR2DVR',IALR%nZ_IA,IALR%vasy,asy_njK, &
                        asy_PO2FBR,asy_auxWP,asy_TDWP,asy_TDWPm)
        call rBasisTrans('FBR2DVR',IALR%nZ_IALR,IALR%vlr,lr_njK, &
                        lr_PO2FBR,lr_auxWP,lr_TDWP,lr_TDWPm)

        call VpotAction()
        call rotAction()

!> Dumping for interaction region
!$OMP parallel do default(shared) private(iPES,ir,ijKi,iZ)
        do iPES = 1, nPES
            do ijKi = 1, int_njK
                do ir = 1, IALR%vint
                    do iZ = 1, IALR%nZ_I
                        int_TDWP(iZ,ir,ijKi,iPES) = int_TDWP(iZ,ir,ijKi,iPES)*int_rFabs(ir)
                        int_TDWPm(iZ,ir,ijKi,iPES) = int_TDWPm(iZ,ir,ijKi,iPES)*int_rFabs(ir)
                    end do
                end do
            end do
        end do
!$OMP end parallel do

!> Trans r back to FBR
        call rBasisTrans('DVR2FBR',IALR%nZ_I,IALR%vint,int_njK, &
                        int_PO2FBR,int_auxWP,int_TDWP,int_TDWPm)
        call rBasisTrans('DVR2FBR',IALR%nZ_IA,IALR%vasy,asy_njK, &
                        asy_PO2FBR,asy_auxWP,asy_TDWP,asy_TDWPm)
        call rBasisTrans('DVR2FBR',IALR%nZ_IALR,IALR%vlr,lr_njK, &
                        lr_PO2FBR,lr_auxWP,lr_TDWP,lr_TDWPm)
    end subroutine totHamAction
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine kinAction()
        implicit none
        integer :: iZ, ir, idjKl, idjKa, idjKi  
        integer :: iPES, K, j, Kmax
        
        Kmax = min(IALR%jasy, initWP%Jtot)

!> r in FBR, Z in DVR now
        do iPES = 1, nPES
            call dcopy(int_njK*IALR%vint*IALR%nZ_I,int_TDWPm(1,1,1,iPES),1,int_auxWP(1,1,1,iPES),1)
            call dcopy(asy_njK*IALR%vasy*IALR%nZ_IA,asy_TDWPm(1,1,1,iPES),1,asy_auxWP(1,1,1,iPES),1)
        end do 
        call dcopy(lr_njK*IALR%vlr*IALR%nZ_IALR,lr_TDWPm(1,1,1,1),1,lr_auxWP(1,1,1,1),1)

!> Copy int_auxWP to asy_auxWP
        do iPES = 1, nPES 
            do K = initWP%Kmin, Kmax
                do j = initWP%jmin, IALR%jasy, initWP%jinc 
                    if (j >= K) then 
                        idjKi = int_seqjK(j,K)
                        idjKa = asy_seqjK(j,K)
                        do ir = 1, IALR%vasy
                            do iZ = 1, IALR%nZ_I
                                asy_auxWP(iZ,ir,idjKa,iPES) = int_auxWP(iZ,ir,idjKi,iPES)
                            end do
                        end do
                    end if
                end do
            end do
        end do

!> Copy asy_auxWP to lr_auxWP
        do K = initWP%Kmin, Kmax 
            if (initWP%j0 >= K) then 
                idjKa  = asy_seqjK(initWP%j0,K)
                idjKl = lr_seqjK(initWP%j0,K)
                do ir = 1, IALR%vlr
                    do iZ = 1, IALR%nZ_IA
                        lr_auxWP(iZ,ir,idjKl,1) = asy_auxWP(iZ,ir,idjKa,initWP%PES0)
                    end do
                end do
            end if
        end do

!> Kinectic action
        !> long range
        do idjkl = 1, lr_njK
            do ir = 1, IALR%vlr
                call dsint(IALR%nZ_IALR, lr_auxWP(1,ir,idjkl,1), lr_WSaveZ)
                do iZ = 1, IALR%nZ_IALR
                    lr_auxWP(iZ,ir,idjkl,1) = lr_kinEigen(iZ,ir)*lr_auxWP(iZ,ir,idjkl,1)
                end do
                call dsint(IALR%nZ_IALR, lr_auxWP(1,ir,idjkl,1), lr_WSaveZ)
            end do
        end do

        !> asymptotic region
        do iPES = 1, nPES
!$OMP parallel do default(shared) private(idjKa,ir,iZ)
            do idjKa = 1, asy_njK
                do ir  = 1, IALR%vasy
                    call dsint(IALR%nZ_IA, asy_auxWP(1,ir,idjKa,iPES), asy_WSaveZ(1,idjKa))
                    do iZ = 1, IALR%nZ_IA
                        asy_auxWP(iZ,ir,idjKa,iPES) = asy_kinEigen(iZ,ir) * asy_auxWP(iZ,ir,idjKa,iPES)
                    end do
                    call dsint(IALR%nZ_IA, asy_auxWP(1,ir,idjKa,iPES), asy_WSaveZ(1,idjKa))
                end do
            end do
!$OMP end parallel do
        end do

        !> interaction region
        do iPES = 1, nPES
!$OMP parallel do default(shared) private(idjKi,ir,iZ)
            do idjKi = 1, int_njK
                do ir  = 1, IALR%vint
                    call dsint(IALR%nZ_I, int_auxWP(1,ir,idjKi,iPES), int_WSaveZ(1,idjKi))
                    do iZ = 1, IALR%nZ_I
                        int_auxWP(iZ,ir,idjKi,iPES) = int_kinEigen(iZ,ir) * int_auxWP(iZ,ir,idjKi,iPES)
                    end do
                    call dsint(IALR%nZ_I, int_auxWP(1,ir,idjKi,iPES), int_WSaveZ(1,idjKi))
                end do
            end do
!$OMP end parallel do
        end do

!> Copy lr_auxWP to asy_auxWP
        do K = initWP%Kmin, Kmax 
            if (initWP%j0 >= K) then 
                idjKa  = asy_seqjK(initWP%j0,K)
                idjKl = lr_seqjK(initWP%j0,K)
                do ir = 1, IALR%vlr
                    do iZ = 1, IALR%nZ_IA
                        asy_auxWP(iZ,ir,idjKa,initWP%PES0) = lr_auxWP(iZ,ir,idjKl,1)
                    end do
                end do
            end if
        end do

!> Copy asy_auxWP to int_auxWP
        do iPES = 1, nPES 
            do K = initWP%Kmin, Kmax
                do j = initWP%jmin, IALR%jasy, initWP%jinc 
                    if (j >= K) then 
                        idjKi = int_seqjK(j,K)
                        idjKa = asy_seqjK(j,K)
                        do ir = 1, IALR%vasy
                            do iZ = 1, IALR%nZ_I
                                int_auxWP(iZ,ir,idjKi,iPES) = asy_auxWP(iZ,ir,idjKa,iPES)
                            end do
                        end do
                    end if
                end do
            end do
        end do

!> Add auxWP on int and asyWP
        do iPES = 1, nPES
            call daxpy(int_njK*IALR%vint*IALR%nZ_I,1.0_f8,int_auxWP(1,1,1,iPES),1,int_TDWP(1,1,1,iPES),1)
            call daxpy(asy_njK*IALR%vasy*IALR%nZ_IA,1.0_f8,asy_auxWP(1,1,1,iPES),1,asy_TDWP(1,1,1,iPES),1)
        end do
        call daxpy(lr_njK*IALR%vlr*IALR%nZ_IALR,1.0_f8,lr_auxWP(1,1,1,1),1,lr_TDWP(1,1,1,1),1)

    end subroutine kinAction
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine CPAction()
        implicit none
        integer :: iPES, idjK, jdjK, idjKa, jdjKa
        integer :: j, jstart, K, Kp, Kmax, iZ, ir

!> r in FBR, Z in DVR
        if (nKBlock == 1) then 
            !> For nK == 1
            !> interaction region
            do iPES = 1, nPES 
!$OMP parallel do default(shared) private(idjK,ir,iZ)
                do idjK = 1, int_njK
                    do ir  = 1, IALR%vint
                        do iZ = 1, IALR%nZ_I
                            int_TDWP(iZ,ir,idjK,iPES) = int_TDWP(iZ,ir,idjK,iPES)+int_CPK1(iZ,idjK)*int_TDWPm(iZ,ir,idjK,iPES)
                        end do
                    end do
                end do
!$OMP end parallel do
            end do
            !> asymptotic region
            do iPES = 1, nPES
!$OMP parallel do default(shared) private(idjKa,ir,iZ)
                do idjK = 1, asy_njK
                    do ir  = 1, IALR%vasy
                        do iZ = IALR%nZ_I+1, IALR%nZ_IA
                            asy_TDWP(iZ,ir,idjK,iPES) = asy_TDWP(iZ,ir,idjK,iPES)+asy_CPK1(iZ,idjK)*asy_TDWPm(iZ,ir,idjK,iPES)
                        end do
                    end do
                end do
!$OMP end parallel do
            end do
            !> long range region
!$OMP parallel do default(shared) private(idjK,ir,iZ)
            do idjK = 1, lr_njK
                do ir = 1, IALR%vlr
                    do iZ = IALR%nZ_IA+1, IALR%nZ_IALR
                        lr_TDWP(iZ,ir,idjK,1) = lr_TDWP(iZ,ir,idjK,1)+lr_CPK1(iZ,idjK)*lr_TDWPm(iZ,ir,idjK,1)
                    end do
                end do 
            end do
!$OMP end parallel do
        else
            !> For nK > 1 - full K-coupling using IA_CPMat and lr_CPMat
            !> Use int_seqjK to directly index all K values for same j,
            !> since jKPair is ordered K-first (K=0 block, then K=1 block, ...),
            !> so same-j entries with different K are NOT adjacent.
            !> interaction region
            Kmax = min(IALR%jint, initWP%Jtot)
            do iPES = 1, nPES
!$OMP parallel do default(shared) private(idjK,jdjK,K,Kp,j,ir,iZ)
                do idjK = 1, int_njK
                    j = int_jKPair(idjK,1)
                    K = int_jKPair(idjK,2)
                    do ir  = 1, IALR%vint
                        do iZ = 1, IALR%nZ_I
                            !> Sum over all K' for same j using int_seqjK
                            do Kp = initWP%Kmin, min(Kmax, j)
                                jdjK = int_seqjK(j,Kp)
                                int_TDWP(iZ,ir,idjK,iPES) = int_TDWP(iZ,ir,idjK,iPES) + &
                                    IA_CPMat(iZ,idjK,jdjK) * int_TDWPm(iZ,ir,jdjK,iPES)
                            end do
                        end do
                    end do
                end do
!$OMP end parallel do
            end do

            !> asymptotic region
            Kmax = min(IALR%jasy, initWP%Jtot)
            do iPES = 1, nPES
!$OMP parallel do default(shared) private(idjKa,jdjKa,idjK,jdjK,K,Kp,j,ir,iZ)
                do idjKa = 1, asy_njK
                    j = asy_jKPair(idjKa,1)
                    K = asy_jKPair(idjKa,2)
                    !> Map to int index for IA_CPMat
                    idjK = int_seqjK(j,K)
                    do ir  = 1, IALR%vasy
                        do iZ = IALR%nZ_I+1, IALR%nZ_IA
                            !> Sum over all K' for same j using asy_seqjK
                            do Kp = initWP%Kmin, min(Kmax, j)
                                jdjKa = asy_seqjK(j,Kp)
                                jdjK = int_seqjK(j,Kp)
                                asy_TDWP(iZ,ir,idjKa,iPES) = asy_TDWP(iZ,ir,idjKa,iPES) + &
                                    IA_CPMat(iZ,idjK,jdjK) * asy_TDWPm(iZ,ir,jdjKa,iPES)
                            end do
                        end do
                    end do
                end do
!$OMP end parallel do
            end do

            !> long range region (only iZ > nZ_IA)
            Kmax = min(initWP%j0, initWP%Jtot)
            do Kp = initWP%Kmin, Kmax
                do K = initWP%Kmin, Kmax
                    idjK = lr_seqjK(initWP%j0,K)
                    jdjK = lr_seqjK(initWP%j0,Kp)
                    do ir = 1, IALR%vlr
                        do iZ = IALR%nZ_IA+1, IALR%nZ_IALR
                            lr_TDWP(iZ,ir,idjK,1) = lr_TDWP(iZ,ir,idjK,1) + &
                                                    lr_CPMat(iZ,idjK,jdjK) * &
                                                    lr_TDWPm(iZ,ir,jdjK,1)
                        end do
                    end do
                end do
            end do
        end if

    end subroutine CPAction
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine rotAction()
        implicit none
        integer :: iPES, iZ, ir, idjK 

!> r in DVR, Z in DVR 
        !> interaction region
        do iPES = 1, nPES
!$OMP parallel do default(shared) private(idjK,ir,iZ)
            do idjK = 1, int_njK
                do ir  = 1, IALR%vint
!$OMP SIMD
                    do iZ = 1, IALR%nZ_I
                        int_TDWP(iZ,ir,idjK,iPES) = int_rotMat(ir,idjK)*int_TDWPm(iZ,ir,idjK,iPES)+int_TDWP(iZ,ir,idjK,iPES)
                    end do
                end do
            end do
!$OMP end parallel do
        end do

        !> asymptotic region
        do iPES = 1, nPES
!$OMP parallel do default(shared) private(idjK,ir,iZ)
            do idjK = 1, asy_njK
                do ir  = 1, IALR%vasy
!$OMP SIMD
                    do iZ = IALR%nZ_I+1, IALR%nZ_IA
                        asy_TDWP(iZ,ir,idjK,iPES) = asy_rotMat(ir,idjK)*asy_TDWPm(iZ,ir,idjK,iPES)+asy_TDWP(iZ,ir,idjK,iPES)
                    end do
                end do
            end do
!$OMP end parallel do
        end do

        !> long range region
!$OMP parallel do default(shared) private(idjK,ir,iZ)
        do idjK = 1, lr_njK
            do ir = 1, IALR%vlr
!$OMP SIMD
                do iZ = IALR%nZ_IA+1, IALR%nZ_IALR
                    lr_TDWP(iZ,ir,idjK,1) = lr_rotMat(ir,idjK)*lr_TDWPm(iZ,ir,idjK,1)+lr_TDWP(iZ,ir,idjK,1)
                end do
            end do 
        end do
!$OMP end parallel do
    end subroutine rotAction
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine VpotAction()
        implicit none
        integer :: iPES, jPES, iZ, ir, ith
        integer :: K, Kmax_int, Kmax_asy, Kmax_lr
        integer :: nti_int, ntk_int, nti_asy, ntk_asy, nti_lr, ntk_lr
        integer :: j
        real(f8) :: vtmp

!> r in DVR, Z in DVR
!> Loop over K blocks for angular FBR-DVR transformation
!> Different K blocks use different associated Legendre functions P_j^K,
!> which are NOT orthogonal to each other under the same Gauss-Legendre quadrature.
!> Therefore, the FBR-DVR transformation must be done per K block.
        Kmax_int = min(IALR%jint, initWP%Jtot)
        Kmax_asy = min(IALR%jasy, initWP%Jtot)
        Kmax_lr  = min(initWP%j0, initWP%Jtot)

        nti_int = 0
        nti_asy = 0
        nti_lr  = 0

        do K = initWP%Kmin, Kmax_int
            !> Count basis functions for this K block in each region
            ntk_int = 0
            do j = initWP%jmin, IALR%jint, initWP%jinc
                if (j >= K) ntk_int = ntk_int + 1
            end do

            ntk_asy = 0
            if (K <= Kmax_asy) then
                do j = initWP%jmin, IALR%jasy, initWP%jinc
                    if (j >= K) ntk_asy = ntk_asy + 1
                end do
            end if

            ntk_lr = 0
            if (K <= Kmax_lr .and. initWP%j0 >= K) ntk_lr = 1

            !> ---- Interaction region ----
            !> Forward transform: FBR to angular DVR for this K block
            do iPES = 1, nPES
                call dgemm('N','T',IALR%nZ_I*IALR%vint,&
                            IALR%int_nA, ntk_int,&
                            1.0_f8,int_TDWPm(1,1,nti_int+1,iPES),IALR%nZ_I*IALR%vint, &
                            int_YMat(1,nti_int+1),IALR%int_nA,&
                            0.0_f8,int_auxWP(1,1,1,iPES),IALR%nZ_I*IALR%vint)
            end do
            !> Apply potential in DVR
            if (nPES == 1) then
!$OMP parallel  do default(shared) private(ith,ir,iZ)
                do ith = 1, IALR%int_nA
                    do ir = 1, IALR%vint
                        do iZ = 1, IALR%nZ_I
                            int_auxWP(iZ,ir,ith,1) = int_auxWP(iZ,ir,ith,1)*int_Vdia(1,1,iZ,ir,ith)
                        end do
                    end do
                end do
!$OMP end parallel do
            else
!$OMP parallel  do default(shared) private(ith,ir,iZ,iPES,jPES,vtmp)
                do ith = 1, IALR%int_nA
                    do ir  = 1, IALR%vint
                        do iZ  = 1, IALR%nZ_I
                            do iPES = 1, nPES
                                vtmp = 0.0_f8
                                do jPES = 1, nPES
                                    vtmp = vtmp + int_Vdia(jPES,iPES,iZ,ir,ith)*int_auxWP(iZ,ir,ith,jPES)
                                end do
                                int_Vtmp(iPES,ith) = vtmp
                            end do
                            do iPES = 1, nPES
                                int_auxWP(iZ,ir,ith,iPES) = int_Vtmp(iPES,ith)
                            end do
                        end do
                    end do
                end do
!$OMP end parallel do
            end if
            !> Back transform: angular DVR to FBR, accumulate into TDWP
            do iPES = 1, nPES
                call dgemm('N','N', IALR%nZ_I*IALR%vint, &
                ntk_int, IALR%int_nA, &
                1.0_f8, int_auxWP(1,1,1,iPES), IALR%nZ_I * IALR%vint, &
                int_YMat(1,nti_int+1), IALR%int_nA, &
                1.0_f8, int_TDWP(1,1,nti_int+1,iPES), IALR%nZ_I * IALR%vint )
            end do

            !> ---- Asymptotic region ----
            if (ntk_asy > 0) then
                !> Forward transform
                do iPES = 1, nPES
                    call dgemm('N','T',IALR%nZ_IA*IALR%vasy,&
                                IALR%asy_nA, ntk_asy,&
                                1.0_f8,asy_TDWPm(1,1,nti_asy+1,iPES),IALR%nZ_IA*IALR%vasy, &
                                asy_YMat(1,nti_asy+1),IALR%asy_nA,&
                                0.0_f8,asy_auxWP(1,1,1,iPES),IALR%nZ_IA*IALR%vasy)
                end do
                !> Apply potential
                if (nPES == 1) then
!$OMP parallel  do default(shared) private(ith,ir,iZ)
                    do ith = 1, IALR%asy_nA
                        do ir = 1, IALR%vasy
                            do iZ = IALR%nZ_I+1, IALR%nZ_IA
                                asy_auxWP(iZ,ir,ith,1) = asy_auxWP(iZ,ir,ith,1)*asy_Vdia(1,1,iZ-IALR%nZ_I,ir,ith)
                            end do
                        end do
                    end do
!$OMP end parallel do
                else
!$OMP parallel  do default(shared) private(ith,ir,iZ,iPES,jPES,vtmp)
                    do ith = 1, IALR%asy_nA
                        do ir  = 1, IALR%vasy
                            do iZ  = IALR%nZ_I+1, IALR%nZ_IA
                                do iPES = 1, nPES
                                    vtmp = 0.0_f8
                                    do jPES = 1, nPES
                                        vtmp = vtmp + asy_Vdia(jPES,iPES,iZ-IALR%nZ_I,ir,ith)*asy_auxWP(iZ,ir,ith,jPES)
                                    end do
                                    asy_Vtmp(iPES,ith) = vtmp
                                end do
                                do iPES = 1, nPES
                                    asy_auxWP(iZ,ir,ith,iPES) = asy_Vtmp(iPES,ith)
                                end do
                            end do
                        end do
                    end do
!$OMP end parallel do
                end if
                !> Back transform
                do iPES = 1, nPES
                    call dgemm('N','N', IALR%nZ_IA*IALR%vasy, &
                    ntk_asy, IALR%asy_nA, &
                    1.0_f8, asy_auxWP(1,1,1,iPES), IALR%nZ_IA * IALR%vasy, &
                    asy_YMat(1,nti_asy+1), IALR%asy_nA, &
                    1.0_f8, asy_TDWP(1,1,nti_asy+1,iPES), IALR%nZ_IA * IALR%vasy )
                end do
            end if

            !> ---- Long-range region ----
            if (ntk_lr > 0) then
                !> Forward transform
                call dgemm('N','T',IALR%nZ_IALR*IALR%vlr,&
                            IALR%lr_nA, ntk_lr,&
                            1.0_f8,lr_TDWPm(1,1,nti_lr+1,1),IALR%nZ_IALR*IALR%vlr, &
                            lr_YMat(1,nti_lr+1),IALR%lr_nA,&
                            0.0_f8,lr_auxWP(1,1,1,1),IALR%nZ_IALR*IALR%vlr)
                !> Apply potential
                do ith = 1, IALR%lr_nA
                    do ir = 1, IALR%vlr
                        do iZ  = IALR%nZ_IA+1, IALR%nZ_IALR
                            lr_auxWP(iZ,ir,ith,1) = lr_Vdia(iZ-IALR%nZ_IA,ir,ith)*lr_auxWP(iZ,ir,ith,1)
                        end do
                    end do
                end do
                !> Back transform
                call dgemm('N','N', IALR%nZ_IALR*IALR%vlr, &
                ntk_lr, IALR%lr_nA, &
                1.0_f8, lr_auxWP(1,1,1,1), IALR%nZ_IALR * IALR%vlr, &
                lr_YMat(1,nti_lr+1), IALR%lr_nA, &
                1.0_f8, lr_TDWP(1,1,nti_lr+1,1), IALR%nZ_IALR * IALR%vlr )
            end if

            nti_int = nti_int + ntk_int
            nti_asy = nti_asy + ntk_asy
            nti_lr  = nti_lr  + ntk_lr
        end do

    end subroutine VpotAction
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine HamiltonianScale()
        implicit none
        real(f8) :: Hmax, Hmin
        real(f8) :: Ev0j0
        real(f8) :: dtMax, dtMin, dtAvg
        integer :: iEtot, iPES, ith, iZ, ir

        nKBlock = min(initWP%Jtot, IALR%jasy) - initWP%Kmin + 1

        call setKinMat()
        call setRotMat()
        call setCPMat()

        call setPot_IALR()

        VMax = max(maxval(int_Vdia), maxval(asy_Vdia), maxval(lr_Vdia))
        VMin = min(minval(int_Vdia), minval(asy_Vdia), minval(lr_Vdia))

        Hmax = TZMax + TrMax + rotMax + CPMax + VMax + 0.1_f8 * ev2au
        Hmin = TZMin + TrMin + rotMin + CPMin + VMin - 0.1_f8 * ev2au

        Hplus = (Hmax + Hmin) / 2.0_f8; Hminus = (Hmax - Hmin) / 2.0_f8

!> Scale Hamiltonian in [-1,1]
        !> kinZ 
        int_ZkinMat = int_ZkinMat / Hminus; asy_ZkinMat = asy_ZkinMat / Hminus; lr_ZkinMat = lr_ZkinMat / Hminus
        !> kinr
        int_rKinMat = int_rKinMat / Hminus; asy_rKinMat = asy_rKinMat / Hminus; lr_rKinMat = lr_rKinMat / Hminus
        !> rot
        int_rotMat = int_rotMat / Hminus; asy_rotMat = asy_rotMat / Hminus; lr_rotMat = lr_rotMat / Hminus
        !> CP
        if (nKBlock == 1) then 
            int_CPK1 = int_CPK1 / Hminus; asy_CPK1 = asy_CPK1 / Hminus; lr_CPK1 = lr_CPK1 / Hminus
        else
            IA_CPMat = IA_CPMat / Hminus; lr_CPMat = lr_CPMat / Hminus
        end if
        !> Vdia
        do ith = 1, IALR%int_nA
            do ir  = 1, IALR%vint
                do iZ  = 1, IALR%nZ_I
                    do iPES = 1, nPES
                        int_Vdia(iPES,iPES,iZ,ir,ith) = int_Vdia(iPES,iPES,iZ,ir,ith) - Hplus
                    end do
                end do
            end do
        end do
        int_Vdia = int_Vdia / Hminus

        do ith = 1, IALR%asy_nA
            do ir  = 1, IALR%vasy
                do iZ  = 1, IALR%nZ_IA-IALR%nZ_I
                    do iPES = 1, nPES
                        asy_Vdia(iPES,iPES,iZ,ir,ith) = asy_Vdia(iPES,iPES,iZ,ir,ith) - Hplus
                    end do
                end do
            end do
        end do
        asy_Vdia = asy_Vdia / Hminus
        lr_Vdia = (lr_Vdia - Hplus) / Hminus

!> Chebyshev angles
        if (IF_inelastic) then
            Ev0j0 = Evj_BC(initWP%v0,initWP%j0,initWP%PES0)
        else
            Ev0j0 = int_POEig(initWP%v0+1)
        end if
        do iEtot = 1, nEtot
            Etot(iEtot) = Ecol(iEtot) + Ev0j0
            ChebyAngle(iEtot) = dacos((Etot(iEtot) - Hplus) / Hminus)
        end do

        dtMin = 1.0_f8/(Hminus*sin(maxval(ChebyAngle)))
        dtMax = 1.0_f8/(Hminus*sin(minval(ChebyAngle)))
        dtAvg = (dtMax+dtMin)/2.0_f8
        timeStep = dtAvg

        call initDFFT()

        write(outFileUnit,'(1x,a)') '=====> Chebyshev scale information <====='
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a,f12.5,2x,a,f12.5)') 'Hmax(eV) = ', Hmax*au2ev, 'Hmin = ', Hmin*au2ev
        write(outFileUnit,'(1x,a,f12.5,2x,a,f12.5)') 'H+(eV) = ', Hplus*au2ev, 'H-(eV) = ', Hminus*au2ev
        write(outFileUnit,'(1x,a,f12.5)') 'Z Kinetic Max(eV) = ', TZMax*au2ev
        write(outFileUnit,'(1x,a,f12.5,a,f12.5)') 'r Kinetic Max(eV) = ', TrMax*au2ev, ', r Kinetic Min(eV) = ', TrMin*au2ev
        write(outFileUnit,'(1x,a,f12.5)') 'Rotational Max(eV) = ', rotMax*au2ev
        write(outFileUnit,'(1x,a,f12.5)') 'CPE Max(eV) = ', CPMax*au2ev
        write(outFileUnit,'(1x,a,f12.5,a,f12.5)') 'Vpot Max(eV) = ', VMax*au2ev, ', Vpot Min(eV) = ', VMin*au2ev
        write(outFileUnit,'(1x,a,f12.5)') 'TimeStep Max = ', dtMax
        write(outFileUnit,'(1x,a,f12.5)') 'TimeStep Min = ', dtMin 
        write(outFileUnit,'(1x,a,f12.5)') 'TimeStep Avg = ', dtAvg
        write(outFileUnit,*) ''

    end subroutine HamiltonianScale
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setKinMat()
        implicit none
        integer :: iZ, ir 

!> Kinetic energy along Z
        allocate(int_ZkinMat(IALR%nZ_I))
        allocate(asy_ZkinMat(IALR%nZ_IA))
        allocate(lr_ZkinMat(IALR%nZ_IALR))

        do iZ = 1, IALR%nZ_IALR
            lr_ZkinMat(iZ) = (iZ*pi/(IALR%Z_range(2)-IALR%Z_range(1)))**2/2.0_f8/massTot
        end do
        do iZ = 1, IALR%nZ_IA
            asy_ZkinMat(iZ) = (iZ*pi/(Z_IALR(IALR%nZ_IA+1)-IALR%Z_range(1)))**2/2.0_f8/massTot
        end do
        do iZ = 1, IALR%nZ_I
            int_ZkinMat(iZ) = (iZ*pi/(Z_IALR(IALR%nZ_I+1)-IALR%Z_range(1)))**2/2.0_f8/massTot
        end do

        TZMin = 0.0_f8
        TZMax = (pi/(Z_IALR(2)-Z_IALR(1)))**2/2.0_f8/massTot 

!> Kinetic energy along r
        allocate(int_rKinMat(IALR%vint))
        allocate(asy_rKinMat(IALR%vasy))
        allocate(lr_rKinMat(IALR%vlr))

        int_rKinMat(:) = int_POEig(1:IALR%vint)
        asy_rKinMat(:) = int_POEig(1:IALR%vasy)
        lr_rKinMat(:) = int_POEig(1:IALR%vlr)

        TrMin = int_POEig(1)
        TrMax = int_POEig(IALR%vint)

    end subroutine setKinMat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setRotMat()
        implicit none
        integer :: ijK, ir, qn_j 
        real(f8) :: jEigen 

        allocate(int_rotMat(IALR%vint,int_njK))
        allocate(asy_rotMat(IALR%vasy,asy_njK))
        allocate(lr_rotMat(IALR%vlr,lr_njK))

!> long range 
        do ijK = 1, lr_njK
            do ir = 1, IALR%vlr 
                qn_j = lr_jKPair(ijK,1)
                jEigen = (qn_j*(qn_j+1.0_f8)) / (2.0_f8*massBC*r_LR(ir)**2)
                lr_rotMat(ir,ijK) = min(jEigen, TMaxCut)
            end do
        end do
!> asymptotic region
        do ijK = 1, asy_njK
            do ir  = 1, IALR%vasy
                qn_j = asy_jKPair(ijK,1)
                jEigen = (qn_j*(qn_j+1.0_f8)) / (2.0_f8*massBC*r_Asy(ir)**2)
                asy_rotMat(ir,ijK) = min(jEigen, TMaxCut)
            end do
        end do
!> interaction region
        do ijK = 1, int_njK
            do ir  = 1, IALR%vint
                qn_j = int_jKPair(ijK,1)
                jEigen = (qn_j*(qn_j+1.0_f8)) / (2.0_f8*massBC*r_Int(ir)**2)
                int_rotMat(ir,ijK) = min(jEigen, TMaxCut)
            end do
        end do

    rotMin = 0.0_f8
    rotMax = maxval(int_rotMat)

    end subroutine setRotMat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setCPMat()
        implicit none
        integer :: iZ, K, Kp, qn_j, qn_jp, idjK, jdjK, Kmax
        integer :: nlo, nhi, ndo 
        integer :: idxj, idxi, lwork, info 
        real(f8) :: term1, term2, term3, delta, tmpE
        real(f8) :: Zfact(IALR%nZ_IALR), phase, lreal, leigen
        real(f8), allocatable :: CPM(:,:), UMat(:,:,:)
        real(f8), allocatable :: Bvj(:,:), CPE(:), work(:)

        allocate(IA_CPMat(IALR%nZ_IA,int_njK,int_njK))
        allocate(lr_CPMat(IALR%nZ_IALR,lr_njK,lr_njK))
        allocate(CPM(int_njK,int_njK))
        allocate(UMat(IALR%nZ_IALR,int_njK,int_njK))
        allocate(int_lo(int_njK))

        CPM = 0.0_f8

        do iZ = 1, IALR%nZ_IALR
            Zfact(iZ) = 1.0_f8/(2.0_f8*massTot*Z_IALR(iZ)**2)
        end do

        if (nKBlock == 1) then 
            allocate(int_CPK1(IALR%nZ_I,int_njK), asy_CPK1(IALR%nZ_IA,asy_njK), lr_CPK1(IALR%nZ_IALR,lr_njK))
            !> Long range region
            do idjK = 1, lr_njK
                do iZ = 1, IALR%nZ_IALR
                    qn_j = lr_jKPair(idjK,1)
                    leigen = initWP%Jtot*(initWP%Jtot+1.0_f8) + qn_j*(qn_j+1.0_f8) &
                            - 2.0_f8*initWP%Kmin*initWP%Kmin
                    lr_CPK1(iZ,idjK) = min(Zfact(iZ)*leigen, TMaxCut)  
                end do
            end do
            !> Asymptotic region
            do idjK = 1, asy_njK
                do iZ = 1, IALR%nZ_IA
                    qn_j = asy_jKPair(idjK,1)
                    leigen = initWP%Jtot*(initWP%Jtot+1.0_f8) + qn_j*(qn_j+1.0_f8) &
                            - 2.0_f8*initWP%Kmin*initWP%Kmin
                    asy_CPK1(iZ,idjK) = min(Zfact(iZ)*leigen, TMaxCut)  
                end do
            end do
            !> Interaction region
            do idjK = 1, int_njK    
                do iZ = 1, IALR%nZ_I
                    qn_j = int_jKPair(idjK,1)
                    leigen = initWP%Jtot*(initWP%Jtot+1.0_f8) + qn_j*(qn_j+1.0_f8) &
                            - 2.0_f8*initWP%Kmin*initWP%Kmin
                    int_CPK1(iZ,idjK) = min(Zfact(iZ)*leigen, TMaxCut)  
                end do
            end do
            CPMin = 0.0_f8
            CPMax = maxval(int_CPK1)

        else
!> Get CPM in whole IALR region — built per j block using int_seqjK for correct indexing
!> Reference: cold_dia_cheby/cheby.f90 set_cpm, sub.f90 cp0
            CPM = 0.0_f8
            Kmax = min(initWP%Jtot, IALR%jint)
            do qn_j = initWP%jmin, IALR%jint, initWP%jinc
                do K = initWP%Kmin, min(Kmax, qn_j)
                    idjK = int_seqjK(qn_j, K)
                    do Kp = initWP%Kmin, min(Kmax, qn_j)
                        jdjK = int_seqjK(qn_j, Kp)

                        term1 = merge((initWP%Jtot*(initWP%Jtot+1.0_f8)+qn_j*(qn_j+1.0_f8)-2.0_f8*K*K), &
                                    0.0_f8, K == Kp)

                        delta = merge(dsqrt(2.0_f8), 1.0_f8, K == 0)
                        term2 = merge(delta*lambdaPlus(initWP%Jtot,K)*lambdaPlus(qn_j,K), &
                                    0.0_f8, Kp == K+1)
                        delta = merge(dsqrt(2.0_f8), 1.0_f8, K == 1)
                        term3 = merge(delta*lambdaMinus(initWP%Jtot,K)*lambdaMinus(qn_j,K), &
                                    0.0_f8, Kp == K-1)
                        CPM(idjK, jdjK) = term1 - term2 - term3
                    end do
                end do
            end do

!> Diagonalize CPM per j block, apply Z-dependent cutoff, reconstruct
!> Use int_seqjK to gather the (j,K) indices for each j block
            CPMax = -100.0_f8
            CPMin = 0.0_f8
            UMat = 0.0_f8

            do qn_j = initWP%jmin, IALR%jint, initWP%jinc
                !> Count K values for this j
                ndo = 0
                do K = initWP%Kmin, min(Kmax, qn_j)
                    ndo = ndo + 1
                end do
                if (ndo == 0) cycle

                allocate(Bvj(ndo,ndo), CPE(ndo))
                !> Gather the K-block sub-matrix from CPM
                idxi = 0
                do K = initWP%Kmin, min(Kmax, qn_j)
                    idxi = idxi + 1
                    idxj = 0
                    idjK = int_seqjK(qn_j, K)
                    do Kp = initWP%Kmin, min(Kmax, qn_j)
                        idxj = idxj + 1
                        jdjK = int_seqjK(qn_j, Kp)
                        Bvj(idxi, idxj) = CPM(idjK, jdjK)
                    end do
                end do

                !> Diagonalize
                allocate(work(1))
                call dsyev('V', 'U', ndo, Bvj, ndo, CPE, work, -1, info)
                lwork = int(work(1))
                deallocate(work)
                allocate(work(lwork))
                call dsyev('V', 'U', ndo, Bvj, ndo, CPE, work, lwork, info)
                if (info /= 0) then
                    write(outFileUnit,*) "Error in DSYEV of int Bvj, info =", info, " j=", qn_j
                    write(outFileUnit,*) "POSITION: propagation.f90, subroutine setCPMat()"
                    stop
                end if
                deallocate(work)

                !> Phase convention: (-1)^(j+Kmax_j)*eigvec(ndo,:) positive
                nhi = min(initWP%Jtot, qn_j)
                phase = (-1.0_f8)**(qn_j+nhi)
                do idxj = 1, ndo
                    if (Bvj(ndo,idxj)*phase < 0.0_f8) then
                        Bvj(:,idxj) = -Bvj(:,idxj)
                    end if
                end do

                !> Reconstruct UMat with Z-dependent eigenvalues and cutoff
                !> UMat(iZ, i, j) = sum_k R(i,k) * min(l_k*(l_k+1)/2μZ², TMaxCut) * R(j,k)
                do iZ = 1, IALR%nZ_IALR
                    idxi = 0
                    do K = initWP%Kmin, min(Kmax, qn_j)
                        idxi = idxi + 1
                        idjK = int_seqjK(qn_j, K)
                        idxj = 0
                        do Kp = initWP%Kmin, min(Kmax, qn_j)
                            idxj = idxj + 1
                            jdjK = int_seqjK(qn_j, Kp)
                            delta = 0.0_f8
                            do nhi = 1, ndo
                                lreal = dabs(dsqrt(CPE(nhi)+0.25_f8)-0.5_f8)
                                leigen = lreal*(lreal+1.0_f8)
                                tmpE = min(leigen*Zfact(iZ), TMaxCut)
                                CPMax = max(CPMax, tmpE)
                                delta = delta + Bvj(idxi,nhi)*tmpE*Bvj(idxj,nhi)
                            end do
                            UMat(iZ, idjK, jdjK) = delta
                        end do
                    end do
                end do

                deallocate(Bvj, CPE)
            end do

!> Rebuild CPM in IA
            IA_CPMat(1:IALR%nZ_IA,:,:) = UMat(1:IALR%nZ_IA,:,:)
!> Rebuild CPM in LR (UMat uses int_seqjK indexing, lr_CPMat uses lr_seqjK)
            Kmax = min(initWP%Jtot, initWP%j0)
            do K = initWP%Kmin, Kmax
                do Kp = initWP%Kmin, Kmax
                    lr_CPMat(:,lr_seqjK(initWP%j0,K),lr_seqjK(initWP%j0,Kp)) = &
                        UMat(:,int_seqjK(initWP%j0,K),int_seqjK(initWP%j0,Kp))
                end do
            end do

            deallocate(CPM,UMat)
        end if

    end subroutine setCPMat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    real(f8) function lambdaPlus(J,K) result(res)
        implicit none
        integer, intent(in) :: J, K

        res = dsqrt(J*(J+1.0_f8) - K*(K+1.0_f8))
        return

    end function lambdaPlus
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    real(f8) function lambdaMinus(J,K) result(res)
        implicit none
        integer, intent(in) :: J, K

        res = dsqrt(J*(J+1.0_f8) - K*(K-1.0_f8))
        return

    end function lambdaMinus
!> ------------------------------------------------------------------------------------------------------------------ <!

end module m_Prop
