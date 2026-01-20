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
!> Auxiliary wave-packet arrays for propagation
    real(f8), allocatable, private :: int_auxWP(:,:,:,:)
    real(f8), allocatable, private :: asy_auxWP(:,:,:,:)
    real(f8), allocatable, private :: lr_auxWP(:,:)

    public
    private :: setKinMat, setRotMat, setCPMat
    private :: kinAction, rotAction, CPAction, VpotAction
    private :: lambdaPlus, lambdaMinus

contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine propProcess()
        implicit none
        real(f8) :: lr_norm, asy_norm, int_norm, tot_norm
        integer :: int_njMax, asy_njMax 
        integer :: iStep, iPrint

        int_njMax = max(int_njK, IALR%int_nA)
        asy_njMax = max(asy_njK, IALR%asy_nA)
        allocate(int_auxWP(IALR%nZ_I, IALR%vint, int_njMax, nPES))
        allocate(asy_auxWP(IALR%nZ_IA, IALR%vasy, asy_njMax, nPES))
        allocate(lr_auxWP(IALR%nZ_IALR, lr_njK))

!> k = 0, |\phi>_0 = |initWP>
        lr_TDWPm = lr_TDWP; lr_TDWP = 0.0_f8
        asy_TDWPm = asy_TDWP; asy_TDWP = 0.0_f8
        int_TDWPm = int_TDWP; int_TDWP = 0.0_f8

!> Propagation loop
        iPrint = 0
        iStep = 0
        write(outFileUnit,'(1x,a)') '=====> Chebyshev propagation <====='
        write(outFileUnit,'(1x,a)') ''
        lr_norm = sum(abs(lr_TDWPm(IALR%nZ_IA+1:IALR%nZ_IALR,:))**2)
        asy_norm = sum(abs(asy_TDWPm(IALR%nZ_I+1:IALR%nZ_IA,:,:,:))**2)
        int_norm = sum(abs(int_TDWPm(1:IALR%nZ_I,:,:,:))**2)
        tot_norm = lr_norm + asy_norm + int_norm
        write(outFileUnit,'(1x,a,i8,a,f12.5,a,f12.5,a,f12.5,a,f12.5,a,f12.5)') 'Step = ',iStep, &
                                                                        ' , intNorm = ', int_norm, &
                                                                        ' , asyNorm = ', asy_norm, &
                                                                        ' , lrNorm = ', lr_norm, &
                                                                        ' , totNorm = ', tot_norm
        do iStep = 1, timeTot
            call ChebyshevRecursion(iStep)
            iPrint = iPrint + 1
            if (iPrint == timePrint) then 
                lr_norm = sum(abs(lr_TDWPm(IALR%nZ_IA+1:IALR%nZ_IALR,:))**2)
                asy_norm = sum(abs(asy_TDWPm(IALR%nZ_I+1:IALR%nZ_IA,:,:,:))**2)
                int_norm = sum(abs(int_TDWPm(1:IALR%nZ_I,:,:,:))**2)
                tot_norm = lr_norm + asy_norm + int_norm
                write(outFileUnit,'(1x,a,i8,a,f12.5,a,f12.5,a,f12.5,a,f12.5,a,f12.5)') 'Step = ',iStep, &
                                                                                ' , intNorm = ', int_norm, &
                                                                                ' , asyNorm = ', asy_norm, &
                                                                                ' , lrNorm = ', lr_norm, &
                                                                                ' , totNorm = ', tot_norm
                iPrint = 0
            end if
        end do
        write(outFileUnit,'(1x,a)') 'Done!'

    end subroutine propProcess
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine ChebyshevRecursion(iStep)
        implicit none
        integer, intent(in) :: iStep
        real(f8) :: WFtmp, fact, VabsTmp 
        integer :: iPES, ir, iZ, ijK, j 

!> Hamiltonian action
        call totHamAction()

!> Recursion relation
        if (iStep == 1) then
            fact = 1.0_f8
        else
            fact = 0.5_f8
        end if
!> Long range region
        do ijK = 1, lr_njK
            do iZ = IALR%nZ_IA+1, IALR%nZ_IALR
                WFtmp = fact*lr_TDWPm(iZ,ijK)
                lr_TDWPm(iZ,ijK) = 2.0_f8*lr_TDWP(iZ,ijK)*lr_ZFabs(iZ)
                lr_TDWP(iZ,ijK) = -WFtmp*lr_ZFabs(iZ)
            end do
        end do
!> Asymptotic region
        do iPES = 1, nPES
!$OMP parallel do default(shared) private(ijK,ir,iZ,j,VabsTmp,WFtmp)
            do ijK = 1, asy_njK
                do ir = 1, IALR%vasy
                    do iZ = IALR%nZ_I+1, IALR%nZ_IA
                        j = asy_jKPair(ijK,1)
                        if (iPES == initWP%PES0 .and. j == initWP%j0 .and. ir == initWP%v0+1) then
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
!$OMP parallel do default(shared) private(ijK,ir,iZ,j,VabsTmp,WFtmp)
            do ijK = 1, int_njK
                do ir = 1, IALR%vint
                    do iZ = 1, IALR%nZ_I
                        j = int_jKPair(ijK,1)
                        if (iPES == initWP%PES0 .and. j == initWP%j0 .and. ir == initWP%v0+1) then
                            VabsTmp = 1.0_f8
                        else if (r_Int(ir) <= r_Asy(IALR%vasy)) then
                            VabsTmp = 1.0_f8
                        else
                            VabsTmp = int_ZFabs(iZ)
                        end if
                        WFtmp = fact*int_TDWPm(iZ,ir,ijK,iPES)
                        int_TDWPm(iZ,ir,ijK,iPES) = 2.0_f8*int_TDWP(iZ,ir,ijK,iPES)*VabsTmp
                        int_TDWP(iZ,ir,ijK,iPES) = -WFtmp*VabsTmp
                    end do
                end do
            end do
!$OMP end parallel do
        end do
        
    end subroutine ChebyshevRecursion
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine totHamAction()
        implicit none
        integer :: iPES, ir
        
!>  r in FBR and Z in DVR now 
        call kinAction()
        call CPAction()

!> Trans r to DVR 
        call rBasisTrans('FBR2DVR',IALR%nZ_I,IALR%vint,int_njK, &
                        int_PO2FBR,int_auxWP,int_TDWP,int_TDWPm)
        call rBasisTrans('FBR2DVR',IALR%nZ_IA,IALR%vasy,asy_njK, &
                        asy_PO2FBR,asy_auxWP,asy_TDWP,asy_TDWPm)

        call VpotAction()
        call rotAction()

!> Dumping for interaction region
!$OMP parallel do default(shared) private(iPES,ir)
        do iPES = 1, nPES
            do ir = 1, IALR%vint
                int_TDWP(:,ir,:,iPES) = int_TDWP(:,ir,:,iPES)*int_rFabs(ir)
                int_TDWPm(:,ir,:,iPES) = int_TDWPm(:,ir,:,iPES)*int_rFabs(ir)
            end do
        end do
!$OMP end parallel do

!> Trans r back to FBR
        call rBasisTrans('DVR2FBR',IALR%nZ_I,IALR%vint,int_njK, &
                        int_PO2FBR,int_auxWP,int_TDWP,int_TDWPm)
        call rBasisTrans('DVR2FBR',IALR%nZ_IA,IALR%vasy,asy_njK, &
                        asy_PO2FBR,asy_auxWP,asy_TDWP,asy_TDWPm)

    end subroutine totHamAction
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine kinAction()
        implicit none
        integer :: iZ, ir, idjKl, idjKa, idjKi  
        integer :: iPES, K, j, Kmax

!> r in FBR, Z in DVR now
        do iPES = 1, nPES
            call dcopy(int_njK*IALR%vint*IALR%nZ_I,int_TDWPm(1,1,1,iPES),1,int_auxWP(1,1,1,iPES),1)
            call dcopy(asy_njK*IALR%vasy*IALR%nZ_IA,asy_TDWPm(1,1,1,iPES),1,asy_auxWP(1,1,1,iPES),1)
        end do 
        call dcopy(lr_njK*IALR%nZ_IALR,lr_TDWPm,1,lr_auxWP,1)

!> Copy int_auxWP to asy_auxWP
        do iPES = 1, nPES 
            Kmax = min(IALR%jasy, initWP%Jtot)
            do K = initWP%Kmin, Kmax 
                do j = initWP%jmin, IALR%jasy, initWP%jinc 
                    if (j >= K) then 
                        idjKi = int_seqjK(j,K)
                        idjKa = asy_seqjK(j,K)
                        asy_auxWP(1:IALR%nZ_I,1:IALR%vasy,idjKa,iPES) = int_auxWP(1:IALR%nZ_I,1:IALR%vasy,idjKi,iPES)
                    end if
                end do
            end do
        end do

!> Copy asy_auxWP to lr_auxWP
        Kmax = min(initWP%j0, initWP%Jtot)
        do K = initWP%Kmin, Kmax 
            if (initWP%j0 >= K) then 
                idjKa  = asy_seqjK(initWP%j0,K)
                idjKl = lr_seqjK(initWP%j0,K)
                lr_auxWP(1:IALR%nZ_IA,idjKl) = asy_auxWP(1:IALR%nZ_IA,initWP%v0+1,idjKa,initWP%PES0)
            end if
        end do

!> Kinectic action
        !> long range
        do idjkl = 1, lr_njK
            call dsint(IALR%nZ_IALR, lr_auxWP(1,idjkl), lr_WSaveZ)
            do iZ = 1, IALR%nZ_IALR
                lr_auxWP(iZ,idjkl) = lr_kinEigen(iZ) * lr_auxWP(iZ,idjkl)
            end do
            call dsint(IALR%nZ_IALR, lr_auxWP(1,idjkl), lr_WSaveZ)
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
        Kmax = min(initWP%j0, initWP%Jtot)
        do K = initWP%Kmin, Kmax 
            if (initWP%j0 >= K) then 
                idjKa  = asy_seqjK(initWP%j0,K)
                idjKl = lr_seqjK(initWP%j0,K)
                asy_auxWP(1:IALR%nZ_IA,initWP%v0+1,idjKa,initWP%PES0) = lr_auxWP(1:IALR%nZ_IA,idjKl)
            end if
        end do
!> Copy asy_auxWP to int_auxWP
        do iPES = 1, nPES 
            Kmax = min(IALR%jasy, initWP%Jtot)
            do K = initWP%Kmin, Kmax 
                do j = initWP%jmin, IALR%jasy, initWP%jinc 
                    if (j >= K) then 
                        idjKi = int_seqjK(j,K)
                        idjKa = asy_seqjK(j,K)
                        int_auxWP(1:IALR%nZ_I,1:IALR%vasy,idjKi,iPES) = asy_auxWP(1:IALR%nZ_I,1:IALR%vasy,idjKa,iPES)
                    end if
                end do
            end do
        end do

!> Add auxWP on int and asyWP
        do iPES = 1, nPES
            call daxpy(int_njK*IALR%vint*IALR%nZ_I,1.0_f8,int_auxWP(1,1,1,iPES),1,int_TDWP(1,1,1,iPES),1)
            call daxpy(asy_njK*IALR%vasy*IALR%nZ_IA,1.0_f8,asy_auxWP(1,1,1,iPES),1,asy_TDWP(1,1,1,iPES),1)
        end do
        call daxpy(lr_njK*IALR%nZ_IALR,1.0_f8,lr_auxWP,1,lr_TDWP,1)

    end subroutine kinAction
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine CPAction()
        implicit none
        integer :: iPES, idjK, jdjK, idjKa, jdjKa
        integer :: j, jstart, K, Kp, Kmax, iZ

!> r in FBR, Z in DVR
        !> interaction region and asymptotic region
        do iPES = 1, nPES
!$OMP parallel do default(shared) private(iZ,idjK,jdjK,j,K,Kp,idjKa,jdjKa)
            do idjK = 1, int_njK
                j = int_jKPair(idjK,1)
                K = int_jKPair(idjK,2)
                do jdjK = 1, int_njK
                    if (int_jKPair(jdjK,1) /= j) cycle
                    Kp = int_jKPair(jdjK,2)
                    do iZ = 1, IALR%nZ_I
                        int_TDWP(iZ,:,idjK,iPES) = int_TDWP(iZ,:,idjK,iPES) + &
                                                    IA_CPMat(iZ,idjK,jdjK) * &
                                                    int_TDWPm(iZ,:,jdjK,iPES)
                    end do
                    if (j <= IALR%jasy) then 
                        idjKa = asy_seqjK(j,K)
                        jdjKa = asy_seqjK(j,Kp)
                        do iZ = IALR%nZ_I+1, IALR%nZ_IA
                            asy_TDWP(iZ,:,idjKa,iPES) = asy_TDWP(iZ,:,idjKa,iPES) + &
                                                        IA_CPMat(iZ,idjKa,jdjKa) * &
                                                        asy_TDWPm(iZ,:,jdjKa,iPES)
                        end do
                    end if
                end do 
            end do
!$OMP end parallel do
        end do 

        !> long range region
        Kmax = min(initWP%j0, initWP%Jtot)
        do Kp = initWP%Kmin, Kmax
            do K = initWP%Kmin, Kmax
                idjK = lr_seqjK(initWP%j0,K)
                jdjK = lr_seqjK(initWP%j0,Kp)
                lr_TDWP(1:IALR%nZ_IALR,idjK) = lr_TDWP(1:IALR%nZ_IALR,idjK) + &
                                                lr_CPMat(1:IALR%nZ_IALR,idjK,jdjK) * &
                                                lr_TDWPm(1:IALR%nZ_IALR,jdjK)
            end do
        end do

    end subroutine CPAction
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine rotAction()
        implicit none
        integer :: iPES, iZ, ir, idjK 

!> r in DVR, Z in DVR 
        !> interaction region
        do iPES = 1, nPES
!$OMP parallel do default(shared) private(idjK,ir)
            do idjK = 1, int_njK
                do ir  = 1, IALR%vint
                    int_TDWP(:,ir,idjK,iPES) = int_rotMat(ir,idjK)*int_TDWPm(:,ir,idjK,iPES)+int_TDWP(:,ir,idjK,iPES)
                end do
            end do
!$OMP end parallel do
        end do

        !> asymptotic region
        do iPES = 1, nPES
!$OMP parallel do default(shared) private(idjK,ir,iZ)
            do idjK = 1, asy_njK
                do ir  = 1, IALR%vasy
                    do iZ = IALR%nZ_I+1, IALR%nZ_IA
                        asy_TDWP(iZ,ir,idjK,iPES) = asy_rotMat(ir,idjK)*asy_TDWPm(iZ,ir,idjK,iPES)+asy_TDWP(iZ,ir,idjK,iPES)
                    end do
                end do
            end do
!$OMP end parallel do
        end do

        !> long range region
        do idjK = 1, lr_njK
            do iZ = IALR%nZ_IA+1, IALR%nZ_IALR
                lr_TDWP(iZ,idjK) = lr_rotMat(idjK)*lr_TDWPm(iZ,idjK)+lr_TDWP(iZ,idjK)
            end do
        end do

    end subroutine rotAction
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine VpotAction()
        implicit none
        integer :: iPES, jPES, iZ, ir, ith
        real(f8) :: vtmp
        real(f8), allocatable :: WFtmp(:,:)
        
!> r in DVR, Z in DVR
!> Trans angle FBR to DVR
        do iPES = 1, nPES
            call dgemm('N','T',IALR%nZ_I*IALR%vint,&
                        IALR%int_nA,int_njK,&
                        1.0_f8,int_TDWPm(1,1,1,iPES),IALR%nZ_I*IALR%vint, &
                        int_YMat(1,1),IALR%int_nA,&
                        0.0_f8,int_auxWP(1,1,1,iPES),IALR%nZ_I*IALR%vint)
            call dgemm('N','T',IALR%nZ_IA*IALR%vasy,&
                        IALR%asy_nA,asy_njK,&
                        1.0_f8,asy_TDWPm(1,1,1,iPES),IALR%nZ_IA*IALR%vasy, &
                        asy_YMat(1,1),IALR%asy_nA,&
                        0.0_f8,asy_auxWP(1,1,1,iPES),IALR%nZ_IA*IALR%vasy)
        end do

        !> interaction region
        allocate(WFtmp(nPES,IALR%int_nA))
!$OMP parallel  do default(shared) private(ith,ir,iZ,iPES,jPES,vtmp)
        do ith = 1, IALR%int_nA
            do ir  = 1, IALR%vint
                do iZ  = 1, IALR%nZ_I
                    do iPES = 1, nPES
                        vtmp = 0.0_f8
                        do jPES = 1, nPES
                            vtmp = vtmp + int_Vdia(jPES,iPES,iZ,ir,ith)*int_auxWP(iZ,ir,ith,jPES)
                        end do
                        WFtmp(iPES,ith) = vtmp
                    end do 
                    do iPES = 1, nPES
                        int_auxWP(iZ,ir,ith,iPES) = WFtmp(iPES,ith)
                    end do
                end do
            end do
        end do
!$OMP end parallel do
!> Trans back to FBR and add to WP
        do iPES = 1, nPES
            call dgemm('N','N', IALR%nZ_I*IALR%vint, & 
            int_njK, IALR%int_nA, & 
            1.0_f8, int_auxWP(1,1,1,iPES), IALR%nZ_I * IALR%vint, &
            int_YMat(1,1), IALR%int_nA, &
            1.0_f8, int_TDWP(1,1,1,iPES), IALR%nZ_I * IALR%vint )
        end do
        deallocate(WFtmp)

        !> asymptotic region
        allocate(WFtmp(nPES,IALR%asy_nA))
!$OMP parallel  do default(shared) private(ith,ir,iZ,iPES,jPES,vtmp)
        do ith = 1, IALR%asy_nA
            do ir  = 1, IALR%vasy
                do iZ  = IALR%nZ_I+1, IALR%nZ_IA
                    do iPES = 1, nPES
                        vtmp = 0.0_f8
                        do jPES = 1, nPES
                            vtmp = vtmp + asy_Vdia(jPES,iPES,iZ-IALR%nZ_I,ir,ith)*asy_auxWP(iZ,ir,ith,jPES)
                        end do
                        WFtmp(iPES,ith) = vtmp
                    end do
                    do iPES = 1, nPES
                        asy_auxWP(iZ,ir,ith,iPES) = WFtmp(iPES,ith)
                    end do
                end do
            end do
        end do
!$OMP end parallel do
!> Trans back to FBR and add to WP
        do iPES = 1, nPES
            call dgemm('N','N', IALR%nZ_IA*IALR%vasy, & 
            asy_njK, IALR%asy_nA, & 
            1.0_f8, asy_auxWP(1,1,1,iPES), IALR%nZ_IA * IALR%vasy, &
            asy_YMat(1,1), IALR%asy_nA, &
            1.0_f8, asy_TDWP(1,1,1,iPES), IALR%nZ_IA * IALR%vasy )
        end do
        deallocate(WFtmp)

        !> long range region
        do ith = 1, lr_njK
            do iZ  = IALR%nZ_IA+1, IALR%nZ_IALR
                lr_TDWP(iZ,ith) = lr_Vdia(iZ-IALR%nZ_IA)*lr_TDWPm(iZ,ith)+lr_TDWP(iZ,ith)
            end do
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
        IA_CPMat = IA_CPMat / Hminus; lr_CPMat = lr_CPMat / Hminus
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
        Ev0j0 = int_POEig(initWP%v0+1)
        do iEtot = 1, nEtot
            Etot(iEtot) = Ecol(iEtot) + Ev0j0
            ChebyAngle(iEtot) = acos((Etot(iEtot) - Hminus) / Hminus)
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
            asy_ZkinMat(iZ) = (iZ*pi/(Z_IALR(IALR%nZ_IA+1-range(1))))**2/2.0_f8/massTot
        end do
        do iZ = 1, IALR%nZ_I
            int_ZkinMat(iZ) = (iZ*pi/(Z_IALR(IALR%nZ_I+1-range(1))))**2/2.0_f8/massTot
        end do

        TZMin = 0.0_f8
        TZMax = (pi/(Z_IALR(2)-Z_IALR(1)))**2/2.0_f8/massTot 

!> Kinetic energy along r
        allocate(int_rKinMat(IALR%vint))
        allocate(asy_rKinMat(IALR%vasy))

        int_rKinMat(:) = int_POEig(1:IALR%vint)
        asy_rKinMat(:) = int_POEig(1:IALR%vasy)
        lr_rKinMat = int_POEig(initWP%v0+1)

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
        allocate(lr_rotMat(lr_njK))

!> long range 
        do ijK = 1, lr_njK
            qn_j = lr_jKPair(ijK,1)
            jEigen = (qn_j*(qn_j+1.0_f8)) / (2.0_f8*massBC*r_LR(1)**2)
            lr_rotMat(ijK) = min(jEigen, TMaxCut)
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

        do iZ = 1, IALR%nZ_IALR
            Zfact(iZ) = 1.0_f8/(2.0_f8*massTot*Z_IALR(iZ)**2)
        end do

!> Get CPM in whole IALR region
        do idjK = 1, int_njK
            qn_j = int_jKPair(idjK,1)
            K = int_jKPair(idjK,2)
            do jdjK = 1, int_njK
                qn_jp = int_jKPair(jdjK,1)
                Kp = int_jKPair(jdjK,2)
                if (qn_j /= qn_jp) cycle
                
                term1 = merge((initWP%Jtot*(initWP%Jtot+1.0_f8)+qn_j*(qn_j+1.0_f8)-2.0_f8*K*K), &
                            0.0_f8, K == Kp)

                delta = merge(dsqrt(2.0_f8), 1.0_f8, K == 0)
                term2 = merge(delta*lambdaPlus(initWP%Jtot,K)*lambdaPlus(qn_j,K), &
                            0.0_f8, Kp == K+1)
                delta = merge(dsqrt(2.0_f8), 1.0_f8, K == 1)
                term3 = merge(delta*lambdaMinus(initWP%Jtot,K)*lambdaMinus(qn_j,K), &
                            0.0_f8, Kp == K-1)
                CPM(idjK, jdjK) = term1 - term2 - term3
                CPM(jdjK, idjK) = CPM(idjK, jdjK)
            end do
        end do

!> The phase is determined by CPM(:,jdjK) and (-1)^(j+Kmax)CPM(:,jdjK) is set to be positive and others are adjusted accordingly
!> Ref: J. Chem. Phys. 158, 054801 (2023)
!> nlo: lower index, nhi: higher index, ndo: numbers of dimension
        nlo = 0
        do while (nlo /= int_njK)
            if (nlo > int_njK) then 
                write(outFileUnit,*) "Error in nlo!"
                write(outFileUnit,*) "POSITION: propagation.f90, subroutine setCPMat()"
            end if 
            nlo = nlo + 1
            qn_j = int_jKPair(nlo,1)
            nhi = nlo + 1
            do while (nhi <= int_njK)
                qn_jp = int_jKPair(nhi,1)
                if (qn_j /= qn_jp) exit
                nhi = nhi + 1
            end do
            nhi = nhi - 1
            ndo = nhi - nlo + 1
            allocate(Bvj(ndo,ndo),CPE(ndo))

            do idxj = 1, ndo
                do idxi = 1, ndo
                    Bvj(idxi,idxj) = CPM(nlo+idxi-1,nlo+idxj-1)
                end do 
            end do

            allocate(work(1))
            call dsyev('V', 'U', ndo, Bvj, ndo, CPE, work, -1, info)
            lwork = int(work(1))
            deallocate(work)

            allocate(work(lwork))
            call dsyev('V', 'U', ndo, Bvj, ndo, CPE, work, lwork, info)
            if (info /= 0) then
                write(outFileUnit,*) "Error in DSYEV of int Bvj, info =", info
                write(outFileUnit,*) "POSITION: propagation.f90, subroutine setCPMat()"
                stop
            end if
            deallocate(work)
            
            Kmax = min(initWP%Jtot, qn_j)
            phase = (-1.0_f8)**(qn_j+Kmax)
            do idxj = 1, ndo 
                if (Bvj(ndo,idxj)*phase >= 0.0_f8) then 
                    do idxi = 1, ndo
                        CPM(nlo+idxi-1,nlo+idxj-1) = Bvj(idxi,idxj)
                    end do 
                else 
                    do idxi = 1, ndo
                        CPM(nlo+idxi-1,nlo+idxj-1) = -Bvj(idxi,idxj)
                    end do
                end if
                int_lo(nlo+idxj-1) = dabs(dsqrt(CPE(idxj)+0.25_f8)-0.5_f8)
            end do

            nlo = nhi
            deallocate(Bvj,CPE)
        end do
!> Cut the CPE if Z is too small
        CPMax = -100.0_f8
        CPMin = 0.0_f8
        do iZ = 1, IALR%nZ_IALR 
            do idjK = 1, int_njK
                do jdjK = 1, int_njK
                    delta = 0.0_f8
                    do K = 1, int_njK
                        lreal = int_lo(K)
                        leigen = dsqrt(lreal*(lreal+1.0_f8))
                        tmpE = min(leigen*Zfact(iZ), TMaxCut)
                        CPMax = max(CPMax, tmpE)
                        delta = delta + CPM(jdjK,K)*tmpE*CPM(idjK,K)
                    end do 
                    UMat(iZ,jdjK,idjK) = delta
                end do 
            end do 
        end do
!> Rebuild CPM in IA 
        IA_CPMat(1:IALR%nZ_IA,:,:) = UMat(1:IALR%nZ_IA,:,:)
!> Rebuild CPM in LR
        Kmax = min(initWP%Jtot, initWP%j0)
        do K = initWP%Kmin, Kmax
            do Kp = initWP%Kmin, Kmax
                idjK = lr_seqjK(initWP%j0,K)
                jdjK = lr_seqjK(initWP%j0,Kp)
                lr_CPMat(:,idjK,jdjK) = UMat(:,idjK,jdjK)
            end do
        end do

        deallocate(CPM,UMat)

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