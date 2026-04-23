module m_Match
    use m_MachinaBasic, only : f8, c8
    use m_gPara
    use m_Basis, only : rBasisTrans
    use m_Potent
    use m_InterCoorRCB, only : Evj_AB, Evj_AC
    implicit none

    real(f8), private, allocatable :: AtDMat_BC(:,:,:,:)

    public
    private :: setAtDMat_BC

contains

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setAtDMat_BC()
        implicit none
        integer :: iPES, ith, ir 
        real(f8) :: bond(3), AtD(nPES,nPES), Vadia(nPES)
        real(f8) :: r, th, Z

        allocate(AtDMat_BC(nPES,nPES,IALR%vasy,IALR%asy_nA))
        do ith = 1, IALR%asy_nA
            do ir = 1, IALR%vasy
                Z = Z_IALR(iInePos)
                r = r_Asy(ir)
                th = acos(asy_ANode(ith))
                call Jacobi2Bond(Z, r, th, atomMass(2), atomMass(3), bond)
                call diagDiaVmat(bond, AtD, Vadia)
                AtDMat_BC(:,:,ir,ith) = AtD
            end do
        end do

    end subroutine setAtDMat_BC
!> ------------------------------------------------------------------------------------------------------------------ <!
        
!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine extractSmat(type)
        implicit none
        character(len=*), intent(in) :: type
        complex(c8), allocatable :: adiaTIDWF(:,:,:)
        complex(c8), allocatable :: KbTIDWF(:,:,:)
        complex(c8), allocatable :: ovlpWF_BF(:,:,:), ovlpWF_SF(:,:,:)
        complex(c8), allocatable :: Smat_SF(:,:,:)
        complex(c8) :: WFtmp, af
        real(f8) :: eps, sinTheta, Ef, kf, rbZ, rbZP, rbZU, rbZUP
        integer :: iEtot, v, j, ijK, iIC, Kr, Kp, K, Kpp, ijKp, l
        integer :: Krmax, Kpmax, iPES, jPES, idjK

        Krmax = min(initWP%Jtot, IALR%jint)

        if (type == 'A+BC->C+AB') then
            Kpmax = min(initWP%Jtot, channel1%jmax)
            allocate(Smat_AB(0:channel1%vmax,channel1%njK,nPES,nEtot))
            allocate(adiaTIDWF(channel1%nIC,initWP%Kmin:Krmax,nPES))
            allocate(KbTIDWF(channel1%nIC,channel1%Kmin:Kpmax,nPES))
            allocate(ovlpWF_BF(0:channel1%vmax,channel1%njK,nPES))
            allocate(ovlpWF_SF(0:channel1%vmax,channel1%njK,nPES))

            Smat_AB = imgZero
            adiaTIDWF = imgZero
            KbTIDWF = imgZero

            do iEtot = 1, nEtot
                eps = epsilon(Etot(iEtot))
!> Trans to adiabatic representation
                do Kr = initWP%Kmin, Krmax
                    do iIC = 1, channel1%nIC
                        do iPES = 1, nPES
                            WFtmp = imgZero
                            do jPES = 1, nPES
                                WFtmp = WFtmp + intAB_TIDWF(iEtot,iIC,Kr,jPES)*AtDMat_AB(jPES,iPES,iIC)
                            end do
                            adiaTIDWF(iIC,Kr,iPES) = WFtmp
                        end do
                    end do
                end do
!> Trans Kr to Kp
                do iPES = 1, nPES
                    do iIC = 1, channel1%nIC
                        do Kp = channel1%Kmin, Kpmax
                            WFtmp = imgZero
                            do Kr = initWP%Kmin, Krmax
                                WFtmp = WFtmp + adiaTIDWF(iIC,Kr,iPES)*TKMat_AB(Kr,Kp,iIC)
                            end do
                            KbTIDWF(iIC,Kp,iPES) = WFtmp
                        end do
                    end do
                end do
!> Overlap with product ICWF in BF
                ovlpWF_BF = imgZero
                ovlpWF_SF = imgZero
                do iPES = 1, nPES
                    do v = 0, channel1%vmax
                        do ijK = 1, channel1%njK
                            j = AB_jKPair(ijK,1)
                            K = AB_jKPair(ijK,2)
                            if ((Etot(iEtot)-Evj_AB(v,j,iPES)) > eps) then
                                do iIC = 1, channel1%nIC
                                    ovlpWF_BF(v,ijK,iPES) = ovlpWF_BF(v,ijK,iPES) + &
                                                            ICWFvj_AB(iIC,v,ijK,iPES) * KbTIDWF(iIC,K,iPES)
                                end do
                            end if
                        end do
                    end do
                end do

!> BF to SF 
                do iPES = 1, nPES
                    do v = 0, channel1%vmax
                        do ijK = 1, channel1%njK
                            j = AB_jKPair(ijK,1)
                            K = AB_jKPair(ijK,2)
                            if ((Etot(iEtot)-Evj_AB(v,j,iPES)) > eps) then
                                WFtmp = imgZero
                                do Kp = channel1%Kmin, min(j, initWP%Jtot)
                                    ijKp = AB_seqjK(j,Kp)
                                    if (ijKp > 0) then
                                        WFtmp = WFtmp + BLK_AB(Kp,ijK) * ovlpWF_BF(v,ijKp,iPES)
                                    end if
                                end do
                                ovlpWF_SF(v,ijK,iPES) = WFtmp
                            end if
                        end do
                    end do
                end do

!> Calculate S-matrix in SF
                sinTheta = dsin(ChebyAngle(iEtot))
                allocate(Smat_SF(0:channel1%vmax,channel1%njK,nPES))
                Smat_SF = imgZero
                do iPES = 1, nPES
                    do ijK = 1, channel1%njK
                        j = AB_jKPair(ijK,1)
                        l = qn_l_AB(ijK)
                        do v = 0, channel1%vmax
                            Ef = Etot(iEtot) - Evj_AB(v,j,iPES)
                            if (Ef > eps) then
                                kf = dsqrt(2.0_f8 * channel1%massTot * Ef)
                                call rbesjy(l*1.0_f8, kf*channel1%Zinf, rbZ, rbZP, rbZU, rbZUP)
                                af = dsqrt(kf/(2.0_f8*pi*channel1%massTot)) / (-rbZU + img*rbZ)
                                Smat_SF(v,ijK,iPES) = af / initAM(iEtot)  &
                                                    * ovlpWF_SF(v,ijK,iPES) / (Hminus*sinTheta) 
                            end if
                        end do
                    end do
                end do
!> Store S-matrix in SF
                do iPES = 1, nPES
                    do ijK = 1, channel1%njK
                        j = AB_jKPair(ijK,1)
                        do v = 0, channel1%vmax
                            if ((Etot(iEtot)-Evj_AB(v,j,iPES)) > eps) then
                                Smat_AB(v,ijK,iPES,iEtot) = Smat_SF(v,ijK,iPES)
                            end if
                        end do
                    end do
                end do
                deallocate(Smat_SF)
            end do

            deallocate(adiaTIDWF, KbTIDWF, ovlpWF_BF, ovlpWF_SF)

        else if (type == 'A+BC->B+AC') then
            Kpmax = min(initWP%Jtot, channel2%jmax)
            allocate(Smat_AC(0:channel2%vmax,channel2%njK,nPES,nEtot))
            allocate(adiaTIDWF(channel2%nIC,initWP%Kmin:Krmax,nPES))
            allocate(KbTIDWF(channel2%nIC,channel2%Kmin:Kpmax,nPES))
            allocate(ovlpWF_BF(0:channel2%vmax,channel2%njK,nPES))
            allocate(ovlpWF_SF(0:channel2%vmax,channel2%njK,nPES))

            Smat_AC = imgZero

            do iEtot = 1, nEtot
            eps = epsilon(Etot(iEtot))
!> Trans to adiabatic representation
                adiaTIDWF = imgZero
                do Kr = initWP%Kmin, Krmax
                    do iIC = 1, channel2%nIC
                        do iPES = 1, nPES
                            WFtmp = imgZero
                            do jPES = 1, nPES
                                WFtmp = WFtmp + intAC_TIDWF(iEtot,iIC,Kr,jPES)*AtDMat_AC(jPES,iPES,iIC)
                            end do
                            adiaTIDWF(iIC,Kr,iPES) = WFtmp
                        end do
                    end do
                end do
!> Trans Kr to Kp
                KbTIDWF = imgZero
                do iPES = 1, nPES
                    do iIC = 1, channel2%nIC
                        do Kp = channel2%Kmin, Kpmax
                            WFtmp = imgZero
                            do Kr = initWP%Kmin, Krmax
                                WFtmp = WFtmp + adiaTIDWF(iIC,Kr,iPES)*TKMat_AC(Kr,Kp,iIC)
                            end do
                            KbTIDWF(iIC,Kp,iPES) = WFtmp
                        end do
                    end do
                end do
!> Overlap with product ICWF in BF
                ovlpWF_BF = imgZero
                ovlpWF_SF = imgZero
                do iPES = 1, nPES
                    do v = 0, channel2%vmax
                        do ijK = 1, channel2%njK
                            j = AC_jKPair(ijK,1)
                            K = AC_jKPair(ijK,2)
                            if ((Etot(iEtot)-Evj_AC(v,j,iPES)) > eps) then
                                do iIC = 1, channel2%nIC
                                    ovlpWF_BF(v,ijK,iPES) = ovlpWF_BF(v,ijK,iPES) + &
                                                            ICWFvj_AC(iIC,v,ijK,iPES) * KbTIDWF(iIC,K,iPES)
                                end do
                            end if
                        end do
                    end do
                end do

!> BF to SF (sum over K components; consistent with reference sub_rcb.f90)
                do iPES = 1, nPES
                    do v = 0, channel2%vmax
                        do ijK = 1, channel2%njK
                            j = AC_jKPair(ijK,1)
                            K = AC_jKPair(ijK,2)
                            if ((Etot(iEtot)-Evj_AC(v,j,iPES)) > eps) then
                                WFtmp = imgZero
                                do Kpp = channel2%Kmin, min(j, initWP%Jtot)
                                    ijKp = AC_seqjK(j,Kpp)
                                    if (ijKp > 0) then
                                        WFtmp = WFtmp + BLK_AC(Kpp,ijK) * ovlpWF_BF(v,ijKp,iPES)
                                    end if
                                end do
                                ovlpWF_SF(v,ijK,iPES) = WFtmp
                            end if
                        end do
                    end do
                end do

!> Calculate S-matrix in SF (same formula as reference sub_rcb.f90:294)
                sinTheta = dsin(ChebyAngle(iEtot))
                allocate(Smat_SF(0:channel2%vmax,channel2%njK,nPES))
                Smat_SF = imgZero
                do iPES = 1, nPES
                    do ijK = 1, channel2%njK
                        j = AC_jKPair(ijK,1)
                        l = qn_l_AC(ijK)
                        do v = 0, channel2%vmax
                            Ef = Etot(iEtot) - Evj_AC(v,j,iPES)
                            if (Ef > eps) then
                                kf = dsqrt(2.0_f8 * channel2%massTot * Ef)
                                call rbesjy(l*1.0_f8, kf*channel2%Zinf, rbZ, rbZP, rbZU, rbZUP)
                                Smat_SF(v,ijK,iPES) = dsqrt(kf/(2.0_f8*pi*channel2%massTot)) / initAM(iEtot) &
                                    * ovlpWF_SF(v,ijK,iPES) / (Hminus*sinTheta) &
                                    / (-rbZU + img*rbZ)
                            end if
                        end do
                    end do
                end do
!> Store S-matrix in SF
                do iPES = 1, nPES
                    do ijK = 1, channel2%njK
                        j = AC_jKPair(ijK,1)
                        do v = 0, channel2%vmax
                            if ((Etot(iEtot)-Evj_AC(v,j,iPES)) > eps) then
                                Smat_AC(v,ijK,iPES,iEtot) = Smat_SF(v,ijK,iPES)
                            end if
                        end do
                    end do
                end do
                deallocate(Smat_SF)
            end do

            deallocate(adiaTIDWF, KbTIDWF, ovlpWF_BF, ovlpWF_SF)

        end if

    end subroutine extractSmat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine writeSmat(type)
        implicit none
        character(len=*), intent(in) :: type
        integer :: iEtot, v, j, l, ijK, iPES
        integer :: fileUnit
        character(len=100) :: filename
        real(f8) :: Ef, kf, Sreal, Simg, S2
        complex(c8) :: Stmp

        if (type == 'A+BC->C+AB') then
            do iPES = 1, nPES
                !> Generate filename: AB_J{Jtot}L{l0}iPES{iPES}.Smat
                write(filename, '(a,i0,a,i0,a,i0,a)') 'AB_J', initWP%Jtot, 'L', initWP%l0, 'iPES', iPES, '.Smat'

                open(newunit=fileUnit, file=trim(filename), status='replace', action='write')

                !> Write header information
                write(fileUnit, '(a)') '# S-matrix for A+BC -> C+AB channel'
                write(fileUnit, '(a,3f18.10)') '# Atom masses (au): ', atomMass(1), atomMass(2), atomMass(3)
                write(fileUnit, '(a,f18.10)') '# Total reduced mass A+BC (au): ', massTot
                write(fileUnit, '(a,f18.10)') '# Total reduced mass C+AB (au): ', channel1%massTot
                write(fileUnit, '(a,3i5)') '# Initial state v0, j0, l0: ', initWP%v0, initWP%j0, initWP%l0
                write(fileUnit, '(a,2i8)') '# Jtot, nEtot: ', initWP%Jtot, nEtot
                write(fileUnit, '(a)') '#'
                write(fileUnit, '(a)') '# iEtot  Ecol(eV)  v  j  l  kf(au)  Re(S)  Im(S)  |S|^2'
                write(fileUnit, '(a)') '#'

                !> Write S-matrix elements
                do iEtot = 1, nEtot
                    do ijK = 1, channel1%njK
                        j = AB_jKPair(ijK,1)
                        l = qn_l_AB(ijK)
                        do v = 0, channel1%vmax
                            Ef = Etot(iEtot) - Evj_AB(v,j,iPES)
                            if (Ef > epsilon(Etot(iEtot))) then
                                kf = dsqrt(2.0_f8 * channel1%massTot * Ef)
                                Stmp = Smat_AB(v,ijK,iPES,iEtot)
                                Sreal = Stmp%re
                                Simg = Stmp%im
                                S2 = cdabs(Stmp)**2
                                write(fileUnit, '(i6,2x,f15.8,2x,3i4,2x,4f13.8)') &
                                    iEtot, Ecol(iEtot)*au2ev, v, j, l, kf, Sreal, Simg, S2
                            end if
                        end do
                    end do
                end do

                close(fileUnit)
                write(outFileUnit,'(1x,a,a)') 'AB S-matrix written to file: ', trim(filename)
            end do

        else if (type == 'A+BC->B+AC') then
            do iPES = 1, nPES
                !> Generate filename: AC_J{Jtot}L{l0}iPES{iPES}.Smat
                write(filename, '(a,i0,a,i0,a,i0,a)') 'AC_J', initWP%Jtot, 'L', initWP%l0, 'iPES', iPES, '.Smat'

                open(newunit=fileUnit, file=trim(filename), status='replace', action='write')

                !> Write header information
                write(fileUnit, '(a)') '# S-matrix for A+BC -> B+AC channel'
                write(fileUnit, '(a,3f18.10)') '# Atom masses (au): ', atomMass(1), atomMass(2), atomMass(3)
                write(fileUnit, '(a,f18.10)') '# Total reduced mass A+BC (au): ', massTot
                write(fileUnit, '(a,f18.10)') '# Total reduced mass B+AC (au): ', channel2%massTot
                write(fileUnit, '(a,3i5)') '# Initial state v0, j0, l0: ', initWP%v0, initWP%j0, initWP%l0
                write(fileUnit, '(a,2i8)') '# Jtot, nEtot: ', initWP%Jtot, nEtot
                write(fileUnit, '(a)') '#'
                write(fileUnit, '(a)') '# iEtot  Ecol(eV)  v  j  l  kf(au)  Re(S)  Im(S)  |S|^2'
                write(fileUnit, '(a)') '#'

                !> Write S-matrix elements
                do iEtot = 1, nEtot
                    do ijK = 1, channel2%njK
                        j = AC_jKPair(ijK,1)
                        l = qn_l_AC(ijK)
                        do v = 0, channel2%vmax
                            Ef = Etot(iEtot) - Evj_AC(v,j,iPES)
                            if (Ef > epsilon(Etot(iEtot))) then
                                kf = dsqrt(2.0_f8 * channel2%massTot * Ef)
                                Stmp = Smat_AC(v,ijK,iPES,iEtot)
                                Sreal = Stmp%re
                                Simg = Stmp%im
                                S2 = cdabs(Stmp)**2
                                write(fileUnit, '(i6,2x,f15.8,2x,3i4,2x,4f13.8)') &
                                    iEtot, Ecol(iEtot)*au2ev, v, j, l, kf, Sreal, Simg, S2
                            end if
                        end do
                    end do
                end do

                close(fileUnit)
                write(outFileUnit,'(1x,a,a)') 'AC S-matrix written to file: ', trim(filename)
            end do

        end if

    end subroutine writeSmat
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine setInelastic()
        implicit none
        real(f8), external :: CG
        integer :: ijK, j, K, l, idx, Kmax
        integer :: iPES, v, no, fileUnit, nCh, iCh, jCh
        integer :: vTmp, jTmp, KTmp
        real(f8) :: fact, eTmp
        real(f8) :: EtotMin
        character(len=5) :: ctype, typeTmp
        character(len=100) :: fileName
        integer, allocatable :: chV(:), chJ(:), chK(:)
        real(f8), allocatable :: chE(:)
        character(len=5), allocatable :: chType(:)

        Kmax = min(IALR%jasy, initWP%Jtot)

!> qn_l for BC: map (j,K) to l via parity rule (-1)^(j+l+Jtot) = tpar
        allocate(qn_l_BC(asy_njK))
        do ijK = 1, asy_njK
            j = asy_jKPair(ijK,1)
            K = asy_jKPair(ijK,2)
            idx = 0
            do l = abs(initWP%Jtot - j), initWP%Jtot + j
                if ((-1)**(l+j+initWP%Jtot) == initWP%tpar) then
                    if (idx == K - initWP%Kmin) then
                        qn_l_BC(ijK) = l
                        exit
                    end if
                    idx = idx + 1
                end if
            end do
        end do

!> BLK_BC: BF to SF transformation via CG coefficients
        allocate(BLK_BC(initWP%Kmin:Kmax,asy_njK))
        BLK_BC = 0.0_f8
        do ijK = 1, asy_njK
            j = asy_jKPair(ijK,1)
            l = qn_l_BC(ijK)
            fact = dsqrt((2.0_f8*l+1.0_f8)/(2.0_f8*initWP%Jtot+1.0_f8))
            do K = initWP%Kmin, min(j, Kmax)
                BLK_BC(K,ijK) = fact*CG(j,K,l,0,initWP%Jtot)
                if (K /= 0) BLK_BC(K,ijK) = BLK_BC(K,ijK)*dsqrt(2.0_f8)
            end do
        end do

        write(outFileUnit,'(1x,a)') '=====> Inelastic channel setup <====='
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a,i4)') 'Number of (j,K) pairs for BC: ', asy_njK
        write(outFileUnit,'(1x,a,i4)') 'Number of vibrational states :  ', IALR%vasy
        write(outFileUnit,'(1x,a,i4)') 'Total number of channels     : ', IALR%vasy*asy_njK
        write(outFileUnit,'(1x,a)') ''
        write(outFileUnit,'(1x,a)') ''

!> Channel list in v-j-K representation at minimum collision energy
        EtotMin = Etot(1)
        nCh = IALR%vasy * asy_njK
        allocate(chV(nCh), chJ(nCh), chK(nCh), chE(nCh), chType(nCh))
        do iPES = 1, nPES
            if (nPES == 1) then
                fileName = 'BC_channel.chk'
            else
                write(fileName,'(a,i0,a)') 'BC_channel_iPES', iPES, '.chk'
            end if
            open(newunit=fileUnit, file=trim(fileName), status='replace', action='write')

            no = 0
            do ijK = 1, asy_njK
                j = asy_jKPair(ijK,1)
                K = asy_jKPair(ijK,2)
                do v = 0, IALR%vasy-1
                    no = no + 1
                    chV(no) = v
                    chJ(no) = j
                    chK(no) = K
                    chE(no) = Evj_BC(v,j,iPES)
                    if ((EtotMin - chE(no)) > epsilon(EtotMin)) then
                        chType(no) = 'open '
                    else
                        chType(no) = 'close'
                    end if
                end do
            end do

!> Sort channel list by Evj (ascending)
            do iCh = 1, nCh-1
                do jCh = iCh+1, nCh
                    if (chE(jCh) < chE(iCh)) then
                        eTmp = chE(iCh); chE(iCh) = chE(jCh); chE(jCh) = eTmp
                        vTmp = chV(iCh); chV(iCh) = chV(jCh); chV(jCh) = vTmp
                        jTmp = chJ(iCh); chJ(iCh) = chJ(jCh); chJ(jCh) = jTmp
                        KTmp = chK(iCh); chK(iCh) = chK(jCh); chK(jCh) = KTmp
                        typeTmp = chType(iCh); chType(iCh) = chType(jCh); chType(jCh) = typeTmp
                    end if
                end do
            end do

            write(fileUnit,'(a)') '# ==============================================================='
            write(fileUnit,'(a,f15.8,a)') '# Minimum collision energy: ', Ecol(1)*au2ev, ' eV'
            write(fileUnit,'(a)') '# Channel list sorted by Evj (ascending)'
            write(fileUnit,'(a)') '# ---------------------------------------------------------------'
            write(fileUnit,'(a)') '#   NO    v    j    K      Evj(eV)      Evj(cm-1)    type'
            write(fileUnit,'(a)') '# ---------------------------------------------------------------'
            do iCh = 1, nCh
                ctype = chType(iCh)
                write(fileUnit,'(i6,1x,i4,1x,i4,1x,i4,2x,f13.9,2x,f13.3,2x,a5)') &
                    iCh, chV(iCh), chJ(iCh), chK(iCh), chE(iCh)*au2ev, chE(iCh)*au2cm, ctype
            end do
            write(fileUnit,'(a)') '# ==============================================================='

            close(fileUnit)
            write(outFileUnit,'(1x,a,a)') 'Inelastic channel list written: ', trim(fileName)
        end do
        deallocate(chV, chJ, chK, chE, chType)

    end subroutine setInelastic
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine extractSmat_ine()
        implicit none
        complex(c8), allocatable :: ovlpWF_BF(:,:,:), ovlpWF_SF(:,:,:)
        complex(c8), allocatable :: Smat_SF(:,:,:)
        complex(c8), allocatable :: adiaTIDWF(:,:,:,:)
        complex(c8) :: WFtmp, af
        real(f8) :: eps, sinTheta, Ef, kf, rbZ, rbZP, rbZU, rbZUP
        integer :: iEtot, v, j, ijK, K, Kp, ijKp, l
        integer :: Kmax, iPES, jPES, ir, ith

        Kmax = min(initWP%Jtot, IALR%jasy)

        call setAtDMat_BC()

        allocate(Smat_BC(0:IALR%vasy-1,asy_njK,nPES,nEtot))
        allocate(ovlpWF_BF(0:IALR%vasy-1,asy_njK,nPES))
        allocate(ovlpWF_SF(0:IALR%vasy-1,asy_njK,nPES))
        allocate(adiaTIDWF(IALR%vasy,IALR%asy_nA,initWP%Kmin:Kmax,nPES))

        Smat_BC = imgZero
        adiaTIDWF = imgZero

        do iEtot = 1, nEtot
            eps = epsilon(Etot(iEtot))
!> Trans to adiabatic representation
            adiaTIDWF = imgZero
            do K = initWP%Kmin, Kmax
                do ith = 1, IALR%asy_nA
                    do ir = 1, IALR%vasy
                        do iPES = 1, nPES
                            WFtmp = imgZero
                            do jPES = 1, nPES
                                WFtmp = WFtmp + ine_TIDWF(iEtot,ir,ith,K,jPES)*AtDMat_BC(jPES,iPES,ir,ith)
                            end do
                            adiaTIDWF(ir,ith,K,iPES) = WFtmp
                        end do
                    end do
                end do
            end do 
!> Overlap with BC ro-vib WF in BF
            ovlpWF_BF = imgZero
            ovlpWF_SF = imgZero
            do iPES = 1, nPES
                do v = 0, IALR%vasy-1
                    do ijK = 1, asy_njK
                        j = asy_jKPair(ijK,1)
                        K = asy_jKPair(ijK,2)
                        if ((Etot(iEtot)-Evj_BC(v,j,iPES)) > eps) then
                            do ith = 1, IALR%asy_nA
                                do ir = 1, IALR%vasy
                                    ovlpWF_BF(v,ijK,iPES) = ovlpWF_BF(v,ijK,iPES) + &
                                        adiaTIDWF(ir,ith,K,iPES)*WFvj_BC(ir,ith,v,j,K,iPES)
                                end do
                            end do
                        end if
                    end do
                end do
            end do
                
!> BF to SF
            do iPES = 1, nPES
                do v = 0, IALR%vasy-1
                    do ijK = 1, asy_njK
                        j = asy_jKPair(ijK,1)
                        K = asy_jKPair(ijK,2)
                        if ((Etot(iEtot)-Evj_BC(v,j,iPES)) > eps) then
                            WFtmp = imgZero
                            do Kp = initWP%Kmin, min(j, initWP%Jtot)
                                ijKp = asy_seqjK(j,Kp)
                                if (ijKp > 0) then
                                    WFtmp = WFtmp + BLK_BC(Kp,ijK) * ovlpWF_BF(v,ijKp,iPES)
                                end if
                            end do
                            ovlpWF_SF(v,ijK,iPES) = WFtmp
                        end if
                    end do
                end do
            end do
!> Calculate S-matrix in SF
            sinTheta = dsin(ChebyAngle(iEtot))
            allocate(Smat_SF(0:IALR%vasy-1,asy_njK,nPES))
            Smat_SF = imgZero
            do iPES = 1, nPES
                do ijK = 1, asy_njK
                    j = asy_jKPair(ijK,1)
                    l = qn_l_BC(ijK)
                    do v = 0, IALR%vasy-1
                        Ef = Etot(iEtot) - Evj_BC(v,j,iPES)
                        if (Ef > eps) then
                            kf = dsqrt(2.0_f8 * massTot * Ef)
                            call rbesjy(l*1.0_f8, kf*Z_IALR(iInePos), rbZ, rbZP, rbZU, rbZUP)
                            af = dsqrt(kf/(2.0_f8*pi*massTot)) / (-rbZU + img*rbZ)
                            Smat_SF(v,ijK,iPES) = af / initAM(iEtot) &
                                                * ovlpWF_SF(v,ijK,iPES) / (Hminus*sinTheta)
                        end if
                    end do
                end do
            end do
!> Store S-matrix in SF
            do iPES = 1, nPES
                do ijK = 1, asy_njK
                    j = asy_jKPair(ijK,1)
                    do v = 0, IALR%vasy-1
                        if ((Etot(iEtot)-Evj_BC(v,j,iPES)) > eps) then
                            Smat_BC(v,ijK,iPES,iEtot) = Smat_SF(v,ijK,iPES)
                        end if
                    end do
                end do
            end do
            deallocate(Smat_SF)
        end do

        deallocate(adiaTIDWF, ovlpWF_BF, ovlpWF_SF)


    end subroutine extractSmat_ine
!> ------------------------------------------------------------------------------------------------------------------ <!

!> ------------------------------------------------------------------------------------------------------------------ <!
    subroutine writeSmat_ine()
!> Write inelastic S-matrix (A+BC -> A+BC) to file
        implicit none
        integer :: iEtot, v, j, l, ijK, iPES
        integer :: fileUnit
        character(len=100) :: filename
        real(f8) :: Ef, kf, Sreal, Simg, S2
        complex(c8) :: Stmp

        do iPES = 1, nPES
            write(filename, '(a,i0,a,i0,a,i0,a)') 'BC_J', initWP%Jtot, 'L', initWP%l0, 'iPES', iPES, '.Smat'

            open(newunit=fileUnit, file=trim(filename), status='replace', action='write')

            write(fileUnit, '(a)') '# S-matrix for A+BC -> A+BC inelastic channel'
            write(fileUnit, '(a,3f18.10)') '# Atom masses (au): ', atomMass(1), atomMass(2), atomMass(3)
            write(fileUnit, '(a,f18.10)') '# Total reduced mass A+BC (au): ', massTot
            write(fileUnit, '(a,3i5)') '# Initial state v0, j0, l0: ', initWP%v0, initWP%j0, initWP%l0
            write(fileUnit, '(a,2i8)') '# Jtot, nEtot: ', initWP%Jtot, nEtot
            write(fileUnit, '(a,f18.10)') '# Projection plane Z (au): ', Z_IALR(iInePos)
            write(fileUnit, '(a)') '#'
            write(fileUnit, '(a)') '# iEtot  Ecol(eV)  v  j  l  kf(au)  Re(S)  Im(S)  |S|^2'
            write(fileUnit, '(a)') '#'

            do iEtot = 1, nEtot
                do ijK = 1, asy_njK
                    j = asy_jKPair(ijK,1)
                    l = qn_l_BC(ijK)
                    do v = 0, IALR%vasy-1
                        Ef = Etot(iEtot) - Evj_BC(v,j,iPES)
                        if (Ef > epsilon(Etot(iEtot))) then
                            kf = dsqrt(2.0_f8 * massTot * Ef)
                            Stmp = Smat_BC(v,ijK,iPES,iEtot)
                            Sreal = Stmp%re
                            Simg = Stmp%im
                            S2 = cdabs(Stmp)**2
                            write(fileUnit, '(i6,2x,f15.8,2x,3i4,2x,4f13.8)') &
                                iEtot, Ecol(iEtot)*au2ev, v, j, l, kf, Sreal, Simg, S2
                        end if
                    end do
                end do
            end do

            close(fileUnit)
            write(outFileUnit,'(1x,a,a)') 'BC S-matrix written to file: ', trim(filename)
        end do

    end subroutine writeSmat_ine
!> ------------------------------------------------------------------------------------------------------------------ <!

end module m_Match
