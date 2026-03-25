module m_Match
    use m_MachinaBasic, only : f8, c8
    use m_gPara
    use m_InterCoorRCB, only : Evj_AB, Evj_AC
    implicit none

    public

contains

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

end module m_Match