program IALR_WP_CRWP
    use m_MachinaBasic, only : f8
    use m_gPara
    use m_Potent
    use m_Basis
    use m_InitWP
    use m_Prop
    use m_InterCoorRCB
    use m_Match
    implicit none

!> Test
    call initPara()
    call setBasis()
    call setInitTotWP()
    call setEnergyAM()
    call initAllVabs()
    call initRCB()
    call HamiltonianScale()
    call propProcess()
    call extractSmat('A+BC->C+AB')
    if (nChannel == 2) then
        call extractSmat('A+BC->B+AC')
    end if
    call writeSmat('A+BC->C+AB')
    if (nChannel == 2) then
        call writeSmat('A+BC->B+AC')
    end if


end program IALR_WP_CRWP
