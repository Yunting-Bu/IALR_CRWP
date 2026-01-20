program IALR_WP_CRWP
    use m_MachinaBasic, only : f8 
    use m_gPara
    use m_Potent
    use m_Basis
    use m_InitWP
    use m_Prop
    use m_InterCoorRCB
    implicit none

!> Test
    call initPara()
    call setBasis()
    call setInitTotWP()
    call setEnergyAM()
    call initRCB()
    call HamiltonianScale()
    call initAllVabs()
    call propProcess()


end program IALR_WP_CRWP
