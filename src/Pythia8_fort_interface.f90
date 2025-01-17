!! Pythia8_fort_interface.f90 is a part of the PACIAE event generator.
!! Copyright (C) 2024 PACIAE Group.
!! PACIAE is licensed under the GNU GPL v2 or later, see LICENSE for details.
!! Open source: https://github.com/ArcsaberHep/PACIAE4
!! Author: An-Ke Lei, January 2024 - November 2024.

!> This is the Fortran interface program to link PACIAE (Fortran 77/90) with
!!  PYTHIA 8 (C++).

!!                                               By An-Ke at CCNU on 16/01/2024
!!                                  Last updated by An-Ke at UiO  on 17/01/2025


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>  Interfaces into PYTHIA 8, *pythia = new Pythia. Instantiates a basic
!!    Pythia obeject.
    SUBROUTINE PYINST_PY8

!---Imports mudules.
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

!---Fortran data type.

!***********************************************************************
!   Local variabls.
    INTEGER, PARAMETER :: KSZJ = 80000, KSZJ_PY8 = 300000
    ! INTEGER :: I, J
    ! INTEGER :: frameType, idA, idB
    ! REAL(KIND=8) :: eCM
!***********************************************************************

!***********************************************************************
!---PACIAE.
!-----------------------------------------------------------------------
    INTEGER :: itorw, iikk, kkii
    REAL(KIND=8) :: cp0, cr0
    COMMON/SA34/ itorw, iikk, cp0, cr0, kkii
!-----------------------------------------------------------------------
    INTEGER :: i_mode, i_tune
    INTEGER :: KF_woDecay(100), KF_proj, KF_targ
    REAL(KIND=8) :: win, energy_B, psno, b_min, b_max
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay, &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!***********************************************************************

!***********************************************************************
!---PYTHIA 6.
!-----------------------------------------------------------------------
    INTEGER :: MINT(400)
    REAL(KIND=8) :: VINT(400)
    COMMON/PYINT1/ MINT, VINT
    SAVE /PYINT1/
!-----------------------------------------------------------------------
    INTEGER :: MSTU(200), MSTJ(200)
    REAL(KIND=8) :: PARU(200), PARJ(200)
    COMMON/PYDAT1/ MSTU, PARU, MSTJ, PARJ
    SAVE /PYDAT1/
!-----------------------------------------------------------------------
    INTEGER :: MSTP(200), MSTI(200)
    REAL(KIND=8) :: PARP(200), PARI(200)
    COMMON/PYPARS/ MSTP, PARP, MSTI, PARI
    SAVE /PYPARS/
!-----------------------------------------------------------------------
    INTEGER :: N, NPAD
    INTEGER :: K(KSZJ,5)
    REAL(KIND=8) :: P(KSZJ,5), V(KSZJ,5)
    COMMON/PYJETS/ N, NPAD, K, P, V
    SAVE /PYJETS/
!***********************************************************************

!***********************************************************************
!---PYTHIA 8 information storage for feeding-back.
!-----------------------------------------------------------------------
    INTEGER :: N_PY8, NPAD_PY8
    ! INTEGER, PARAMETER :: KSZJ_PY8 = 300000
    INTEGER :: K_PY8(KSZJ_PY8,8)
    REAL(KIND=8) :: P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
    COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, K_PY8, P_PY8, V_PY8
    SAVE /PYJETS_PY8/
!***********************************************************************

!---C++ data type.

!***********************************************************************
!   Local variabls.
    ! INTEGER(KIND=C_INT) :: frameType_c, idA_c, idB_c
    ! REAL(KIND=C_DOUBLE) :: eCM_c
    INTEGER(KIND=C_INT) :: idStable(100)
    ! INTEGER(KIND=C_INT) :: iFail
!***********************************************************************

!***********************************************************************
    INTEGER(KIND=C_INT) :: MINT_c(400)
    REAL(KIND=C_DOUBLE) :: VINT_c(400)
    COMMON/PYINT1_c/ MINT_c, VINT_c
    SAVE /PYINT1_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: MSTU_c(200), MSTJ_c(200)
    REAL(KIND=C_DOUBLE) :: PARU_c(200), PARJ_c(200)
    COMMON/PYDAT1_c/ MSTU_c, PARU_c, MSTJ_c, PARJ_c
    SAVE /PYDAT1_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: MSTP_c(200), MSTI_c(200)
    REAL(KIND=C_DOUBLE) :: PARP_c(200), PARI_c(200)
    COMMON/PYPARS_c/ MSTP_c, PARP_c, MSTI_c, PARI_c
    SAVE /PYPARS_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: nPY8, nPadPY8
    ! INTEGER(C_INT), PARAMETER :: kSZJ_c = 300000
    INTEGER(KIND=C_INT) :: kPY8(KSZJ_PY8,8)
    REAL(KIND=C_DOUBLE) :: pPY8(KSZJ_PY8,7), vPY8(KSZJ_PY8,5)
    COMMON/PYJETS_c/ nPY8, nPadPY8, kPY8, pPY8, vPY8
    SAVE /PYJETS_c/
!***********************************************************************

!---Pointer.

!***********************************************************************
!---NOTE: DO NOT TOUCH THESES POINTERS IF YOU DON'T KNOW WHAT THEY ARE !!!
!-----------------------------------------------------------------------
!---Pointers to PYTHIA8 obejects.
    TYPE(C_PTR) :: pythia8, pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
    COMMON/PYTHIA8_PTR/ pythia8, &
                        pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
!-----------------------------------------------------------------------
!---Pointers to PACIAE4 obejects.
    TYPE(C_PTR) :: pahooks, paHIhooks, &
                   pahooks_pp, pahooks_pn, pahooks_np, pahooks_nn
    COMMON/PACIAE4_PTR/ pahooks, paHIhooks, &
                        pahooks_pp, pahooks_pn, pahooks_np, pahooks_nn
!***********************************************************************

!---Interfaces to C++.

!***********************************************************************
    INTERFACE
        SUBROUTINE instantiation_PY8( MINT_c, VINT_c, &
                                      MSTU_c, PARU_c, MSTJ_c, PARJ_c, &
                                      MSTP_c, PARP_c, MSTI_c, PARI_c, &
                                      nPY8, kPY8, pPY8, vPY8, idStable, &
                                      pythia8, pahooks, paHIhooks ) &
                                      BIND( C, NAME="instantiation_PY8" )
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND=C_INT) :: idStable(*)
            INTEGER(KIND=C_INT) :: MINT_c(*)
            REAL(KIND=C_DOUBLE) :: VINT_c(*)
            INTEGER(KIND=C_INT) :: MSTU_c(*), MSTJ_c(*)
            REAL(KIND=C_DOUBLE) :: PARU_c(*), PARJ_c(*)
            INTEGER(KIND=C_INT) :: MSTP_c(*), MSTI_c(*)
            REAL(KIND=C_DOUBLE) :: PARP_c(*), PARI_c(*)
            INTEGER(C_INT) :: nPY8
            INTEGER(C_INT) :: kPY8(*)
            REAL(C_DOUBLE) :: pPY8(*), vPY8(*)
            ! INTEGER(KIND=C_INT) :: iFail
            TYPE(C_PTR) :: pythia8
            TYPE(C_PTR) :: pahooks, paHIhooks
        END SUBROUTINE instantiation_PY8
    END INTERFACE
!***********************************************************************


!  Warning: the INTEGER would be converted from integer*8 -> integer*4.
!           But dont't worry, if we use PYTHIA 8. Because the storage of
!           the color flow was redesigned in PYTHIA 8.
    ! frameType_c = frameType
    ! idA_c = idA
    ! idB_c = idB
    ! eCM_c = eCM
    idStable = KF_woDecay

    MINT_c = MINT
    VINT_c = VINT
    MSTU_c = MSTU
    PARU_c = PARU
    MSTJ_c = MSTJ
    PARJ_c = PARJ
    MSTP_c = MSTP
    PARP_c = PARP
    MSTI_c = MSTI
    PARI_c = PARI

!   Accesses to C++ program in Pythia8_cpp_interface.cpp.
    CALL instantiation_PY8( MINT_c, VINT_c, &
                            MSTU_c, PARU_c, MSTJ_c, PARJ_c, &
                            MSTP_c, PARP_c, MSTI_c, PARI_c, &
                            nPY8, kPY8, pPY8, vPY8, idStable, &
                            pythia8, pahooks, paHIhooks )


    RETURN
    END



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>  Interfaces into C++ (PYTHIA 8), delete *pythia. Deletes the Pythia obeject.
    SUBROUTINE PYDELE_PY8

!---Imports mudules.
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

!---C++ data type.

!***********************************************************************
    INTEGER(KIND=C_INT) :: MINT_c(400)
    REAL(KIND=C_DOUBLE) :: VINT_c(400)
    COMMON/PYINT1_c/ MINT_c, VINT_c
    SAVE /PYINT1_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: MSTP_c(200), MSTI_c(200)
    REAL(KIND=C_DOUBLE) :: PARP_c(200), PARI_c(200)
    COMMON/PYPARS_c/ MSTP_c, PARP_c, MSTI_c, PARI_c
    SAVE /PYPARS_c/
!-----------------------------------------------------------------------
!***********************************************************************

!---Pointer.

!***********************************************************************
!---NOTE: DO NOT TOUCH THESES POINTERS IF YOU DON'T KNOW WHAT THEY ARE !!!
!-----------------------------------------------------------------------
!---Pointers to PYTHIA8 obejects.
    TYPE(C_PTR) :: pythia8, pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
    COMMON/PYTHIA8_PTR/ pythia8, &
                        pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
!-----------------------------------------------------------------------
!---Pointers to PACIAE4 obejects.
    TYPE(C_PTR) :: pahooks, paHIhooks, &
                   pahooks_pp, pahooks_pn, pahooks_np, pahooks_nn
    COMMON/PACIAE4_PTR/ pahooks, paHIhooks, &
                        pahooks_pp, pahooks_pn, pahooks_np, pahooks_nn
!***********************************************************************

!---Interfaces to C++.

!***********************************************************************
    INTERFACE
        SUBROUTINE delete_object_PY8( MINT_c, MSTP_c, pythia8, &
                                      pythia8_pp, pythia8_pn, &
                                      pythia8_np, pythia8_nn, &
                                      pahooks, paHIhooks, &
                                      pahooks_pp, pahooks_pn, &
                                      pahooks_np, pahooks_nn ) &
            BIND( C, NAME="delete_object_PY8" )
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND=C_INT) :: MINT_c(*), MSTP_c(*)
            TYPE(C_PTR) :: pythia8, &
                           pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
            TYPE(C_PTR) :: pahooks, paHIhooks, &
                           pahooks_pp, pahooks_pn, pahooks_np, pahooks_nn
        END SUBROUTINE delete_object_PY8
    END INTERFACE
!***********************************************************************


!   Accesses to C++ program in Pythia8_cpp_interface.cpp .
    CALL delete_object_PY8( MINT_c, MSTP_c, pythia8, &
                            pythia8_pp, pythia8_pn, &
                            pythia8_np, pythia8_nn, &
                            pahooks, paHIhooks, &
                            pahooks_pp, pahooks_pn, &
                            pahooks_np, pahooks_nn )


    RETURN
    END



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>  Interfaces into PYTHIA 8, Pythia::init(). Initilizes the PYTHIA 8 generator.
    SUBROUTINE PYINIT_PY8( frameType, idA, idB, eCM )

!---Imports mudules.
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

!---Fortran data type.

!***********************************************************************
!   Local variabls.
    INTEGER, PARAMETER :: KSZJ = 80000, KSZJ_PY8 = 300000
    INTEGER :: I !, J
    INTEGER :: frameType, idA, idB
    REAL(KIND=8) :: eCM
!***********************************************************************

!***********************************************************************
!---PACIAE.
!-----------------------------------------------------------------------
    INTEGER :: itorw, iikk, kkii
    REAL(KIND=8) :: cp0, cr0
    COMMON/SA34/ itorw, iikk, cp0, cr0, kkii
!-----------------------------------------------------------------------
    INTEGER :: i_mode, i_tune
    INTEGER :: KF_woDecay(100), KF_proj, KF_targ
    REAL(KIND=8) :: win, energy_B, psno, b_min, b_max
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay, &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!***********************************************************************

!***********************************************************************
!---PYTHIA 6.
!-----------------------------------------------------------------------
    INTEGER :: MINT(400)
    REAL(KIND=8) :: VINT(400)
    COMMON/PYINT1/ MINT, VINT
    SAVE /PYINT1/
!-----------------------------------------------------------------------
    INTEGER :: MSTU(200), MSTJ(200)
    REAL(KIND=8) :: PARU(200), PARJ(200)
    COMMON/PYDAT1/ MSTU, PARU, MSTJ, PARJ
    SAVE /PYDAT1/
!-----------------------------------------------------------------------
    INTEGER :: MSTP(200), MSTI(200)
    REAL(KIND=8) :: PARP(200), PARI(200)
    COMMON/PYPARS/ MSTP, PARP, MSTI, PARI
    SAVE /PYPARS/
!-----------------------------------------------------------------------
    INTEGER :: N, NPAD
    INTEGER :: K(KSZJ,5)
    REAL(KIND=8) :: P(KSZJ,5), V(KSZJ,5)
    COMMON/PYJETS/ N, NPAD, K, P, V
    SAVE /PYJETS/
!***********************************************************************

!***********************************************************************
!---PYTHIA 8 information storage for feeding-back.
!-----------------------------------------------------------------------
    INTEGER :: N_PY8, NPAD_PY8
    ! INTEGER, PARAMETER :: KSZJ_PY8 = 80000
    INTEGER :: K_PY8(KSZJ_PY8,8)
    REAL(KIND=8) :: P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
    COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, K_PY8, P_PY8, V_PY8
    SAVE /PYJETS_PY8/
!***********************************************************************

!---C++ data type.

!***********************************************************************
!   Local variabls.
    ! INTEGER(KIND=C_INT) :: frameType_c, idA_c, idB_c
    INTEGER(KIND=C_INT) :: frameType_c, idA_c, idB_c
    REAL(KIND=C_DOUBLE) :: eCM_c
    INTEGER(KIND=C_INT) :: idStable(100)
    ! INTEGER(KIND=C_INT) :: iFail
!***********************************************************************

!***********************************************************************
    INTEGER(KIND=C_INT) :: MINT_c(400)
    REAL(KIND=C_DOUBLE) :: VINT_c(400)
    COMMON/PYINT1_c/ MINT_c, VINT_c
    SAVE /PYINT1_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: MSTU_c(200), MSTJ_c(200)
    REAL(KIND=C_DOUBLE) :: PARU_c(200), PARJ_c(200)
    COMMON/PYDAT1_c/ MSTU_c, PARU_c, MSTJ_c, PARJ_c
    SAVE /PYDAT1_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: MSTP_c(200), MSTI_c(200)
    REAL(KIND=C_DOUBLE) :: PARP_c(200), PARI_c(200)
    COMMON/PYPARS_c/ MSTP_c, PARP_c, MSTI_c, PARI_c
    SAVE /PYPARS_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: nPY8, nPadPY8
    ! INTEGER(C_INT), PARAMETER :: kSZJ_c = 300000
    INTEGER(KIND=C_INT) :: kPY8(KSZJ_PY8,8)
    REAL(KIND=C_DOUBLE) :: pPY8(KSZJ_PY8,7), vPY8(KSZJ_PY8,5)
    COMMON/PYJETS_c/ nPY8, nPadPY8, kPY8, pPY8, vPY8
    SAVE /PYJETS_c/
!***********************************************************************

!---Pointer.

!***********************************************************************
!---NOTE: DO NOT TOUCH THESES POINTERS IF YOU DON'T KNOW WHAT THEY ARE !!!
!-----------------------------------------------------------------------
!---Pointers to PYTHIA8 obejects.
    TYPE(C_PTR) :: pythia8, pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
    COMMON/PYTHIA8_PTR/ pythia8, &
                        pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
!-----------------------------------------------------------------------
!---Pointers to PACIAE4 obejects.
    TYPE(C_PTR) :: pahooks, paHIhooks, &
                   pahooks_pp, pahooks_pn, pahooks_np, pahooks_nn
    COMMON/PACIAE4_PTR/ pahooks, paHIhooks, &
                        pahooks_pp, pahooks_pn, pahooks_np, pahooks_nn
!***********************************************************************

!---Interfaces to C++.

!***********************************************************************
    INTERFACE
        SUBROUTINE init_PY8( frameType_c, idA_c, idB_c, eCM_c, &
                             MINT_c, VINT_c, &
                             MSTU_c, PARU_c, MSTJ_c, PARJ_c, &
                             MSTP_c, PARP_c, MSTI_c, PARI_c, &
                             nPY8, kPY8, pPY8, vPY8, &
                             pythia8, &
                             pythia8_pp, pythia8_pn, &
                             pythia8_np, pythia8_nn, &
                             pahooks, paHIhooks ) &
                             BIND( C, NAME="init_PY8" )
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(C_INT), VALUE :: frameType_c, idA_c, idB_c
            REAL(C_DOUBLE), VALUE :: eCM_c
            INTEGER(KIND=C_INT) :: MINT_c(*)
            REAL(KIND=C_DOUBLE) :: VINT_c(*)
            INTEGER(KIND=C_INT) :: MSTU_c(*), MSTJ_c(*)
            REAL(KIND=C_DOUBLE) :: PARU_c(*), PARJ_c(*)
            INTEGER(KIND=C_INT) :: MSTP_c(*), MSTI_c(*)
            REAL(KIND=C_DOUBLE) :: PARP_c(*), PARI_c(*)
            INTEGER(C_INT) :: nPY8
            INTEGER(C_INT) :: kPY8(*)
            REAL(C_DOUBLE) :: pPY8(*), vPY8(*)
            TYPE(C_PTR) :: pythia8, &
                           pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
            TYPE(C_PTR) :: pahooks, paHIhooks
        END SUBROUTINE init_PY8
    END INTERFACE

!  Warning: the INTEGER would be converted from integer*8 -> integer*4.
!           But dont't worry, if we use PYTHIA 8. Because the storage of
!           the color flow was redesigned in PYTHIA 8.
    frameType_c = frameType
    idA_c = idA
    idB_c = idB
    eCM_c = eCM
    idStable = KF_woDecay

    MINT_c = MINT
    VINT_c = VINT
    MSTU_c = MSTU
    PARU_c = PARU
    MSTJ_c = MSTJ
    PARJ_c = PARJ
    MSTP_c = MSTP
    PARP_c = PARP
    MSTI_c = MSTI
    PARI_c = PARI

!   Universal usage of PACIAE, "3MOM" collision.
    nPY8 = 2
    do I=1,5,1
        pPY8(1,I) = P(1,I)
        pPY8(2,I) = P(2,I)
        vPY8(1,I) = V(1,I)
        vPY8(2,I) = V(2,I)
    END DO


!   Accesses to C++ program in Pythia8_cpp_interface.cpp.
    CALL init_PY8( frameType_c, idA_c, idB_c, eCM_c, &
                   MINT_c, VINT_c, &
                   MSTU_c, PARU_c, MSTJ_c, PARJ_c, &
                   MSTP_c, PARP_c, MSTI_c, PARI_c, &
                   nPY8, kPY8, pPY8, vPY8, &
                   pythia8, &
                   pythia8_pp, pythia8_pn, &
                   pythia8_np, pythia8_nn, &
                   pahooks, paHIhooks )


    ! MINT = MINT_c
    ! VINT = VINT_c
    ! MSTU = MSTU_c
    ! PARU = PARU_c
    ! MSTJ = MSTJ_c
    ! PARJ = PARJ_c
    ! MSTP = MSTP_c
    ! PARP = PARP_c
    ! MSTI = MSTI_c
    ! PARI = PARI_c


    RETURN
    END



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>  Interfaces into PYTHIA 8, Pythia::next(). Generates the core processes.
    SUBROUTINE PYEVNT_PY8

!---Imports mudules.
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

!---Fortran data type.

!***********************************************************************
!   Logical functions.
    LOGICAL IS_NUCLEUS
!***********************************************************************

!***********************************************************************
!   Local variabls.
    INTEGER, PARAMETER :: KSZJ = 80000, KSZJ_PY8 = 300000
    INTEGER :: I, J
    ! INTEGER :: frameType, idA, idB
    ! REAL(KIND=8) :: eCM
!***********************************************************************

!***********************************************************************
!---PACIAE.
!-----------------------------------------------------------------------
    INTEGER :: itorw, iikk, kkii
    REAL(KIND=8) :: cp0, cr0
    COMMON/SA34/ itorw, iikk, cp0, cr0, kkii
!-----------------------------------------------------------------------
    INTEGER :: i_mode, i_tune
    INTEGER :: KF_woDecay(100), KF_proj, KF_targ
    REAL(KIND=8) :: win, energy_B, psno, b_min, b_max
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay, &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!***********************************************************************

!***********************************************************************
!---PYTHIA 6.
!-----------------------------------------------------------------------
    INTEGER :: MINT(400)
    REAL(KIND=8) :: VINT(400)
    COMMON/PYINT1/ MINT, VINT
    SAVE /PYINT1/
!-----------------------------------------------------------------------
    INTEGER :: MSTU(200), MSTJ(200)
    REAL(KIND=8) :: PARU(200), PARJ(200)
    COMMON/PYDAT1/ MSTU, PARU, MSTJ, PARJ
    SAVE /PYDAT1/
!-----------------------------------------------------------------------
    INTEGER :: MSTP(200), MSTI(200)
    REAL(KIND=8) :: PARP(200), PARI(200)
    COMMON/PYPARS/ MSTP, PARP, MSTI, PARI
    SAVE /PYPARS/
!-----------------------------------------------------------------------
    INTEGER :: N, NPAD
    INTEGER :: K(KSZJ,5)
    REAL(KIND=8) :: P(KSZJ,5), V(KSZJ,5)
    COMMON/PYJETS/ N, NPAD, K, P, V
    SAVE /PYJETS/
!***********************************************************************

!***********************************************************************
!---PYTHIA 8 information storage for feeding-back.
!-----------------------------------------------------------------------
    INTEGER :: N_PY8, NPAD_PY8
    ! INTEGER, PARAMETER :: KSZJ_PY8 = 80000
    INTEGER :: K_PY8(KSZJ_PY8,8)
    REAL(KIND=8) :: P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
    COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, K_PY8, P_PY8, V_PY8
    SAVE /PYJETS_PY8/
!-----------------------------------------------------------------------
    INTEGER :: NO_JUNC_PY8, KIND_JUNC_PY8(KSZJ_PY8), &
               COL_JUNC_PY8(KSZJ_PY8,3), ENDC_JUNC_PY8(KSZJ_PY8,3), &
               STAT_JUNC_PY8(KSZJ_PY8,3)
    COMMON/PYJUNC_PY8/ NO_JUNC_PY8, KIND_JUNC_PY8, &
                       COL_JUNC_PY8, ENDC_JUNC_PY8, STAT_JUNC_PY8
    SAVE /PYJUNC_PY8/
!***********************************************************************

!---C++ data type.

!***********************************************************************
!   Local variabls.
    ! INTEGER(KIND=C_INT) :: frameType_c, idA_c, idB_c
    ! REAL(KIND=C_DOUBLE) :: eCM_c
    ! INTEGER(KIND=C_INT) :: idStable(100)
    ! INTEGER(KIND=C_INT) :: iFail
!***********************************************************************

!***********************************************************************
    INTEGER(KIND=C_INT) :: MINT_c(400)
    REAL(KIND=C_DOUBLE) :: VINT_c(400)
    COMMON/PYINT1_c/ MINT_c, VINT_c
    SAVE /PYINT1_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: MSTU_c(200), MSTJ_c(200)
    REAL(KIND=C_DOUBLE) :: PARU_c(200), PARJ_c(200)
    COMMON/PYDAT1_c/ MSTU_c, PARU_c, MSTJ_c, PARJ_c
    SAVE /PYDAT1_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: MSTP_c(200), MSTI_c(200)
    REAL(KIND=C_DOUBLE) :: PARP_c(200), PARI_c(200)
    COMMON/PYPARS_c/ MSTP_c, PARP_c, MSTI_c, PARI_c
    SAVE /PYPARS_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: nPY8, nPadPY8
    ! INTEGER(C_INT), PARAMETER :: kSZJ_c = 300000
    INTEGER(KIND=C_INT) :: kPY8(KSZJ_PY8,8)
    REAL(KIND=C_DOUBLE) :: pPY8(KSZJ_PY8,7), vPY8(KSZJ_PY8,5)
    COMMON/PYJETS_c/ nPY8, nPadPY8, kPY8, pPY8, vPY8
    SAVE /PYJETS_c/
!-----------------------------------------------------------------------
    INTEGER(C_INT) :: noJuncPY8, kindJuncPY8(KSZJ_PY8), &
                      colJuncPY8(KSZJ_PY8,3), endcJuncPY8(KSZJ_PY8,3), &
                      statJuncPY8(KSZJ_PY8,3)
    COMMON/PYJUNC_c/ noJuncPY8, kindJuncPY8, colJuncPY8, endcJuncPY8, &
                     statJuncPY8
    SAVE /PYJUNC_c/
!***********************************************************************

!---Pointer.

!***********************************************************************
!---NOTE: DO NOT TOUCH THESES POINTERS IF YOU DON'T KNOW WHAT THEY ARE !!!
!-----------------------------------------------------------------------
!---Pointers to PYTHIA8 obejects.
    TYPE(C_PTR) :: pythia8, pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
    COMMON/PYTHIA8_PTR/ pythia8, &
                        pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
    TYPE(C_PTR) :: localPythia8
!-----------------------------------------------------------------------
!---Pointers to PACIAE4 obejects.
    TYPE(C_PTR) :: pahooks, paHIhooks, &
                   pahooks_pp, pahooks_pn, pahooks_np, pahooks_nn
    COMMON/PACIAE4_PTR/ pahooks, paHIhooks, &
                        pahooks_pp, pahooks_pn, pahooks_np, pahooks_nn
    TYPE(C_PTR) :: localPahooks, localPaHIhooks
!***********************************************************************

!---Interfaces to C++.

!***********************************************************************
    INTERFACE
        SUBROUTINE next_PY8( MINT_c, VINT_c, &
                             MSTU_c, PARU_c, MSTJ_c, PARJ_c, &
                             MSTP_c, PARP_c, MSTI_c, PARI_c, &
                             nPY8, kPY8, pPY8, vPY8, &
                             noJuncPY8, kindJuncPY8, &
                             colJuncPY8, endcJuncPY8, statJuncPY8, &
                             localPythia8, &
                             localPahooks, localPaHIhooks ) &
                             BIND( C, NAME="next_PY8" )
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND=C_INT) :: MINT_c(*)
            REAL(KIND=C_DOUBLE) :: VINT_c(*)
            INTEGER(KIND=C_INT) :: MSTU_c(*), MSTJ_c(*)
            REAL(KIND=C_DOUBLE) :: PARU_c(*), PARJ_c(*)
            INTEGER(KIND=C_INT) :: MSTP_c(*), MSTI_c(*)
            REAL(KIND=C_DOUBLE) :: PARP_c(*), PARI_c(*)
            INTEGER(C_INT) :: nPY8
            INTEGER(C_INT) :: kPY8(*)
            REAL(C_DOUBLE) :: pPY8(*), vPY8(*)
            INTEGER(C_INT) :: noJuncPY8
            INTEGER(C_INT) :: kindJuncPY8(*)
            INTEGER(C_INT) :: colJuncPY8(*), endcJuncPY8(*), statJuncPY8(*)
            TYPE(C_PTR) :: localPythia8
            TYPE(C_PTR) :: localPahooks, localPaHIhooks
        END SUBROUTINE next_PY8
    END INTERFACE
!***********************************************************************


!  Warning: the INTEGER would be converted from integer*8 -> integer*4.
!           But dont't worry, if we use PYTHIA 8. Because the storage of
!           the color flow was redesigned in PYTHIA 8.
    ! frameType_c = frameType
    ! idA_c = idA
    ! idB_c = idB
    ! eCM_c = eCM
    ! idStable = KF_woDecay

    MINT_c = MINT
    VINT_c = VINT
    MSTU_c = MSTU
    PARU_c = PARU
    MSTJ_c = MSTJ
    PARJ_c = PARJ
    MSTP_c = MSTP
    PARP_c = PARP
    MSTI_c = MSTI
    PARI_c = PARI

!   Universal usage of PACIAE, "3MOM" collision.
    nPY8 = 2
    DO I=1,5,1
        pPY8(1,I) = P(1,I)
        pPY8(2,I) = P(2,I)
        vPY8(1,I) = V(1,I)
        vPY8(2,I) = V(2,I)
    END DO

!   Selects the instance.
    localPythia8 = pythia8
    localPahooks = pahooks
    localPaHIhooks = paHIhooks
!   For lA(Al), hA(Ah) and AB; typically pA(Ap) and AA collisions.
    IF( ( IS_NUCLEUS(KF_proj) .OR. IS_NUCLEUS(KF_targ) ) &
        .AND. ( MSTP(191) == 5 .OR. MSTP(191) == 6 ) )THEN
        ! AB; typically AA
        IF( IS_NUCLEUS(KF_proj) .AND. IS_NUCLEUS(KF_targ) )THEN
            IF(      MINT(11) == 2212 .AND. MINT(12) == 2212 )THEN
                 localPythia8 = pythia8_pp
            ELSE IF( MINT(11) == 2212 .AND. MINT(12) == 2112 )THEN
                 localPythia8 = pythia8_pn
            ELSE IF( MINT(11) == 2112 .AND. MINT(12) == 2212 )THEN
                 localPythia8 = pythia8_np
            ELSE IF( MINT(11) == 2112 .AND. MINT(12) == 2112 )THEN
                 localPythia8 = pythia8_nn
            END IF
        ! lA(Al), hA(hA); typically pA(Ap)
        ELSE
            IF(      MINT(11) == KF_proj .AND. MINT(12) == 2212 )THEN
                 localPythia8 = pythia8_pp
            ELSE IF( MINT(11) == KF_proj .AND. MINT(12) == 2112 )THEN
                 localPythia8 = pythia8_pn
            ELSE IF( MINT(11) == 2112 .AND. MINT(12) == KF_targ )THEN
                 localPythia8 = pythia8_np
            ELSE IF( MINT(11) == 2212 .AND. MINT(12) == KF_targ )THEN
                 localPythia8 = pythia8_nn
            END IF
        END IF
    END IF


!   Accesses to C++ program in Pythia8_cpp_interface.cpp.
    CALL next_PY8( MINT_c, VINT_c, &
                   MSTU_c, PARU_c, MSTJ_c, PARJ_c, &
                   MSTP_c, PARP_c, MSTI_c, PARI_c, &
                   nPY8, kPY8, pPY8, vPY8, &
                   noJuncPY8, kindJuncPY8, &
                   colJuncPY8, endcJuncPY8, statJuncPY8, &
                   localPythia8, &
                   localPahooks, localPaHIhooks )


    MINT = MINT_c
    VINT = VINT_c
    MSTU = MSTU_c
    PARU = PARU_c
    MSTJ = MSTJ_c
    PARJ = PARJ_c
    MSTP = MSTP_c
    PARP = PARP_c
    MSTI = MSTI_c
    PARI = PARI_c

!   Note: the first entry from PYTHIA8 is a dummy "particle" event[0],
!         which represents "system", corresponding to kPY8(1,*).
!   Its information:
!    " no        id  name            status     mothers   daughters     colours"
!    "  0        90  (system)           -11     0     0     0     0     0     0"

!   K1: status (KS in PYTHIA6); K2: id (KF in PYTHIA6);
!   K3: mother1; K4: daughter1; K5: daughter2.
!   K6: mother2; K7: color; K8: anti-color.
!   P1: px; P2: py; P3: pz; P4: e; P5: mass.
!   P6: scale; P7: polarization.
!   V1: x; V2: y; V3: z; V4: t; V5: tau, i.e. proper life time.

!   Stores partonic iformation from PYTHIA 8 for the feeding-back to
!    the hadronization of PYTHIA 8.
    ! Warning: from integer*4 -> integer*8.
    N_PY8 = nPY8
    DO J=1,5,1
        DO I=1,nPY8,1
            K_PY8(I,J) = kPY8(I,J)
            P_PY8(I,J) = pPY8(I,J)
            V_PY8(I,J) = vPY8(I,J)
        END DO
    END DO
    DO J=6,7,1
        DO I=1,nPY8,1
            K_PY8(I,J) = kPY8(I,J)
            P_PY8(I,J) = pPY8(I,J)
        END DO
    END DO
    DO I=1,nPY8,1
        K_PY8(I,8) = kPY8(I,8)
    END DO

!   Stores the junction configuration from PYTHIA 8 for the feeding-back to
!    the hadronization of PYTHIA 8.
    ! Warning: from integer*4 -> integer*8.
    NO_JUNC_PY8 = noJuncPY8
    DO I = 1, noJuncPY8, 1
        KIND_JUNC_PY8( I )  = kindJuncPY8( I )
    END DO
    DO J=1,3,1
        DO I = 1, noJuncPY8, 1
            COL_JUNC_PY8(  I, J ) =  colJuncPY8( I, J )
            ENDC_JUNC_PY8( I, J ) = endcJuncPY8( I, J )
            STAT_JUNC_PY8( I, J ) = statJuncPY8( I, J )
        END DO
    END DO


    RETURN
    END



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>  Interface into PYTHIA8, Pythia::forceHadronLevel().
    SUBROUTINE PYEXEC_PY8

!---Imports mudules.
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

!---Fortran data type.

!***********************************************************************
!   Logical functions.
    LOGICAL IS_NUCLEUS
!***********************************************************************

!***********************************************************************
!   Local variabls.
    INTEGER, PARAMETER :: KSZJ = 80000, KSZJ_PY8 = 300000
    INTEGER :: I, J
    ! INTEGER :: frameType, idA, idB
    ! REAL(KIND=8) :: eCM
!***********************************************************************

!***********************************************************************
!---PACIAE.
!-----------------------------------------------------------------------
    INTEGER :: itorw, iikk, kkii
    REAL(KIND=8) :: cp0, cr0
    COMMON/SA34/ itorw, iikk, cp0, cr0, kkii
!-----------------------------------------------------------------------
    INTEGER :: i_mode, i_tune
    INTEGER :: KF_woDecay(100), KF_proj, KF_targ
    REAL(KIND=8) :: win, energy_B, psno, b_min, b_max
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay, &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!***********************************************************************

!***********************************************************************
!---PYTHIA 6.
!-----------------------------------------------------------------------
    INTEGER :: MINT(400)
    REAL(KIND=8) :: VINT(400)
    COMMON/PYINT1/ MINT, VINT
    SAVE /PYINT1/
!-----------------------------------------------------------------------
    INTEGER :: MSTU(200), MSTJ(200)
    REAL(KIND=8) :: PARU(200), PARJ(200)
    COMMON/PYDAT1/ MSTU, PARU, MSTJ, PARJ
    SAVE /PYDAT1/
!-----------------------------------------------------------------------
    INTEGER :: MSTP(200), MSTI(200)
    REAL(KIND=8) :: PARP(200), PARI(200)
    COMMON/PYPARS/ MSTP, PARP, MSTI, PARI
    SAVE /PYPARS/
!-----------------------------------------------------------------------
    INTEGER :: N, NPAD
    INTEGER :: K(KSZJ,5)
    REAL(KIND=8) :: P(KSZJ,5), V(KSZJ,5)
    COMMON/PYJETS/ N, NPAD, K, P, V
    SAVE /PYJETS/
!***********************************************************************

!***********************************************************************
!---PYTHIA 8 information storage for feeding-back.
!-----------------------------------------------------------------------
    INTEGER :: N_PY8, NPAD_PY8
    ! INTEGER, PARAMETER :: KSZJ_PY8 = 80000
    INTEGER :: K_PY8(KSZJ_PY8,8)
    REAL(KIND=8) :: P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
    COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, K_PY8, P_PY8, V_PY8
    SAVE /PYJETS_PY8/
!-----------------------------------------------------------------------
    INTEGER :: NO_JUNC_PY8, KIND_JUNC_PY8(KSZJ_PY8), &
               COL_JUNC_PY8(KSZJ_PY8,3), ENDC_JUNC_PY8(KSZJ_PY8,3), &
               STAT_JUNC_PY8(KSZJ_PY8,3)
    COMMON/PYJUNC_PY8/ NO_JUNC_PY8, KIND_JUNC_PY8, &
                       COL_JUNC_PY8, ENDC_JUNC_PY8, STAT_JUNC_PY8
    SAVE /PYJUNC_PY8/
!***********************************************************************

!---C++ data type.

!***********************************************************************
!   Local variabls.
    ! INTEGER(KIND=C_INT) :: frameType_c, idA_c, idB_c
    ! REAL(KIND=C_DOUBLE) :: eCM_c
    ! INTEGER(KIND=C_INT) :: idStable(100)
    INTEGER(KIND=C_INT) :: iFail
!***********************************************************************

!***********************************************************************
    INTEGER(KIND=C_INT) :: MINT_c(400)
    REAL(KIND=C_DOUBLE) :: VINT_c(400)
    COMMON/PYINT1_c/ MINT_c, VINT_c
    SAVE /PYINT1_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: MSTU_c(200), MSTJ_c(200)
    REAL(KIND=C_DOUBLE) :: PARU_c(200), PARJ_c(200)
    COMMON/PYDAT1_c/ MSTU_c, PARU_c, MSTJ_c, PARJ_c
    SAVE /PYDAT1_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: MSTP_c(200), MSTI_c(200)
    REAL(KIND=C_DOUBLE) :: PARP_c(200), PARI_c(200)
    COMMON/PYPARS_c/ MSTP_c, PARP_c, MSTI_c, PARI_c
    SAVE /PYPARS_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: nPY8, nPadPY8
    ! INTEGER(C_INT), PARAMETER :: kSZJ_c = 300000
    INTEGER(KIND=C_INT) :: kPY8(KSZJ_PY8,8)
    REAL(KIND=C_DOUBLE) :: pPY8(KSZJ_PY8,7), vPY8(KSZJ_PY8,5)
    COMMON/PYJETS_c/ nPY8, nPadPY8, kPY8, pPY8, vPY8
    SAVE /PYJETS_c/
!-----------------------------------------------------------------------
    INTEGER(C_INT) :: noJuncPY8, kindJuncPY8(KSZJ_PY8), &
                      colJuncPY8(KSZJ_PY8,3), endcJuncPY8(KSZJ_PY8,3), &
                      statJuncPY8(KSZJ_PY8,3)
    COMMON/PYJUNC_c/ noJuncPY8, kindJuncPY8, colJuncPY8, endcJuncPY8, &
                     statJuncPY8
    SAVE /PYJUNC_c/
!***********************************************************************

!---Pointer.

!***********************************************************************
!---NOTE: DO NOT TOUCH THESES POINTERS IF YOU DON'T KNOW WHAT THEY ARE !!!
!-----------------------------------------------------------------------
!---Pointers to PYTHIA8 obejects.
    TYPE(C_PTR) :: pythia8, pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
    COMMON/PYTHIA8_PTR/ pythia8, &
                        pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
    TYPE(C_PTR) :: localPythia8
!-----------------------------------------------------------------------
!---Pointers to PACIAE4 obejects.
    TYPE(C_PTR) :: pahooks, paHIhooks, &
                   pahooks_pp, pahooks_pn, pahooks_np, pahooks_nn
    COMMON/PACIAE4_PTR/ pahooks, paHIhooks, &
                        pahooks_pp, pahooks_pn, pahooks_np, pahooks_nn
    TYPE(C_PTR) :: localPahooks, localPaHIhooks
!***********************************************************************

!---Interfaces to C++.

!***********************************************************************
    INTERFACE
        SUBROUTINE forceHadronLevel_PY8( MINT_c, VINT_c, &
                                         MSTU_c, PARU_c, MSTJ_c, PARJ_c, &
                                         MSTP_c, PARP_c, MSTI_c, PARI_c, &
                                         nPY8, kPY8, pPY8, vPY8, &
                                         noJuncPY8, kindJuncPY8, &
                                         colJuncPY8, endcJuncPY8, statJuncPY8, &
                                         iFail, localPythia8, &
                                         localPahooks, localPaHIhooks ) &
                                         BIND( C, NAME="forceHadronLevel_PY8" )
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND=C_INT) :: MINT_c(400)
            REAL(KIND=C_DOUBLE) :: VINT_c(400)
            INTEGER(KIND=C_INT) :: MSTU_c(200), MSTJ_c(200)
            REAL(KIND=C_DOUBLE) :: PARU_c(200), PARJ_c(200)
            INTEGER(KIND=C_INT) :: MSTP_c(200), MSTI_c(200)
            REAL(KIND=C_DOUBLE) :: PARP_c(200), PARI_c(200)
            INTEGER(C_INT) :: nPY8
            INTEGER(C_INT) :: kPY8(*)
            REAL(C_DOUBLE) :: pPY8(*), vPY8(*)
            INTEGER(C_INT) :: noJuncPY8
            INTEGER(C_INT) :: kindJuncPY8(*)
            INTEGER(C_INT) :: colJuncPY8(*), endcJuncPY8(*), statJuncPY8(*)
            INTEGER(KIND=C_INT) :: iFail
            TYPE(C_PTR) :: localPythia8
            TYPE(C_PTR) :: localPahooks, localPaHIhooks
        END SUBROUTINE forceHadronLevel_PY8
    END INTERFACE
!***********************************************************************


!  Warning: the INTEGER would be converted from integer*8 -> integer*4.
!           But dont't worry, if we use PYTHIA 8. Because the storage of
!           the color flow was redesigned in PYTHIA 8.
    ! frameType_c = frameType
    ! idA_c = idA
    ! idB_c = idB
    ! eCM_c = eCM
    ! frameType_c = MINT(111)
    ! idA_c = MINT(11)
    ! idB_c = MINT(12)
    ! eCM_c = VINT(290)
    ! idStable = KF_woDecay
    iFail = kkii

    MINT_c = MINT
    VINT_c = VINT
    MSTU_c = MSTU
    PARU_c = PARU
    MSTJ_c = MSTJ
    PARJ_c = PARJ
    MSTP_c = MSTP
    PARP_c = PARP
    MSTI_c = MSTI
    PARI_c = PARI

!   Note: the first entry from PYTHIA8 is a dummy "particle" event[0],
!         which represents "system", corresponding to kPY8(1,*).
!   Its information:
!    " no        id  name            status     mothers   daughters     colours"
!    "  0        90  (system)           -11     0     0     0     0     0     0"

!   K1: status (KS in PYTHIA6); K2: id (KF in PYTHIA6);
!   K3: mother1; K4: daughter1; K5: daughter2.
!   K6: mother2; K7: color; K8: anti-color.
!   P1: px; P2: py; P3: pz; P4: e; P5: mass.
!   P6: scale; P7: polarization.
!   V1: x; V2: y; V3: z; V4: t; V5: tau, i.e. proper life time.

!   Feeds back stored status, id, color etc. information from /PYJETS_PY8/
!    into PYTHIA 8.
    nPY8 = N_PY8
    DO J=1,5,1
        DO I=1,N_PY8,1
            kPY8(I,J) = K_PY8(I,J)
            pPY8(I,J) = P_PY8(I,J)
            vPY8(I,J) = V_PY8(I,J)
        END DO
    END DO
    DO J=6,7,1
        DO I=1,N_PY8,1
            kPY8(I,J) = K_PY8(I,J)
            pPY8(I,J) = P_PY8(I,J)
        END DO
    END DO
    DO I=1,N_PY8,1
        kPY8(I,8) = K_PY8(I,8)
    END DO

!   Feeds back stored junction configuration, from /PYJUNC_PY8/
!    into PYTHIA 8.
    ! Warning: from integer*8 -> integer*4
    noJuncPY8 = NO_JUNC_PY8
    DO I = 1, NO_JUNC_PY8, 1
        kindJuncPY8( I ) = KIND_JUNC_PY8( I )
    END DO
    DO J=1,3,1
        DO I = 1, NO_JUNC_PY8, 1
            colJuncPY8(  I, J ) =  COL_JUNC_PY8( I, J )
            endcJuncPY8( I, J ) = ENDC_JUNC_PY8( I, J )
            statJuncPY8( I, J ) = STAT_JUNC_PY8( I, J )
        END DO
    END DO

!   Selects the instance.
    localPythia8 = pythia8
    localPahooks = pahooks
    localPaHIhooks = paHIhooks
!   For lA(Al), hA(Ah) and AB; typically pA(Ap) and AA collisions.
    IF( ( IS_NUCLEUS(KF_proj) .OR. IS_NUCLEUS(KF_targ) ) &
        .AND. ( MSTP(191) == 5 .OR. MSTP(191) == 6 ) )THEN
        ! AB; typically AA
        IF( IS_NUCLEUS(KF_proj) .AND. IS_NUCLEUS(KF_targ) )THEN
            IF(      MINT(11) == 2212 .AND. MINT(12) == 2212 )THEN
                 localPythia8 = pythia8_pp
            ELSE IF( MINT(11) == 2212 .AND. MINT(12) == 2112 )THEN
                 localPythia8 = pythia8_pn
            ELSE IF( MINT(11) == 2112 .AND. MINT(12) == 2212 )THEN
                 localPythia8 = pythia8_np
            ELSE IF( MINT(11) == 2112 .AND. MINT(12) == 2112 )THEN
                 localPythia8 = pythia8_nn
            END IF
        ! lA(Al), hA(hA); typically pA(Ap)
        ELSE
            IF(      MINT(11) == KF_proj .AND. MINT(12) == 2212 )THEN
                 localPythia8 = pythia8_pp
            ELSE IF( MINT(11) == KF_proj .AND. MINT(12) == 2112 )THEN
                 localPythia8 = pythia8_pn
            ELSE IF( MINT(11) == 2112 .AND. MINT(12) == KF_targ )THEN
                 localPythia8 = pythia8_np
            ELSE IF( MINT(11) == 2112 .AND. MINT(12) == KF_targ )THEN
                 localPythia8 = pythia8_nn
            END IF
        END IF
    END IF


!   Accesses to C++ program in Pythia8_cpp_interface.cpp .
    CALL forceHadronLevel_PY8( MINT_c, VINT_c, &
                               MSTU_c, PARU_c, MSTJ_c, PARJ_c, &
                               MSTP_c, PARP_c, MSTI_c, PARI_c, &
                               nPY8, kPY8, pPY8, vPY8, &
                               noJuncPY8, kindJuncPY8, &
                               colJuncPY8, endcJuncPY8, statJuncPY8, &
                               iFail, localPythia8, &
                               localPahooks, localPaHIhooks )


!   Hadronization failed.
    kkii = iFail
    if( kkii == 2 ) return

!   Hadronization succeeded.

    MINT = MINT_c
    VINT = VINT_c
    MSTU = MSTU_c
    PARU = PARU_c
    MSTJ = MSTJ_c
    PARJ = PARJ_c
    MSTP = MSTP_c
    PARP = PARP_c
    MSTI = MSTI_c
    PARI = PARI_c

!   Stores hadronic iformation from PYTHIA 8.
    ! Warning: from integer*4 -> integer*8.
    N_PY8 = nPY8
    DO J=1,5,1
        DO I=1,nPY8,1
            K_PY8(I,J) = kPY8(I,J)
            P_PY8(I,J) = pPY8(I,J)
            V_PY8(I,J) = vPY8(I,J)
        END DO
    END DO
    DO J=6,7,1
        DO I=1,nPY8,1
            K_PY8(I,J) = kPY8(I,J)
            P_PY8(I,J) = pPY8(I,J)
        END DO
    END DO
    DO I=1,nPY8,1
        K_PY8(I,8) = kPY8(I,8)
    END DO


    RETURN
    END



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>  Lists event records.
    SUBROUTINE PYLIST_PY8( I_LIST )

!---Imports mudules.
    ! USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

    INTEGER :: I_LIST


    CALL PYLIST( I_LIST )


    RETURN
    END



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!>  Interfaces into C++ (PYTHIA 8), pythia.stat. Prints out cross sections
!!    information.
    SUBROUTINE PYSTAT_PY8(I_STAT)

!---Imports mudules.
    USE, INTRINSIC :: ISO_C_BINDING
    IMPLICIT NONE

!---Fortran data type.

!***********************************************************************
!   Local variabls.
    INTEGER :: I_STAT
!***********************************************************************

!***********************************************************************
!---PYTHIA 6.
!-----------------------------------------------------------------------
    INTEGER :: MINT(400)
    REAL(KIND=8) :: VINT(400)
    COMMON/PYINT1/ MINT, VINT
    SAVE /PYINT1/
!-----------------------------------------------------------------------
    INTEGER :: MSTU(200), MSTJ(200)
    REAL(KIND=8) :: PARU(200), PARJ(200)
    COMMON/PYDAT1/ MSTU, PARU, MSTJ, PARJ
    SAVE /PYDAT1/
!-----------------------------------------------------------------------
    INTEGER :: MSTP(200), MSTI(200)
    REAL(KIND=8) :: PARP(200), PARI(200)
    COMMON/PYPARS/ MSTP, PARP, MSTI, PARI
    SAVE /PYPARS/
!***********************************************************************

!---C++ data type.

!***********************************************************************
    INTEGER(KIND=C_INT) :: MINT_c(400)
    REAL(KIND=C_DOUBLE) :: VINT_c(400)
    COMMON/PYINT1_c/ MINT_c, VINT_c
    SAVE /PYINT1_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: MSTU_c(200), MSTJ_c(200)
    REAL(KIND=C_DOUBLE) :: PARU_c(200), PARJ_c(200)
    COMMON/PYDAT1_c/ MSTU_c, PARU_c, MSTJ_c, PARJ_c
    SAVE /PYDAT1_c/
!-----------------------------------------------------------------------
    INTEGER(KIND=C_INT) :: MSTP_c(200), MSTI_c(200)
    REAL(KIND=C_DOUBLE) :: PARP_c(200), PARI_c(200)
    COMMON/PYPARS_c/ MSTP_c, PARP_c, MSTI_c, PARI_c
    SAVE /PYPARS_c/
!-----------------------------------------------------------------------
!***********************************************************************

!---Pointer.

!***********************************************************************
!---NOTE: DO NOT TOUCH THESES POINTERS IF YOU DON'T KNOW WHAT THEY ARE !!!
!-----------------------------------------------------------------------
!---Pointers to PYTHIA8 obejects.
    TYPE(C_PTR) :: pythia8, pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
    COMMON/PYTHIA8_PTR/ pythia8, &
                        pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
!***********************************************************************

!---Interfaces to C++.

!***********************************************************************
    INTERFACE
        SUBROUTINE stat_PY8( MINT_c, MSTP_c, MSTI_c, &
                             pythia8, &
                             pythia8_pp, pythia8_pn, &
                             pythia8_np, pythia8_nn ) &
            BIND( C, NAME="stat_PY8" )
            USE, INTRINSIC :: ISO_C_BINDING
            INTEGER(KIND=C_INT) :: MINT_c(*), MSTP_c(*), MSTI_c(*)
            TYPE(C_PTR) :: pythia8, &
                           pythia8_pp, pythia8_pn, pythia8_np, pythia8_nn
        END SUBROUTINE stat_PY8
    END INTERFACE
!***********************************************************************

    MINT_c = MINT
    VINT_c = VINT
    MSTU_c = MSTU
    PARU_c = PARU
    MSTJ_c = MSTJ
    PARJ_c = PARJ
    MSTP_c = MSTP
    PARP_c = PARP
    MSTI_c = MSTI
    PARI_c = PARI

!   Accesses to C++ program in Pythia8_cpp_interface.cpp .
    CALL stat_PY8( MINT_c, MSTP_c, MSTI_c, &
                   pythia8, &
                   pythia8_pp, pythia8_pn, &
                   pythia8_np, pythia8_nn )


    RETURN
    END



!*******************************************************************************
!*******************************************************************************
!240124 Lei
!
!       A series of general functions and subroutines for the execution with
!        both PYTHIA 6 and 8.
!
!                                            Written by Anke at CCNU on 01/2024.
!
!===============================================================================
!       Free variables in PYTHIA 6 COMMONs:
!         MINT(355) - MINT(400)
!         VINT(360) - VINT(400)
!         MSTU(163) - MSTU(200)
!         PARU(196) - PARU(200)
!         MSTJ(122) - MSTJ(200)
!         PARJ(196) - PARJ(200)
!         MSTP(186) - MSTP(200)
!         PARP(195) - PARP(200)
!         MSTI(73)  - MSTI(200)
!         PARI(115) - PARI(200)
!
!===============================================================================
!       If PYTHIA 8 or Angantyr mode was activated:
!-------------------------------------------------------------------------------
!         Parameters/values to be input into PYTHIA 8:
!--------------------
!         MINT & VINT
!           MINT(5): iii
!           MINT(11): KF_BEAM of each PAINIT call.
!           MINT(12): KF_TARGET of each PAINIT call.
!           MINT(111): I_FRAME of each PAINIT call.
!           MINT(39): psno
!           MINT(41): ipden
!           MINT(42): itden
!           MINT(355): i_NN_pair
!
!           VINT(139): bp
!           VINT(290): win
!--------------------
!         MSTU & PARU
!           MSTU(21): mstu21   ! Not used temporarily.
!
!           PARU
!--------------------
!         MSTJ & PARJ
!           MSTJ(1): kjp22
!
!           PARJ(21): adj1(34)
!           PARJ(41): adj1(6)
!           PARJ(42): adj1(7)
!--------------------
!         MSTP & PARP
!           MSTP(5): i_tune
!           MSTP(95): i_color_reconnection
!           MSTP(111): mstptj
!           MSTP(186): adj1(12)
!           MSTP(187): i_deex
!           MSTP(191): i_mode
!           MSTP(192): i_seed
!           MSTP(194): adj1(5)
!           MSTP(199): nout
!           MSTP(195): neve
!
!           PARP(31): adj1(10)
!--------------------
!         MSTI & PARI
!           MSTI(11): KF_proj of the entire simulation system
!           MSTI(12): KF_targ of the entire simulation system
!           MSTI(111): ifram of the entire simulation system
!
!           PARI
!-------------------------------------------------------------------------------
!         Parameters/values to be returned from PYTHIA 8:
!--------------------
!         MINT & VINT
!           MINT(385): nCollTot     from Angantyr
!           MINT(386): nCollND      from Angantyr
!           MINT(387): deprecated from PYTHIA 8.313
!           MINT(388): nCollSDP     from Angantyr
!           MINT(389): nCollSDT     from Angantyr
!           MINT(390): nCollDD      from Angantyr
!           MINT(391): nCollCD      from Angantyr
!           MINT(392): nCollEL      from Angantyr
!           MINT(393): nPartProj    from Angantyr
!           MINT(394): nAbsProj     from Angantyr
!           MINT(395): nDiffProj    from Angantyr
!           MINT(396): nElProj      from Angantyr
!           MINT(397): nPartTarg    from Angantyr
!           MINT(398): nAbsTarg     from Angantyr
!           MINT(399): nDiffTarg    from Angantyr
!           MINT(400): nElTarg      from Angantyr
!
!           VINT(98): weightSum
!           VINT(99): 1/weight
!           VINT(100): weight
!           VINT(139): bp (the impact parameter         from Angantyr)
!           VINT(370): the impact parameter angle phi   from Angantyr
!           VINT(371): glauberTot                       from Angantyr
!           VINT(372): glauberND                        from Angantyr
!           VINT(373): glauberINEL                      from Angantyr
!           VINT(374): glauberEL                        from Angantyr
!           VINT(375): glauberDiffP                     from Angantyr
!           VINT(376): glauberDiffT                     from Angantyr
!           VINT(377): glauberDDiff                     from Angantyr
!           VINT(378): glauberBSlope                    from Angantyr
!           VINT(381): glauberTotErr                    from Angantyr
!           VINT(382): glauberNDErr                     from Angantyr
!           VINT(383): glauberINELErr                   from Angantyr
!           VINT(384): glauberELErr                     from Angantyr
!           VINT(385): glauberDiffPErr                  from Angantyr
!           VINT(386): glauberDiffTErr                  from Angantyr
!           VINT(387): glauberDDiffErr                  from Angantyr
!           VINT(388): glauberBSlopeErr                 from Angantyr
!           VINT(396): version number of PYTHIA 8.
!           VINT(397): total px of the collision system
!           VINT(398): total py of the collision system
!           VINT(399): total pz of the collision system
!           VINT(400): total E of  the collision system
!--------------------
!         MSTU & PARU
!           MSTU(23): error counter   ! Not used temporarily.
!           MSTU(30): error counter with junctions   ! Not used temporarily.
!
!           PARU
!--------------------
!         MSTJ & PARJ
!           MSTJ
!
!           PARJ
!--------------------
!         MSTP & PARP
!           MSTP
!
!           PARP
!--------------------
!         MSTI & PARI
!           MSTI
!
!           PARI(1): sigmaGen
!           PARI(2): sigmaGen / weightSum
!           PARI(9): 1/weight
!           PARI(10): weight
!       MINT(360): number of NN pairs. !#TODO(Lei20240220): use or not?
!===============================================================================



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PATUNE
!!  This subroutine was used to set basic parameters and random
!!   number seed for PYTHIA.
!!  Parameters will be set from COMMON blocks of
!!   PYTHIA 6 -> PYINIT_PY8 -> init_PY8 -> PYTHIA 8.
!
!   PY6: PYTHIA 6;
!   PY8: PYTHIA 8;
!   ANG: PYTHIA8/Angantyr;
!   PA: PACIAE;
!   P2011: Perugia 2011 Tune;
!   M2013: Monash 2013 Tune;
!   W2024: PA-Wuhan 2024 Tune (...waiting, by Anke).
!
!   Select basic FSR/ISR/UE parameter set = "tune".
!
!   By default, "Perugia 2011 Tune" (350) for PYTHIA 6,
!    and "Monash 2013 tune" (7 for e+e-, 14 for pp) in PYTHIA 8 /
!    Angantyr, respectively.
!
!   i_tune: = 0, select basic "Perugia 2011 Tune" or "Monash 2013 tune",
!           > 0, select corresponding tune set of "i_tune", i.e. MSTP(5).
!
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    INTEGER PYK,PYCHGE,PYCOMP
    LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_PARTON,IS_QUARK,IS_DIQUARK
    COMMON/PYINT1/MINT(400),VINT(400)
    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    COMMON/PYDATR/MRPY(6),RRPY(100)
    common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
    common/sa10/csnn,cspin,cskn,cspipi,cspsn,cspsm,rcsit,ifram, &
     iabsb,iabsm,i_sigma_AQM,ajpsi,csspn,csspm,csen
    common/sa12/ppsa(5),nchan,nsjp,sjp,ttaup,taujp
    common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
    common/sa24/adj1(40),nnstop,non24,zstop
    common/sa27/itime,kjp22,gtime,astr,akapa(6),parj1,parj2,parj3, &
     parj21,parj4,adiv,gpmax,nnc
    common/sa29/i_color_reconnection,NCR,parp2
    common/sa34/itorw,iikk,cp0,cr0,kkii
    common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
     nap,nat,nzp,nzt,pio
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
    CHARACTER*100 PARAM_PYTHIA6, PACIAE_INPUT


!-------------------------------------------------------------------------------
!-------------------------------- Tune Setting ---------------------------------
!   B- and C-framework via PYTHIA 6.
    IF( i_mode /= 1 )THEN
        IF( .NOT.IS_PYTHIA8(i_mode) .AND. i_tune > 0 )THEN
            CALL PYTUNE( i_tune )
        ELSE
            ! "Perugia 2011 Tune"
            CALL PYTUNE( 350 )
        END IF
!       Note "Monash 2013 Tune" was set automatically in PYTHIA 8.
    END IF
!-------------------------------- Tune Setting ---------------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------   User Parameters Setting   -------------------------
    IF( i_mode /= 1 )THEN   ! If i_mode
        ! For PYTHIA 6. Using "PYGIVE" prompts the changes.
        IF( .NOT.IS_PYTHIA8(I_MODE) )THEN
            ! K factor for the hard process and the MPI (uderlying events).
            ! MSTP(33) = 1
            CALL PYGIVE( "MSTP(33)=1" )
            ! PARP(31) = adj1(10)
            WRITE( PACIAE_INPUT, "(F11.5)" ) adj1(10)
            PARAM_PYTHIA6 = "PARP(31)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
            CALL PYGIVE( PARAM_PYTHIA6 )
            ! Parameters in Lund string fragmentation function.
            ! PARJ(41) = adj1(6)
            WRITE( PACIAE_INPUT, "(F11.5)" ) adj1(6)
            PARAM_PYTHIA6 = "PARP(41)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
            CALL PYGIVE( PARAM_PYTHIA6 )
            ! PARJ(42) = adj1(7)
            WRITE( PACIAE_INPUT, "(F11.5)" ) adj1(7)
            PARAM_PYTHIA6 = "PARP(42)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
            CALL PYGIVE( PARAM_PYTHIA6 )
            ! PARJ(21) = adj1(34)
            WRITE( PACIAE_INPUT, "(F11.5)" ) adj1(34)
            PARAM_PYTHIA6 = "PARJ(21)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
            CALL PYGIVE( PARAM_PYTHIA6 )
            ! Keeps conservation in the independent fragmaentation of PY6.
            ! MSTJ(3) = 4
            CALL PYGIVE( "MSTJ(3)=4" )
            ! Selection of CR models.
            ! MSTP(95) = i_color_reconnection
            WRITE( PACIAE_INPUT, "(I5)" ) i_color_reconnection
            PARAM_PYTHIA6 = "MSTP(95)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
            CALL PYGIVE( PARAM_PYTHIA6 )
            ! Set PARP(85)=0D0 to avoid "double count" CR when use old PYEVNT.
            IF( MSTP(81) <= 1 .AND. MSTP(95) > 0 ) &
                CALL PYGIVE( "PARP(85)=0D0" )

        ! For PYTHIA 8.
        ELSE IF( IS_PYTHIA8(i_mode) )THEN
            MSTP(5)   = i_tune
            ! Default MPI-CR.
            MSTP(95)  = 0
            MSTP(191) = i_mode
            ! MRPY(1) is the random number seed, will be set in PASEED.
            ! MSTP(192) = MRPY(1)
            ! nPDF.
            MSTP(194) = INT( adj1(5) )
            ! Without nPDF.
            IF( INT( adj1(5) ) ==  0 ) MSTP(194) = -1
            ! Only Isospin effect.
            IF( INT( adj1(5) ) == -1 ) MSTP(194) = 0
            ! MSTP(195) =
            ! MSTP(196) =
            ! MSTP(197) =
            ! MSTP(198) =
            MSTP(199) = nout
            MSTP(200) = neve
            MINT(39)  = INT( psno )
            MINT(1)   = nchan
            MINT(41)  = ipden
            MINT(42)  = itden
            MINT(355) = 1
            MSTI(11)  = KF_proj
            MSTI(12)  = KF_targ
            MSTI(111) = ifram
            ! K factor for the hard process and the MPI (uderlying events).
            PARP(31) = adj1(10)
            ! Parameters in Lund string fragmentation function.
            PARJ(41) = adj1(6)
            PARJ(42) = adj1(7)
            PARJ(21) = adj1(34)
            ! Chooses fragmentation subclass for PYTHIA 8 and PYTHIA8/Angantyr.
            MSTJ(1) = kjp22
            ! Selection of CR models.
            MSTP(95) = 0
            IF( i_color_reconnection == 0 ) MSTP(95) = -1
            ! 1, 2, 3, 4.
            IF( i_color_reconnection >= 15 ) &
            MSTP(95) = i_color_reconnection - 14
            ! For the string-melting.
            MSTP(186) = INT( adj1(12) )
            MSTP(187) = i_deex
        END IF
    END IF   ! IF i_mode
!-------------------------   User Parameters Setting   -------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!------------------------   Random Number Seed Giving   ------------------------
    ! In Pythia8_fort_interface.f90
    CALL PASEED
!------------------------   Random Number Seed Giving   ------------------------
!-------------------------------------------------------------------------------


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PASEED
!!  Generate random number seed for PYTHIA 6 and 8.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    INTEGER PYK,PYCHGE,PYCOMP
    LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_PARTON,IS_QUARK,IS_DIQUARK
    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    COMMON/PYDATR/MRPY(6),RRPY(100)
    common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
    common/sa18/i_deex,n_deex_step,i_pT_coal,i_pT_endpoint,a_FF,aPS_c,aPS_b
    common/sa24/adj1(40),nnstop,non24,zstop
    common/sa29/i_color_reconnection,NCR,parp2
    common/sa34/itorw,iikk,cp0,cr0,kkii
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!   MRPY(1) is the seed of PYTHIA random number generator. (D=20240116)
!   For the intrinsic subroutine DATE_AND_TIME call.
    DIMENSION n_current_date_and_time(8)
    ! character*4 c_date_and_time(8)
!   Sets default seed as 20240116.
    DATA i_seed /20240116/
    CHARACTER*100 PARAM_PYTHIA6, PACIAE_INPUT


!   Note it is an array. One doesn't need to initialize it using do-loop.
    n_current_date_and_time = 0

    i_seed_mode = INT( adj1(26) )
!   Sets default seed as 20240116.
    IF( i_seed_mode <= 0 )THEN
        ! MRPY(1) = i_seed
        WRITE( PACIAE_INPUT, "(I20)" ) i_seed
        PARAM_PYTHIA6 = "MRPY(1)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
        CALL PYGIVE( PARAM_PYTHIA6 )
!   The initial seed of the random number generator will be given by
!    real-time clock from computer.
    ELSE IF( i_seed_mode == 1 )THEN
!   Gets the current real-date_and_time from the computer.
        ! Fortran intrinsic function.
        CALL DATE_AND_TIME( VALUES = n_current_date_and_time )
!   Only use hour-to-milliseconds because of the limitation of INTEGER type.
!            n_current_date_and_time(1)   ! year   (CCYY)
!            n_current_date_and_time(2)   ! month  (1-12)
!            n_current_date_and_time(3)   ! day    (1-31)
!            n_current_date_and_time(4)   ! The time difference, in minutes,
!                                         !  with respect to UTC.
!            n_current_date_and_time(5)   ! hour   (1-23)
!            n_current_date_and_time(6)   ! minute (1-59)
!            n_current_date_and_time(7)   ! second (0-60)
!            n_current_date_and_time(8)   ! milliseconds (0-999)
        i_seed = n_current_date_and_time(5) * 10000000 + &
                 n_current_date_and_time(6) * 100000   + &
                 n_current_date_and_time(7) * 1000     + &
                 n_current_date_and_time(8) * 1
!   The seed will be 9-digit HH-MM-SS-ms-ms-ms.
        ! MRPY(1) = i_seed
        WRITE( PACIAE_INPUT, "(I20)" ) i_seed
        PARAM_PYTHIA6 = "MRPY(1)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
        CALL PYGIVE( PARAM_PYTHIA6 )
    ELSE IF( i_seed_mode > 1 )THEN
        ! MRPY(1) = i_seed_mode
        WRITE( PACIAE_INPUT, "(I20)" ) i_seed_mode
        PARAM_PYTHIA6 = "MRPY(1)=" // TRIM(ADJUSTL( PACIAE_INPUT ))
        CALL PYGIVE( PARAM_PYTHIA6 )
    END IF
!   If PY8/Angantyr mode was activated, use MSTP(192)/MRPY(1) as seed in PY8.
    IF( IS_PYTHIA8(i_mode) )THEN
        MSTP(192) = MRPY(1)
        ! MSTP(192) = i_seed_mode
    END IF


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PAINST
!!  Instantiates Pythia objects for PYTHIA 8, i.e. allocates memory on heap.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    INTEGER PYK,PYCHGE,PYCOMP
    LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_PARTON,IS_QUARK,IS_DIQUARK
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


    ! In Pythia8_fort_interface.f90.
    IF( IS_PYTHIA8(i_mode) ) CALL PYINST_PY8


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PADELE
!!  Deletes the Pythia object from PYTHIA 8, i.e. releases memory from heap.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    INTEGER PYK,PYCHGE,PYCOMP
    LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_PARTON,IS_QUARK,IS_DIQUARK
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


    ! In Pythia8_fort_interface.f90.
    IF( IS_PYTHIA8(i_mode) ) CALL PYDELE_PY8


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PAINIT(FRAME,BEAM,TARGET,WIN)
!!  Identification of frame, beam, target and initilization for
!!   both PYTHIA 6 and 8.
!!  Note that "FRAME", "BEAM" and "TARGET" are CHARACTER type.
!!  "WIN" is DOUBLE PRECISION.
!
!   Part of the statement was borrowed from "SUBROUTINE PYINBM".
!
!   This subroutine can accept:
!    (1) For PYTHIA 6 mode, the beam and targert of the hadrons, leptons,
!        reggeon and pomeron defined in the following "CHCDE" before "kl0".
!    (2) For PYTHIA 8 and PYTHIA8/Angantyr modes, the beam and targert of
!        the particles defined in the following "CHCDE",
!        including ions "197Au", "208Pb", etc. If one wants to simulate
!        the user-defiend ion collisions, please input CHCDE(72) : "myIon"
!        and change the following KCDE(72) as corresponding PDG code.
!
!   The standard PDG (hyper)nucleus codes are 10-digit:
!           +-10LZZZAAAI.
!   +- denotes (hyper)nucleus and anti-one.
!   L = n_Lambda, total number of strange quarks (Lambda hyperons).
!   Z = n_p, total charges (protons).
!   A = n_p + n_n + n_Lambda, total baryon number.
!       n_n, total number of neutrons.
!   I = isomer level, = 0: groud state,
!                     > 0: excitations. m, n, p, q = 1, 2, 3, 4.
!   Cf. PDG "Monte Carlo particle numbering scheme" for more details.
!   Here for normal ions, we use
!       id = 100ZZZAAA0 = 1,000,000,000 + Z*10,000 + A*10 + 0.
!   For example, deuteron D/2H: 1000010020; 1000020040 = 4He;
!                        197Au: 1000791970; 208Pb: 1000822080; etc.
!
!   It is worth mentioning that, there is another "overload" "KF" version
!    named "PAINIT_KF". It is more flexible for PYTHIA 8.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    INTEGER PYK,PYCHGE,PYCOMP
    LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_PARTON,IS_QUARK,IS_DIQUARK
    PARAMETER (KSZJ=80000)
    COMMON/PYINT1/MINT(400),VINT(400)
    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
    common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
    common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
     nap,nat,nzp,nzt,pio
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win_o, energy_B, psno, b_min, b_max
!   Local arrays and character variables.
    CHARACTER*(*) FRAME,BEAM,TARGET
    CHARACTER CHFRAM*12,CHBEAM*12,CHTARG*12,CHCOM(3)*12,CHALP(2)*26, &
    CHIDNT(3)*12,CHTEMP*12,CHCDE(73)*12,CHNAME*16
    DIMENSION LEN(3), KCDE(73)
    DATA CHALP/'abcdefghijklmnopqrstuvwxyz', &
               'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
    DATA CHCDE/    "e-          ","e+          ","nu_e        ", &
    "nu_ebar     ","mu-         ","mu+         ","nu_mu       ", &
    "nu_mubar    ","tau-        ","tau+        ","nu_tau      ", &
    "nu_taubar   ","pi+         ","pi-         ","n0          ", &
    "nbar0       ","p+          ","pbar-       ","gamma       ", &
    "lambda0     ","sigma-      ","sigma0      ","sigma+      ", &
    "xi-         ","xi0         ","omega-      ","pi0         ", &
    "reggeon     ","pomeron     ","gamma/e-    ","gamma/e+    ", &
    "gamma/mu-   ","gamma/mu+   ","gamma/tau-  ","gamma/tau+  ", &
    "k+          ","k-          ","ks0         ","kl0         ", &
    "2H          ","2Hbar       ","4He         ","4Hebar      ", &
    "6Li         ","6Libar      ","9Be         ","9Bebar      ", &
    "12C         ","12Cbar      ","14N         ","14Nbar      ", &
    "16O         ","16Obar      ","27Al        ","27Albar     ", &
    "40Ar        ","40Arbar     ","56Fe        ","56Febar     ", &
    "63Cu        ","63Cubar     ","84Kr        ","84Krbar     ", &
    "107Ag       ","107Agbar    ","129Xe       ","129Xebar    ", &
    "197Au       ","197Aubar    ","208Pb       ","208Pbbar    ", &
    "myIon       ","null        "/
    DATA KCDE/11,-11,12,-12,13,-13,14,-14,15,-15,16,-16, &
              211,-211,2112,-2112,2212,-2212,22,3122,3112,3212,3222, &
              3312,3322,3334,111,110,990,6*22,321,-321,310,130, &
              1000010020, -1000010020, 1000020040, -1000020040, &
              1000030060, -1000030060, 1000040090, -1000040090, &
              1000060120, -1000060120, 1000070140, -1000070140, &
              1000080160, -1000080160, 1000130270, -1000130270, &
              1000180400, -1000180400, 1000260560, -1000260560, &
              1000290630, -1000290630, 1000360840, -1000360840, &
              1000471070, -1000471070, 1000541290, -1000541290, &
              1000791970, -1000791970, 1000822080, -1000822080, &
              1009999990, 0/
!   Local arrays and character variables.
    INTEGER frameType
!   KF_GAMMA has not been used to PYTHIA 8 yet.
    INTEGER KF_BEAM(2), KF_GAMMA(2)
    CHARACTER*100 PARAM_PYTHIA6
    DATA I_CALL /0/


    I_CALL = I_CALL + 1


!-------------------------------------------------------------------------------
!-------------------------   PYTHIA 6 Initialization   -------------------------
!   Reads the extra configuration for PYTHIA 6.
    OPEN( 11, FILE="pythia6_extra.cfg", STATUS="UNKNOWN" )
    DO WHILE( .TRUE. )
        READ( 11, "(A)", END=100 ) PARAM_PYTHIA6
        CALL PYGIVE( PARAM_PYTHIA6 )
    END DO
100 CLOSE(11)
    IF( .NOT.IS_PYTHIA8(i_mode) )THEN
!   Prohibits writing of information on variable values changed by a PYGIVE
!     and a PYTUNE.
        IF( I_CALL >= 1 ) CALL PYGIVE( "MSTU(13)=0" )

!   Core initialization of PYTHIA 6.
        CALL PYINIT( FRAME, BEAM, TARGET, WIN )
!-------------------------   PYTHIA 6 Initialization   -------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!--------------------   PYTHIA 8 / Angantyr Initialization   -------------------
    ELSE IF( IS_PYTHIA8(i_mode) )THEN

!----------------------------   Beam Identifying   -----------------------------
!   Identifies beam and target particles and frame of process.
        CHFRAM = TRIM( ADJUSTL(FRAME) )  // ' '
        CHBEAM = TRIM( ADJUSTL(BEAM) )   // ' '
        CHTARG = TRIM( ADJUSTL(TARGET) ) // ' '

!   Converts character variables to lowercase and find their length.
        CHCOM(1) = CHFRAM
        CHCOM(2) = CHBEAM
        CHCOM(3) = CHTARG
        DO I=1,3
            LEN(I) = 12
            DO LL=12,1,-1
                IF( LEN(I) == LL .AND. CHCOM(I)(LL:LL) == ' ' ) &
                LEN(I) = LL-1
                DO LA=1,26
                    IF( CHCOM(I)(LL:LL) == CHALP(2)(LA:LA) ) &
                    CHCOM(I)(LL:LL) = CHALP(1)(LA:LA)
                END DO
            END DO
            CHIDNT(I) = CHCOM(I)
!   Fixes up bar, underscore and charge in particle name (if needed).
            DO LL=1,10
                IF( CHIDNT(I)(LL:LL) == '~' )THEN
                    CHTEMP = CHIDNT(I)
                    CHIDNT(I) = CHTEMP(1:LL-1) // 'bar' // &
                                CHTEMP(LL+1:10) // '  '
                END IF
            END DO
            IF(CHIDNT(I)(1:2) == 'nu'.AND.CHIDNT(I)(3:3) /= '_')THEN
                CHTEMP=CHIDNT(I)
                CHIDNT(I)='nu_'//CHTEMP(3:7)
            ELSE IF( CHIDNT(I)(1:2) == 'n ' )THEN
                CHIDNT(I)(1:3)='n0 '
            ELSE IF( CHIDNT(I)(1:4) == 'nbar' )THEN
                CHIDNT(I)(1:5)='nbar0'
            ELSE IF( CHIDNT(I)(1:2) == 'p ' )THEN
                CHIDNT(I)(1:3)='p+ '
            ELSE IF( CHIDNT(I)(1:4) == 'pbar'.OR. &
                    CHIDNT(I)(1:2) == 'p-' )THEN
                CHIDNT(I)(1:5)='pbar-'
            ELSE IF( CHIDNT(I)(1:6) == 'lambda' )THEN
                CHIDNT(I)(7:7)='0'
            ELSE IF( CHIDNT(I)(1:3) == 'reg' )THEN
                CHIDNT(I)(1:7)='reggeon'
            ELSE IF( CHIDNT(I)(1:3) == 'pom' )THEN
                CHIDNT(I)(1:7)='pomeron'
            END IF
        END DO

!   Identifies incoming beam and target particles.
!   Translates to KF codes.
        DO I=1,2
            DO J=1,73
                IF( CHIDNT(I+1) == CHCDE(J) ) KF_BEAM(I) = KCDE(J)
            END DO
            KF_GAMMA(I) = 0
            IF( KF_BEAM(I) == 22 .AND. CHIDNT(I+1)(6:6) == '/' )THEN
                CHTEMP = CHIDNT(I+1)(7:12)//' '
                DO J=1,12
                    IF( CHTEMP == CHCDE(J) ) KF_GAMMA(I) = KCDE(J)
                END DO
            END IF
        END DO
        IF( KF_BEAM(1) == 0) WRITE(MSTU(11),5000) CHBEAM(1:LEN(2))
        IF( KF_BEAM(2) == 0) WRITE(MSTU(11),5100) CHTARG(1:LEN(3))
        IF( KF_BEAM(1) == 0 .OR. KF_BEAM(2) == 0) STOP

!   Identifies the choice of frame for PYTHIA 8. Defalut "CMS", i.e. "1".
        i_frame_PY8 = 1
!   Events defined in the CM frame.
        IF( CHCOM(1)(1:2) == 'cm'  )THEN
            i_frame_PY8 = 1
!   Events defined in fixed target frame.
        ELSE IF( CHCOM(1)(1:3) == 'fix'  )THEN
            i_frame_PY8 = 2
!   Frame defined by user three-vectors.
        ELSE IF(CHCOM(1)(1:1) == '3' )THEN
            i_frame_PY8 = 3
!   Frame defined by user four-vectors. There is no "4MOM" in PYTHIA 8.
        ELSE IF(CHCOM(1)(1:1) == '4' )THEN
            i_frame_PY8 = 3
!   Frame defined by user five-vectors. There is no "5MOM" in PYTHIA 8.
        ELSE IF(CHCOM(1)(1:1) == '5' )THEN
            i_frame_PY8 = 3
        END IF

!   Initializes the error counters.
        MSTU(23) = 0
        MSTU(30) = 0
!   Records the frame of the collision of each "PAINIT" call.
        MINT(111) = i_frame_PY8
!   Records the frame of the entire simulation system.
        ! MSTI(111) = iframe
!   Records the collision energy.
        VINT(290) = WIN
!   Records the type of the collision.
        MINT(41) = ipden
        MINT(42) = itden
!   Records the number of events.
        MINT(5) = iii
!   Specifies whether or not to hadronize.
        ! MSTP(111) = mstptj
!   Specifies the b parameter sampling method for NA(AN) and AA, if Angantyr.
        ! MINT(39) = INT( psno )
!   Records id/KF codes of projectile and target of each "PAINIT" call.
        MINT(11) = KF_BEAM(1)
        MINT(12) = KF_BEAM(2)
!   Records id/KF codes of proj. and targ. of the entire simulation system.
        ! MSTI(11) = KF_proj
        ! MSTI(12) = KF_targ

!   Formats for initialization and error information.
 5000   FORMAT( 1X,"Error: unrecognized beam particle in PAINIT ",A, &
                "D0"/1X, "Execution stopped!" )
 5100   FORMAT( 1X,"Error: unrecognized target particle in PAINIT ", &
                A, "D0"/1X, "Execution stopped!" )
!----------------------------   Beam Identifying   -----------------------------

!-----------------------------   Initialization   ------------------------------
        frameType = i_frame_PY8
        idA = KF_BEAM(1)
        idB = KF_BEAM(2)
        eCM = WIN
        ! In Pythia8_fort_interface.f90.
        CALL PYINIT_PY8( frameType, idA, idB, eCM )
!-----------------------------   Initialization   ------------------------------

    END IF
!--------------------   PYTHIA 8 / Angantyr Initialization   -------------------
!-------------------------------------------------------------------------------


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PAINIT_KF(I_FRAME,KF_BEAM,KF_TARGET,WIN)
!!  The "overload" "KF" version of "PAINIT".
!!  Identification of frame, beam, target and initilization for
!!   both PYTHIA 6 and 8.
!!  Note that "I_FRAME", "KF_BEAM" and "KF_TARGET" are INTEGER type.
!!  "WIN" is DOUBLE PRECISION.
!
!   This subroutine is more flexible than "PAINIT".
!
!   This subroutine can accept:
!    (1) For PYTHIA 6 mode, particles defined in "PYINIT".
!             +-11 = e-+;   +-12 = nu_e,   nu_ebar;
!             +-13 = mu-+;  +-14 = nu_mu,  nu_mubar;
!             +-15 = tau-+; +-16 = nu_tau, nu_taubar;
!               22 = gamma (on-shell);
!                ? = gamma/e-+, gamma/mu-+, gamma/tau-+
!            +-211 = pi+-; 111 = pi0;
!            +-321 = K+-;  310 = KS0; 120 = KL0;
!           +-2212 = p+-;   +-2112 = n0, n0bar;
!             3122 = Lambda0; 3112 = Sigma-; 3212 = Sigma0; 3222 = Sigma+;
!             3312 = Xi-;     3322 = Xi0;    3334 = Omega-;
!              110 = reggeon;  990 = pomeron.
!    (2) For PYTHIA 8 mode, almost all of subatomic particles.
!    (3) For PYTHIA8/Angantyr mode, hadrons and nuclei.
!         nuclei: 1000010020 = D/2H;  1000020040 = 4He;
!                 1000791970 = 197Au; 1000822080 = 208Pb; etc.
!
!   The standard PDG (hyper)nucleus codes (id/KF code) are 10-digit:
!           +-10LZZZAAAI.
!   +- denotes (hyper)nucleus and anti-one.
!   L = n_Lambda, total number of strange quarks (Lambda hyperons).
!   Z = n_p, total charges (protons).
!   A = n_p + n_n + n_Lambda, total baryon number.
!       n_n, total number of neutrons.
!   I = isomer level, = 0: groud state,
!                     > 0: excitations. m, n, p, q = 1, 2, 3, 4.
!   Cf. PDG "Monte Carlo particle numbering scheme" for more details.
!   Here for normal ions, we use
!       id = 100ZZZAAA0 = 1,000,000,000 + Z*10,000 + A*10 + 0.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    INTEGER PYK,PYCHGE,PYCOMP
    LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_PARTON,IS_QUARK,IS_DIQUARK
    PARAMETER (KSZJ=80000)
    COMMON/PYINT1/MINT(400),VINT(400)
    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
    common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
    common/sa30/vneump,vneumt,mstptj
    common/syspar/ipden,itden,suppm,suptm,suppc,suptc,r0p,r0t, &
     nap,nat,nzp,nzt,pio
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win_o, energy_B, psno, b_min, b_max
!   Local arrays and character variables.
    INTEGER frameType
    CHARACTER FRAME*16, BEAM*16, TARGET*16, NAME_X*16
    CHARACTER*100 PARAM_PYTHIA6
    DATA I_CALL /0/


    I_CALL = I_CALL + 1


!-------------------------------------------------------------------------------
!-------------------------   PYTHIA 6 Initialization   -------------------------
!   Reads the extra configuration for PYTHIA 6.
    OPEN( 11, FILE="pythia6_extra.cfg", STATUS="UNKNOWN" )
    DO WHILE( .TRUE. )
        READ( 11, "(A)", END=100 ) PARAM_PYTHIA6
        CALL PYGIVE( PARAM_PYTHIA6 )
    END DO
100 CLOSE(11)
    IF( .NOT.IS_PYTHIA8(i_mode) )THEN
!   Sets the name of collision frame. Defalut "CMS".
        FRAME = "CMS"
        IF( I_FRAME == 0 )THEN
            FRAME = "NONE"
        ELSE IF(  I_FRAME == 1 )THEN
            FRAME = "CMS"
        ELSE IF(  I_FRAME == 2 )THEN
            FRAME = "FIXT"
        ELSE IF(  I_FRAME == 3 )THEN
            FRAME = "3MOM"
        ELSE IF(  I_FRAME == 4 )THEN
            FRAME = "4MOM"
        ELSE IF(  I_FRAME == 5 )THEN
            FRAME = "5MOM"
        ELSE IF(  I_FRAME == 11 )THEN
            FRAME = "USER"
        ELSE IF(  I_FRAME == 12 )THEN
            FRAME = "USER_LHA"
        ELSE
            FRAME = "ERROR_FRAME"
            WRITE(*,*) "Wrong frame in PAINIT_KF, I_FRAME =",I_FRAME
            STOP
        END IF
!   Gets names of beam and target.
        CALL PYNAME( KF_BEAM, BEAM )
        CALL PYNAME( KF_TARGET, TARGET )

!   Prohibits writing of information on variable values changed by a PYGIVE
!     and a PYTUNE.
        IF( I_CALL >= 1 ) CALL PYGIVE( "MSTU(13)=0" )

!   Core initialization of PYTHIA 6.
        CALL PYINIT( FRAME, BEAM, TARGET, WIN )
!-------------------------   PYTHIA 6 Initialization   -------------------------
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!--------------------   PYTHIA 8 / Angantyr Initialization   -------------------
    ELSE IF( IS_PYTHIA8(i_mode) )THEN
!   Sets the name of collision frame. 1 = "CMS", 2 = "back-to-back"...
        i_frame_PY8 = I_FRAME
!   Initializes the error counters.
        MSTU(23) = 0
        MSTU(30) = 0
!   Records the frame of the collision of each "PAINIT" call.
        MINT(111) = i_frame_PY8
!   Records the frame of the entire simulation system.
        ! MSTI(111) = iframe
!   Records the collision energy.
        VINT(290) = WIN
!   Records the type of the collision.
        MINT(41) = ipden
        MINT(42) = itden
!   Records the number of events.
        MINT(5) = iii
!   Specifies whether or not to hadronize.
        ! MSTP(111) = mstptj
!   Specifies the b parameter sampling method for NA(AN) and AA, if Angantyr.
        ! MINT(39) = INT( psno )
!   Records id/KF codes of projectile and target of each "PAINIT_KF" call.
        MINT(11) = KF_BEAM
        MINT(12) = KF_TARGET
!   Records id/KF codes of proj. and targ. of the entire simulation system.
        ! MSTI(11) = KF_proj
        ! MSTI(12) = KF_targ

!   Core initialization and execution by PYTHIA 8.
        frameType = i_frame_PY8
        idA = KF_BEAM
        idB = KF_TARGET
        eCM = WIN
        ! In Pythia8_fort_interface.f90.
        CALL PYINIT_PY8( frameType, idA, idB, eCM )
    END IF
!--------------------   PYTHIA 8 / Angantyr Initialization   -------------------
!-------------------------------------------------------------------------------


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PAEVNT
!!  Generates core processes via PYTHIA 6 or 8.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_PARTON,IS_QUARK,IS_DIQUARK
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


!   Via PYTHIA 6.
!   Calling PYEVNW is the default (PYEVNW will be called inside PYEVNT).
    IF( .NOT.IS_PYTHIA8(i_mode) )THEN
        ! In p_40.f
        CALL PYEVNT
!   Via PYTHIA 8 and PYTHIA8/Angantyr.
    ELSE IF( IS_PYTHIA8(i_mode) )THEN
        ! In Pythia8_fort_interface.f90
        CALL PYEVNT_PY8
    END IF


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PAEXEC
!!  Hadronizaion via PYTHIA 6 or 8.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_PARTON,IS_QUARK,IS_DIQUARK
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


!   Via PYTHIA 6.
    IF( .NOT.IS_PYTHIA8(i_mode) )THEN
        ! In p_40.f
        CALL PYEXEC
!   Via PYTHIA 8 and PYTHIA8/Angantyr.
    ELSE IF( IS_PYTHIA8(i_mode) )THEN
        ! In Pythia8_fort_interface.f90.
        CALL PYEXEC_PY8
    END IF


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PAEDIT( MEDIT )
!!  Performs global manipulations on the event record, in particular
!!   to exclude unexisted unstable or undetectable partons/particles.
!!  In PYTHIA 6, existed, stable and detectable entries in "PYJETS"
!!   are those with 0 < KS < 11, while in PYTHIA 8, KS > 0.
!!  "MEDIT" is dummy now, reserved for future expansion.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_PARTON,IS_QUARK,IS_DIQUARK
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


!   Use PYEDIT(1) for PYTHIA 6 and PYEDIT(11) for PYTHIA 8 safely.
    IF( .NOT.IS_PYTHIA8(i_mode) )THEN
        CALL PYEDIT(1)
    ELSE IF( IS_PYTHIA8(i_mode) )THEN
        CALL PYEDIT(11)   ! Note here!
        CALL PYEDIT(12)   ! Note here!
    END IF


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PARECO( I_NN_PAIR )
!!  Records the original contents of each collision after PYTHIA execution.
!!  I.e. appends the information of one NN collision from /PYJETS/
!!   (PYTHIA 6), or /PYJETS_PY8/ & /PYJUNC_PY8/ (PYTHIA 8), etc. into
!!   /PYJETS_AA/, /PYJUNC_AA/ , etc.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    INTEGER PYK,PYCHGE,PYCOMP
    LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_PARTON,IS_QUARK,IS_DIQUARK
    PARAMETER (KSZJ=80000,KSZJ_PY8=300000)
    COMMON/PYINT1/MINT(400),VINT(400)
    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
    common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!   Arrays of particle information of one NN pair for PYTHIA 8.
    COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, &
           K_PY8(KSZJ_PY8,8), P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
!   Arrays of particle information of a AA collision (total NN pairs).
    COMMON/PYJETS_AA/ N_AA(KSZJ_PY8), N_TOT_AA, N_COLL_NN_AA, &
           K_AA(KSZJ_PY8,8), P_AA(KSZJ_PY8,7), V_AA(KSZJ_PY8,5)
!   Arrays of junction configuration of one NN pair for PYTHIA 8.
    INTEGER NO_JUNC_PY8, KIND_JUNC_PY8(KSZJ_PY8), &
            COL_JUNC_PY8(KSZJ_PY8,3), ENDC_JUNC_PY8(KSZJ_PY8,3), &
            STAT_JUNC_PY8(KSZJ_PY8,3)
    COMMON/PYJUNC_PY8/ NO_JUNC_PY8, KIND_JUNC_PY8, &
                       COL_JUNC_PY8, ENDC_JUNC_PY8, STAT_JUNC_PY8
!   Arrays of junction configuration of a AA collision (total NN pairs) for PY8.
    INTEGER NO_JUNC_AA(KSZJ_PY8), N_TOT_JUNC_AA, KIND_JUNC_AA(KSZJ_PY8), &
            COL_JUNC_AA(KSZJ_PY8,3), ENDC_JUNC_AA(KSZJ_PY8,3), &
            STAT_JUNC_AA(KSZJ_PY8,3)
    COMMON/PYJUNC_AA/ NO_JUNC_AA, N_TOT_JUNC_AA, KIND_JUNC_AA, &
                      COL_JUNC_AA, ENDC_JUNC_AA, STAT_JUNC_AA
    DATA I_NN_PAIR_PREVIOUS / 0 /
    DATA N_PREVIOUS / 0 /
    DATA N_JUNC_PREVIOUS / 0 /


!   Checks if a new event was generated.
    IF( I_NN_PAIR <= 1 )THEN
        N_TOT_AA      = 0
        N_TOT_JUNC_AA = 0
        N_COLL_NN_AA  = 0
        I_NN_PAIR_PREVIOUS = 0
        N_PREVIOUS         = 0
        N_JUNC_PREVIOUS    = 0
    END IF

!   Checks if the NN paire was the one regenerated (due to some errors).
!   Deduction due to the regenerated event.
    IF( I_NN_PAIR == I_NN_PAIR_PREVIOUS )THEN
        N_TOT_AA      = N_TOT_AA      - N_PREVIOUS
        N_TOT_JUNC_AA = N_TOT_JUNC_AA - N_JUNC_PREVIOUS
    END IF

    N_COLL_NN_AA = I_NN_PAIR
    ! Reserved for the re-generation check.
    I_NN_PAIR_PREVIOUS = I_NN_PAIR
    N_PREVIOUS = N

!   PYTHIA 6 mode.
    IF( .NOT.IS_PYTHIA8(i_mode) )THEN
        N_AA( N_COLL_NN_AA ) = N
!   "PYJETS" (superposition) -> "PYJETS_AA"
        DO J=1,5,1
            DO I=1,N,1
                K_AA( N_TOT_AA + I , J ) = K( I, J )
                P_AA( N_TOT_AA + I , J ) = P( I, J )
                V_AA( N_TOT_AA + I , J ) = V( I, J )
            END DO
        END DO
        N_TOT_AA = N_TOT_AA + N
!   PYTHIA 8 mode.
    ELSE IF( IS_PYTHIA8(i_mode) )THEN
!   "PYJETS_PY8" (superposition) -> "PYJETS_AA"
        N_AA( N_COLL_NN_AA ) = N_PY8
        ! Index 1-5.
        DO J=1,5,1
            DO I=1,N_PY8,1
                K_AA( N_TOT_AA + I , J ) = K_PY8( I, J )
                P_AA( N_TOT_AA + I , J ) = P_PY8( I, J )
                V_AA( N_TOT_AA + I , J ) = V_PY8( I, J )
            END DO
        END DO
        ! Index 6 and 7.
        DO J=6,7,1
            DO I=1,N_PY8,1
                K_AA( N_TOT_AA + I , J ) = K_PY8( I, J )
                P_AA( N_TOT_AA + I , J ) = P_PY8( I, J )
            END DO
        END DO
        ! Index 8.
        DO I=1,N_PY8,1
            K_AA( N_TOT_AA + I , 8 ) = K_PY8( I, 8 )
        END DO
!       "PYJUNC_PY8" (superposition) -> "PYJUNC_AA"
        NO_JUNC_AA( N_COLL_NN_AA ) = NO_JUNC_PY8
        N_JUNC = NO_JUNC_PY8
        DO I=1,N_JUNC,1
            KIND_JUNC_AA( N_TOT_JUNC_AA + I ) = KIND_JUNC_PY8( I )
            DO J=1,3,1
                COL_JUNC_AA(  N_TOT_JUNC_AA + I, J ) =  COL_JUNC_PY8( I, J )
                ENDC_JUNC_AA( N_TOT_JUNC_AA + I, J ) = ENDC_JUNC_PY8( I, J )
                STAT_JUNC_AA( N_TOT_JUNC_AA + I, J ) = STAT_JUNC_PY8( I, J )
            END DO
        END DO
        N_TOT_AA      = N_TOT_AA      + N_PY8
        N_TOT_JUNC_AA = N_TOT_JUNC_AA + N_JUNC
        ! Reserved for the re-generation check.
        N_PREVIOUS      = N_PY8
        N_JUNC_PREVIOUS = N_JUNC
    END IF


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PAFEED( I_NN,I_BEGIN,I_END, SUM_CHARGE_AND_MOMENTUM )
!!  Feeds the contents of one NN collision from /PYJETS_AA/ into
!!   /PYJETS/, /PYJETS_PY8/, /PYJUNC_PY8/, etc.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    INTEGER PYK,PYCHGE,PYCOMP
    LOGICAL IS_PYTHIA8,IS_EXIST,IS_NUCLEUS,IS_PARTON,IS_QUARK,IS_DIQUARK
    PARAMETER (KSZJ=80000,KSZJ_PY8=300000)
    COMMON/PYINT1/MINT(400),VINT(400)
    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
    common/sa1/kjp21,non1,bp,iii,neve,nout,nosc
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
!   Arrays of particle information of one NN pair for PYTHIA 8.
    COMMON/PYJETS_PY8/ N_PY8, NPAD_PY8, &
           K_PY8(KSZJ_PY8,8), P_PY8(KSZJ_PY8,7), V_PY8(KSZJ_PY8,5)
!   Arrays of particle information of a AA collision (total NN pairs).
    COMMON/PYJETS_AA/ N_AA(KSZJ_PY8), N_TOT_AA, N_COLL_NN_AA, &
                      K_AA(KSZJ_PY8,8), P_AA(KSZJ_PY8,7), V_AA(KSZJ_PY8,5)
!   Arrays of junction configuration of one NN pair for PYTHIA 8.
    INTEGER NO_JUNC_PY8, KIND_JUNC_PY8(KSZJ_PY8), &
            COL_JUNC_PY8(KSZJ_PY8,3), ENDC_JUNC_PY8(KSZJ_PY8,3), &
            STAT_JUNC_PY8(KSZJ_PY8,3)
    COMMON/PYJUNC_PY8/ NO_JUNC_PY8, KIND_JUNC_PY8, &
                       COL_JUNC_PY8, ENDC_JUNC_PY8, STAT_JUNC_PY8
!   Arrays of junction configuration of a AA collision (total NN pairs) for PY8.
    INTEGER NO_JUNC_AA(KSZJ_PY8), N_TOT_JUNC_AA, KIND_JUNC_AA(KSZJ_PY8), &
            COL_JUNC_AA(KSZJ_PY8,3), ENDC_JUNC_AA(KSZJ_PY8,3), &
            STAT_JUNC_AA(KSZJ_PY8,3)
    COMMON/PYJUNC_AA/ NO_JUNC_AA, N_TOT_JUNC_AA, KIND_JUNC_AA, &
                      COL_JUNC_AA, ENDC_JUNC_AA, STAT_JUNC_AA
    DIMENSION SUM_CHARGE_AND_MOMENTUM( 6 )
    DATA I_END_JUNC   / 0 /


!   Checks if a new event was generated.
    IF( I_NN == 1 )THEN
        I_BEGIN_JUNC = 0
        I_END_JUNC   = 0
    end if

!   PYTHIA 6 mode. ("i_mode == 3" is enough.)
    IF( .NOT.IS_PYTHIA8(i_mode) )THEN
!       Strings of one NN pair in "PYJETS_AA" -> "PYJETS".
        I_PYJ = 0
        DO I1 = I_BEGIN, I_END, 1
            I_PYJ = I_PYJ + 1
            DO I2=1,5,1
                K( I_PYJ, I2 ) = K_AA( I1, I2 )
                P( I_PYJ, I2 ) = P_AA( I1, I2 )
                V( I_PYJ, I2 ) = V_AA( I1, I2 )
                KS = K_AA( I1, 1 )
                KF = K_AA( I1, 2 )
                IF( KS >= 1 .AND. KS <= 10 ) &
                SUM_CHARGE_AND_MOMENTUM( I2 ) &
                    = SUM_CHARGE_AND_MOMENTUM( I2 ) + P_AA( I1, I2 )
            END DO
            IF( KS >= 1 .AND. KS <= 10 )THEN
                SUM_CHARGE_AND_MOMENTUM( 5 ) &
                    = SUM_CHARGE_AND_MOMENTUM( 4 )**2 &
                    - SUM_CHARGE_AND_MOMENTUM( 1 )**2 &
                    - SUM_CHARGE_AND_MOMENTUM( 2 )**2 &
                    - SUM_CHARGE_AND_MOMENTUM( 3 )**2
                SUM_CHARGE_AND_MOMENTUM( 6 ) &
                    = SUM_CHARGE_AND_MOMENTUM( 6 ) + PYCHGE( KF )
            END IF
        END DO
        ! N = i_PYJ   ! It should be.
!   PYTHIA 8 mode. ("i_mode == 6/9" are enough.)
    ELSE IF( IS_PYTHIA8(i_mode) )THEN
!       Strings of one NN pair in "PYJETS_AA" -> "PYJETS_PY8".
!       Index 1-5.
        I_PYJ = 0
        DO I1 = I_BEGIN, I_END, 1
            I_PYJ = I_PYJ + 1
            DO I2=1,5,1
                K_PY8( I_PYJ, I2 ) = K_AA( I1, I2 )
                P_PY8( I_PYJ, I2 ) = P_AA( I1, I2 )
                V_PY8( I_PYJ, I2 ) = V_AA( I1, I2 )
                KS = K_AA( I1, 1 )
                KF = K_AA( I1, 2 )
                IF( KS >= 1 ) &
                SUM_CHARGE_AND_MOMENTUM( I2 ) &
                    = SUM_CHARGE_AND_MOMENTUM( I2 ) + P_AA( I1, I2 )
            END DO
            IF( KS >= 1 )THEN
                SUM_CHARGE_AND_MOMENTUM( 5 ) &
                    = SUM_CHARGE_AND_MOMENTUM( 4 )**2 &
                    - SUM_CHARGE_AND_MOMENTUM( 1 )**2 &
                    - SUM_CHARGE_AND_MOMENTUM( 2 )**2  &
                    - SUM_CHARGE_AND_MOMENTUM( 3 )**2
                SUM_CHARGE_AND_MOMENTUM( 6 ) &
                    = SUM_CHARGE_AND_MOMENTUM( 6 ) + PYCHGE( KF )
            END IF
        END DO
!       Index 6-7.
        I_PYJ = 0
        DO I1 = I_BEGIN, I_END, 1
            I_PYJ = I_PYJ + 1
            DO I2=6,7,1
                K_PY8( I_PYJ, I2 ) = K_AA( I1, I2 )
                P_PY8( I_PYJ, I2 ) = P_AA( I1, I2 )
            END DO
        END DO
!       Index 8.
        I_PYJ = 0
        DO I1 = I_BEGIN, I_END, 1
            I_PYJ = I_PYJ + 1
            K_PY8( I_PYJ, 8 ) = K_AA( I1, 8 )
        END DO
        ! N_PY8 = I_PYJ   ! It should be.
!       The junct. configuration of one NN pair in "PYJUNC_AA" -> "PYJUNC_PY8".
        ! Updaes I_BEGIN_JUNC.
        I_BEGIN_JUNC     = I_END_JUNC + 1
        N_JUNC           = NO_JUNC_AA( I_NN )
        NO_JUNC_PY8 = NO_JUNC_AA( I_NN )
        ! Updaes I_END_JUNC.
        I_END_JUNC       = I_END_JUNC + N_JUNC
        I_JUNC = 0
        DO I = I_BEGIN_JUNC, I_END_JUNC, 1
            I_JUNC = I_JUNC + 1
            KIND_JUNC_PY8( I_JUNC ) = KIND_JUNC_AA( I )
            DO J=1,3,1
                COL_JUNC_PY8(  I_JUNC, J ) =  COL_JUNC_AA( I, J )
                ENDC_JUNC_PY8( I_JUNC, J ) = ENDC_JUNC_AA( I, J )
                STAT_JUNC_PY8( I_JUNC, J ) = STAT_JUNC_AA( I, J )
            END DO
        END DO
    END IF


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    LOGICAL FUNCTION IS_PYTHIA8( i_mode )
!!  Determines whether the mode is PYTHIA 8 (Angantyr).
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i_mode


    IS_PYTHIA8 = i_mode == 5 .OR. i_mode == 6 .OR. i_mode == 8 .OR. i_mode == 9


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    LOGICAL FUNCTION IS_EXIST( KS, i_mode )
!!  Determines whether an entry (KS, i.e. K(I,1) ) exists.
!   In PYTHIA 6, an ordinary particles with KS between 1 and 10
!    exist currently, while other KS present the historical entries.
!    Special attention should be paid to the "junction": KS=41,42
!    exists, KS=51,52 dose not.
!   In PYTHIA 8, the entries with KS > 0 exists currently.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: KS, i_mode


    IS_EXIST = .NOT. ( KS <= 0 &
                .OR. ( (i_mode == 2 .OR. i_mode == 3 .OR. i_mode == 4) &
                .AND. (KS > 10 .AND. KS /= 41 .AND. KS /= 42) ) )


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    LOGICAL FUNCTION IS_NUCLEUS( KF )
!!  Determines whether a particle (KF code) is a nucleus.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: KF


    IS_NUCLEUS = ABS(KF) > 1000000000


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    LOGICAL FUNCTION IS_PARTON( KF )
!!  Determines whether a particle (KF code) is a parton (quark, gluon, diquark).
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: KF


    IS_PARTON = ( ABS(KF) >= 1 .AND. ABS(KF) <= 8 ) .OR. KF == 21 &
            .OR. ( ( ABS(KF) >= 1103 .AND. ABS(KF) <= 5503 ) &
            .AND. MOD( ABS(KF) / 10, 10 ) == 0 )


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    LOGICAL FUNCTION IS_QUARK( KF )
!!  Determines whether a particle (KF code) is a quark/antiquark.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: KF


    IS_QUARK = ABS(KF) >= 1 .AND. ABS(KF) <= 8


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    LOGICAL FUNCTION IS_DIQUARK( KF )
!!  Determines whether a particle (KF code) is a diquark.
!   Includes light and heavy diquarks ABS(KF): 1103 ~ 5503.
!   The tens position of the diquark KF is 0, differing from that
!    of the baryon.
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: KF


    IS_DIQUARK = ( ABS(KF) >= 1103 .AND. ABS(KF) <= 5503 ) &
            .AND. MOD( ABS(KF) / 10, 10 ) == 0


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    FUNCTION PAPYP(I,J)
!!  Provides various real-valued event related data in /PYJETS/.
!   Rewritten from the function PYP of PYTHIA 6. See PYP in the PYTHIA 6 manual.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    INTEGER PYCHGE
    LOGICAL IS_EXIST
    PARAMETER (KSZJ=80000)
    COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
    SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
!   For the simulation control.
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
    DIMENSION PSUM(4)


!   Set default value. For I = 0 sum of momenta or charges,
!    or invariant mass of system.
    PAPYP = 0D0
    IF( I < 0 .OR. I > MSTU(4) .OR. J <= 0 )THEN
    ELSE IF( I == 0 .AND. J <= 4 )THEN
        DO I1=1,N,1
            KS = K(I1,1)
            KF = K(I1,2)
            IF( .NOT.IS_EXIST(KS,i_mode) .OR. KF == 88 ) CYCLE
            PAPYP = PAPYP + P(I1,J)
        END DO
    ELSE IF( I == 0 .AND. J == 5 )THEN
        DO J1=1,4,1
            PSUM(J1) = 0D0
            DO I1=1,N,1
                KS = K(I1,1)
                KF = K(I1,2)
                IF( IS_EXIST(KS,i_mode) .AND. KF /= 88 ) &
                    PSUM(J1) = PSUM(J1) + P(I1,J1)
            END DO
        END DO
        PAPYP = SQRT( MAX(0D0, PSUM(4)**2 -PSUM(1)**2 -PSUM(2)**2 -PSUM(3)**2) )
    ELSE IF( I == 0 .AND. J == 6 )THEN
        DO I1=1,N
            KS = K(I1,1)
            KF = K(I1,2)
            IF( IS_EXIST(KS,i_mode) .AND. KF /= 88 ) &
                PAPYP = PAPYP + PYCHGE(K(I1,2)) / 3D0
        END DO
    ELSE IF( I == 0 )THEN

!   Direct readout of P matrix.
    ELSE IF( J <= 5 )THEN
        PAPYP = P(I,J)

!   Charge, total momentum, transverse momentum, transverse mass.
    ELSE IF( J <= 12 )THEN
        IF( J == 6 ) PAPYP = PYCHGE(K(I,2)) / 3D0
        IF( J == 7  .OR. J == 8  ) PAPYP = P(I,1)**2 + P(I,2)**2 + P(I,3)**2
        IF( J == 9  .OR. J == 10 ) PAPYP = P(I,1)**2 + P(I,2)**2
        IF( J == 11 .OR. J == 12 ) PAPYP = P(I,5)**2 + P(I,1)**2 + P(I,2)**2
        IF( J == 8  .OR. J == 10 .OR. J == 12 ) PAPYP = SQRT(PAPYP)

!   Theta and phi angle in radians or degrees.
    ELSE IF( J <= 16 )THEN
        IF( J <= 14 ) PAPYP = PYANGL( P(I,3), SQRT( P(I,1)**2 + P(I,2)**2) )
        IF( J >= 15 ) PAPYP = PYANGL( P(I,1), P(I,2) )
        IF( J == 14 .OR. J == 16 ) PAPYP = PAPYP * 180D0 / PARU(1)

!   True rapidity, rapidity with pion mass, pseudorapidity.
    ELSE IF( J <= 19 )THEN
        PMR = 0D0
        IF( J == 17 ) PMR = P(I,5)
        IF( J == 18 ) PMR = PYMASS(211)
        PR  = MAX( 1D-20, PMR**2 + P(I,1)**2 + P(I,2)**2 )
        PAPYP = SIGN( LOG( MIN( ( SQRT(PR + P(I,3)**2) + ABS(P(I,3)) ) &
            / SQRT(PR), 1D20 ) ), P(I,3) )

!   Energy and momentum fractions (only to be used in CM frame).
    ELSE IF( J <= 25 )THEN
        IF( J == 20 ) PAPYP = 2D0 &
                            * SQRT(P(I,1)**2 +P(I,2)**2 +P(I,3)**2)/PARU(21)
        IF( J == 21 ) PAPYP = 2D0 * P(I,3) / PARU(21)
        IF( J == 22 ) PAPYP = 2D0 * SQRT( P(I,1)**2 + P(I,2)**2 ) / PARU(21)
        IF( J == 23 ) PAPYP = 2D0 * P(I,4) / PARU(21)
        IF( J == 24 ) PAPYP = ( P(I,4) + P(I,3) ) / PARU(21)
        IF( J == 25 ) PAPYP = ( P(I,4) - P(I,3) ) / PARU(21)
    END IF


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    FUNCTION IPAPYK(I,J)
!!  Provides various integer-valued event related data in /PYJETS/.
!   Rewritten from the function PYK of PYTHIA 6. See PYK in the PYTHIA 6 manual.
!TODO(Lei20241109): not so good for PYTHIA 8.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    INTEGER PYCHGE,PYCOMP
    LOGICAL IS_EXIST
    PARAMETER (KSZJ=80000)
    COMMON/PYJETS/N,NPAD,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
    SAVE /PYJETS/,/PYDAT1/,/PYDAT2/
!   For the simulation control.
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max


!   Default value. For I=0 number of entries, number of stable entries
!    or 3 times total charge.
    IPAPYK = 0
    IF( I < 0 .OR. I > MSTU(4) .OR. J <= 0  )THEN
    ELSE IF( I == 0 .AND. J == 1 )THEN
        IPAPYK = N
    ELSE IF( I == 0 .AND. (J == 2 .OR. J == 6) )THEN
        DO I1=1,N
            KS = K(I1,1)
            KF = K(I1,2)
            IF( .NOT.IS_EXIST(KS,i_mode) .OR. KF == 88 ) cycle
            IF( J == 2 ) IPAPYK = IPAPYK + 1
            IF( J == 6 ) IPAPYK = IPAPYK + PYCHGE( K(I1,2) )
        END DO
    ELSE IF( I == 0  )THEN

!   For I > 0 direct readout of K matrix or charge.
    ELSE IF( J <= 5 )THEN
        IPAPYK = K(I,J)
    ELSE IF( J == 6 )THEN
        IPAPYK = PYCHGE( K(I,2) )

!   Status (existing/fragmented/decayed), parton/hadron separation.
    ELSE IF( J <= 8 )THEN
        KS = K(I,1)
        KF = K(I,2)
        IF( IS_EXIST(KS,i_mode) .AND. KF /= 88 ) IPAPYK = 1
        IF( J == 8 ) IPAPYK = IPAPYK * K(I,2)
    ELSE IF( J <= 12 )THEN
        KFA = ABS( K(I,2) )
        KC = PYCOMP(KFA)
        KQ = 0
        IF( KC /= 0 ) KQ = KCHG(KC,2)
        IF( J == 9  .AND. KC /= 0 .AND. KQ /= 0 ) IPAPYK = K(I,2)
        IF( J == 10 .AND. KC /= 0 .AND. KQ == 0 ) IPAPYK = K(I,2)
        IF( J == 11 ) IPAPYK = KC
        IF( J == 12 ) IPAPYK = KQ * SIGN( 1, K(I,2) )

!   Heaviest flavour in hadron/diquark.
    ELSE IF( J == 13 )THEN
        KFA = ABS( K(I,2) )
        IPAPYK = MOD( KFA/100, 10 ) * (-1)**MOD( KFA / 100, 10 )
        IF(KFA < 10) IPAPYK = KFA
        IF( MOD( KFA / 1000, 10 ) /= 0 ) IPAPYK = MOD( KFA / 1000, 10 )
        IPAPYK = IPAPYK * ISIGN( 1, K(I,2) )

!   Particle history: generation, ancestor, rank.
    ELSE IF( J <= 15 )THEN
        I2 = I
        I1 = I
  110   IPAPYK = IPAPYK + 1
        I2 = I1
        I1 = K(I1,3)
        IF( I1 > 0 )THEN
            IF( K(I1,1) > 0 .AND. K(I1,1) <= 20 ) GOTO 110
        END IF
        IF( J == 15 ) IPAPYK = I2
    ELSE IF( J == 16 )THEN
        KFA = ABS( K(I,2) )
        IF( K(I,1) <= 20 .AND. ( ( KFA >= 11 .AND. KFA <= 20 ) .OR. KFA == 22 &
            .OR. ( KFA > 100 .AND. MOD( KFA / 10, 10 ) /= 0 ) ) )THEN
            I1 = I
    120     I2 = I1
            I1 = K(I1,3)
            IF( I1 > 0 )THEN
                KFAM = ABS( K(I1,2) )
                ILP = 1
                IF( KFAM /= 0 .AND. KFAM <= 10 ) ILP = 0
                IF( KFAM == 21 .OR. KFAM == 91 .OR. KFAM == 92 .OR. &
                    KFAM == 93 ) ILP = 0
                IF( KFAM > 100 .AND. MOD( KFAM / 10, 10 ) == 0 ) ILP = 0
                IF( ILP == 1 ) GOTO 120
            END IF
            IF( K(I1,1) == 12 )THEN
                DO I3=I1+1,I2
                    IF( K(I3,3) == K(I2,3) .AND. K(I3,2) /= 91    &
                        .AND. K(I3,2) /= 92 .AND. K(I3,2) /= 93 ) &
                        IPAPYK = IPAPYK + 1
                END DO
            ELSE
                I3 = I2
140             IPAPYK = IPAPYK + 1
                I3 = I3 + 1
                IF( I3 < N .AND. K(I3,3) == K(I2,3) ) GOTO 140
            END IF
        END IF

!   Particle coming from collapsing jet system or not.
    ELSE IF( J == 17 )THEN
        I1 = I
  150   IPAPYK = IPAPYK + 1
        I3 = I1
        I1 = K(I1,3)
        I0 = MAX(1,I1)
        KC = PYCOMP( K(I0,2) )
        IF( I1 == 0 .OR. K(I0,1) <= 0 .OR. K(I0,1) > 20 .OR. KC == 0 )THEN
            IF( IPAPYK == 1 ) IPAPYK = -1
            IF( IPAPYK > 1  ) IPAPYK = 0
            RETURN
        END IF
        IF( KCHG(KC,2) == 0 ) GOTO 150
        IF( K(I1,1) /= 12 ) IPAPYK = 0
        IF( K(I1,1) /= 12 ) RETURN
        I2 = I1
  160   I2 = I2 + 1
        IF( I2 < N .AND. K(I2,1) /= 11 ) GOTO 160
        K3M = K( I3-1, 3 )
        IF( K3M >= I1 .AND. K3M <= I2 ) IPAPYK = 0
        K3P = K( I3+1, 3 )
        IF( I3 < N .AND. K3P >= I1 .AND. K3P <= I2 ) IPAPYK = 0

!   Number of decay products. Colour flow.
    ELSE IF( J == 18 )THEN
        IF( K(I,1) == 11 .OR. K(I,1) == 12 ) IPAPYK = MAX( 0, K(I,5)-K(I,4)+1 )
        IF( K(I,4) == 0  .OR. K(I,5) == 0  ) IPAPYK = 0
    ELSE IF( J <= 22 )THEN
        IF( K(I,1) /= 3.AND.K(I,1) /= 13.AND.K(I,1) /= 14) RETURN
        IF( J == 19 ) IPAPYK = MOD( K(I,4) / MSTU(5), MSTU(5) )
        IF( J == 20 ) IPAPYK = MOD( K(I,5) / MSTU(5), MSTU(5) )
        IF( J == 21 ) IPAPYK = MOD( K(I,4), MSTU(5) )
        IF( J == 22 ) IPAPYK = MOD( K(I,5), MSTU(5) )
    ELSE
    END IF


    RETURN
    END



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    SUBROUTINE PAXTOT( KF1_IN, KF2_IN, ECM_IN, XTOT_OUT )
!!  Parametrizes total, elastic and diffractive cross-sections
!!   for different energies and beams. Donnachie-Landshoff for
!!   total and Schuler-Sjostrand for elastic and diffractive.
!
!   Process code IPROC:
!    =  1 : p + p;
!    =  2 : pbar + p;
!    =  3 : pi+ + p;
!    =  4 : pi- + p;
!    =  5 : pi0 + p;
!    =  6 : phi + p;
!    =  7 : J/psi + p;
!    = 11 : rho + rho;
!    = 12 : rho + phi;
!    = 13 : rho + J/psi;
!    = 14 : phi + phi;
!    = 15 : phi + J/psi;
!    = 16 : J/psi + J/psi;
!
!   This routine is borrowed from PYXTOT of PYTHIA 6.428 with
!    some modifications.
!
!   KF1_IN, KF2_IN: incoming particles, e.g. 2212 - proton.
!   ECM_IN: incoming CM energy in GeV.
!   XTOT_OUT: the total cross section in mb to be returned.
!
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
!   Local arrays.
    DIMENSION NPROC(30),XPAR(30),YPAR(30),IHADA(20),IHADB(20), &
    PMHAD(4),BHAD(4),BETP(4),IFITSD(20),IFITDD(20),CEFFS(10,8), &
    CEFFD(10,9),SIGTMP(6,0:5)
!   Common constants.
    DATA EPS/0.0808D0/, ETA/-0.4525D0/, ALP/0.25D0/, CRES/2D0/, &
    PMRC/1.062D0/, SMP/0.880D0/, FACEL/0.0511D0/, FACSD/0.0336D0/, &
    FACDD/0.0084D0/
!   Number of multiple processes to be evaluated (= 0 : undefined).
    DATA NPROC/7*1,3*0,6*1,4*0,4*3,2*6,4*0/
!   X and Y parameters of sigmatot = X * s**epsilon + Y * s**(-eta).
    DATA XPAR/2*21.70D0,3*13.63D0,10.01D0,0.970D0,3*0D0, &
    8.56D0,6.29D0,0.609D0,4.62D0,0.447D0,0.0434D0,4*0D0, &
    0.0677D0,0.0534D0,0.0425D0,0.0335D0,2.11D-4,1.31D-4,4*0D0/
    DATA YPAR/ &
    56.08D0,98.39D0,27.56D0,36.02D0,31.79D0,-1.51D0,-0.146D0,3*0D0, &
    13.08D0,-0.62D0,-0.060D0,0.030D0,-0.0028D0,0.00028D0,4*0D0, &
    0.129D0,0.115D0,0.081D0,0.072D0,2.15D-4,1.70D-4,4*0D0/
!   Beam and target hadron class:
!   = 1 : p/n ; = 2 : pi/rho/omega; = 3 : phi; = 4 : J/psi.
    DATA IHADA/2*1,3*2,3,4,3*0,3*2,2*3,4,4*0/
    DATA IHADB/7*1,3*0,2,3,4,3,2*4,4*0/
!   Characteristic class masses, slope parameters, beta = sqrt(X).
    DATA PMHAD/0.938D0,0.770D0,1.020D0,3.097D0/
    DATA BHAD/2.3D0,1.4D0,1.4D0,0.23D0/
    DATA BETP/4.658D0,2.926D0,2.149D0,0.208D0/
!   Fitting constants used in parametrizations of diffractive results.
    DATA IFITSD/2*1,3*2,3,4,3*0,5,6,7,8,9,10,4*0/
    DATA IFITDD/2*1,3*2,3,4,3*0,5,6,7,8,9,10,4*0/
    DATA ((CEFFS(J1,J2),J2=1,8),J1=1,10)/ &
    0.213D0, 0.0D0, -0.47D0, 150D0, 0.213D0, 0.0D0, -0.47D0, 150D0, &
    0.213D0, 0.0D0, -0.47D0, 150D0, 0.267D0, 0.0D0, -0.47D0, 100D0, &
    0.213D0, 0.0D0, -0.47D0, 150D0, 0.232D0, 0.0D0, -0.47D0, 110D0, &
    0.213D0, 7.0D0, -0.55D0, 800D0, 0.115D0, 0.0D0, -0.47D0, 110D0, &
    0.267D0, 0.0D0, -0.46D0,  75D0, 0.267D0, 0.0D0, -0.46D0,  75D0, &
    0.232D0, 0.0D0, -0.46D0,  85D0, 0.267D0, 0.0D0, -0.48D0, 100D0, &
    0.115D0, 0.0D0, -0.50D0,  90D0, 0.267D0, 6.0D0, -0.56D0, 420D0, &
    0.232D0, 0.0D0, -0.48D0, 110D0, 0.232D0, 0.0D0, -0.48D0, 110D0, &
    0.115D0, 0.0D0, -0.52D0, 120D0, 0.232D0, 6.0D0, -0.56D0, 470D0, &
    0.115D0, 5.5D0, -0.58D0, 570D0, 0.115D0, 5.5D0, -0.58D0, 570D0/
    DATA ((CEFFD(J1,J2),J2=1,9),J1=1,10)/ &
    3.11D0,  -7.34D0,  9.71D0, 0.068D0, -0.42D0,  1.31D0, &
    -1.37D0,  35.0D0,   118D0,  3.11D0, -7.10D0,  10.6D0, &
    0.073D0, -0.41D0,  1.17D0, -1.41D0,  31.6D0,    95D0, &
    3.12D0,  -7.43D0,  9.21D0, 0.067D0, -0.44D0,  1.41D0, &
    -1.35D0,  36.5D0,   132D0,  3.13D0, -8.18D0, -4.20D0, &
    0.056D0, -0.71D0,  3.12D0, -1.12D0,  55.2D0,  1298D0, &
    3.11D0,  -6.90D0,  11.4D0, 0.078D0, -0.40D0,  1.05D0, &
    -1.40D0,  28.4D0,    78D0,  3.11D0, -7.13D0,  10.0D0, &
    0.071D0, -0.41D0,  1.23D0, -1.34D0,  33.1D0,   105D0, &
    3.12D0,  -7.90D0, -1.49D0, 0.054D0, -0.64D0,  2.72D0, &
    -1.13D0,  53.1D0,   995D0,  3.11D0, -7.39D0,  8.22D0, &
    0.065D0, -0.44D0,  1.45D0, -1.36D0,  38.1D0,   148D0, &
    3.18D0,  -8.95D0, -3.37D0, 0.057D0, -0.76D0,  3.32D0, &
    -1.12D0,  55.6D0,  1472D0,  4.18D0, -29.2D0,  56.2D0, &
    0.074D0, -1.36D0,  6.67D0, -1.14D0, 116.2D0,  6532D0/
    DATA PARU101 /0.00729735D0/
    DATA PARP102 /0.232D0/
    DATA PARP104 /0.8D0/ ! /1.0D0/
    DIMENSION SIGT(6)


    XTOT_OUT = 40D0
    SIGT = 0D0
    XTOT = 0D0
    XEL  = 0D0
    XDIFF_XB  = 0D0
    XDIFF_AX  = 0D0
    XDIFF_AXB = 0D0
    XINEL_ND  = 0D0
    XINEL_NSD = 0D0
    XINEL_NDD = 0D0

!   Orders flavours of incoming particles: KF1 < KF2.
    IF( IABS(KF1_IN)  <=  IABS(KF2_IN)  )THEN
        KF1  = IABS( KF1_IN )
        KF2  = IABS( KF2_IN )
        IORD = 1
    ELSE
        KF1  = IABS( KF2_IN )
        KF2  = IABS( KF1_IN )
        IORD = 2
    END IF
    ISGN12 = ISIGN( 1, KF1_IN*KF2_IN )

!   Finds process number (for lookup tables).
    IPROC = 30
    IF( KF1 > 1000 )THEN
        IPROC = 1
        IF( ISGN12 < 0 ) IPROC = 2
    ELSE IF( KF1 > 100 .AND. KF2 > 1000 )THEN
        IPROC = 3
        IF( ISGN12 < 0 ) IPROC = 4
        IF( KF1 == 111 ) IPROC = 5
    ELSE IF( KF1 > 100 )THEN
        IPROC = 11
    END IF

!   Number of multiple processes to be stored; beam/target side.
    NPR = NPROC(IPROC)

!   Do not do any more for user-set or undefined cross-sections.
    IF(NPR /= 1)THEN
        WRITE(MSTU(11),*) "PAXTOT: cross section for this process " &
                // "not yet implemented, KF1, KF2:", KF1_IN, KF2_IN
        RETURN
    END IF

!   Parameters. Combinations of the energy. (Duplication?)
    AEM  = PARU101
    PMTH = PARP102
    S    = ECM_IN**2
    SRT  = ECM_IN
    SEPS = S**EPS
    SETA = S**ETA
    SLOG = LOG(S)

!   Loops over multiple processes (for VDM).
    I = 1
    IPR = IPROC

!   Evaluates hadron species, mass, slope contribution and fit number.
    IHA = IHADA(IPR)
    IHB = IHADB(IPR)
    PMA = PMHAD(IHA)
    PMB = PMHAD(IHB)
    BHA = BHAD(IHA)
    BHB = BHAD(IHB)
    ISD = IFITSD(IPR)
    IDD = IFITDD(IPR)

!   Skips if energy too low relative to masses.
    SIGTMP = 0D0
    IF( SRT  <  (PMA + PMB + PARP104) )THEN
        WRITE(MSTU(11),*)
        WRITE(MSTU(11),*) "PAXTOT warning: energy too low relative" &
                // " to masses. KF1, KF2, ECM, (MASSES + E_threshold): ", &
                                KF1_IN, KF2_IN, SRT, PMA + PMB + PARP104
        WRITE(MSTU(11),*) "A total NN cross section of 40 mb will " &
                       // "be used. Or please input it via ""para1_1""."
        WRITE(MSTU(11),*)
        RETURN
    END IF

!   Total cross-section. Elastic slope parameter and cross-section.
    SIGT(1) = XPAR(IPR)*SEPS + YPAR(IPR)*SETA
    BEL = 2D0*BHA + 2D0*BHB + 4D0*SEPS - 4.2D0
    SIGT(2) = FACEL*SIGT(1)**2 / BEL

!   Diffractive scattering A + B -> X + B.
    BSD  = 2D0*BHB
    SQML = (PMA + PMTH)**2
    SQMU = S*CEFFS(ISD,1) + CEFFS(ISD,2)
    SUM1 = LOG( ( BSD + 2D0*ALP*LOG(S/SQML) ) / &
                ( BSD + 2D0*ALP*LOG(S/SQMU) ) ) / (2D0*ALP)
    BXB  = CEFFS(ISD,3) + CEFFS(ISD,4) / S
    SUM2 = CRES*LOG( 1D0 + ( (PMA+PMRC) / (PMA+PMTH) )**2 ) / &
            (BSD + 2D0*ALP*LOG( S/( (PMA+PMTH)*(PMA+PMRC) ) ) + BXB)
    SIGT(3) = FACSD*XPAR(IPR)*BETP(IHB) * MAX(0D0,SUM1+SUM2)

!   Diffractive scattering A + B -> A + X.
    BSD  = 2D0*BHA
    SQML = (PMB + PMTH)**2
    SQMU = S*CEFFS(ISD,5) + CEFFS(ISD,6)
    SUM1 = LOG( ( BSD + 2D0*ALP*LOG(S/SQML) ) / &
                ( BSD + 2D0*ALP*LOG(S/SQMU) ) ) / (2D0*ALP)
    BAX  = CEFFS(ISD,7) + CEFFS(ISD,8) / S
    SUM2 = CRES*LOG( 1D0 + ( (PMB+PMRC) / (PMB+PMTH) )**2) / &
            (BSD + 2D0*ALP*LOG( S/( (PMB+PMTH)*(PMB+PMRC) ) ) + BAX)
    SIGT(4) = FACSD*XPAR(IPR)*BETP(IHA) * MAX(0D0,SUM1+SUM2)

!   OrderS single diffractive correctly.
    IF( IORD == 2 )THEN
        SIGSAV      = SIGT(3)
        SIGT(3) = SIGT(4)
        SIGT(4) = SIGSAV
    END IF

!   Double diffractive scattering A + B -> X1 + X2.
    YEFF=LOG(S*SMP/((PMA+PMTH)*(PMB+PMTH))**2)
    DEFF=CEFFD(IDD,1)+CEFFD(IDD,2)/SLOG+CEFFD(IDD,3)/SLOG**2
    SUM1=(DEFF+YEFF*(LOG(MAX(1D-10,YEFF/DEFF))-1D0))/(2D0*ALP)
    IF( YEFF <= 0 ) SUM1 = 0D0
    SQMU=S*(CEFFD(IDD,4)+CEFFD(IDD,5)/SLOG+CEFFD(IDD,6)/SLOG**2)
    SLUP=LOG(MAX(1.1D0,S/(ALP*(PMA+PMTH)**2*(PMB+PMTH)*(PMB+PMRC))))
    SLDN=LOG(MAX(1.1D0,S/(ALP*SQMU*(PMB+PMTH)*(PMB+PMRC))))
    SUM2=CRES*LOG(1D0+((PMB+PMRC)/(PMB+PMTH))**2)*LOG(SLUP/SLDN)/ &
    (2D0*ALP)
    SLUP=LOG(MAX(1.1D0,S/(ALP*(PMB+PMTH)**2*(PMA+PMTH)*(PMA+PMRC))))
    SLDN=LOG(MAX(1.1D0,S/(ALP*SQMU*(PMA+PMTH)*(PMA+PMRC))))
    SUM3=CRES*LOG(1D0+((PMA+PMRC)/(PMA+PMTH))**2)*LOG(SLUP/SLDN)/ &
    (2D0*ALP)
    BXX=CEFFD(IDD,7)+CEFFD(IDD,8)/SRT+CEFFD(IDD,9)/S
    SLRR=LOG(S/(ALP*(PMA+PMTH)*(PMA+PMRC)*(PMB+PMTH)*(PMB+PMRC)))
    SUM4=CRES**2*LOG(1D0+((PMA+PMRC)/(PMA+PMTH))**2)* &
    LOG(1D0+((PMB+PMRC)/(PMB+PMTH))**2)/MAX(0.1D0,2D0*ALP*SLRR+BXX)
    SIGT(5)=FACDD*XPAR(IPR)*MAX(0D0,SUM1+SUM2+SUM3+SUM4)

!   Non-diffractive by unitarity.
    SIGT(6) = SIGT(1) - SIGT(2) - SIGT(3) - SIGT(4) - SIGT(5)

!   Total cross section.
    XTOT      = SIGT(1)
!   Elastic cross section.
    XEL       = SIGT(2)
!   Single diffractive cross section (AB -> XB).
    XDIFF_XB  = SIGT(3)
!   Single diffractive cross section (AB -> AX).
    XDIFF_AX  = SIGT(4)
!   Double diffractive cross section.
    XDIFF_AXB = SIGT(5)
!   The inelastic, non-diffractive cross section.
    XINEL_ND  = SIGT(6)
!   The inelastic, non-single-diffractive cross section.
    XINEL_NSD = SIGT(6) + SIGT(5)
!   Non-double diffractive inelastic cross section.
    XINEL_NDD = SIGT(6) + SIGT(3) + SIGT(4)

    XTOT_OUT = XTOT


    RETURN
    END



!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    SUBROUTINE PASTAT(I_STAT,I_CALL)
!!  Prints out information about cross-sections.
!!#TODO(Lei20240220): not work for PYTHIA 8 mode.
!
!   The blocks used for cross section are /PYINT5/ and /PYINT5_S/
!    MSUB(I) stored switch of subprocesses.
!    NGEN(I,1): the number of times the differential cross section has
!               been evaluated for subprocess ISUB. NGEN(0,1) is the sum.
!    NGEN(I,3): the number of times an event of subprocess
!               type ISUB is generated. NGEN(0,3) is the sum of these.
!    XSEC(I,3): the estimated integrated cross section for subproces ISUB.
!    XSEC(0,3): the estimated the estimated total cross for section for
!               all subprocesses included. (mb)
!    PROC(I): character strings for the different possible subprocesses.
!    PROC(0): denotes all processes.
!
!   I_STAT: controlls different usage
!       =-2, initializes the global variables.
!       =-1, counts single-event cross sections.
!       = 0, counts total-event cross sections.
!       = 1, prints the cross sections.
!
!   I_CALL: the number of PASTAT call.
!
!   Double precision and integer declarations.
    IMPLICIT DOUBLE PRECISION(A-H, O-Z)
    IMPLICIT INTEGER(I-N)
    LOGICAL IS_PYTHIA8
    COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
    COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
    COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
    COMMON/PYINT1/MINT(400),VINT(400)
    COMMON/PYINT5/NGENPD,NGEN(0:500,3),XSEC(0:500,3)
    COMMON/PYINT6/PROC(0:500)
    CHARACTER PROC*28
    SAVE /PYDAT1/,/PYSUBS/,/PYINT1/,/PYINT5/,/PYINT6/
!   Local arrays, character variables and data.
    CHARACTER PROGA(6)*28,DISGA(2)*28, &
    PROGG9(13)*28,PROGG4(4)*28,PROGG2(2)*28,PROGP4(4)*28
    DATA PROGA/ &
    'VMD/hadron * VMD            ','VMD/hadron * direct         ', &
    'VMD/hadron * anomalous      ','direct * direct             ', &
    'direct * anomalous          ','anomalous * anomalous       '/
    DATA DISGA/'e * VMD','e * anomalous'/
    DATA PROGG9/ &
    'direct * direct             ','direct * VMD                ', &
    'direct * anomalous          ','VMD * direct                ', &
    'VMD * VMD                   ','VMD * anomalous             ', &
    'anomalous * direct          ','anomalous * VMD             ', &
    'anomalous * anomalous       ','DIS * VMD                   ', &
    'DIS * anomalous             ','VMD * DIS                   ', &
    'anomalous * DIS             '/
    DATA PROGG4/ &
    'direct * direct             ','direct * resolved           ', &
    'resolved * direct           ','resolved * resolved         '/
    DATA PROGG2/ &
    'direct * hadron             ','resolved * hadron           '/
    DATA PROGP4/ &
    'VMD * hadron                ','direct * hadron             ', &
    'anomalous * hadron          ','DIS * hadron                '/
    COMMON/SA1/KJP21,NON1,BP,III,NEVE,NOUT,NOSC
    COMMON/SYSPAR/IPDEN,ITDEN,SUPPM,SUPTM,SUPPC,SUPTC,R0P,R0T, &
        NAP,NAT,NZP,NZT,PIO
    COMMON/PYINT5_1/NGENPD_1,NGEN_1(0:500,3),XSEC_1(0:500,3)
    COMMON/PYDAT1_1/MSTU_1(200),PARU_1(200),MSTJ_1(200),PARJ_1(200)
    COMMON/PYINT5_S/NGENPD_S,NGEN_S(0:500,3),XSEC_S(0:500,3)
    COMMON/PYDAT1_S/MSTU_S(200),PARU_S(200),MSTJ_S(200),PARJ_S(200)
    COMMON/SA1_PY8/ i_mode, i_tune, KF_woDecay(100), &
           KF_proj, KF_targ, win, energy_B, psno, b_min, b_max
    DIMENSION NGEN_O(0:500,3),XSEC_O(0:500,3),MSTU_O(200)


!   Prints out the PYTHIA 8 information (not available now).
    IF( i_mode == 1 ) RETURN
    IF( IS_PYTHIA8(i_mode) )THEN
        ! In Pythia8_fort_interface.f90
        IF( I_STAT == 1 ) CALL PYSTAT_PY8( 1 )
        RETURN
    END IF

!   Initializes global variables.
    IF(I_STAT == -2 )THEN
        NGEN_S = 0
        XSEC_S = 0D0
        MSTU_S = 0
        PARU_S = 0D0
        MSTJ_S = 0
        PARJ_S = 0D0
        RETURN
!   Counts single-event cross sections.
    ELSE IF( I_STAT == -1 )THEN
        IF(I_CALL <= 1 )THEN
            NGEN_1 = 0
            XSEC_1 = 0D0
            MSTU_1 = 0
            PARU_1 = 0D0
            MSTJ_1 = 0
            PARJ_1 = 0D0
            IF(I_CALL == 0 )THEN
                NGEN = 0
                XSEC = 0D0
            END IF
        END IF
        NGEN_1 = NGEN_1 + NGEN
        XSEC_1 = XSEC_1 + XSEC
        RETURN
!   Counts total-event cross sections and numbers of errors and warnings.
    ELSE IF( I_STAT == 0 )THEN
        NGEN_S = NGEN_S + NGEN_1
        XSEC_S = XSEC_S + XSEC_1
        MSTU_S(23)  = MSTU_S(23)  + MSTU(23)
        MSTU_S(27)  = MSTU_S(27)  + MSTU(27)
        MSTU_S(30)  = MSTU_S(30)  + MSTU(30)
        RETURN
    ELSE IF( I_STAT == 1 )THEN
        write(MSTU(11),"(/)")
        write(MSTU(11),*) "All of events finished, iii=", iii
        write(MSTU(11),"(/)")
        write(MSTU(11),*) "************************* PACIAE " // &
        "Statistics (PASTAT) **************************"
        write(MSTU(11),*)
    END IF
    IF(I_STAT /= 1 .OR. III /= NEVE) RETURN
    NGEN_O = NGEN
    XSEC_O = XSEC
    MSTU_O = MSTU
    NGEN = NGEN_S / III
    XSEC = XSEC_S / DBLE(III)
    MSTU(23) = MSTU_S(23)
    MSTU(27) = MSTU_S(27)
    MSTU(30) = MSTU_S(30)

!   Cross-sections.
    IF(MINT(121) > 1) CALL PYSAVE(5,0)
    WRITE(MSTU(11),5000)
    WRITE(MSTU(11),5100)
    WRITE(MSTU(11),5200) 0,PROC(0),NGEN(0,3),NGEN(0,1),XSEC(0,3)
    DO 100 I=1,500
        IF(MSUB(I) /= 1) GOTO 100
        WRITE(MSTU(11),5200) I,PROC(I),NGEN(I,3),NGEN(I,1),XSEC(I,3)
100   CONTINUE
    IF(MINT(121) > 1 )THEN
        WRITE(MSTU(11),5300)
        DO 110 IGA=1,MINT(121)
            CALL PYSAVE(3,IGA)
            IF(MINT(121) == 2.AND.MSTP(14) == 10 )THEN
                WRITE(MSTU(11),5200) IGA,DISGA(IGA),NGEN(0,3),NGEN(0,1), &
                XSEC(0,3)
            ELSE IF( MINT(121) == 9.OR.MINT(121) == 13 )THEN
                WRITE(MSTU(11),5200) IGA,PROGG9(IGA),NGEN(0,3),NGEN(0,1), &
                XSEC(0,3)
            ELSE IF( MINT(121) == 4.AND.MSTP(14) == 30 )THEN
                WRITE(MSTU(11),5200) IGA,PROGP4(IGA),NGEN(0,3),NGEN(0,1), &
                XSEC(0,3)
            ELSE IF( MINT(121) == 4 )THEN
                WRITE(MSTU(11),5200) IGA,PROGG4(IGA),NGEN(0,3),NGEN(0,1), &
                XSEC(0,3)
            ELSE IF( MINT(121) == 2 )THEN
                WRITE(MSTU(11),5200) IGA,PROGG2(IGA),NGEN(0,3),NGEN(0,1), &
                XSEC(0,3)
            ELSE
                WRITE(MSTU(11),5200) IGA,PROGA(IGA),NGEN(0,3),NGEN(0,1), &
                XSEC(0,3)
            END IF
110       CONTINUE
        CALL PYSAVE(5,0)
    END IF
    WRITE(MSTU(11),5400) MSTU(23),MSTU(30),MSTU(27), &
    1D0-DBLE(NGEN(0,3))/MAX(1D0,DBLE(NGEN(0,2)))

    NGEN = NGEN_O
    XSEC = XSEC_O
    MSTU = MSTU_O

!   Formats for printouts.
! 5000 FORMAT('1',9('*'),1X,'PYSTAT:  Statistics on Number of ', &
!             'Events and Cross-sections',1X,9('*'))
5000 FORMAT('1',5('*'),1X,'PASTAT: Statistics on Number of ', &
            'NN Events and Total Cross-sections',1X,6('*'))
5100 FORMAT(/1X,78('=')/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',12X, &
            'Subprocess',12X,'I',6X,'Number of points',6X,'I',4X,'Sigma',3X, &
            'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',34('-'),'I',28('-'), &
            'I',4X,'(mb)',4X,'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,'I',1X, &
            'N:o',1X,'Type',25X,'I',4X,'Generated',9X,'Tried',1X,'I',12X, &
            'I'/1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')/1X,'I',34X,'I',28X, &
            'I',12X,'I')
5200 FORMAT(1X,'I',1X,I3,1X,A28,1X,'I',1X,I12,1X,I13,1X,'I',1X,1P, &
            D10.3,1X,'I')
5300 FORMAT(1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')/ &
            1X,'I',34X,'I',28X,'I',12X,'I')
5400 FORMAT(1X,'I',34X,'I',28X,'I',12X,'I'/1X,78('=')// &
            1X,'********* Total number of errors, excluding junctions =', &
            1X,I8,' *************'/ &
            1X,'********* Total number of errors, including junctions =', &
            1X,I8,' *************'/ &
            1X,'********* Total number of warnings =                   ', &
            1X,I8,' *************'/ &
            1X,'********* Fraction of events that fail fragmentation ', &
            'cuts =',1X,F8.5,' *********'/)


    RETURN
    END
!
!                                            Written by Anke at CCNU on 01/2024.
!
!       A series of general functions and subroutines for the execution with
!        both PYTHIA 6 and 8.
!
!240124 Lei
!*******************************************************************************
!*******************************************************************************