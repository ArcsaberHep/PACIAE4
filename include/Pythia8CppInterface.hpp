#pragma once
// Pythia8CppInterface.hpp is a part of the PACIAE event generator.
// Copyright (C) 2024 PACIAE Group.
// PACIAE is licensed under the GNU GPL v2 or later, see LICENSE for details.
// Open source: https://github.com/ArcsaberHep/PACIAE4
// Author: An-Ke Lei, January 2024 - January 2024.

// This is the C++ interface header file to link PACIAE (Fortran 77/90) with
//   PYTHIA 8 (C++).

//                                               By An-Ke at CCNU on 16/01/2024
//                                  Last updated by An-Ke at UiO  on 17/01/2025

// PYTHIA 8 header files.
#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8Plugins/ProgressLog.h"

// PACIAE header files.
#include "PaciaeUserHooks.hpp"

// PYTHIA namespace.
using namespace Pythia8;

// PACIAE namespace.
using namespace Paciae4;

//--------------------------------------------------------------------------

extern "C" {
    void instantiation_PY8( int (&MINT_c)[400], double (&VINT_c)[400],
                            int (&MSTU_c)[200], double (&PARU_c)[200],
                            int (&MSTJ_c)[200], double (&PARJ_c)[200],
                            int (&MSTP_c)[200], double (&PARP_c)[200],
                            int (&MSTI_c)[200], double (&PARI_c)[200],
                            int& nPY8, int (&kPY8)[8][300000],
                            double (&pPY8)[7][300000],
                            double (&vPY8)[5][300000],
                            int idStable[100],
                            Pythia** pythia,
                            std::shared_ptr< PaciaeUserHooks >**
                                             paciaeUserHooks,
                            std::shared_ptr< PaciaeHIUserHooks >**
                                             paciaeHIUserHooks );

//--------------------------------------------------------------------------

    void delete_object_PY8( int (&MINT_c)[400], int (&MSTP_c)[200],
                            Pythia** pythia,
                            Pythia** pythia_pp, Pythia** pythia_pn,
                            Pythia** pythia_np, Pythia** pythia_nn,
                            std::shared_ptr< PaciaeUserHooks >**
                                             paciaeUserHooks,
                            std::shared_ptr< PaciaeHIUserHooks >**
                                             paciaeHIUserHooks,
                            std::shared_ptr< PaciaeUserHooks >**
                                             paciaeUserHooks_pp,
                            std::shared_ptr< PaciaeUserHooks >**
                                             paciaeUserHooks_pn,
                            std::shared_ptr< PaciaeUserHooks >**
                                             paciaeUserHooks_np,
                            std::shared_ptr< PaciaeUserHooks >**
                                             paciaeUserHooks_nn );

//--------------------------------------------------------------------------

    void init_PY8( int frameType, int idA, int idB, double eCM,
                   int (&MINT_c)[400], double (&VINT_c)[400],
                   int (&MSTU_c)[200], double (&PARU_c)[200],
                   int (&MSTJ_c)[200], double (&PARJ_c)[200],
                   int (&MSTP_c)[200], double (&PARP_c)[200],
                   int (&MSTI_c)[200], double (&PARI_c)[200],
                   int& nPY8, int (&kPY8)[8][300000],
                   double (&pPY8)[7][300000],
                   double (&vPY8)[5][300000],
                   Pythia** pythia,
                   Pythia** pythia_pp, Pythia** pythia_pn,
                   Pythia** pythia_np, Pythia** pythia_nn,
                   std::shared_ptr< PaciaeUserHooks >** paciaeUserHooks,
                   std::shared_ptr< PaciaeHIUserHooks >** paciaeHIUserHooks );

//--------------------------------------------------------------------------
    void next_PY8( int (&MINT_c)[400], double (&VINT_c)[400],
                   int (&MSTU_c)[200], double (&PARU_c)[200],
                   int (&MSTJ_c)[200], double (&PARJ_c)[200],
                   int (&MSTP_c)[200], double (&PARP_c)[200],
                   int (&MSTI_c)[200], double (&PARI_c)[200],
                   int& nPY8, int (&kPY8)[8][300000],
                   double (&pPY8)[7][300000],
                   double (&vPY8)[5][300000],
                   int& noJuncPY8,
                   int (&kindJuncPY8)[300000],
                   int (&colJuncPY8)[3][300000],
                   int (&endcJuncPY8)[3][300000],
                   int (&statJuncPY8)[3][300000],
                   Pythia** pythia,
                   std::shared_ptr< PaciaeUserHooks >** paciaeUserHooks,
                   std::shared_ptr< PaciaeHIUserHooks >** paciaeHIUserHooks );

//--------------------------------------------------------------------------

    void forceHadronLevel_PY8( int (&MINT_c)[400], double (&VINT_c)[400],
                               int (&MSTU_c)[200], double (&PARU_c)[200],
                               int (&MSTJ_c)[200], double (&PARJ_c)[200],
                               int (&MSTP_c)[200], double (&PARP_c)[200],
                               int (&MSTI_c)[200], double (&PARI_c)[200],
                               int& nPY8, int (&kPY8)[8][300000],
                               double (&pPY8)[7][300000],
                               double (&vPY8)[5][300000],
                               int& noJuncPY8,
                               int (&kindJuncPY8)[300000],
                               int (&colJuncPY8)[3][300000],
                               int (&endcJuncPY8)[3][300000],
                               int (&statJuncPY8)[3][300000],
                               int (&iFail),
                               Pythia** pythia,
                               std::shared_ptr< PaciaeUserHooks >**
                                                               paciaeUserHooks,
                               std::shared_ptr< PaciaeHIUserHooks >**
                                                            paciaeHIUserHooks );

//--------------------------------------------------------------------------

    void decayAll_PY8( int& nPY8, int (&kPY8)[8][300000],
                       double (&pPY8)[7][300000],
                       double (&vPY8)[5][300000],
                       Pythia** pythia );

//--------------------------------------------------------------------------

    void stat_PY8( int (&MINT_c)[400], int (&MSTP_c)[200], int (&MSTI_c)[200],
                   Pythia** pythia,
                   Pythia** pythia_pp, Pythia** pythia_pn,
                   Pythia** pythia_np, Pythia** pythia_nn );

//--------------------------------------------------------------------------

    void debug_PY8( );
/*  void debug_PY8( int frameType, int idA, int idB, double eCM,
                    int (&MINT_c)[400], double (&VINT_c)[400],
                    int (&MSTU_c)[200], double (&PARU_c)[200],
                    int (&MSTJ_c)[200], double (&PARJ_c)[200],
                    int (&MSTP_c)[200], double (&PARP_c)[200],
                    int (&MSTI_c)[200], double (&PARI_c)[200],
                    int& nPY8, int (&kPY8)[8][300000],
                    double (&pPY8)[7][300000],
                    double (&vPY8)[5][300000],
                    int idStable[100], int (&iFail), */
}