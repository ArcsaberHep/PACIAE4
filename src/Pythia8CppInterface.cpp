// Pythia8CppInterface.cpp is a part of the PACIAE event generator.
// Copyright (C) 2024 PACIAE Group.
// PACIAE is licensed under the GNU GPL v2 or later, see LICENCE for details.
// Open source: https://github.com/ArcsaberHep/PACIAE4
// Author: An-Ke Lei, January 2024 - November 2024.

// This is the C++ interface program to link PACIAE (Fortran 77/90) with
//   PYTHIA 8 (C++).

//                                               By An-Ke at CCNU on 16/01/2024
//                                  Last updated by An-Ke at UiO  on 11/11/2024

// PYTHIA 8 header files.
#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8Plugins/ProgressLog.h"

// PACIAE header files.
#include "Pythia8CppInterface.hpp"
#include "PaciaeUserHooks.hpp"

// PYTHIA namespace.
using namespace Pythia8;

// PACIAE namespace.
using namespace Paciae4;

//--------------------------------------------------------------------------

// This function will instantiate Pythia objects including Userhooks,
//  i.e. allocates memory on the heap.
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
                                             paciaeHIUserHooks ) {

// Specifies output target file, i.e. "cout" redirection.
    std::streambuf* coutBuff = std::cout.rdbuf();
    std::ofstream outputFile;
    outputFile.open( "pythia8_report.log" );
    std::streambuf* fileBuff = outputFile.rdbuf();
    std::cout.rdbuf( fileBuff );

// Generator instantiation. On heap, not stack!
    *pythia = new Pythia();
    // Shorthand for some public members of pythia (also static ones).
    Settings& settings = (*pythia) -> settings;

// Note that in Fortran, an array index begins from (1), while in C++ it
//  is from [0] . So e.g. MSTP_c[4] here represents MSTP(5) in the
//  Fortran-based PYTHIA 6.

// Sets the parameter tune.
    int& iTune = MSTP_c[4];
    if( iTune > 0 ) {
        // 14: Monash 2013.
        settings.mode( "Tune:pp", iTune );
        // 7: Monash 2013. Tune::ee would be set in the Tune::pp ?
        // settings.mode( "Tune:ee",  iTune );
    }

// K factor of the hard process and the multiparton interaction (MPI).
    double& Kfactor = PARP_c[30];
    // settings.parm( "SigmaProcess:Kfactor", Kfactor );  // Should we?
    settings.parm( "MultipartonInteractions:Kfactor", Kfactor );

// Parameters in Lund string fragmentation function.
    double& aLund = PARJ_c[40];
    settings.parm("StringZ:aLund", aLund );
    double& bLund = PARJ_c[41];
    settings.parm("StringZ:bLund", bLund );
    double& sigma = PARJ_c[20];
    settings.parm("StringPT:sigma", sigma );

// Switches on/off Rope Hadronization in PYTHIA 8.
    int& iRope = MSTJ_c[0];
    if( iRope >= 6 && iRope <= 8 ) {
        (*pythia) -> readString( "Ropewalk:RopeHadronization = on" );
        (*pythia) -> readString( "Ropewalk:doFlavour = on" );
        (*pythia) -> readString( "Ropewalk:doShoving = on" );
        if( iRope == 6 )
            (*pythia) -> readString( "Ropewalk:doShoving = off" );
        if( iRope == 7 )
            (*pythia) -> readString( "Ropewalk:doFlavour = off" );
        // Requires parton coordinate information.
        (*pythia) -> readString( "PartonVertex:setVertex = on" );
        (*pythia) -> readString( "Fragmentation:setVertices = on" );
    }
    else if( iRope == 9 ) {
        (*pythia) -> readString( "StringPT:thermalModel = on" );
    }

// Selection of CR models.
    int& CRmode = MSTP_c[94];
    if( CRmode == -1 ) {
        settings.flag( "ColourReconnection:reconnect", 0 );
    }
    else {
        settings.mode( "ColourReconnection:mode", CRmode );
        if( CRmode == 1 ) {
            // Necessary option!
            settings.mode( "BeamRemnants:remnantMode", 1 );
        }
        else if( CRmode == 3 || CRmode == 4 ) {
            // Necessary option!
            settings.flag( "ColourReconnection:forceResonance", 1 );
        }
    }

//  Primary long-lived particle definition in ALICE-PUBLIC-2017-005.
    // settings.flag( "ParticleDecays:limitTau0", 1 );
    // settings.parm( "ParticleDecays:tau0Max", 10.0 );  // mm/c, i.e. 1 cm/c.
    (*pythia) -> readString( "211:mayDecay = false"   );   // pi+
    (*pythia) -> readString( "-211:mayDecay = false"  );   // pi-
    (*pythia) -> readString( "130:mayDecay = false"   );   // K0_L
    (*pythia) -> readString( "310:mayDecay = false"   );   // K0_S
    (*pythia) -> readString( "311:mayDecay = false"   );   // K0
    (*pythia) -> readString( "-311:mayDecay = false"  );   // Kbar0
    (*pythia) -> readString( "321:mayDecay = false"   );   // K+
    (*pythia) -> readString( "-321:mayDecay = false"  );   // K-
    (*pythia) -> readString( "2212:mayDecay = false"  );   // p+
    (*pythia) -> readString( "-2212:mayDecay = false" );   // pbar-
    (*pythia) -> readString( "2112:mayDecay = false"  );   // n0
    (*pythia) -> readString( "-2112:mayDecay = false" );   // nbar0
    (*pythia) -> readString( "3122:mayDecay = false"  );   // Lambda0
    (*pythia) -> readString( "-3122:mayDecay = false" );   // Lambdabar0
    (*pythia) -> readString( "3112:mayDecay = false"  );   // Sigma-
    (*pythia) -> readString( "-3112:mayDecay = false" );   // Sigmabar+
    (*pythia) -> readString( "3222:mayDecay = false"  );   // Sigma+
    (*pythia) -> readString( "-3222:mayDecay = false" );   // Sigmabar-
    (*pythia) -> readString( "3312:mayDecay = false"  );   // Xi-
    (*pythia) -> readString( "-3312:mayDecay = false" );   // Xibar+
    (*pythia) -> readString( "3322:mayDecay = false"  );   // Xi0
    (*pythia) -> readString( "-3322:mayDecay = false" );   // Xibar0
    (*pythia) -> readString( "3334:mayDecay = false"  );   // Omega-
    (*pythia) -> readString( "-3334:mayDecay = false" );   // Omegabar+
    (*pythia) -> readString( "11:mayDecay = false"    );   // e-
    (*pythia) -> readString( "-11:mayDecay = false"   );   // e+
    (*pythia) -> readString( "13:mayDecay = false"    );   // mu-
    (*pythia) -> readString( "-13:mayDecay = false"   );   // mu+
    (*pythia) -> readString( "22:mayDecay = false"    );   // gamma
//  mayDecay(id, false) means the particle id will not decay (stable).
    for( int iParticle = 0; iParticle < 100; ++iParticle ) {
        int& id = idStable[iParticle];
        if( id != 0 ) (*pythia) -> particleData.mayDecay(id, false);
    }

// Selects subprocesses.
    int& iProcess = MINT_c[0];
    switch (iProcess) {
        // Inelastic (INEL)
        case 0 :
            (*pythia) -> readString( "SoftQCD:inelastic = on" );
            break;
        // Non-Single Difractive (NSD)
        case 1 :
            (*pythia) -> readString( "SoftQCD:nonDiffractive = on" );
            (*pythia) -> readString( "SoftQCD:doubleDiffractive = on" );
            (*pythia) -> readString( "SoftQCD:centralDiffractive = on" );
            break;
        // Drell-Yan
        case 2 :
            (*pythia) -> readString( "WeakSingleBoson:ffbar2gmZ = on" );
            break;
        // J/psi (color-singlet via NRQCD )
        case 3 :
            (*pythia) -> readString( "Charmonium:" \
                                     "gg2ccbar(3S1)[3S1(1)]g = on,off" );
            (*pythia) -> readString( "Charmonium:" \
                                     "gg2ccbar(3S1)[3S1(1)]gm = on,off" );
            (*pythia) -> readString( "Charmonium:" \
                                     "gg2doubleccbar(3S1)[3S1(1)] " \
                                     "= on,on,off" );
            (*pythia) -> readString( "Charmonium:" \
                                     "qqbar2doubleccbar(3S1)[3S1(1)] " \
                                     "= on,on,off" );
            break;
        // Heavy-flavor production
        case 4 :
            (*pythia) -> readString( "HardQCD:hardccbar = on" );
            (*pythia) -> readString( "HardQCD:hardbbbar = on" );
            break;
        // Direct photon
        case 5 :
            (*pythia) -> readString( "PromptPhoton:all = on" );
            break;
        // Soft QCD
        case 6 :
            (*pythia) -> readString( "SoftQCD:all = on" );
            break;
        // W+/- production
        case 7 :
            (*pythia) -> readString( "WeakSingleBoson:ffbar2W = on" );
            break;
        // Hard QCD
        case 8 :
            (*pythia) -> readString( "HardQCD:all = on" );
            break;
        // Z0 production
        case 9 :
            (*pythia) -> readString( "WeakSingleBoson:ffbar2gmZ = on" );
            break;
        // Inelastic nondiffrative. (Minimum Bias)
        case 10 :
            (*pythia) -> readString( "SoftQCD:nonDiffractive = on" );
            break;
        case -61 :
            (*pythia) -> readString( "Charmonium:all = on" );
            break;
        case -62 :
            (*pythia) -> readString( "Bottomonium:all = on" );
            break;
        case -63 :
            (*pythia) -> readString( "Onia:all = on" );
            break;
        default :
            std::cout << "\nYou have chosen an undefined nchan = " << iProcess
                      << ".\nPlease specify subprocesses in " \
                         "\"pythia8_extra.cfg\".\n"
                      << std::endl;
    }
    // For e+e-. e+e- -> gamma*/Z -> qqbar (e+e- -> W+W- -> qqbar ?).
    if( std::abs(MSTI_c[10]) == 11 && std::abs(MSTI_c[11]) == 11
        && MSTI_c[10]*MSTI_c[11] < 0 ) {
        // Hard process, with hadronic Z decays.
        (*pythia) -> readString( "WeakSingleBoson:ffbar2gmZ = on" );
        // Switches off all Z0 decays and then switches back on those to quarks.
        (*pythia) -> readString( "23:onMode = off");
        (*pythia) -> readString( "23:onIfAny = 1 2 3 4 5");
        // Hard process, with hadronic W decays.
        // (*pythia) -> readString( "WeakDoubleBoson:ffbar2WW = on" );
        // Switches off all W+- decays then switches back on those to quarks.
        // (*pythia) -> readString( "24:onMode = off" );
        // (*pythia) -> readString( "24:onIfAny = 1 2 3 4 5" );
    }

// Switches of Angantyr machinery.
    int& iExecMode = MSTP_c[190];
    if( iExecMode == 8 || iExecMode == 9 ) {
        // Switches on Angantyr machinery for pA and AA, even NN collisions !
        (*pythia) -> readString( "HeavyIon:mode = 2" );
        // Gives partons transverse-coordinates x and y ( z = t = 0 ).
        (*pythia) -> readString( "PartonVertex:setVertex = on" );
        // Gives fragmented hadrons four-coordinates x, y, z and t.
        (*pythia) -> readString( "Fragmentation:setVertices = on" );
        // Defines nuclei, just giving names.
        if( !(*pythia)->particleData.isParticle( std::abs(MSTI_c[10]) ) ) {
            (*pythia) -> particleData.addParticle( std::abs(MSTI_c[10]),
                                                   "myIon", "myIonbar" );
        }
        if( !(*pythia)->particleData.isParticle( std::abs(MSTI_c[11]) ) ) {
            (*pythia) -> particleData.addParticle( std::abs(MSTI_c[11]),
                                                   "myIon2", "myIonbar2" );
        }
    }

// Sets up user hooks.
    *paciaeUserHooks = new std::shared_ptr<  PaciaeUserHooks >(
                           std::make_shared< PaciaeUserHooks >() );
    (*pythia) -> setUserHooksPtr( **paciaeUserHooks );
    // Sets up heavy-ion hooks.
    int bMode = 1;
    if( MINT_c[38] == 3 ) bMode = 0;
    *paciaeHIUserHooks = new std::shared_ptr<  PaciaeHIUserHooks >(
                             std::make_shared< PaciaeHIUserHooks >( bMode ) );
    (*pythia) -> setHIHooks( **paciaeHIUserHooks );

    // Recovers the std::cout .
    std::cout.rdbuf( coutBuff );
    outputFile.close();

    return;
    }

//--------------------------------------------------------------------------

//  This function will delete the Pythia object,
//   i.e. releases memory from the heap.
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
                                             paciaeUserHooks_nn ) {

// Specifies output target file, i.e. "cout" redirection.
    std::streambuf* coutBuff = std::cout.rdbuf();
    std::ofstream outputFile;
    outputFile.open( "pythia8_report.log", ios::out | ios::app );
    std::streambuf* fileBuff = outputFile.rdbuf();
    std::cout.rdbuf( fileBuff );

// Deletes instances.
    delete *pythia;
    delete *pythia_pp;
    delete *pythia_pn;
    delete *pythia_np;
    delete *pythia_nn;
    delete *paciaeUserHooks;
    delete *paciaeHIUserHooks;
    delete *paciaeUserHooks_pp;
    delete *paciaeUserHooks_pn;
    delete *paciaeUserHooks_np;
    delete *paciaeUserHooks_nn;

// Recovers the std::cout .
    std::cout.rdbuf( coutBuff );
    outputFile.close();

    return;
    }

//--------------------------------------------------------------------------

//  This function will initialize the collision system for PYTHIA 8.
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
                   std::shared_ptr< PaciaeHIUserHooks >** paciaeHIUserHooks ) {

// Specifies output target file, i.e. "cout" redirection.
    std::streambuf* coutBuff = std::cout.rdbuf();
    std::ofstream outputFile;
    outputFile.open( "pythia8_report.log", ios::out | ios::app );
    std::streambuf* fileBuff = outputFile.rdbuf();
    std::cout.rdbuf( fileBuff );

// Shorthand for some public members of pythia (also static ones).
// The pythia object was passed in from the external calling.
    Settings& settings = (*pythia) -> settings;
    ParticleData& particleData = (*pythia) -> particleData;

// Sets collision frame, id, and energy.
    settings.mode( "Beams:frameType", frameType );
    settings.mode( "Beams:idA", idA );
    settings.mode( "Beams:idB", idB );
    if( frameType == 2 ) {
        double& eA = pPY8[3][0];
        double& eB = pPY8[3][1];
        settings.parm( "Beams:eA", eA );
        settings.parm( "Beams:eB", eB );
    }
    else if( frameType == 3 ) {
        double& pxA = pPY8[0][0];
        double& pyA = pPY8[1][0];
        double& pzA = pPY8[2][0];
        double& pxB = pPY8[0][1];
        double& pyB = pPY8[1][1];
        double& pzB = pPY8[2][1];
        settings.parm( "Beams:pxA", pxA );
        settings.parm( "Beams:pyA", pyA );
        settings.parm( "Beams:pzA", pzA );
        settings.parm( "Beams:pxB", pxB );
        settings.parm( "Beams:pyB", pyB );
        settings.parm( "Beams:pzB", pzB );
    }
    else {
        settings.parm( "Beams:eCM", eCM );
    }

// Nuclear parton distribution functions.
    int& useHardNPDF = MSTP_c[193];
    if( useHardNPDF != -1 ) {
        // Projectile.
        (*pythia) -> readString( "PDF:useHardNPDFA = on" );
        settings.mode( "PDF:nPDFSetA", useHardNPDF );
        settings.mode( "PDF:nPDFBeamA", std::abs(MSTI_c[10]) );
        // Target.
        (*pythia) -> readString( "PDF:useHardNPDFB = on" );
        settings.mode( "PDF:nPDFSetB", useHardNPDF );
        settings.mode( "PDF:nPDFBeamB", std::abs(MSTI_c[11]) );
    }

// Sets the seed of the random number generator.
    settings.flag( "Random:setSeed", 1 );
    int& seed = MSTP_c[191];
    settings.mode( "Random:seed",  seed );

// Switches on/off PYTHIA 8 reports and event listing.
    // Switches on automatic event reports.
    (*pythia) -> readString( "Init:showProcesses = on" );
    (*pythia) -> readString( "Init:showMultipartonInteractions = on" );
    (*pythia) -> readString( "Init:showChangedSettings = on" );
    (*pythia) -> readString( "Init:showAllSettings = off" );
    (*pythia) -> readString( "Init:showChangedparticleData = off" );
    (*pythia) -> readString( "Init:showChangedResonanceData = off" );
    (*pythia) -> readString( "Init:showAllparticleData = off" );
    (*pythia) -> readString( "Init:showOneparticleData = 0" );
    // Switch off automatic event listing.
    (*pythia) -> readString( "Next:numberCount = 0" );
    (*pythia) -> readString( "Next:numberShowLHA = 0" );
    (*pythia) -> readString( "Next:numberShowInfo = 0" );
    (*pythia) -> readString( "Next:numberShowProcess = 0" );
    (*pythia) -> readString( "Next:numberShowEvent = 0" );
    // Other settings
    // Event-generation settings
    (*pythia) -> readString( "Next:showScaleAndVertex = off" );
    (*pythia) -> readString( "Next:showMothersAndDaughters = off" );
    // Statistics reports.
    (*pythia) -> readString( "Stat:showProcessLevel = on" );
    (*pythia) -> readString( "Stat:showPartonLevel = on" );
    (*pythia) -> readString( "Stat:showErrors = on" );
    (*pythia) -> readString( "Main:numberOfEvents = 1" );
    // Main-program settings
    (*pythia) -> readString( "Main:numberOfTriedEvents = 0" );
    (*pythia) -> readString( "Main:numberOfSelectedEvents = 0" );
    (*pythia) -> readString( "Main:numberOfAcceptedEvents = 0" );
    (*pythia) -> readString( "Main:timesAllowErrors = 0" );

// Key requirement: Switches on/off hadronization and decay.
    int& doHadronLevel = MSTP_c[110];
    // Does the string fragmentation indeed for the string-melting.
    int& hadronizationModel = MSTP_c[185];
    int& isStringMelting = MSTP_c[186];
    if( hadronizationModel == 1 && isStringMelting > 99 ) {
        doHadronLevel = 1;
        // Forbids particles decay. Retains primary hadrons.
        (*pythia) -> readString( "HadronLevel:Decay = off" );
    }
    settings.flag( "HadronLevel:all", doHadronLevel );
    // Avoid the standard scrutiny of mother/daughter relations.
    // But note that other event checks are done below.
    (*pythia) -> readString( "Check:event = off" );

// Reads extra setting from the external file.
    (*pythia) -> readFile("pythia8_extra.cfg");

// Instantiates 4 PYTHIA8 objects and does initialization for lA(Al), hA(Ah)
//   AB collisions. Typically pA(Ap) and AA. On heap, not stack!
    if( (  std::abs(MSTI_c[10])  > 1000000000 
        || std::abs(MSTI_c[11])  > 1000000000 )
        && ( MSTP_c[190] == 5 || MSTP_c[190] == 6 ) ) {
        // Multiple generator instantiations.
        *pythia_pp = new Pythia( settings, particleData, false );
        *pythia_pn = new Pythia( settings, particleData, false );
        *pythia_np = new Pythia( settings, particleData, false );
        *pythia_nn = new Pythia( settings, particleData, false );
        // Sets up UserHooks.
        (*pythia_pp) -> setUserHooksPtr( **paciaeUserHooks );
        (*pythia_pn) -> setUserHooksPtr( **paciaeUserHooks );
        (*pythia_np) -> setUserHooksPtr( **paciaeUserHooks );
        (*pythia_nn) -> setUserHooksPtr( **paciaeUserHooks );
        // Specifies beam and target particles.
        int idLeptonHadron = 2212;
        // lA, hA; pA
        if(    std::abs(MSTI_c[10]) <= 1000000000
            && std::abs(MSTI_c[11]) > 1000000000 ) {
            idLeptonHadron = MSTI_c[10];
        }
        // Al, Ah; Ap
        else if( std::abs(MSTI_c[10]) > 1000000000
              && std::abs(MSTI_c[11]) <= 1000000000 ) {
            idLeptonHadron = MSTI_c[11];
        }
        // lp, hp; typically pp
        (*pythia_pp) -> settings.mode( "Beams:idA", idLeptonHadron );
        (*pythia_pp) -> settings.mode( "Beams:idA",
                                       std::copysign( 2212, MSTI_c[11] ) );
        // ln, hn; typically pn
        (*pythia_pn) -> settings.mode( "Beams:idA", idLeptonHadron );
        (*pythia_pn) -> settings.mode( "Beams:idB",
                                       std::copysign( 2112, MSTI_c[11] ) );
        // nl, nh; typically np
        (*pythia_np) -> settings.mode( "Beams:idA",
                                       std::copysign( 2112, MSTI_c[10] ) );
        (*pythia_np) -> settings.mode( "Beams:idB", idLeptonHadron );
        // pl, ph
        (*pythia_nn) -> settings.mode( "Beams:idA",
                                       std::copysign( 2212, MSTI_c[10] ) );
        (*pythia_nn) -> settings.mode( "Beams:idB", idLeptonHadron );
        // nn for AA and AB
        if(    std::abs(MSTI_c[10]) > 1000000000
            && std::abs(MSTI_c[11]) > 1000000000 ) {
            (*pythia_nn) -> settings.mode( "Beams:idA",
                                            std::copysign( 2112, MSTI_c[10] ) );
            (*pythia_nn) -> settings.mode( "Beams:idB",
                                            std::copysign( 2112, MSTI_c[11] ) );
        }

        // Initialization.
        string equal_sign(119,'=');

        std::cout << "\n " << equal_sign
                  << "\n           p + p Instance Initialization"
                  << "\n " << equal_sign << std::endl;
        (*pythia_pp) -> init();

        std::cout << "\n " << equal_sign
                  << "\n           p + n Instance Initialization"
                  << "\n " << equal_sign << std::endl;
        (*pythia_pn) -> init();

        std::cout << "\n " << equal_sign
                  << "\n           n + p Instance Initialization"
                  << "\n " << equal_sign << std::endl;
        (*pythia_np) -> init();

        std::cout << "\n " << equal_sign
                  << "\n           n + n Instance Initialization"
                  << "\n " << equal_sign << std::endl;
        (*pythia_nn) -> init();
    }

// Initialization.
    string equal_sign(119,'=');
    std::cout << "\n " << equal_sign
              << "\n           Basic Instance Initialization (Default p + p )"
              << "\n " << equal_sign << std::endl;

    (*pythia) -> init();

    std::cout << "\n " << equal_sign
              << "\n           End Instance Initialization"
              << "\n " << equal_sign << std::endl;

// Recovers the std::cout .
    std::cout.rdbuf( coutBuff );
    outputFile.close();

    return;
    }

//--------------------------------------------------------------------------

//  This function will generate the core processes event via PYTHIA 8.
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
                   std::shared_ptr< PaciaeHIUserHooks >** paciaeHIUserHooks ) {

// Specifies output target file, i.e. "cout" redirection.
    std::streambuf* coutBuff = std::cout.rdbuf();
    std::ofstream outputFile;
    outputFile.open( "pythia8_report.log", ios::out | ios::app );
    std::streambuf* fileBuff = outputFile.rdbuf();
    std::cout.rdbuf( fileBuff );
    // Recovers cout at the end of function.
    // std::cout.rdbuf( coutBuff );

// Shorthand for some public members of pythia (also static ones).
// The pythia object was passed in from the external calling.
    Event& event = (*pythia) -> event;
    const Info& info = (*pythia) -> info;
    Settings& settings = (*pythia) -> settings;

// Sets the impact parameter generator if Angantyr pA/Ap or AA.
    int& bSampleMode = MINT_c[38];
    int& iExecMode = MSTP_c[190];
    if( (  std::abs(MSTI_c[10]) > 1000000000
        || std::abs(MSTI_c[11]) > 1000000000 )
        && bSampleMode != 3 && (iExecMode == 8 || iExecMode == 9) ) {
        double bp = VINT_c[138];
        (**paciaeHIUserHooks) -> getBGeneratorPtr()
                              -> setImpactParameter( bp, 0.0, 1.0 );
    }

// Event generation.
    int& iProcess = MINT_c[0];
    bool doVetoProcess = false;
    do {
        while ( !( (*pythia) -> next() ) ) ;
        // For the vetoing of the Angantyr machinery.
        if( iExecMode == 8 || iExecMode == 9 ) {
            doVetoProcess = true;
            switch (iProcess) {
                // Inelastic (INEL)
                case 0 :
                    if( info.hiInfo->nCollTot() > info.hiInfo->nCollEL() )
                        doVetoProcess = false;
                    break;
                // Non-Single Difractive (NSD)
                case 1 :
                    if( info.hiInfo->nCollND() > 0 ||
                        info.hiInfo->nCollDD() > 0 ||
                        info.hiInfo->nCollCD() > 0 ) doVetoProcess = false;
                    break;
                // Inelastic nondiffrative. (Minimum Bias)
                case 10 :
                    if( info.hiInfo->nCollND() > 0 ) doVetoProcess = false;
                    break;
                default :
                    doVetoProcess = false;
            }
        }
    } while ( doVetoProcess );

// Store the parton configuration for feeding back to PACIAE.
    nPY8 = 0;
    // nPY8 = event.size();
    for ( int iEntry = 0; iEntry < event.size(); ++iEntry ) {
        nPY8 += 1;
        // Status, id, mothers, daughters, color and anti-color.
        kPY8[0][iEntry] = event[iEntry].status() ;
        kPY8[1][iEntry] = event[iEntry].id() ;
        kPY8[2][iEntry] = event[iEntry].mother1() ;
        kPY8[3][iEntry] = event[iEntry].daughter1() ;
        kPY8[4][iEntry] = event[iEntry].daughter2() ;
        kPY8[5][iEntry] = event[iEntry].mother2() ;   // Note here!
        kPY8[6][iEntry] = event[iEntry].col() ;   // Color.
        kPY8[7][iEntry] = event[iEntry].acol() ;   // Anti-color.
        // 4-momentum, mass, scale and polarization.
        pPY8[0][iEntry] = event[iEntry].px() ;
        pPY8[1][iEntry] = event[iEntry].py() ;
        pPY8[2][iEntry] = event[iEntry].pz() ;
        pPY8[3][iEntry] = event[iEntry].e() ;
        pPY8[4][iEntry] = event[iEntry].m() ;
        pPY8[5][iEntry] = event[iEntry].scale() ;
        pPY8[6][iEntry] = event[iEntry].pol() ;
        // Vertex 4-coordinate and proper life time.
        // Converts units from mm, mm/c to fm, fm/c, respectively.
        vPY8[0][iEntry] = event[iEntry].xProd() * MM2FM ;
        vPY8[1][iEntry] = event[iEntry].yProd() * MM2FM ;
        vPY8[2][iEntry] = event[iEntry].zProd() * MM2FM ;
        vPY8[3][iEntry] = event[iEntry].tProd() * MM2FM ;
        vPY8[4][iEntry] = event[iEntry].tau()   * MM2FM ;
    }
    // Store the junction configuration for feeding back to PACIAE.
    noJuncPY8 = 0;
    // noJuncPY8 = event.sizeJunction();
    for ( int iJun = 0; iJun < event.sizeJunction(); ++iJun ) {
        noJuncPY8 += 1;
        kindJuncPY8[iJun] = event.kindJunction(iJun) ;
        for( int j =0; j < 3; ++j ) {
            colJuncPY8[j][iJun]  = event.colJunction(iJun,j) ;
            endcJuncPY8[j][iJun] = event.endColJunction(iJun,j) ;
            statJuncPY8[j][iJun] = event.statusJunction(iJun,j) ;
        }
    }
    // The information of the collision system.
    VINT_c[396] = event[0].px();
    VINT_c[397] = event[0].py();
    VINT_c[398] = event[0].pz();
    VINT_c[399] = event[0].e();
    // The weight of the event.
    // There are some bugs in PYTHIA8/Anganty pA/AA + HardQCD!
    VINT_c[97] = info.weightSum();
    VINT_c[98] = 1. / info.weight();
    PARI_c[8]  = 1. / info.weight();
    VINT_c[99] = info.weight();
    PARI_c[9]  = info.weight();
    // Totak cross sections generated.
    PARI_c[0]  = info.sigmaGen();
    // Normallized factor of cross section.
    PARI_c[1]  = info.sigmaGen() / info.weightSum();
    // Impact parameter in MPI. Not physical unit like fm, but only 
    //  rescaled so that the average should be unity for minimum-bias 
    //  events (meaning less than that for events with hard processes).
    VINT_c[138] = info.bMPI();
    // More information from Angantyr.
    if( iExecMode == 8 || iExecMode == 9 ) {
        // This is the real impact parameter in fm.
        VINT_c[138] = info.hiInfo -> b();
        VINT_c[391] = info.hiInfo -> phi();
        VINT_c[392] = info.hiInfo -> sigmaTot();
        VINT_c[393] = info.hiInfo -> sigmaTotErr();
        VINT_c[394] = info.hiInfo -> sigmaND();
        VINT_c[395] = info.hiInfo -> sigmaNDErr();
        MINT_c[384] = info.hiInfo -> nCollTot();
        MINT_c[385] = info.hiInfo -> nCollND();
        MINT_c[386] = info.hiInfo -> nCollNDTot();
        MINT_c[387] = info.hiInfo -> nCollSDP();
        MINT_c[388] = info.hiInfo -> nCollSDT();
        MINT_c[389] = info.hiInfo -> nCollDD();
        MINT_c[390] = info.hiInfo -> nCollCD();
        MINT_c[391] = info.hiInfo -> nCollEL();
        MINT_c[392] = info.hiInfo -> nPartProj();
        MINT_c[393] = info.hiInfo -> nAbsProj();
        MINT_c[394] = info.hiInfo -> nDiffProj();
        MINT_c[395] = info.hiInfo -> nElProj();
        MINT_c[396] = info.hiInfo -> nPartTarg();
        MINT_c[397] = info.hiInfo -> nAbsTarg();
        MINT_c[398] = info.hiInfo -> nDiffTarg();
        MINT_c[399] = info.hiInfo -> nElTarg();
    }

    // Recovers the std::cout .
    std::cout.rdbuf( coutBuff );
    outputFile.close();

    return;
    }

//--------------------------------------------------------------------------

//  This function will hadronize the partonic state to hadrons via PYTHIA 8.
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
                                                paciaeHIUserHooks ) {

// Specifies output target file, i.e. "cout" redirection.
    std::streambuf* coutBuff = std::cout.rdbuf();
    std::ofstream outputFile;
    outputFile.open( "pythia8_report.log", ios::out | ios::app );
    std::streambuf* fileBuff = outputFile.rdbuf();
    std::cout.rdbuf( fileBuff );
    // Recovers cout at the end of function.
    // std::cout.rdbuf( coutBuff );

// Shorthand for some public members of pythia (also static ones).
// The pythia object was passed in from the external calling.
    Event& event = (*pythia) -> event;
    // Settings& settings = (*pythia) -> settings;
    // const Info& info = (*pythia) -> info;

// Resets the event record for filling the new one.
    event.clear();   // With following "iParton = 0"
    // event.free();   // With following "iParton = 0"
    // event.reset();   // With following "iParton = 1"
    // Fills parton-level configuration from PACIAE.
    for( int iParton = 0; iParton < nPY8; ++iParton ) {
        int& id        = kPY8[1][iParton] ;
        int& status    = kPY8[0][iParton] ;
        // if( status < 0 && status != -11 ) status = -21;
        // if( status > 0 ) status =  23;
        int& mother1   = kPY8[2][iParton] ;
        int& mother2   = kPY8[5][iParton] ;   // Note here!
        int& daughter1 = kPY8[3][iParton] ;
        int& daughter2 = kPY8[4][iParton] ;
        int& col       = kPY8[6][iParton] ;
        int& acol      = kPY8[7][iParton] ;
        double& px     = pPY8[0][iParton] ;
        double& py     = pPY8[1][iParton] ;
        double& pz     = pPY8[2][iParton] ;
        double& e      = pPY8[3][iParton] ;
        double& m      = pPY8[4][iParton] ;
        double& scale  = pPY8[5][iParton] ;
        double& pol    = pPY8[6][iParton] ;
        int iEntry    = event.append( id, status, mother1, mother2,
                                      daughter1, daughter2, col, acol,
                                      px, py, pz, e, m, scale, pol );
        // Sets 4-coordinate for PACIAE.
        // Converts units from fm, fm/c to mm, mm/c.
        double x   = vPY8[0][iParton] * FM2MM ;
        double y   = vPY8[1][iParton] * FM2MM ;
        double z   = vPY8[2][iParton] * FM2MM ;
        double t   = vPY8[3][iParton] * FM2MM ;
        double tau = vPY8[4][iParton] * FM2MM;   // It will be set by id?
        event[iEntry].vProd(x, y, z, t) ;
        event[iEntry].tau(tau) ;
    }

// Resets the junction record for filling the new one.
    event.clearJunctions();
    // Fills the junction configuration from PACIAE.
    for( int iJun = 0; iJun < noJuncPY8; ++iJun ) {
        int& kind    = kindJuncPY8[iJun] ;
        int& col0    = colJuncPY8[0][iJun] ;
        int& col1    = colJuncPY8[1][iJun] ;
        int& col2    = colJuncPY8[2][iJun] ;
        int iEntry = event.appendJunction( kind, col0, col1, col2  ) ;
        // int iEntry = event.appendJunction( kindJuncPY8[iJun],
        //                                    colJuncPY8[0][iJun],
        //                                    colJuncPY8[1][iJun],
        //                                    colJuncPY8[2][iJun] ) ;
        for( int j = 0; j < 3; ++j ) {
            event.endColJunction( iEntry, j, endcJuncPY8[j][iJun] ) ;
            event.statusJunction( iEntry, j, statJuncPY8[j][iJun] ) ;
        }
    }

// Hadronization. Uses "false", because we have appended the junction
//  configuration.
    if( !(*pythia) -> forceHadronLevel(false) ) {
    // if( ! ( (*pythia) -> forceHadronLevel(true) ) ) {
        // Hadronization failed.
        std::cout << "\n PACIAE event: " << MINT_c[4]
                  << " Hadronization failed! Gives up this event."
                  << std::endl;
        iFail = 2;
        return;
    }

// Store the particle information for returning back to PACIAE.
    nPY8 = 0;
    for ( int iParticle = 0; iParticle < event.size(); ++iParticle ) {
        nPY8 += 1;
        // Status, id, mothers, daughters, color and anti-color.
        kPY8[0][iParticle] = event[iParticle].status() ;
        kPY8[1][iParticle] = event[iParticle].id() ;
        kPY8[2][iParticle] = event[iParticle].mother1() ;
        kPY8[3][iParticle] = event[iParticle].daughter1() ;
        kPY8[4][iParticle] = event[iParticle].daughter2() ;
        kPY8[5][iParticle] = event[iParticle].mother2() ;   // Note here!
        kPY8[6][iParticle] = event[iParticle].col() ;   // Color.
        kPY8[7][iParticle] = event[iParticle].acol() ;   // Anti-color.
        // 4-momentum, mass, scale and polarization.
        pPY8[0][iParticle] = event[iParticle].px() ;
        pPY8[1][iParticle] = event[iParticle].py() ;
        pPY8[2][iParticle] = event[iParticle].pz() ;
        pPY8[3][iParticle] = event[iParticle].e() ;
        pPY8[4][iParticle] = event[iParticle].m() ;
        pPY8[5][iParticle] = event[iParticle].scale() ;
        pPY8[6][iParticle] = event[iParticle].pol() ;
        // Vertex 4-coordinate and proper life time.
        // Converts units from mm, mm/c to fm, fm/c.
        vPY8[0][iParticle] = event[iParticle].xProd() * MM2FM ;
        vPY8[1][iParticle] = event[iParticle].yProd() * MM2FM ;
        vPY8[2][iParticle] = event[iParticle].zProd() * MM2FM ;
        vPY8[3][iParticle] = event[iParticle].tProd() * MM2FM ;
        vPY8[4][iParticle] = event[iParticle].tau()   * MM2FM ;
    }

    // Recovers the std::cout .
    std::cout.rdbuf( coutBuff );
    outputFile.close();

    return;
    }

//--------------------------------------------------------------------------

//  This function will print out cross sections from PYTHIA 8.
    void stat_PY8( int (&MINT_c)[400], int (&MSTP_c)[200], int (&MSTI_c)[200],
                   Pythia** pythia,
                   Pythia** pythia_pp, Pythia** pythia_pn,
                   Pythia** pythia_np, Pythia** pythia_nn ) {

// Specifies output target file, i.e. "cout" redirection.
    std::streambuf* coutBuff = std::cout.rdbuf();
    std::ofstream outputFile;
    outputFile.open( "pythia8_report.log", ios::out | ios::app );
    std::streambuf* fileBuff = outputFile.rdbuf();
    std::cout.rdbuf( fileBuff );

    string equal_sign(119,'=');

    if( (  std::abs(MSTI_c[10]) > 1000000000
        || std::abs(MSTI_c[11]) > 1000000000 )
        && ( MSTP_c[190] == 5 || MSTP_c[190] == 6 ) ) {
        std::cout << "\n\n " << equal_sign
                    << "\n           p + p Instance Statistics"
                    << "\n " << equal_sign << std::endl;
        (*pythia_pp) -> stat();

        std::cout << "\n " << equal_sign
                    << "\n           p + n Instance Statistics"
                    << "\n " << equal_sign << std::endl;
        (*pythia_pn) -> stat();

        std::cout << "\n " << equal_sign
                    << "\n           n + p Instance Statistics"
                    << "\n " << equal_sign << std::endl;
        (*pythia_np) -> stat();

        std::cout << "\n " << equal_sign
                    << "\n           n + n Instance Statistics"
                    << "\n " << equal_sign << std::endl;
        (*pythia_nn) -> stat();
    }

    std::cout << "\n " << equal_sign
                << "\n           Basic Instance Statistics"
                << "\n " << equal_sign << std::endl;

    (*pythia) -> stat();

    std::cout << "\n " << equal_sign
              << "\n           End Instance Statistics"
              << "\n " << equal_sign << std::endl;

    // Recovers the std::cout .
    std::cout.rdbuf( coutBuff );
    outputFile.close();

    return;
    }

//--------------------------------------------------------------------------

//  This is a debug function to test the Pythia object.
    void debug_PY8( ) {
/*  void debug_PY8( int frameType, int idA, int idB, double eCM,
                    int (&MINT_c)[400], double (&VINT_c)[400],
                    int (&MSTU_c)[200], double (&PARU_c)[200],
                    int (&MSTJ_c)[200], double (&PARJ_c)[200],
                    int (&MSTP_c)[200], double (&PARP_c)[200],
                    int (&MSTI_c)[200], double (&PARI_c)[200],
                    int& nPY8, int (&kPY8)[8][300000],
                    double (&pPY8)[7][300000],
                    double (&vPY8)[5][300000],
                    int idStable[100], int (&iFail),
                    Pythia** pythia ) { */

// Specifies output target file, i.e. "cout" redirection.
    std::streambuf* coutBuff = std::cout.rdbuf();
    std::ofstream outputFile;
    outputFile.open( "pythia8_report.log", ios::out | ios::app );
    std::streambuf* fileBuff = outputFile.rdbuf();
    std::cout.rdbuf( fileBuff );
    // Recovers cout at the end of function.
    // std::cout.rdbuf( coutBuff );

// Shorthand for some public members of pythia (also static ones).
    Pythia* pythia = new Pythia();
    Event& event = pythia -> event;
    // Settings& settings = (*pythia) -> settings;

// Lists events.
    std::cout << "\n Before delete in debug_PY8, event:" << std::endl;
    event.list(true);
    // Also lists junctions.
    event.listJunctions();
    // Prints statistics.
    pythia -> stat();


// Release objects avoiding memory segment-fault errors.
    delete pythia;

    // Recovers the std::cout .
    std::cout.rdbuf( coutBuff );
    outputFile.close();

    return;
    }