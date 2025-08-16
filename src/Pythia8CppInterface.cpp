// Pythia8CppInterface.cpp is a part of the PACIAE event generator.
// Copyright (C) 2025 PACIAE Group.
// PACIAE is licensed under the GNU GPL v2 or later, see LICENSE for details.
// Open source: https://github.com/ArcsaberHep/PACIAE4
// Author: An-Ke Lei, January 2024 - August 2025.

// This is the C++ interface program to link PACIAE (Fortran 77/90) with
//   PYTHIA 8 (C++).

//                                               By An-Ke at CCNU on 16/01/2024
//                                  Last updated by An-Ke at CCNU on 16/08/2025

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
    ParticleData& particleData = (*pythia) -> particleData;

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
    int& idA = MSTI_c[10];
    int& idB = MSTI_c[11];
    int absIdA = std::abs( MSTI_c[10] );
    int absIdB = std::abs( MSTI_c[11] );
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
            (*pythia) -> readString( "WeakSingleBoson:ffbar2ffbar(s:gmZ) = on");
            (*pythia) -> readString( "23:onMode = off" );
            (*pythia) -> readString( "23:onIfAny = 11 12 13 14 15 16" );
            (*pythia) -> readString( "PhaseSpace:mHatMin = 0.001" );
            break;
        // J/psi (color-singlet via NRQCD)
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
        // c-cbar production
        case 4 :
            (*pythia) -> readString( "HardQCD:hardccbar = on" );
            (*pythia) -> readString( "PhaseSpace:mHatMin = 0.001" );
            break;
        // Prompt photon
        case 5 :
            (*pythia) -> readString( "PromptPhoton:all = on" );
            (*pythia) -> readString( "PhaseSpace:mHatMin = 0.001" );
            break;
        // Soft QCD
        case 6 :
            (*pythia) -> readString( "SoftQCD:all = on" );
            break;
        // Single W+/- production
        case 7 :
            (*pythia) -> readString( "WeakSingleBoson:ffbar2W = on" );
            break;
        // Hard QCD
        case 8 :
            (*pythia) -> readString( "HardQCD:all = on" );
            (*pythia) -> readString( "PhaseSpace:mHatMin = 0.001" );
            break;
        // Single Z0 production
        case 9 :
            (*pythia) -> readString( "WeakSingleBoson:ffbar2gmZ = on" );
            (*pythia) -> readString( "WeakZ0:gmZmode = 2" );
            break;
        // Default.
        case 10 :
            // For hh/hA/AB, inelastic nondiffrative (minimum-bias, MB).
            if( ( particleData.isHadron( idA ) && particleData.isHadron( idB ) )
                || ( particleData.isHadron( idA ) && absIdB > 1000000000 )
                || ( particleData.isHadron( idB ) && absIdA > 1000000000 )
                || ( absIdA > 1000000000 && absIdB > 1000000000 ) ) {
                (*pythia) -> readString( "SoftQCD:nonDiffractive = on" );
            }
            // For lh/lA, deep inel. scatterings(DIS). t-channel boson exchange.
            else if( ( particleData.isLepton(idA) && particleData.isHadron(idB))
                || ( particleData.isLepton(idB) && particleData.isHadron(idA) )
                || ( particleData.isLepton(idA) && absIdB > 1000000000 )
                || ( particleData.isLepton(idB) && absIdA > 1000000000 ) ) {
                (*pythia) -> readString( "WeakBosonExchange:all = on" );
                (*pythia) -> readString( "PhaseSpace:mHatMin = 0.001" );
            }
            // For ll, ll -> gamma*/Z/W -> quark(s) + anti-quark(s).
            else if( particleData.isLepton( idA )
                && particleData.isLepton( idB ) ) {
                // Hadronic gamma*/Z/W decays. Switches off all decays and then
                //  switches back on those to quarks.
                (*pythia) -> readString( "23:onMode = off" );
                (*pythia) -> readString( "23:onIfAny = 1 2 3 4 5" );
                (*pythia) -> readString( "24:onMode = off" );
                (*pythia) -> readString( "24:onIfAny = 1 2 3 4 5" );
                (*pythia) -> readString( "PhaseSpace:mHatMin = 1e-20" );
                // For l+l- annihilation, e.g. e+e-.
                if( ( idA + idB ) == 0 ) {
                    (*pythia) -> readString( 
                                    "WeakSingleBoson:ffbar2ffbar(s:gmZ) = on" );
                    (*pythia) -> readString( 
                                    "WeakDoubleBoson:ffbar2gmZgmZ = on" );
                    (*pythia) -> readString( 
                                    "WeakDoubleBoson:ffbar2WW = on" );
                }
                // l1 + l2bar -> W+- -> q1 + q2bar.
                else {
                    (*pythia) -> readString( "WeakSingleBoson:ffbar2W = on" );
                }
            }
            // Inelastic non-diffractive.
            else {
                (*pythia) -> readString( "SoftQCD:nonDiffractive = on" );
            }
            break;
        default :
            std::cout << "\n PACIAE Warning: You have chosen an undefined " \
                         "nchan = " << iProcess
                      << ".\n Please specify subprocesses in " \
                         "\"pythia8_extra.cfg\".\n"
                      << std::endl;
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

    int absIdA = std::abs( MSTI_c[10] );
    int absIdB = std::abs( MSTI_c[11] );

// Nuclear parton distribution functions.
    int useHardNPDF = MSTP_c[193];
    if( useHardNPDF != -1 ) {
        // Preset 208Pb EPS09 LO nPDF.
        int setNPDFIdA = 100822080;
        int setNPDFIdB = 100822080;
        if( useHardNPDF < -1 || useHardNPDF > 3 ) {
            std::cout << "\n PACIAE Warning: You are using the preset 208Pb " \
                        "EPS09 LO nPDF.\n To set other nPDF, modify " \
                        "\"adj1(5)\" in \"usu.dat\" of PACIAE,\n  and place " \
                        "the appropriate nPDF grid file into\n  " \
                        "your-PYTHIA8-directory/share/Pythia8/pdfdata/ .\n"
                      << std::endl;
            useHardNPDF = 1;
        }
        else {
            setNPDFIdA = absIdA;
            setNPDFIdB = absIdB;
        }
        // Projectile.
        if( absIdA > 99 ) {
            (*pythia) -> readString( "PDF:useHardNPDFA = on" );
            settings.mode( "PDF:nPDFSetA", useHardNPDF );
            settings.mode( "PDF:nPDFBeamA", setNPDFIdA );
        }
        // Target.
        if( absIdB > 99 ) {
            (*pythia) -> readString( "PDF:useHardNPDFB = on" );
            settings.mode( "PDF:nPDFSetB", useHardNPDF );
            settings.mode( "PDF:nPDFBeamB", setNPDFIdB );
        }
    }

// Sets the seed of the random number generator.
    settings.flag( "Random:setSeed", 1 );
    int& seed = MSTP_c[191];
    settings.mode( "Random:seed", seed );

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
    (*pythia) -> readFile( "pythia8_extra.cfg" );

// Instantiates 4 PYTHIA8 objects and does initialization for lA(Al), hA(Ah)
//   AB collisions. Typically pA(Ap) and AA. On heap, not stack!
    if( ( absIdA > 1000000000 || absIdB > 1000000000 )
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
        if( absIdA <= 1000000000 && absIdB > 1000000000 ) {
            idLeptonHadron = MSTI_c[10];
        }
        // Al, Ah; Ap
        else if( absIdA > 1000000000 && absIdB <= 1000000000 ) {
            idLeptonHadron = MSTI_c[11];
        }
        // lp, hp; typically pp
        (*pythia_pp) -> settings.mode( "Beams:idA", idLeptonHadron );
        (*pythia_pp) -> settings.mode( "Beams:idB",
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
        if( absIdA > 1000000000 && absIdB > 1000000000 ) {
            (*pythia_nn) -> settings.mode( "Beams:idA",
                                            std::copysign( 2112, MSTI_c[10] ) );
            (*pythia_nn) -> settings.mode( "Beams:idB",
                                            std::copysign( 2112, MSTI_c[11] ) );
        }

        // Nuclear parton distribution functions.
        // Lepton projectile/target.
        if( absIdA < 100 || absIdB < 100 ) {
            (*pythia_pp) -> readString( "PDF:useHardNPDFA = off" );
            (*pythia_pn) -> readString( "PDF:useHardNPDFA = off" );
            (*pythia_np) -> readString( "PDF:useHardNPDFB = off" );
            (*pythia_nn) -> readString( "PDF:useHardNPDFB = off" );
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
    (*pythia) -> readString( "Init:showChangedparticleData = on" );

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

// Shorthand for some public members of pythia (also static ones).
// The pythia object was passed in from the external calling.
    Event& event = (*pythia) -> event;
    const Info& info = (*pythia) -> info;

// Sets the impact parameter generator if Angantyr pA/Ap or AA.
    int& bSampleMode = MINT_c[38];
    int& iExecMode = MSTP_c[190];
    if( (  std::abs( MSTI_c[10] ) > 1000000000
        || std::abs( MSTI_c[11] ) > 1000000000 )
        && bSampleMode != 3 && ( iExecMode == 8 || iExecMode == 9 ) ) {
        double bp = VINT_c[138];
        double bpMin = VINT_c[359];
        double bpMax = VINT_c[360];
        double phi = 0.0;
        int phiSwitch = 0;
        if( bSampleMode == 4 || bSampleMode == 5 || bSampleMode == 6 ) {
            phi = 2.0 * M_PI * ( (*pythia) -> rndm ).flat();
            phiSwitch = 1;
        }
        (**paciaeHIUserHooks) -> getBGeneratorPtr() -> setImpactParameter( bp,
            bpMin, bpMax, phi, phiSwitch, 1.0 );
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
                    // Old PYTHIA version.
                    #if PYTHIA_VERSION_INTEGER < 8313
                        if( info.hiInfo->nCollTot() > info.hiInfo->nCollEL() )
                            doVetoProcess = false;
                    // From PYTHIA 8.313 onwards.
                    #else
                        if( info.hiInfo->nCollTot() > info.hiInfo->nCollEl() )
                            doVetoProcess = false;
                    #endif
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
    // PYTHIA 8 version number.
    VINT_c[395] = PYTHIA_VERSION;
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
        VINT_c[369] = info.hiInfo -> phi();
        // Number of nucleons.
        MINT_c[384] = info.hiInfo -> nCollTot();
        MINT_c[385] = info.hiInfo -> nCollND();
        // MINT_c[386], see below. Note that there is no "nCollNDTot" from 8313.
        MINT_c[387] = info.hiInfo -> nCollSDP();
        MINT_c[388] = info.hiInfo -> nCollSDT();
        MINT_c[389] = info.hiInfo -> nCollDD();
        MINT_c[390] = info.hiInfo -> nCollCD();
        // MINT_c[391], see below. Note the uppercase "EL" and lowercase "El".
        MINT_c[392] = info.hiInfo -> nPartProj();
        MINT_c[393] = info.hiInfo -> nAbsProj();
        MINT_c[394] = info.hiInfo -> nDiffProj();
        MINT_c[395] = info.hiInfo -> nElProj();
        MINT_c[396] = info.hiInfo -> nPartTarg();
        MINT_c[397] = info.hiInfo -> nAbsTarg();
        MINT_c[398] = info.hiInfo -> nDiffTarg();
        MINT_c[399] = info.hiInfo -> nElTarg();
        // Old PYTHIA version.
        #if PYTHIA_VERSION_INTEGER < 8313
            MINT_c[386] = info.hiInfo -> nCollNDTot();
            MINT_c[391] = info.hiInfo -> nCollEL();
            // Cross sections.
            VINT_c[370] = info.hiInfo -> sigmaTot();
            VINT_c[371] = info.hiInfo -> sigmaND();
            // Errors.
            VINT_c[380] = info.hiInfo -> sigmaTotErr();
            VINT_c[381] = info.hiInfo -> sigmaNDErr();
        // From PYTHIA 8.313 onwards.
        #else
            MINT_c[391] = info.hiInfo -> nCollEl();
            // Cross sections.
            VINT_c[370] = info.hiInfo -> glauberTot();
            VINT_c[371] = info.hiInfo -> glauberND();
            VINT_c[372] = info.hiInfo -> glauberINEL();
            VINT_c[373] = info.hiInfo -> glauberEL();
            VINT_c[374] = info.hiInfo -> glauberDiffP();
            VINT_c[375] = info.hiInfo -> glauberDiffT();
            VINT_c[376] = info.hiInfo -> glauberDDiff();
            VINT_c[377] = info.hiInfo -> glauberBSlope();
            // Errors.
            VINT_c[380] = info.hiInfo -> glauberTotErr();
            VINT_c[381] = info.hiInfo -> glauberNDErr();
            VINT_c[382] = info.hiInfo -> glauberINELErr();
            VINT_c[383] = info.hiInfo -> glauberELErr();
            VINT_c[384] = info.hiInfo -> glauberDiffPErr();
            VINT_c[385] = info.hiInfo -> glauberDiffTErr();
            VINT_c[386] = info.hiInfo -> glauberDDiffErr();
            VINT_c[387] = info.hiInfo -> glauberBSlopeErr();
        #endif
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

// Shorthand for some public members of pythia (also static ones).
// The pythia object was passed in from the external calling.
    Event& event = (*pythia) -> event;

// Resets the event record for filling the new one.
    event.clear();
    // Fills parton-level configuration from PACIAE.
    for( int iParton = 0; iParton < nPY8; ++iParton ) {
        int& id        = kPY8[1][iParton] ;
        int& status    = kPY8[0][iParton] ;
        int& mother1   = kPY8[2][iParton] ;
        int& mother2   = kPY8[5][iParton] ;
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
        double tau = vPY8[4][iParton] * FM2MM ;
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
        kPY8[5][iParticle] = event[iParticle].mother2() ;
        kPY8[6][iParticle] = event[iParticle].col() ;
        kPY8[7][iParticle] = event[iParticle].acol() ;
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

//  This function will perform sequential decays via PYTHIA 8.
    void decayAll_PY8( int& nPY8, int (&kPY8)[8][300000],
                       double (&pPY8)[7][300000],
                       double (&vPY8)[5][300000],
                       Pythia** pythia ) {

// Specifies output target file, i.e. "cout" redirection.
    std::streambuf* coutBuff = std::cout.rdbuf();
    std::ofstream outputFile;
    outputFile.open( "pythia8_report.log", ios::out | ios::app );
    std::streambuf* fileBuff = outputFile.rdbuf();
    std::cout.rdbuf( fileBuff );

// Shorthand for some public members of pythia (also static ones).
// The pythia object was passed in from the external calling.
    Event& event = (*pythia) -> event;
    Rndm&   rndm = (*pythia) -> rndm;

// Resets the event record for filling the new one.
    if( std::abs( kPY8[0][0] ) == 11 ) {
        event.reset();
    }
    else {
        event.clear();
    }
    for( int iParticle = 0; iParticle < nPY8; ++iParticle ) {
        int& id        = kPY8[1][iParticle] ;
        int& status    = kPY8[0][iParticle] ;
        int& mother1   = kPY8[2][iParticle] ;
        int& mother2   = kPY8[5][iParticle] ;
        int& daughter1 = kPY8[3][iParticle] ;
        int& daughter2 = kPY8[4][iParticle] ;
        int& col       = kPY8[6][iParticle] ;
        int& acol      = kPY8[7][iParticle] ;
        double& px     = pPY8[0][iParticle] ;
        double& py     = pPY8[1][iParticle] ;
        double& pz     = pPY8[2][iParticle] ;
        double& e      = pPY8[3][iParticle] ;
        double& m      = pPY8[4][iParticle] ;
        double& scale  = pPY8[5][iParticle] ;
        double& pol    = pPY8[6][iParticle] ;
        int iEntry    = event.append( id, status, mother1, mother2,
                                      daughter1, daughter2, col, acol,
                                      px, py, pz, e, m, scale, pol );
        // Converts units from fm, fm/c to mm, mm/c.
        double x   = vPY8[0][iParticle] * FM2MM ;
        double y   = vPY8[1][iParticle] * FM2MM ;
        double z   = vPY8[2][iParticle] * FM2MM ;
        double t   = vPY8[3][iParticle] * FM2MM ;
        event[iEntry].vProd(x, y, z, t) ;
        // Generates lifetime, to give decay away from primary vertex.
        double tau = event[iEntry].tau0() * rndm.exp() ;
        event[iEntry].tau(tau) ;
    }

// Resets the junction record.
    event.clearJunctions();

// Decays of unstable particles.
    (*pythia) -> forceHadronLevel(false) ;

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
        kPY8[5][iParticle] = event[iParticle].mother2() ;
        kPY8[6][iParticle] = event[iParticle].col() ;
        kPY8[7][iParticle] = event[iParticle].acol() ;
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