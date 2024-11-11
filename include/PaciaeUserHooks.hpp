#pragma once
// PaciaeUserHooks.hpp is a part of the PACIAE event generator.
// Copyright (C) 2024 PACIAE Group.
// PACIAE is licensed under the GNU GPL v2 or later, see LICENCE for details.
// Open source: https://github.com/ArcsaberHep/PACIAE4
// Author: An-Ke Lei, September 2024 - November 2024.

// Header file to allow user access to program at different stages.

//                                               By An-Ke at UiO  on 14/09/2024
//                                  Last updated by An-Ke at UiO  on 11/11/2024

// PYTHIA 8 header files.
#include "Pythia8/Pythia.h"
#include "Pythia8/HeavyIons.h"
#include "Pythia8Plugins/ProgressLog.h"

// PACIAE header files.
// #include "***"

// PYTHIA namespace.
using namespace Pythia8;

namespace Paciae4 {

//==========================================================================

//  PACIAE user hook class derived from PYTHIA 8.
class PaciaeUserHooks : public UserHooks {

public:
    PaciaeUserHooks() {}
    ~PaciaeUserHooks() {}

//--------------------------------------------------------------------------

    // Allows changing fragmentation parameters.
    // virtual bool canChangeFragPar() override {
    //     return true;
    // }

//--------------------------------------------------------------------------

    // Does change fragmentation parameters.
    // Input: flavPtr, zPtr, pTPtr, idEnd, m2Had, iParton and posEnd (or
    // negEnd).
    // virtual bool doChangeFragPar( StringFlav*, StringZ*, StringPT*, int,
    //     double, vector<int>, const StringEnd* ) override {
    //     return true;
    // }

};

//==========================================================================

//  PACIAE impact parameter generator class derived from PYTHIA 8.
class PaciaeImpactParameterGenerator : public ImpactParameterGenerator {

public:
    PaciaeImpactParameterGenerator() : b(0.0), phi(0.0), weight(1.0) {}
    ~PaciaeImpactParameterGenerator() {}

//--------------------------------------------------------------------------

//  Sets b, phi and weight parameters via externally passed ones.
    void setImpactParameter( double bIn, double phiIn, double weightIn ) {
        b = bIn;
        phi = phiIn;
        weight = weightIn;
        return;
    }

//--------------------------------------------------------------------------

//  Method for generating impact parameters.
    virtual Vec4 generate( double & weightOut ) const override {
        weightOut = weight;
        return Vec4( b*cos(phi), b*sin(phi), 0.0, 0.0 );
    }

private:
    double b;
    double phi;
    double weight;

};

//==========================================================================

//  PACIAE heavy-ion user hook class derived from PYTHIA 8.
class PaciaeHIUserHooks : public HIUserHooks {

public:
    PaciaeHIUserHooks(int bModeIn) : bMode(bModeIn) {
        bGeneratorPtr = std::make_shared<PaciaeImpactParameterGenerator>();
    }
    ~PaciaeHIUserHooks() {}

//--------------------------------------------------------------------------

    // Overrides hasImpactParameterGenerator and returns true.
    // Allows the external impact parameter generator for Angantyr.
    virtual bool hasImpactParameterGenerator() const override {
        return bMode;
    }
    // Overrides impactParameterGenerator and returns the pointer of the
    //   user-defined impact parameter generator.
    virtual std::shared_ptr<ImpactParameterGenerator> impactParameterGenerator()
    const override {
        return bGeneratorPtr;
    }

//--------------------------------------------------------------------------

    // Interface to access the generator.
    std::shared_ptr<PaciaeImpactParameterGenerator> getBGeneratorPtr() {
        return bGeneratorPtr;
    }

private:
    int bMode;
    std::shared_ptr<PaciaeImpactParameterGenerator> bGeneratorPtr;

};

//==========================================================================

} // end namespace Paciae4