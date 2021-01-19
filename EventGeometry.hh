// EventGeometry Package
//
//  Questions/comments? pkomiske@mit.edu
//
//  Copyright (c) 2019-2021
//  Patrick T. Komiske III, Eric M. Metodiev, Jesse Thaler
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef EVENTGEOMETRY_HH
#define EVENTGEOMETRY_HH

#include "fastjet/PseudoJet.hh"

// include Wasserstein package in the proper namespace
#define BEGIN_EMD_NAMESPACE FASTJET_BEGIN_NAMESPACE namespace contrib { namespace emd {
#define END_EMD_NAMESPACE } } FASTJET_END_NAMESPACE
#define EMDNAMESPACE fastjet::contrib::emd

#include "wasserstein/internal/Event.hh"
#include "wasserstein/internal/PairwiseDistance.hh"
#include "wasserstein/CorrelationDimension.hh"
// EMD.hh incuded at the end, to have access to FastJet classes defined here

BEGIN_EMD_NAMESPACE

////////////////////////////////////////////////////////////////////////////////
// FastJetEvent - an event consisting of FastJet PseudoJets and weights
////////////////////////////////////////////////////////////////////////////////

// FastJet events will derive from this, for checking types later
class FastJetEventBase {};

// PW : a particle weight class (see below)
template<class PW>
struct FastJetEvent : public EventBase<std::vector<PseudoJet>, std::vector<double>>,
                      public FastJetEventBase {
  typedef PW ParticleWeight;
  typedef std::vector<PseudoJet> ParticleCollection;
  typedef std::vector<double> WeightCollection;

  // constructor from PseudoJet, possibly with constituents
  FastJetEvent(const PseudoJet & pj) :
    EventBase<ParticleCollection, WeightCollection>(pj.has_constituents() ? 
                                                    pj.constituents() : 
                                                    ParticleCollection{pj}),
    axis_(pj)
  {}

  // constructor from vector of PseudoJets
  FastJetEvent(const ParticleCollection & pjs) :
    EventBase<ParticleCollection, WeightCollection>(pjs)
  {}

  FastJetEvent() {}

  // name of event
  static std::string name() {
    std::ostringstream oss;
    oss << "FastJetEvent<" << ParticleWeight::name() << '>';
    return oss.str();
  }

  // determine weights
  void ensure_weights() {
    if (!has_weights_) {
      weights_.reserve(particles_.size());
      for (const PseudoJet & pj : particles_) {
        weights_.push_back(ParticleWeight::weight(pj));
        total_weight_ += weights_.back();
      }
      has_weights_ = true;
    }
  }

  // access/set PseudoJet
  PseudoJet & axis() { return axis_; }

private:

  // hold original PseudoJet if given one
  PseudoJet axis_;

}; // FastJetEvent

////////////////////////////////////////////////////////////////////////////////
// FastJetParticleWeight - extracts weight from a PseudoJet
////////////////////////////////////////////////////////////////////////////////

// base class to use for checking types later
struct FastJetParticleWeight {
  typedef double Value;
};

// use pT as weight, most typical choice for hadronic colliders
struct TransverseMomentum : FastJetParticleWeight {
  static std::string name() { return "TransverseMomentum"; }
  static Value weight(const PseudoJet & pj) { return pj.pt(); }
  static void set_weight(PseudoJet & pj, Value w) {
    pj.reset_momentum_PtYPhiM(w, pj.rap(), pj.phi(), pj.m());
  }
};

// use ET as weight, typical for hadronic colliders if mass is relevant
struct TransverseEnergy : FastJetParticleWeight {
  static std::string name() { return "TransverseEnergy"; }
  static Value weight(const PseudoJet & pj) { return pj.Et(); }
  static void set_weight(PseudoJet & pj, Value w) {
    Value pt2(w*w - pj.m2()), pt(pt2 > 0 ? std::sqrt(pt2) : -std::sqrt(-pt2));
    pj.reset_momentum_PtYPhiM(pt, pj.rap(), pj.phi(), pj.m());
  }
};

// use |p3| as weight, typical of e+e- colliders treating pjs as massless
struct Momentum : FastJetParticleWeight {
  static std::string name() { return "Momentum"; }
  static Value weight(const PseudoJet & pj) { return pj.modp(); }
  static void set_weight(PseudoJet & pj, Value w) {
    Value e2(w*w + pj.m2()), e(e2 > 0 ? std::sqrt(e2) : -std::sqrt(-e2));
    pj.reset_momentum(pj.px(), pj.py(), pj.pz(), e);
  }
};

// use E as weight, typical of e+e- colliders
struct Energy : FastJetParticleWeight {
  static std::string name() { return "Energy"; }
  static Value weight(const PseudoJet & pj) { return pj.E(); }
  static void set_weight(PseudoJet & pj, Value w) {
    pj.reset_momentum(pj.px(), pj.py(), pj.pz(), w);
  }
};

////////////////////////////////////////////////////////////////////////////////
// FastJet-specific pairwise distances
////////////////////////////////////////////////////////////////////////////////

// Hadronic Delta_R measure with proper checking for phi
struct DeltaR : public PairwiseDistanceBase<DeltaR, std::vector<PseudoJet>, double> {
  typedef PseudoJet Particle;
  typedef double Value;

  DeltaR(Value R, Value beta) :
    PairwiseDistanceBase<DeltaR, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "DeltaR"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    Value dphiabs(std::fabs(p0.phi() - p1.phi()));
    Value dy(p0.rap() - p1.rap()), dphi(dphiabs > PI ? TWOPI - dphiabs : dphiabs);
    return dy*dy + dphi*dphi;
  }
}; // DeltaR

// Massless dot product measure normalized with transverse momenta
struct HadronicDot : public PairwiseDistanceBase<HadronicDot, std::vector<PseudoJet>, double> {
  typedef PseudoJet Particle;
  typedef double Value;

  HadronicDot(Value R, Value beta) :
    PairwiseDistanceBase<HadronicDot, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "HadronicDot"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    Value d(2*(p0.E()*p1.E() - p0.px()*p1.px() - p0.py()*p1.py() - p0.pz()*p1.pz())/(p0.pt()*p1.pt()));
    return (d > 0 ? d : 0);
  }  
}; // HadronicDot

// Massless dot product measure normalized by total momenta
struct EEDot : public PairwiseDistanceBase<EEDot, std::vector<PseudoJet>, double> {
  typedef PseudoJet Particle;
  typedef double Value;

  EEDot(Value R, Value beta) :
    PairwiseDistanceBase<EEDot, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "EEDot"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    Value d(2 - 2*(p0.px()*p1.px() + p0.py()*p1.py() + p0.pz()*p1.pz())/(p0.modp()*p1.modp()));
    return (d > 0 ? d : 0);
  }
}; // EEDot

// Massive dot product measure normalized with transverse energies
struct HadronicDotMassive : public PairwiseDistanceBase<HadronicDotMassive, std::vector<PseudoJet>, double> {
  typedef PseudoJet Particle;
  typedef double Value;

  HadronicDotMassive(Value R, Value beta) :
    PairwiseDistanceBase<HadronicDotMassive, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "HadronicDotMassive"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    Value d(2*(p0.E()*p1.E() - p0.px()*p1.px() - p0.py()*p1.py() - p0.pz()*p1.pz())/(p0.Et()*p1.Et())
             - p0.m2()/p0.Et2() - p1.m2()/p1.Et2());
    return (d > 0 ? d : 0);
  }
}; // HadronicDotMassive

// Massive dot product measure normalized with energies
struct EEDotMassive : public PairwiseDistanceBase<EEDotMassive, std::vector<PseudoJet>, double> {
  typedef PseudoJet Particle;
  typedef double Value;

  EEDotMassive(Value R, Value beta) :
    PairwiseDistanceBase<EEDotMassive, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "EEDotMassive"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    Value d(2 - 2*(p0.px()*p1.px() + p0.py()*p1.py() + p0.pz()*p1.pz())/(p0.E()*p1.E())
               - p0.m2()/(p0.E()*p0.E()) - p1.m2()/(p1.E()*p1.E()));
    return (d > 0 ? d : 0);
  }
}; // EEDotMassive

// Arc length between momentum vectors
struct EEArcLength : public PairwiseDistanceBase<EEArcLength, std::vector<PseudoJet>, double> {
  typedef PseudoJet Particle;
  typedef double Value;

  EEArcLength(Value R, Value beta) :
    PairwiseDistanceBase<EEArcLength, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "EEArcLength"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    Value dot((p0.px()*p1.px() + p0.py()*p1.py() + p0.pz()*p1.pz())/(p0.modp()*p1.modp()));
    return (dot > 1 ? 0 : (dot < -1 ? pi : std::acos(dot)));
  }
}; // EEArcLength

// Arc length between momentum vectors, normalized by the energy
struct EEArcLengthMassive : public PairwiseDistanceBase<EEArcLengthMassive, std::vector<PseudoJet>, double> {
  typedef PseudoJet Particle;
  typedef double Value;

  EEArcLengthMassive(Value R, Value beta) :
    PairwiseDistanceBase<EEArcLengthMassive, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "EEArcLengthMassive"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    Value dot((p0.px()*p1.px() + p0.py()*p1.py() + p0.pz()*p1.pz())/(p0.E()*p1.E()));
    return (dot > 1 ? 0 : (dot < -1 ? pi : std::acos(dot)));
  }
}; // EEArcLengthMassive

////////////////////////////////////////////////////////////////////////////////
// FastJet-specific Preprocessors
////////////////////////////////////////////////////////////////////////////////

// center all the particles in a vector according to a given rapidity and azimuth
template<class ParticleWeight>
void center_event(FastJetEvent<ParticleWeight> & event,
                  typename ParticleWeight::Value rap,
                  typename ParticleWeight::Value phi) {

  PseudoJet & axis(event.axis());
  axis.reset_momentum_PtYPhiM(axis.pt(), axis.rap() - rap, phi_fix(axis.phi(), phi) - phi, axis.m());
  for (PseudoJet & pj: event.particles())
    pj.reset_momentum_PtYPhiM(pj.pt(), pj.rap() - rap, phi_fix(pj.phi(), phi) - phi, pj.m());
}

////////////////////////////////////////////////////////////////////////////////
// CenterEScheme - center all the particles according to their Escheme axis
////////////////////////////////////////////////////////////////////////////////

template<class EMD>
class CenterEScheme : public Preprocessor<typename EMD::Self> {
public:
  typedef typename EMD::Event Event;

  static_assert(std::is_base_of<FastJetEventBase, Event>::value,
                "CenterEScheme works only with FastJet events.");

  std::string description() const { return "Center according to E-scheme axis"; }
  Event & operator()(Event & event) const {

    // set pj to Escheme axis if it isn't already
  #ifdef __FASTJET_JETDEFINITION_HH__
    if (!event.axis().has_valid_cs() || event.axis().validated_cs()->jet_def().recombination_scheme() != E_scheme)
  #endif
    {
      event.axis().reset_momentum_PtYPhiM(0, 0, 0, 0);
      for (const PseudoJet & pj : event.particles())
        event.axis() += pj;
    }

    // center the particles
    center_event(event, event.axis().rap(), event.axis().phi());

    return event;
  }

}; // CenterEScheme

////////////////////////////////////////////////////////////////////////////////
// CenterPtCentroid - center all the particles according to their pT centroid
////////////////////////////////////////////////////////////////////////////////

template<class EMD>
class CenterPtCentroid : public Preprocessor<typename EMD::Self> {
public:
  typedef typename EMD::Event Event;
  typedef typename EMD::ValuePublic Value;

  static_assert(std::is_base_of<FastJetEventBase, Event>::value,
                "CenterPtCentroid works only with FastJet events.");

  std::string description() const { return "Center according to pT centroid"; }
  Event & operator()(Event & event) const {

    // determine pt centroid
    Value pttot(0), y(0), phi(0);
    for (const PseudoJet & pj : event.particles()) {
      Value pt(pj.pt());
      pttot += pt;
      y += pt * pj.rap();
      phi += pt * phi_fix(pj.phi(), event.particles()[0].phi());
    }
    y /= pttot;
    phi /= pttot;

    // set PtCentroid as axis
    event.axis().reset_momentum_PtYPhiM(pttot, y, phi, 0);

    // center the particles
    center_event(event, y, phi);

    return event;
  }

}; // CenterPtCentroid

////////////////////////////////////////////////////////////////////////////////
// CenterWeightedCentroid - definition of the case specialized for fastjet
////////////////////////////////////////////////////////////////////////////////

template<class EMD>
FastJetEvent<typename EMD::ParticleWeight> & 
CenterWeightedCentroid<EMD>::center(FastJetEvent<typename EMD::ParticleWeight> & event) const {
  event.ensure_weights();

  const ParticleCollection & ps(event.particles());
  const WeightCollection & ws(event.weights());

  // determine weighted centroid
  Value x(0), y(0);
  for (std::size_t i = 0; i < ps.size(); i++) {
    x += ws[i] * ps[i].rap();
    y += ws[i] * phi_fix(ps[i].phi(), ps[0].phi());
  }
  x /= event.total_weight();
  y /= event.total_weight();

  // set PtCentroid as pj of the event
  event.axis().reset_momentum_PtYPhiM(event.total_weight(), x, y, 0);

  // center the particles
  center_event(event, x, y);

  return event;
}

// mask out particles farther than a certain distance from the pseudojet of the event
// since this uses the pseudojet, make sure it is properly set first
template<class EMD>
class MaskCircleRapPhi : public Preprocessor<typename EMD::Self> {
public:
  typedef typename EMD::Event Event;
  typedef typename EMD::ValuePublic Value;

  MaskCircleRapPhi(Value R) : R_(R), R2_(R*R) {}
  std::string description() const { return "Mask particles farther than " + std::to_string(R_) + " from axis"; }
  Event & operator()(Event & event) const {

    std::vector<PseudoJet> & ps(event.particles());

    // get indices of particles to remove
    std::vector<std::size_t> inds;
    for (std::size_t i = 0; i < ps.size(); i++)
      if (EMD::PairwiseDistance::plain_distance(event.axis(), ps[i]) > R2_)
        inds.push_back(i);

    // remove particles and weights if they exist
    if (!inds.empty()) {
      std::reverse(inds.begin(), inds.end());

      for (std::size_t i : inds)
        ps.erase(ps.begin() + i);

      if (event.has_weights()) {
        for (std::size_t i : inds) {
          event.total_weight() -= event.weights()[i];
          event.weights().erase(event.weights().begin() + i);
        }
      }
    }

    return event;
  }

private:

  Value R_, R2_;

}; // MaskCircleRapPhi

END_EMD_NAMESPACE

#include "wasserstein/EMD.hh"

BEGIN_EMD_NAMESPACE

////////////////////////////////////////////////////////////////////////////////
// Aliases
////////////////////////////////////////////////////////////////////////////////

typedef EMD<TransverseMomentum, DeltaR> EMDTransverseMomentumDeltaR;
typedef EMD<TransverseMomentum, HadronicDot> EMDTransverseMomentumHadronicDot;
typedef EMD<TransverseMomentum, HadronicDotMassive> EMDTransverseMomentumHadronicDotMassive;
typedef EMD<TransverseEnergy, DeltaR> EMDTransverseEnergyDeltaR;
typedef EMD<TransverseEnergy, HadronicDot> EMDTransverseEnergyHadronicDot;
typedef EMD<TransverseEnergy, HadronicDotMassive> EMDTransverseEnergyHadronicDotMassive;
typedef EMD<Momentum, EEDot> EMDMomentumEEDot;
typedef EMD<Momentum, EEDotMassive> EMDMomentumEEDotMassive;
typedef EMD<Momentum, EEArcLength> EMDMomentumEEArcLength;
typedef EMD<Momentum, EEArcLengthMassive> EMDMomentumEEArcLengthMassive;
typedef EMD<Energy, EEDot> EMDEnergyEEDot;
typedef EMD<Energy, EEDotMassive> EMDEnergyEEDotMassive;
typedef EMD<Energy, EEArcLength> EMDEnergyEEArcLength;
typedef EMD<Energy, EEArcLengthMassive> EMDEnergyEEArcLengthMassive;

typedef PairwiseEMD<EMD<TransverseMomentum, DeltaR>> PairwiseEMDTransverseMomentumDeltaR;
typedef PairwiseEMD<EMD<TransverseMomentum, HadronicDot>> PairwiseEMDTransverseMomentumHadronicDot;
typedef PairwiseEMD<EMD<TransverseMomentum, HadronicDotMassive>> PairwiseEMDTransverseMomentumHadronicDotMassive;
typedef PairwiseEMD<EMD<TransverseEnergy, DeltaR>> PairwiseEMDTransverseEnergyDeltaR;
typedef PairwiseEMD<EMD<TransverseEnergy, HadronicDot>> PairwiseEMDTransverseEnergyHadronicDot;
typedef PairwiseEMD<EMD<TransverseEnergy, HadronicDotMassive>> PairwiseEMDTransverseEnergyHadronicDotMassive;
typedef PairwiseEMD<EMD<Momentum, EEDot>> PairwiseEMDMomentumEEDot;
typedef PairwiseEMD<EMD<Momentum, EEDotMassive>> PairwiseEMDMomentumEEDotMassive;
typedef PairwiseEMD<EMD<Momentum, EEArcLength>> PairwiseEMDMomentumEEArcLength;
typedef PairwiseEMD<EMD<Momentum, EEArcLengthMassive>> PairwiseEMDMomentumEEArcLengthMassive;
typedef PairwiseEMD<EMD<Energy, EEDot>> PairwiseEMDEnergyEEDot;
typedef PairwiseEMD<EMD<Energy, EEDotMassive>> PairwiseEMDEnergyEEDotMassive;
typedef PairwiseEMD<EMD<Energy, EEArcLength>> PairwiseEMDEnergyEEArcLength;
typedef PairwiseEMD<EMD<Energy, EEArcLengthMassive>> PairwiseEMDEnergyEEArcLengthMassive;

END_EMD_NAMESPACE

#endif // EVENTGEOMETRY_HH