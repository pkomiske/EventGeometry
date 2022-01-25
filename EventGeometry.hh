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

// handle using PyFJCore for PseudoJet
#ifdef EVENTGEOMETRY_USE_PYFJCORE
# include "pyfjcore/fjcore.hh"
#elif !defined(__FASTJET_PSEUDOJET_HH__) && !defined(__FJCORE__)
# include "fastjet/PseudoJet.hh"
#endif

// include Wasserstein package in the proper namespace
#ifndef BEGIN_WASSERSTEIN_NAMESPACE
# define BEGIN_WASSERSTEIN_NAMESPACE namespace fastjet { namespace contrib { namespace eventgeometry {
# define END_WASSERSTEIN_NAMESPACE } } }
# define WASSERSTEIN_NAMESPACE fastjet::contrib::eventgeometry
#endif

// macro for controling template visibility
#ifndef EVENTGEOMETRY_TEMPLATE_VISIBILITY
# define EVENTGEOMETRY_TEMPLATE_VISIBILITY extern
#endif

// handle Wasserstein's template visibility
#ifndef WASSERSTEIN_TEMPLATE_VISIBILITY
# define WASSERSTEIN_TEMPLATE_VISIBILITY EVENTGEOMETRY_TEMPLATE_VISIBILITY
#endif

// macros for declaring templated types
#define EVENTGEOMETRY_TEMPLATE(...) EVENTGEOMETRY_TEMPLATE_VISIBILITY template __VA_ARGS__;
#define EVENTGEOMETRY_TEMPLATE_CLASS(...) EVENTGEOMETRY_TEMPLATE_VISIBILITY template class __VA_ARGS__;
#define EVENTGEOMETRY_TEMPLATE_STRUCT(...) EVENTGEOMETRY_TEMPLATE_VISIBILITY template struct __VA_ARGS__;

#define EVENTGEOMETRY_PREPROCESSOR_TEMPLATE(...) \
  EVENTGEOMETRY_TEMPLATE_CLASS(CenterEScheme<__VA_ARGS__>) \
  EVENTGEOMETRY_TEMPLATE_CLASS(CenterPtCentroid<__VA_ARGS__>) \
  EVENTGEOMETRY_TEMPLATE_CLASS(CenterWeightedCentroid<__VA_ARGS__>) \
  EVENTGEOMETRY_TEMPLATE_CLASS(MaskCircle<__VA_ARGS__>)

// macro for declaring specific EMD template
#define EVENTGEOMETRY_EMD_TEMPLATE(PW, PD) \
  EVENTGEOMETRY_TEMPLATE_CLASS(EMD<double, PW, PD>) \
  EVENTGEOMETRY_TEMPLATE_CLASS(PairwiseEMD<EMD<double, PW, PD>>) \
  EVENTGEOMETRY_PREPROCESSOR_TEMPLATE(EMD<double, PW, PD>) \
  EVENTGEOMETRY_TEMPLATE(double EMD<double, PW, PD>::operator()(const std::vector<PseudoJet> &, const std::vector<PseudoJet> &)) \
  EVENTGEOMETRY_TEMPLATE(double EMD<double, PW, PD>::operator()(const PseudoJet &, const PseudoJet &)) \
  EVENTGEOMETRY_TEMPLATE(double EMD<double, PW, PD>::operator()(const std::vector<PseudoJet> &, const PseudoJet &)) \
  EVENTGEOMETRY_TEMPLATE(double EMD<double, PW, PD>::operator()(const PseudoJet &, const std::vector<PseudoJet> &)) \
  using EMD##PW##PD = EMDFloat64<PW, PD>; \
  using PairwiseEMD##PW##PD = PairwiseEMD<EMD##PW##PD>;

#define EVENTGEOMETRY_EVENT_TEMPLATE(E) \
  EVENTGEOMETRY_TEMPLATE_STRUCT(E) \
  EVENTGEOMETRY_TEMPLATE_STRUCT(FastJetEvent<E>)

#define EVENTGEOMETRY_PAIRWISEDISTANCE_TEMPLATE(PD) \
  EVENTGEOMETRY_TEMPLATE_CLASS(PD) \
  EVENTGEOMETRY_TEMPLATE_CLASS(PairwiseDistanceBase<PD, std::vector<PseudoJet>, double>)

#define EVENTGEOMETRY_PAIRWISEDISTANCE_NONREDUCED_TEMPLATE(PD) \
  EVENTGEOMETRY_PAIRWISEDISTANCE_TEMPLATE(PD) \
  EVENTGEOMETRY_TEMPLATE_CLASS(PairwiseDistanceBaseNonReduced<PD, std::vector<PseudoJet>, double>)  

// declare all EMD templates
#define EVENTGEOMETRY_TEMPLATES \
  WASSERSTEIN_TEMPLATES \
  EVENTGEOMETRY_TEMPLATE_STRUCT(EventBase<std::vector<double>, std::vector<PseudoJet>>) \
  EVENTGEOMETRY_EVENT_TEMPLATE(TransverseMomentum<double>) \
  EVENTGEOMETRY_EVENT_TEMPLATE(TransverseEnergy<double>) \
  EVENTGEOMETRY_EVENT_TEMPLATE(Momentum<double>) \
  EVENTGEOMETRY_EVENT_TEMPLATE(Energy<double>) \
  EVENTGEOMETRY_PAIRWISEDISTANCE_TEMPLATE(DeltaR<double>) \
  EVENTGEOMETRY_PAIRWISEDISTANCE_TEMPLATE(HadronicDot<double>) \
  EVENTGEOMETRY_PAIRWISEDISTANCE_TEMPLATE(HadronicDotMassive<double>) \
  EVENTGEOMETRY_PAIRWISEDISTANCE_TEMPLATE(EEDot<double>) \
  EVENTGEOMETRY_PAIRWISEDISTANCE_TEMPLATE(EEDotMassless<double>) \
  EVENTGEOMETRY_PAIRWISEDISTANCE_NONREDUCED_TEMPLATE(EEArcLength<double>) \
  EVENTGEOMETRY_PAIRWISEDISTANCE_NONREDUCED_TEMPLATE(EEArcLengthMassive<double>) \
  EVENTGEOMETRY_PAIRWISEDISTANCE_NONREDUCED_TEMPLATE(CosPhi<double>) \
  EVENTGEOMETRY_EMD_TEMPLATE(TransverseMomentum, DeltaR) \
  EVENTGEOMETRY_EMD_TEMPLATE(TransverseMomentum, HadronicDot) \
  EVENTGEOMETRY_EMD_TEMPLATE(TransverseMomentum, HadronicDotMassive) \
  EVENTGEOMETRY_EMD_TEMPLATE(TransverseEnergy, DeltaR) \
  EVENTGEOMETRY_EMD_TEMPLATE(TransverseEnergy, HadronicDot) \
  EVENTGEOMETRY_EMD_TEMPLATE(TransverseEnergy, HadronicDotMassive) \
  EVENTGEOMETRY_EMD_TEMPLATE(Momentum, EEDot) \
  EVENTGEOMETRY_EMD_TEMPLATE(Momentum, EEDotMassless) \
  EVENTGEOMETRY_EMD_TEMPLATE(Momentum, EEArcLength) \
  EVENTGEOMETRY_EMD_TEMPLATE(Momentum, EEArcLengthMassive) \
  EVENTGEOMETRY_EMD_TEMPLATE(Momentum, CosPhi) \
  EVENTGEOMETRY_EMD_TEMPLATE(Energy, EEDot) \
  EVENTGEOMETRY_EMD_TEMPLATE(Energy, EEDotMassless) \
  EVENTGEOMETRY_EMD_TEMPLATE(Energy, EEArcLength) \
  EVENTGEOMETRY_EMD_TEMPLATE(Energy, EEArcLengthMassive) \
  EVENTGEOMETRY_EMD_TEMPLATE(Energy, CosPhi)

// Wasserstein options
#define WASSERSTEIN_FASTJET
#define WASSERSTEIN_NO_FLOAT32

// the Wasserstein library (with FastJet compatibility)
#include "wasserstein/Wasserstein.hh"


BEGIN_WASSERSTEIN_NAMESPACE

////////////////////////////////////////////////////////////////////////////////
// FastJetEvent - an event consisting of FastJet PseudoJets and weights
////////////////////////////////////////////////////////////////////////////////

// FastJet events will derive from this, for checking types later
struct FastJetEventBase {};

// PW : a particle weight class (see below)
template<class _ParticleWeight>
struct FastJetEvent : public EventBase<std::vector<typename _ParticleWeight::value_type>, std::vector<PseudoJet>>,
                      public FastJetEventBase {

  typedef _ParticleWeight ParticleWeight;
  typedef typename ParticleWeight::value_type value_type;
  typedef std::vector<PseudoJet> ParticleCollection;
  typedef std::vector<value_type> WeightCollection;

  // constructor from PseudoJet, possibly with constituents
  FastJetEvent(const PseudoJet & pj, value_type event_weight = 1) :
    EventBase<WeightCollection, ParticleCollection>(
        pj.has_constituents() ?
        static_cast<const std::vector<PseudoJet> &>(pj.constituents()) :
        ParticleCollection{pj},
      event_weight)
  {
    axis().reset_momentum(pj);
  }

  // constructor from vector of PseudoJets
  FastJetEvent(const ParticleCollection & pjs, value_type event_weight = 1) :
    EventBase<WeightCollection, ParticleCollection>(pjs, event_weight)
  {}

  FastJetEvent() = default;
  virtual ~FastJetEvent() = default;

  // this ensures that the PseudoJets stored in the events can be accessed in a thread-safe manner
  virtual ParticleCollection preprocess_particles(const ParticleCollection & pjs) const {
    std::vector<PseudoJet> pjs_out(pjs.size());
    for (std::size_t i = 0; i < pjs.size(); i++)
      pjs_out[i].reset_momentum(pjs[i]);

    return pjs_out;
  }

  // name of event
  static std::string name() {
    std::ostringstream oss;
    oss << "FastJetEvent<" << ParticleWeight::name() << '>';
    return oss.str();
  }

  // determine weights
  void ensure_weights() {
    if (!this->has_weights()) {
      this->weights().reserve(this->particles().size());
      for (const PseudoJet & pj : this->particles()) {
        this->weights().push_back(ParticleWeight::weight(pj));
        this->total_weight() += this->weights().back();
      }
      this->has_weights_ = true;
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
struct FastJetParticleWeight {};

// use pT as weight, most typical choice for hadronic colliders
template<typename Value>
struct TransverseMomentum : FastJetParticleWeight {
  typedef Value value_type;
  static std::string name() { return "TransverseMomentum"; }
  static Value weight(const PseudoJet & pj) { return pj.pt(); }
  static void set_weight(PseudoJet & pj, double w) {
    pj.reset_momentum_PtYPhiM(w, pj.rap(), pj.phi(), pj.m());
  }
};

// use ET as weight, typical for hadronic colliders if mass is relevant
template<typename Value>
struct TransverseEnergy : FastJetParticleWeight {
  typedef Value value_type;
  static std::string name() { return "TransverseEnergy"; }
  static Value weight(const PseudoJet & pj) { return pj.Et(); }
  static void set_weight(PseudoJet & pj, double w) {
    double pt2(w*w - pj.m2()), pt(pt2 > 0 ? std::sqrt(pt2) : -std::sqrt(-pt2));
    pj.reset_momentum_PtYPhiM(pt, pj.rap(), pj.phi(), pj.m());
  }
};

// use |p3| as weight, typical of e+e- colliders treating pjs as massless
template<typename Value>
struct Momentum : FastJetParticleWeight {
  typedef Value value_type;
  static std::string name() { return "Momentum"; }
  static Value weight(const PseudoJet & pj) { return pj.modp(); }
  static void set_weight(PseudoJet & pj, double w) {
    double e2(w*w + pj.m2()), e(e2 > 0 ? std::sqrt(e2) : -std::sqrt(-e2));
    pj.reset_momentum(pj.px(), pj.py(), pj.pz(), e);
  }
};

// use E as weight, typical of e+e- colliders
template<typename Value>
struct Energy : FastJetParticleWeight {
  typedef Value value_type;
  static std::string name() { return "Energy"; }
  static Value weight(const PseudoJet & pj) { return pj.E(); }
  static void set_weight(PseudoJet & pj, double w) {
    pj.reset_momentum(pj.px(), pj.py(), pj.pz(), w);
  }
};


////////////////////////////////////////////////////////////////////////////////
// FastJet-specific pairwise distances
////////////////////////////////////////////////////////////////////////////////

// for those distances which don't naturally need a square root operation (i.e. the arc length ones)
template <class PairwiseDistance, class ParticleCollection, typename Value>
class PairwiseDistanceBaseNonReduced : public PairwiseDistanceBase<PairwiseDistance, ParticleCollection, Value> {
public:

  typedef typename ParticleCollection::const_iterator ParticleIterator;
  using PairwiseDistanceBase<PairwiseDistance, ParticleCollection, Value>::PairwiseDistanceBase;

  // returns the distance divided by R, all to beta power
  Value distance(const ParticleIterator & p0, const ParticleIterator & p1) const {
    Value pd(PairwiseDistance::plain_distance_(p0, p1));
    
    if (this->beta() == 1.0)
      return pd/this->R();

    if (this->beta() == 2.0)
      return pd*pd/(this->R()*this->R());

    return std::pow(pd/this->R(), this->beta());
  }

protected:

  ~PairwiseDistanceBaseNonReduced() = default;

}; // PairwiseDistanceBaseNonReduced

// Hadronic Delta_R measure with proper checking for phi
template<typename Value>
class DeltaR : public PairwiseDistanceBase<DeltaR<Value>, std::vector<PseudoJet>, Value> {
public:
  typedef PseudoJet Particle;

  DeltaR(Value R, Value beta) :
    PairwiseDistanceBase<DeltaR<Value>, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "DeltaR"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    double dphiabs(std::fabs(p0.phi() - p1.phi()));
    double dy(p0.rap() - p1.rap()), dphi(dphiabs > PI ? TWOPI - dphiabs : dphiabs);
    return dy*dy + dphi*dphi;
  }
}; // DeltaR

// Dot product measure normalized with transverse momenta
template<typename Value>
class HadronicDot : public PairwiseDistanceBase<HadronicDot<Value>, std::vector<PseudoJet>, Value> {
public:
  typedef PseudoJet Particle;

  HadronicDot(Value R, Value beta) :
    PairwiseDistanceBase<HadronicDot<Value>, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "HadronicDot"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    return std::max(2*fastjet::dot_product(p0, p1) / std::sqrt(p0.pt2()*p1.pt2()), 0.0);
  }  
}; // HadronicDot

// Dot product measure normalized by energy
template<typename Value>
class EEDot : public PairwiseDistanceBase<EEDot<Value>, std::vector<PseudoJet>, Value> {
public:
  typedef PseudoJet Particle;

  EEDot(Value R, Value beta) :
    PairwiseDistanceBase<EEDot<Value>, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "EEDot"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    return std::max(2*fastjet::dot_product(p0, p1) / (p0.E()*p1.E()), 0.0);
  }
}; // EEDot

// Dot product measure normalized with transverse energies
template<typename Value>
class HadronicDotMassive : public PairwiseDistanceBase<HadronicDotMassive<Value>, std::vector<PseudoJet>, Value> {
public:
  typedef PseudoJet Particle;

  HadronicDotMassive(Value R, Value beta) :
    PairwiseDistanceBase<HadronicDotMassive<Value>, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "HadronicDotMassive"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    return std::max(2*fastjet::dot_product(p0, p1) / std::sqrt(p0.Et2()*p1.Et2()), 0.0);
  }
}; // HadronicDotMassive

// Dot product measure normalized with momenta
template<typename Value>
class EEDotMassless : public PairwiseDistanceBase<EEDotMassless<Value>, std::vector<PseudoJet>, Value> {
public:
  typedef PseudoJet Particle;

  EEDotMassless(Value R, Value beta) :
    PairwiseDistanceBase<EEDotMassless<Value>, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "EEDotMassless"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    return std::max(2*fastjet::dot_product(p0, p1) / std::sqrt(p0.modp2()*p1.modp2()), 0.0);
  }
}; // EEDotMassless

// Arc length between momentum vectors
template<typename Value>
class EEArcLength : public PairwiseDistanceBaseNonReduced<EEArcLength<Value>, std::vector<PseudoJet>, Value> {
public:
  typedef PseudoJet Particle;
  typedef std::vector<PseudoJet>::const_iterator ParticleIterator;

  EEArcLength(Value R, Value beta) :
    PairwiseDistanceBaseNonReduced<EEArcLength<Value>, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "EEArcLength"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    return fastjet::theta(p0, p1);
  }
}; // EEArcLength

// Arc length between momentum vectors, normalized by the energy
template<typename Value>
class EEArcLengthMassive : public PairwiseDistanceBaseNonReduced<EEArcLengthMassive<Value>, std::vector<PseudoJet>, Value> {
public:
  typedef PseudoJet Particle;
  typedef std::vector<PseudoJet>::const_iterator ParticleIterator;

  EEArcLengthMassive(Value R, Value beta) :
    PairwiseDistanceBaseNonReduced<EEArcLengthMassive<Value>, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "EEArcLengthMassive"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    return std::acos(std::min(1.0, std::max(-1.0, (p0.px()*p1.px() + p0.py()*p1.py() + p0.pz()*p1.pz())/(p0.E()*p1.E()))));
  }
}; // EEArcLengthMassive

// Cosine of azimuthal angle difference, for e.g. ring-like isotropy (2004.06125)
template<typename Value>
class CosPhi : public PairwiseDistanceBaseNonReduced<CosPhi<Value>, std::vector<PseudoJet>, Value> {
public:
  typedef PseudoJet Particle;
  typedef std::vector<PseudoJet>::const_iterator ParticleIterator;

  CosPhi(Value R, Value beta) :
    PairwiseDistanceBaseNonReduced<CosPhi<Value>, std::vector<PseudoJet>, double>(R, beta)
  {}
  static std::string name() { return "CosPhi"; }
  static Value plain_distance(const PseudoJet & p0, const PseudoJet & p1) {
    double dphiabs(std::fabs(p0.phi() - p1.phi()));
    double dphi(dphiabs > PI ? TWOPI - dphiabs : dphiabs);
    return (1-std::cos(dphi));
  }
}; // CosPhi

////////////////////////////////////////////////////////////////////////////////
// FastJet-specific Preprocessors
////////////////////////////////////////////////////////////////////////////////

inline double phi_fix(double phi, double ref_phi) {
  double diff(phi - ref_phi);
  if (diff > PI) phi -= TWOPI;
  else if (diff < -PI) phi += TWOPI;
  return phi; 
}

// center all the particles in a vector according to a given rapidity and azimuth
template<class ParticleWeight>
void center_event(FastJetEvent<ParticleWeight> & event,
                  typename ParticleWeight::value_type rap,
                  typename ParticleWeight::value_type phi) {

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

    // calculate E scheme axis
    event.axis().reset_momentum_PtYPhiM(0, 0, 0, 0);
    for (const PseudoJet & pj : event.particles())
      event.axis() += pj;

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
  typedef typename EMD::value_type value_type;

  static_assert(std::is_base_of<FastJetEventBase, Event>::value,
                "CenterPtCentroid works only with FastJet events.");

  std::string description() const { return "Center according to pT centroid"; }
  Event & operator()(Event & event) const {

    // determine pt centroid
    value_type pttot(0), y(0), phi(0);
    for (const PseudoJet & pj : event.particles()) {
      value_type pt(pj.pt());
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

  typedef typename EMD::value_type value_type;

  event.ensure_weights();
  const ParticleCollection & ps(event.particles());
  const WeightCollection & ws(event.weights());

  // determine weighted centroid
  value_type x(0), y(0);
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


////////////////////////////////////////////////////////////////////////////////
// MaskCircle - mask out particles farther than a certain distance from the axis of the event
////////////////////////////////////////////////////////////////////////////////

template<class EMD>
class MaskCircle : public Preprocessor<typename EMD::Self> {
public:
  typedef typename EMD::Event Event;
  typedef typename EMD::value_type value_type;

  MaskCircle(double R) : R_(R), R2_(R*R) {}
  std::string description() const {
    std::ostringstream oss;
    oss << "Mask particles farther than " << R_ << " from axis";
    return oss.str();
  }

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
        if (event.total_weight() < 0)
          event.total_weight() = 0;
      }
    }

    return event;
  }

private:

  value_type R_, R2_;

}; // MaskCircle

// declares templates, by default extern (so that library must be linked)
#ifdef DECLARE_EVENTGEOMETRY_TEMPLATES
  EVENTGEOMETRY_TEMPLATES
#endif

END_WASSERSTEIN_NAMESPACE

#endif // EVENTGEOMETRY_HH
