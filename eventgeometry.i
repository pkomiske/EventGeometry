// -*- C++ -*-
%define EVENTGEOMETRY_DOCSTRING
"# EventGeometry FastJet Contrib
"
%enddef

%module("docstring"=EVENTGEOMETRY_DOCSTRING, "threads"=1) eventgeometry
%nothreadallow;

// C++ standard library wrappers
%include <std_vector.i>

// this makes SWIG aware of the types contained in the main fastjet library
// but does not generate new wrappers for them here
#ifdef FASTJET_PREFIX
%import FASTJET_PREFIX/share/fastjet/pyinterface/fastjet.i
#else
%import PyFJCore/pyfjcore/swig/pyfjcore.i
%init %{
fastjet::set_pseudojet_format(fastjet::PseudoJetRepresentation::ptyphim);
%}
#endif

// converts fastjet::Error into a FastJetError Python exception
FASTJET_ERRORS_AS_PYTHON_EXCEPTIONS(eventgeometry)

// turn off exception handling for now, since fastjet::Error is not thrown here
%exception;

// include headers in source file
%{
#ifndef SWIG
#define SWIG
#endif

#define WASSERSTEIN_FASTJET
#include "EventGeometry.hh"
%}

// include EMD wrappers from Wasserstein package
#define WASSERSTEIN_NAMESPACE fastjet::contrib::eventgeometry
#define BEGIN_WASSERSTEIN_NAMESPACE namespace fastjet { namespace contrib { namespace eventgeometry {
#define END_WASSERSTEIN_NAMESPACE } } }
#define WASSERSTEIN_NO_FLOAT32
%include "wasserstein/swig/wasserstein_common.i"

// this needs to come after numpy has been included by wasserstein_common.i
%{
#include "pyfjcore/PyFJCoreExtensions.hh"
%}

%define EVENTGEOMETRY_ADD_EXPLICIT_PREPROCESSORS
void preprocess_CenterEScheme() { $self->preprocess<WASSERSTEIN_NAMESPACE::CenterEScheme>(); }
//void preprocess_CenterWeightedCentroid() { $self->preprocess<CenterWeightedCentroid>(); }
void preprocess_CenterPtCentroid() { $self->preprocess<WASSERSTEIN_NAMESPACE::CenterPtCentroid>(); }
void preprocess_MaskCircleRapPhi(double R) { $self->preprocess<WASSERSTEIN_NAMESPACE::MaskCircleRapPhi>(R); }
%enddef

%ignore WASSERSTEIN_NAMESPACE::FastJetEventBase;
%ignore WASSERSTEIN_NAMESPACE::FastJetParticleWeight;

%template(vectorPseudoJetContainer) std::vector<fastjet::PseudoJetContainer>;
%template(vectorVectorPseudoJetr) std::vector<std::vector<fastjet::PseudoJet>>;

%include "EventGeometry.hh"

namespace WASSERSTEIN_NAMESPACE {

  // extend [Pairwise]EMD to explicitly support PseudoJet arguments and printing
  %extend EMD {
    EVENTGEOMETRY_ADD_EXPLICIT_PREPROCESSORS

    double operator()(const PseudoJetContainer & pjs0, const PseudoJetContainer & pjs1) {
      return (*$self)(pjs0, pjs1);
    }
    double operator()(const std::vector<fastjet::PseudoJet> & pjs0, const std::vector<fastjet::PseudoJet> & pjs1) {
      return (*$self)(pjs0, pjs1);
    }
    double operator()(const fastjet::PseudoJet & pj0, const fastjet::PseudoJet & pj1) {
      return (*$self)(pj0, pj1);
    }

    /*double operator()(const std::vector<fastjet::PseudoJet> & pjs0, const fastjet::PseudoJet & pj1) {
      return (*$self)(pjs0, pj1);
    }
    double operator()(const fastjet::PseudoJet & pj0, const std::vector<fastjet::PseudoJet> & pjs1) {
      return (*$self)(pj0, pjs1);
    }*/
  }

  %extend PairwiseEMD {
    EVENTGEOMETRY_ADD_EXPLICIT_PREPROCESSORS

    // add a single event to the PairwiseEMD object
    void _add_event(const PseudoJetContainer & pjc, double event_weight = 1) {
      $self->events().emplace_back(pjc, event_weight);
      $self->preprocess_back_event();
    }
    void _add_event(const std::vector<PseudoJet> & pjs, double event_weight = 1) {
      $self->events().emplace_back(pjs, event_weight);
      $self->preprocess_back_event();
    }
    void _add_event(const PseudoJet & pj, double event_weight = 1) {
      $self->events().emplace_back(pj, event_weight);
      $self->preprocess_back_event();
    }
  }

  // instantiate EMD and PairwiseEMD templates
  %define EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Weight, Distance)
    %template(EMD##Weight##Distance) EMD<double, Weight, Distance>;
    %template(PairwiseEMD##Weight##Distance) PairwiseEMD<EMD<double, Weight, Distance>, double>;
  %enddef

  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(TransverseMomentum, DeltaR)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(TransverseMomentum, HadronicDot)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(TransverseMomentum, HadronicDotMassive)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(TransverseEnergy, DeltaR)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(TransverseEnergy, HadronicDot)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(TransverseEnergy, HadronicDotMassive)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Momentum, EEDot)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Momentum, EEDotMassive)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Momentum, EEArcLength)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Momentum, EEArcLengthMassive)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Energy, EEDot)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Energy, EEDotMassive)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Energy, EEArcLength)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Energy, EEArcLengthMassive)
}

// add convenience functions for accessing templated EMD, PairwiseEMD, and ApolloniusGroomer classes
%pythoncode %{

from fastjet import FastJetError

__version__ = '1.0.0a2'

# function for storing events in a pairwise_emd object
def _store_events(pairwise_emd, events, event_weights, gdim, mask):

    if gdim or mask:
        raise ValueError('`gdim` and `mask` are not supported')

    for event, event_weight in zip(events, event_weights):
        pairwise_emd._add_event(event, event_weight)

def EMD(*args, weight='TransverseMomentum', pairwise_distance='DeltaR', **kwargs):

    if weight == 'TransverseMomentum':
        if pairwise_distance == 'DeltaR':
            return EMDTransverseMomentumDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return EMDTransverseMomentumHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return EMDTransverseMomentumHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    elif weight == 'TransverseEnergy':
        if pairwise_distance == 'DeltaR':
            return EMDTransverseEnergyDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return EMDTransverseEnergyHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return EMDTransverseEnergyHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    elif weight == 'Energy':
        if pairwise_distance == 'EEDot':
            return EMDEnergyEEDot(*args, **kwargs)
        elif pairwise_distance == 'EEDotMassive':
            return EMDEnergyEEDotMassive(*args, **kwargs)
        elif pairwise_distance == 'EEArcLength':
            return EMDEnergyEEArcLength(*args, **kwargs)
        elif pairwise_distance == 'EEArcLengthMassive':
            return EMDEnergyEEArcLengthMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    elif weight == 'Momentum':
        if pairwise_distance == 'EEDot':
            return EMDMomentumEEDot(*args, **kwargs)
        elif pairwise_distance == 'EEDotMassive':
            return EMDMomentumEEDotMassive(*args, **kwargs)
        elif pairwise_distance == 'EEArcLength':
            return EMDMomentumEEArcLength(*args, **kwargs)
        elif pairwise_distance == 'EEArcLengthMassive':
            return EMDMomentumEEArcLengthMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    else:
        raise TypeError('weight `{}` not allowed'.format(weight))

def PairwiseEMD(*args, weight='TransverseMomentum', pairwise_distance='DeltaR', **kwargs):

    if weight == 'TransverseMomentum':
        if pairwise_distance == 'DeltaR':
            return PairwiseEMDTransverseMomentumDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return PairwiseEMDTransverseMomentumHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return PairwiseEMDTransverseMomentumHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    elif weight == 'TransverseEnergy':
        if pairwise_distance == 'DeltaR':
            return PairwiseEMDTransverseEnergyDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return PairwiseEMDTransverseEnergyHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return PairwiseEMDTransverseEnergyHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    elif weight == 'Energy':
        if pairwise_distance == 'EEDot':
            return PairwiseEMDEnergyEEDot(*args, **kwargs)
        elif pairwise_distance == 'EEDotMassive':
            return PairwiseEMDEnergyEEDotMassive(*args, **kwargs)
        elif pairwise_distance == 'EEArcLength':
            return PairwiseEMDEnergyEEArcLength(*args, **kwargs)
        elif pairwise_distance == 'EEArcLengthMassive':
            return PairwiseEMDEnergyEEArcLengthMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    elif weight == 'Momentum':
        if pairwise_distance == 'EEDot':
            return PairwiseEMDMomentumEEDot(*args, **kwargs)
        elif pairwise_distance == 'EEDotMassive':
            return PairwiseEMDMomentumEEDotMassive(*args, **kwargs)
        elif pairwise_distance == 'EEArcLength':
            return PairwiseEMDMomentumEEArcLength(*args, **kwargs)
        elif pairwise_distance == 'EEArcLengthMassive':
            return PairwiseEMDMomentumEEArcLengthMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    else:
        raise TypeError('weight `{}` not allowed'.format(weight))
%}