// -*- C++ -*-
%define EVENTGEOMETRY_DOCSTRING
"# EventGeometry FastJet Contrib
"
%enddef

%module("docstring"=EVENTGEOMETRY_DOCSTRING, "threads"=1) eventgeometry
%nothreadallow;

#define EVENTGEOMETRY_NAMESPACE fastjet::contrib::eventgeometry

// C++ standard library wrappers
%include <std_vector.i>

// ensure FASTJET_PREFIX is always defined even if not used
#ifndef FASTJET_PREFIX
# define FASTJET_PREFIX /usr/local
#endif

// import either pyfjcore or fastjet
#ifdef EVENTGEOMETRY_USE_PYFJCORE
  %import PyFJCore/pyfjcore/swig/pyfjcore.i
  %pythoncode %{
    from pyfjcore import FastJetError
  %}
#else
  %import FASTJET_PREFIX/share/fastjet/pyinterface/fastjet.i
  %pythoncode %{
    from fastjet import FastJetError
  %}
#endif

// converts fastjet::Error into a FastJetError Python exception
FASTJET_ERRORS_AS_PYTHON_EXCEPTIONS(eventgeometry)

// turn off exception handling for now, since fastjet::Error is not thrown here
%exception;

// include headers in source file
%{
#ifndef SWIG
# define SWIG
#endif

#include "EventGeometry.hh"
%}

// set Wasserstein namespace appropriately for this package
#define WASSERSTEIN_NAMESPACE EVENTGEOMETRY_NAMESPACE
#define BEGIN_WASSERSTEIN_NAMESPACE namespace fastjet { namespace contrib { namespace eventgeometry {
#define END_WASSERSTEIN_NAMESPACE } } }

// include EMD wrappers from Wasserstein package, skipping single precision classes
#define WASSERSTEIN_NO_FLOAT32
%include "wasserstein/swig/wasserstein_common.i"

// this needs to come after numpy has been included by wasserstein_common.i
%{
#include "pyfjcore/PyFJCoreExtensions.hh"
%}

// this needs to be early so that the symbols are loaded from the dynamic library
%pythonbegin %{
  import pyfjcore
%}

%define EVENTGEOMETRY_ADD_EXPLICIT_PREPROCESSORS
void preprocess_CenterEScheme() { $self->preprocess<EVENTGEOMETRY_NAMESPACE::CenterEScheme>(); }
void preprocess_CenterPtCentroid() { $self->preprocess<EVENTGEOMETRY_NAMESPACE::CenterPtCentroid>(); }
void preprocess_MaskCircle(double R) { $self->preprocess<EVENTGEOMETRY_NAMESPACE::MaskCircle>(R); }
%enddef

%ignore EVENTGEOMETRY_NAMESPACE::FastJetEventBase;
%ignore EVENTGEOMETRY_NAMESPACE::FastJetParticleWeight;

%template(vectorPseudoJetContainer) std::vector<fastjet::PseudoJetContainer>;
%template(vectorVectorPseudoJet) std::vector<std::vector<fastjet::PseudoJet>>;

%include "EventGeometry.hh"

namespace WASSERSTEIN_NAMESPACE {

  // extend [Pairwise]EMD to explicitly support PseudoJet arguments and printing
  %extend EMD {
    EVENTGEOMETRY_ADD_EXPLICIT_PREPROCESSORS

    double operator()(const PseudoJetContainer & pjc0, const PseudoJetContainer & pjc1) {
      return $self->operator()(pjc0.as_vector(), pjc1.as_vector());
    }
    double operator()(const std::vector<PseudoJet> & pjs0, const std::vector<PseudoJet> & pjs1) {
      return $self->operator()(pjs0, pjs1);
    }
    double operator()(const PseudoJetContainer & pjc0, const std::vector<PseudoJet> & pjs1) {
      return $self->operator()(pjc0.as_vector(), pjs1);
    }
    double operator()(const std::vector<PseudoJet> & pjs0, const PseudoJetContainer & pjc1) {
      return $self->operator()(pjs0, pjc1.as_vector());
    }
    double operator()(const PseudoJet & pj0, const PseudoJet & pj1) {
      return $self->operator()(pj0, pj1);
    }
  }

  %extend PairwiseEMD {
    EVENTGEOMETRY_ADD_EXPLICIT_PREPROCESSORS

    // add a single event to the PairwiseEMD object
    void _add_event(const PseudoJetContainer & pjc, double event_weight = 1) {
      $self->events().emplace_back(pjc.as_vector(), event_weight);
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
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Momentum, EEDotMassless)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Momentum, EEArcLength)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Momentum, EEArcLengthMassive)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Energy, EEDot)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Energy, EEDotMassless)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Energy, EEArcLength)
  EVENTGEOMETRY_EMD_SWIG_TEMPLATE(Energy, EEArcLengthMassive)
}

// needed by wasserstein_common.i to add events to PairwiseEMDBase
%pythoncode %{
def _store_events(pairwise_emd, events, event_weights, gdim, mask):

    if gdim or mask:
        raise ValueError('`gdim` and `mask` are not supported')

    for event, event_weight in zip(events, event_weights):
        pairwise_emd._add_event(event, event_weight)
%}

// macro to avoid duplicating code for EMD and PairwiseEMD
%define EVENTGEOMETRY_PYTHON_EMD_FACTORY(Class)
%pythoncode %{
def Class(*args, weight='TransverseMomentum', pairwise_distance='DeltaR', **kwargs):
    if weight == 'TransverseMomentum':
        if pairwise_distance == 'DeltaR':
            return Class##TransverseMomentumDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return Class##TransverseMomentumHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return Class##TransverseMomentumHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    elif weight == 'TransverseEnergy':
        if pairwise_distance == 'DeltaR':
            return Class##TransverseEnergyDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return Class##TransverseEnergyHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return Class##TransverseEnergyHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    elif weight == 'Energy':
        if pairwise_distance == 'EEDot':
            return Class##EnergyEEDot(*args, **kwargs)
        elif pairwise_distance == 'EEDotMassless':
            return Class##EnergyEEDotMassless(*args, **kwargs)
        elif pairwise_distance == 'EEArcLength':
            return Class##EnergyEEArcLength(*args, **kwargs)
        elif pairwise_distance == 'EEArcLengthMassive':
            return Class##EnergyEEArcLengthMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    elif weight == 'Momentum':
        if pairwise_distance == 'EEDot':
            return Class##MomentumEEDot(*args, **kwargs)
        elif pairwise_distance == 'EEDotMassless':
            return Class##MomentumEEDotMassless(*args, **kwargs)
        elif pairwise_distance == 'EEArcLength':
            return Class##MomentumEEArcLength(*args, **kwargs)
        elif pairwise_distance == 'EEArcLengthMassive':
            return Class##MomentumEEArcLengthMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not allowed'.format(pairwise_distance))

    else:
        raise TypeError('weight `{}` not allowed'.format(weight))
%}
%enddef

EVENTGEOMETRY_PYTHON_EMD_FACTORY(EMD)
EVENTGEOMETRY_PYTHON_EMD_FACTORY(PairwiseEMD)
