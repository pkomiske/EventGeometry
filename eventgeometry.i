%define EVENTGEOMETRY_DOCSTRING
"# EventGeometry FastJet Contrib
"
%enddef

%module("docstring"=EVENTGEOMETRY_DOCSTRING, "threads"=1) eventgeometry
%nothreadallow;

#define EMDNAMESPACE fastjet::contrib::emd

// C++ standard library wrappers
%include <std_list.i>

// this makes SWIG aware of the types contained in the main fastjet library
// but does not generate new wrappers for them here
%import FASTJET_PREFIX/share/fastjet/pyinterface/fastjet.i

// converts fastjet::Error into a FastJetError Python exception
FASTJET_ERRORS_AS_PYTHON_EXCEPTIONS(eventgeometry)

// turn off exception handling for now, since fastjet::Error is not thrown here
%exception;

// include headers in source file
%{
#ifndef SWIG
#define SWIG
#endif

#include "EventGeometry.hh"
using namespace fastjet::contrib::emd;
%}

// include EMD wrappers from Wasserstein package
#define BEGIN_EMD_NAMESPACE namespace fastjet { namespace contrib { namespace emd {
#define END_EMD_NAMESPACE } } }
%include "wasserstein/swig/wasserstein_common.i"

%define ADD_EXPLICIT_PREPROCESSING(Class)
Class & preprocess_CenterEScheme() { return $self->preprocess<CenterEScheme>(); }
Class & preprocess_CenterWeightedCentroid() { return $self->preprocess<CenterWeightedCentroid>(); }
Class & preprocess_CenterPtCentroid() { return $self->preprocess<CenterPtCentroid>(); }
Class & preprocess_MaskCircleRapPhi(double R) { return $self->preprocess<MaskCircleRapPhi>(R); }
%enddef

// extend [Pairwise]EMD to explicitly support PseudoJet arguments and printing
%extend fastjet::contrib::emd::EMD {
  ADD_EXPLICIT_PREPROCESSING(EMD)

  Value operator()(const fastjet::PseudoJet & pj0, const fastjet::PseudoJet & pj1) {
    return (*$self)(pj0, pj1);
  }
  Value operator()(const std::vector<fastjet::PseudoJet> & pjs0, const fastjet::PseudoJet & pj1) {
    return (*$self)(pjs0, pj1);
  }
  Value operator()(const fastjet::PseudoJet & pj0, const std::vector<fastjet::PseudoJet> & pjs1) {
    return (*$self)(pj0, pjs1);
  }
  Value operator()(const std::vector<fastjet::PseudoJet> & pjs0, const std::vector<fastjet::PseudoJet> & pjs1) {
    return (*$self)(pjs0, pjs1);
  }
}

%extend fastjet::contrib::emd::PairwiseEMD {
  ADD_EXPLICIT_PREPROCESSING(PairwiseEMD)

  void operator()(const std::vector<fastjet::PseudoJet> & evs) {
    (*$self)(evs);
  }
  void operator()(const std::vector<std::vector<fastjet::PseudoJet>> & evs) {
    (*$self)(evs);
  }
  void operator()(const std::vector<fastjet::PseudoJet> & evs0, const std::vector<fastjet::PseudoJet> & evs1) {
    (*$self)(evs0, evs1);
  }
  void operator()(const std::vector<std::vector<fastjet::PseudoJet>> & evs0, const std::vector<fastjet::PseudoJet> & evs1) {
    (*$self)(evs0, evs1);
  }
  void operator()(const std::vector<fastjet::PseudoJet> & evs0, const std::vector<std::vector<fastjet::PseudoJet>> & evs1) {
    (*$self)(evs0, evs1);
  }
  void operator()(const std::vector<std::vector<fastjet::PseudoJet>> & evs0, const std::vector<std::vector<fastjet::PseudoJet>> & evs1) {
    (*$self)(evs0, evs1);
  }
}

// instantiate hadronic and ee [Pairwise]EMD templates
namespace fastjet {
  namespace contrib {
    namespace emd {
      %template(EMDTransverseMomentumDeltaR) EMD<TransverseMomentum, DeltaR>;
      %template(EMDTransverseMomentumHadronicDot) EMD<TransverseMomentum, HadronicDot>;
      %template(EMDTransverseMomentumHadronicDotMassive) EMD<TransverseMomentum, HadronicDotMassive>;
      %template(EMDTransverseEnergyDeltaR) EMD<TransverseEnergy, DeltaR>;
      %template(EMDTransverseEnergyHadronicDot) EMD<TransverseEnergy, HadronicDot>;
      %template(EMDTransverseEnergyHadronicDotMassive) EMD<TransverseEnergy, HadronicDotMassive>;
      %template(EMDMomentumEEDot) EMD<Momentum, EEDot>;
      %template(EMDMomentumEEDotMassive) EMD<Momentum, EEDotMassive>;
      %template(EMDMomentumEEArcLength) EMD<Momentum, EEArcLength>;
      %template(EMDMomentumEEArcLengthMassive) EMD<Momentum, EEArcLengthMassive>;
      %template(EMDEnergyEEDot) EMD<Energy, EEDot>;
      %template(EMDEnergyEEDotMassive) EMD<Energy, EEDotMassive>;
      %template(EMDEnergyEEArcLength) EMD<Energy, EEArcLength>;
      %template(EMDEnergyEEArcLengthMassive) EMD<Energy, EEArcLengthMassive>;

      %template(PairwiseEMDTransverseMomentumDeltaR) PairwiseEMD<EMD<TransverseMomentum, DeltaR>>;
      %template(PairwiseEMDTransverseMomentumHadronicDot) PairwiseEMD<EMD<TransverseMomentum, HadronicDot>>;
      %template(PairwiseEMDTransverseMomentumHadronicDotMassive) PairwiseEMD<EMD<TransverseMomentum, HadronicDotMassive>>;
      %template(PairwiseEMDTransverseEnergyDeltaR) PairwiseEMD<EMD<TransverseEnergy, DeltaR>>;
      %template(PairwiseEMDTransverseEnergyHadronicDot) PairwiseEMD<EMD<TransverseEnergy, HadronicDot>>;
      %template(PairwiseEMDTransverseEnergyHadronicDotMassive) PairwiseEMD<EMD<TransverseEnergy, HadronicDotMassive>>;
      %template(PairwiseEMDMomentumEEDot) PairwiseEMD<EMD<Momentum, EEDot>>;
      %template(PairwiseEMDMomentumEEDotMassive) PairwiseEMD<EMD<Momentum, EEDotMassive>>;
      %template(PairwiseEMDMomentumEEArcLength) PairwiseEMD<EMD<Momentum, EEArcLength>>;
      %template(PairwiseEMDMomentumEEArcLengthMassive) PairwiseEMD<EMD<Momentum, EEArcLengthMassive>>;
      %template(PairwiseEMDEnergyEEDot) PairwiseEMD<EMD<Energy, EEDot>>;
      %template(PairwiseEMDEnergyEEDotMassive) PairwiseEMD<EMD<Energy, EEDotMassive>>;
      %template(PairwiseEMDEnergyEEArcLength) PairwiseEMD<EMD<Energy, EEArcLength>>;
      %template(PairwiseEMDEnergyEEArcLengthMassive) PairwiseEMD<EMD<Energy, EEArcLengthMassive>>;
    }
  }
}

// add convenience functions for accessing templated EMD, PairwiseEMD, and ApolloniusGroomer classes
%pythoncode %{

from fastjet import FastJetError

__version__ = '1.0.0a1'

def EMD(*args, weight='TransverseMomentum', pairwise_distance='DeltaR', **kwargs):

    if weight == 'TransverseMomentum':
        if pairwise_distance == 'DeltaR':
            return EMDTransverseMomentumDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return EMDTransverseMomentumHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return EMDTransverseMomentumHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

    elif weight == 'TransverseEnergy':
        if pairwise_distance == 'DeltaR':
            return EMDTransverseEnergyDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return EMDTransverseEnergyHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return EMDTransverseEnergyHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

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
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

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
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

    else:
        raise TypeError('weight `{}` not recognized'.format(weight))

def PairwiseEMD(*args, weight='TransverseMomentum', pairwise_distance='DeltaR', **kwargs):

    if weight == 'TransverseMomentum':
        if pairwise_distance == 'DeltaR':
            return PairwiseEMDTransverseMomentumDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return PairwiseEMDTransverseMomentumHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return PairwiseEMDTransverseMomentumHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

    elif weight == 'TransverseEnergy':
        if pairwise_distance == 'DeltaR':
            return PairwiseEMDTransverseEnergyDeltaR(*args, **kwargs)
        elif pairwise_distance == 'HadronicDot':
            return PairwiseEMDTransverseEnergyHadronicDot(*args, **kwargs)
        elif pairwise_distance == 'HadronicDotMassive':
            return PairwiseEMDTransverseEnergyHadronicDotMassive(*args, **kwargs)
        else:
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

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
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

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
            raise TypeError('pairwise distance `{}` not recognized'.format(pairwise_distance))

    else:
        raise TypeError('weight `{}` not recognized'.format(weight))
%}