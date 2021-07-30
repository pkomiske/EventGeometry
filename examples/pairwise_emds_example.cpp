//------------------------------------------------------------------------
// This file is part of Wasserstein, a C++ library with a Python wrapper
// that computes the Wasserstein/EMD distance. If you use it for academic
// research, please cite or acknowledge the following works:
//
//   - Komiske, Metodiev, Thaler (2019) arXiv:1902.02346
//       https://doi.org/10.1103/PhysRevLett.123.041801
//   - Komiske, Metodiev, Thaler (2020) arXiv:2004.04159
//       https://doi.org/10.1007/JHEP07%282020%29006
//   - Boneel, van de Panne, Paris, Heidrich (2011)
//       https://doi.org/10.1145/2070781.2024192
//   - LEMON graph library https://lemon.cs.elte.hu/trac/lemon
//
// Copyright (C) 2019-2021 Patrick T. Komiske III
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//------------------------------------------------------------------------

// EventGeometry library
//  - this includes fastjet/PseudoJet.hh by default
//  - this includes the Wasserstein library (header only)
#include "EventGeometry.hh"

// classes and functions for reading/preparing events from EnergyFlow npz files
#include "Wasserstein/examples/include/ExampleUtils.hh"

// without this line, `contrib::eventgeometry` should be prefixed with `fastjet::`
using namespace fastjet;

// `EMDFloat64` uses `double` for the floating-point type
// first template parameter is an Event type, second is a PairwiseDistance type
using EMD = contrib::eventgeometry::EMDFloat64<contrib::eventgeometry::TransverseMomentum,
                                                contrib::eventgeometry::DeltaR>;

// `PairwiseEMD` is a templated class accepting a fully qualified `EMD` type
using PairwiseEMD = contrib::eventgeometry::PairwiseEMD<EMD>;

// empty angle brackets use the default floating-point type `double`
using CorrelationDimension = contrib::eventgeometry::CorrelationDimension<>;

// this function is useful for getting events as a vector of PseudoJets
// `Particle` here is an internal type used to read in events from the npz files
std::vector<fastjet::PseudoJet> convert2pjs(const std::vector<Particle> & particles) {
  std::vector<fastjet::PseudoJet> pjs;
  pjs.reserve(particles.size());

  for (const Particle & particle : particles)
    pjs.push_back(fastjet::PtYPhiM(particle.pt, particle.y, particle.phi));

  return pjs;
}

int main(int argc, char** argv) {

  // load events
  EventProducer * evp(load_events(argc, argv));
  if (evp == nullptr)
    return 1;

  // demonstrate calculating pairwise EMDs
  PairwiseEMD pairwise_emd_obj(0.4, 1.0, false);

  // preprocess events to center
  pairwise_emd_obj.preprocess<contrib::eventgeometry::CenterWeightedCentroid>();

  // print description
  std::cout << pairwise_emd_obj.description() << std::endl;

  // get vector of events
  std::vector<std::vector<fastjet::PseudoJet>> events;

  // loop over events and compute the EMD between each successive pair
  evp->reset();
  for (int i = 0; i < 1000 && evp->next(); i++)
    events.push_back(convert2pjs(evp->particles()));

  // run computation
  pairwise_emd_obj(events);

  // get max and min EMD value
  const std::vector<double> & emds(pairwise_emd_obj.emds(true));
  std::cout << "Min. EMD - " << *std::min_element(emds.begin(), emds.end()) << '\n'
            << "Max. EMD - " << *std::max_element(emds.begin(), emds.end()) << '\n'
            << '\n';

  // setup correlation dimension
  CorrelationDimension corrdim(50, 10., 250.);
  pairwise_emd_obj.set_external_emd_handler(corrdim);

  // rerun computation
  pairwise_emd_obj(events);

  // print out correlation dimensions
  auto corrdims(corrdim.corrdims());
  auto corrdim_bins(corrdim.corrdim_bins());
  std::cout << "\nEMD         Corr. Dim.  Error\n" << std::left;
  for (unsigned i = 0; i < corrdims.first.size(); i++)
    std::cout << std::setw(12) << corrdim_bins[i]
              << std::setw(12) << corrdims.first[i]
              << std::setw(12) << corrdims.second[i]
              << '\n';

  return 0;
}
