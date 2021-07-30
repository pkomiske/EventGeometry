//------------------------------------------------------------------------
// This file is part of EventGeometry, a C++ library with a Python wrapper
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

  // demonstrate calculating single EMDs
  EMD emd_obj(0.4, 1.0, true);

  // preprocess events to center
  emd_obj.preprocess<contrib::eventgeometry::CenterWeightedCentroid>();

  // print description
  std::cout << emd_obj.description() << std::endl;

  // container to hold emd values
  std::vector<double> emds;

  // loop over events and compute the EMD between each successive pair
  evp->reset();
  while (true) {

    // get first event
    if (!evp->next()) break;
    auto event0(convert2pjs(evp->particles()));

    // get second event
    if (!evp->next()) break;
    auto event1(convert2pjs(evp->particles()));

    // compute emd and add it to vector
    emds.push_back(emd_obj(event0, event1));
  }

  // get max and min EMD value
  std::cout << '\n'
            << emds.size() << " EMDs computed\n"
            << "Min. EMD - " << *std::min_element(emds.begin(), emds.end()) << '\n'
            << "Max. EMD - " << *std::max_element(emds.begin(), emds.end()) << '\n'
            << '\n';

  return 0;
}
