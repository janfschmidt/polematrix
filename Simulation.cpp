#include "Simulation.hpp"

void Simulation::setModel()
{
  config->updateSimToolSettings(*lattice);

  auto& palattice = config->getSimToolInstance();
  lattice.reset( new pal::AccLattice(palattice) );
  orbit.reset( new pal::FunctionOfPos<pal::AccPair>(palattice) );
  orbit->simToolClosedOrbit( palattice );

  if (config->gammaMode() == GammaMode::simtool
      || config->gammaMode() == GammaMode::simtool_plus_linear
      || config->gammaMode() == GammaMode::simtool_no_interpolation
      || config->gammaMode() == GammaMode::linear) {
    // no model setup needed
  }
  else {
    config->autocomplete(*lattice);
  }
}
