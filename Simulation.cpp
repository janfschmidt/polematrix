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


SingleParticleSimulation::SingleParticleSimulation(unsigned int id, const std::shared_ptr<Configuration> c)
  : particleId(id), config(c)
{
  switch (config->trajectoryMode()) {
  case TrajectoryMode::simtool:
    trajectory.reset( new SimtoolTrajectory(particleId, config) );
    break;
  default:
    trajectory.reset( new Orbit(particleId, config) );
  }
}


void SingleParticleSimulation::setModel(std::shared_ptr<const pal::AccLattice> l, std::shared_ptr<const pal::FunctionOfPos<pal::AccPair>> o)
{
  lattice = l;
  orbit = o;
  trajectory->setOrbit(orbit);
}
