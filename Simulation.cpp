/* Simulation base classes
 * used by Tracking/TrackingTask and ResStrengths/ParticleResStrengths
 *
 * Copyright (C) 2017 Jan Felix Schmidt <janschmidt@mailbox.org>
 *   
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Simulation.hpp"



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


std::string SingleParticleSimulation::getProgressBar(unsigned int barWidth) const
{
  std::stringstream bar;

  bar << particleId << ":";

  if (barWidth!=0) {
    unsigned int steps = (1.0 * barWidth * getProgress()) + 0.5;
    unsigned int i=0;
    bar << "[";
    for (; i<steps; i++)
      bar << "=";
    for (; i<barWidth; i++)
      bar << " ";
    bar << "]";
  }
  bar << std::fixed<<std::setprecision(0)<<std::setw(2)<<std::setfill('0')<< getProgress()*100. << "%";
  
  return bar.str();
}
