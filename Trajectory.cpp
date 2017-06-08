/* Trajectory Classes
 * It represents the trajectory of a single particle, which can be
 * an Elegant/MadX trajectory or just a closed orbit
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

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "Trajectory.hpp"
#include "debug.hpp"


void Trajectory::init()
{
  if (!initDone) {
    initImplementation();
    initDone = true;
    polematrix::debug(__PRETTY_FUNCTION__, "init done");
  }
}


SimtoolTrajectory::SimtoolTrajectory(unsigned int id, const std::shared_ptr<Configuration> c)
  : Trajectory(id,c), simtoolTrajectory(config->getSimToolInstance(), gsl_interp_akima) {}


void SimtoolTrajectory::initImplementation()
{
  // simtool: sdds import thread safe since SDDSToolKit-devel-3.3.1-2
  simtoolTrajectory.simToolTrajectory( config->getSimToolInstance(), particleId+1 );
}


void SimtoolTrajectory::saveSimtoolData()
{
  if ( !config->saveGamma(particleId) )
    return;
  
  simtoolTrajectory.info.add("polematrix particle ID", particleId);
  std::stringstream file;
  file <<std::setw(4)<<std::setfill('0')<< particleId << ".dat";
  simtoolTrajectory.print( (config->outpath()/"trajectorySimtool_").string() + file.str() );
}


Oscillation::Oscillation(unsigned int id, const std::shared_ptr<Configuration> c)
  : Trajectory(id,c), beta(config->getSimToolInstance()) {}

pal::AccPair Oscillation::get(const double& pos)
{
  double s = orbit->posInTurn(pos);
  // oscillation amplitude from (single particle) emittance & beta function
  // oscillation frequency from tune (see initImplementation())
  pal::AccPair b = beta.interp(s);
  pal::AccPair phase = freq * pos + phase0;
  pal::AccPair traj;
  traj.x = std::sqrt(emittance.x * b.x) * std::cos(phase.x);
  traj.z = std::sqrt(emittance.z * b.z) * std::cos(phase.z);
  return orbit->interp(s) + traj;
}


void Oscillation::initImplementation()
{
  // init twiss functions beta & phase
  std::string betaX, betaZ;
  if (config->getSimToolInstance().tool == pal::madx) {
    betaX = "BETX";
    betaZ = "BETY";
  }
  else if (config->getSimToolInstance().tool == pal::elegant) {
    betaX = "betax";
    betaZ = "betay";
  }
  beta.readTwissColumn(config->getSimToolInstance(), betaX, betaZ);
  if (particleId == 0) {
  std::cout << "* " << beta.size() << " beta function sampling points read" << std::endl
	    << "  from " << config->getSimToolInstance().twiss() << std::endl;
  }

  // init single particle emittance (gaussian distribution)
  //      & start phase (uniform distribution)
  boost::random::normal_distribution<> exDistr(0.0, config->emittance().x);
  boost::random::normal_distribution<> ezDistr(0.0, config->emittance().z);
  boost::random::uniform_real_distribution<> phase0Distr(0.0, 2*M_PI);

  boost::random::mt11213b rng(config->seed()+particleId);
  emittance.x = std::fabs( exDistr(rng) );
  emittance.z = std::fabs( ezDistr(rng) );
  phase0.x = phase0Distr(rng);
  phase0.z = phase0Distr(rng);

  // init oscillation frequency from tune
  freq = config->tune() * 2*M_PI / orbit->circumference();
}
