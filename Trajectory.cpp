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

#include "Trajectory.hpp"
#include "debug.hpp"

SimtoolTrajectory::SimtoolTrajectory(unsigned int id, const std::shared_ptr<Configuration> c)
  : Trajectory(id,c), simtoolTrajectory(config->getSimToolInstance(), gsl_interp_akima), initDone(false) {}


void SimtoolTrajectory::init()
{
  if (!initDone) {
    // simtool: sdds import thread safe since SDDSToolKit-devel-3.3.1-2
    simtoolTrajectory.simToolTrajectory( config->getSimToolInstance(), particleId+1 );
    initDone = true;
    polematrix::debug(__PRETTY_FUNCTION__, "init done");
  }
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
