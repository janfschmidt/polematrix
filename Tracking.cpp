/* Tracking Class
 * manages the tracking of multiple spins via a thread pool.
 * It contains a queue of TrackingTasks (each representing a single spin tracking)
 * and configuration.
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

#include <iostream>
#include "Tracking.hpp"
#include "version.hpp"



void Tracking::start()
{ 
  if (!modelReady())
    throw TrackError("Cannot start tracking, if model is not specified (Lattice, Orbit).");

  if (config->t_stop() <= config->t_start()) {
    std::stringstream msg;
    msg << "Cannot track backwards: t_stop=" << config->t_stop() << " < t_start=" << config->t_start();
    throw TrackError(msg.str());
  }

  // fill queue
  for (unsigned int i=0; i<config->nParticles(); i++) {
    queue.emplace_back( TrackingTask(i,config) );
  }
  // set iterator to begin of queue
  queueIt = queue.begin();

  // write current config to file
  config->save( config->confOutFile().string() );

  std::cout << "Start tracking "<<config->nParticles()<<" Spins..." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  //start threads (incl. progress bars)
  startThreads();

  waitForThreads();

  // finished: calc. time & error output
  auto stop = std::chrono::high_resolution_clock::now();
  auto secs = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout << std::endl
	    << "-----------------------------------------------------------------" << std::endl;
  if (errors.size()>0) {
    std::cout << "ERRORS occured during tracking!" << std::endl;
  }
  std::cout << "Tracking "<<numSuccessful()<< " Spins done. Tracking took ";
  std::cout << secs.count() << " s = "<< int(secs.count()/60.+0.5) << " min." << std::endl;
  std::cout << "Thanks for using polematrix " << polemversion() << std::endl;
  std::cout << printErrors();
  std::cout << "-----------------------------------------------------------------" << std::endl;

  if (numSuccessful() > 0)
    calcPolarization();
}




//calculate polarization: average over all successfully tracked spin vectors for each time step
void Tracking::calcPolarization()
{
  unsigned int i=0;
  for (; i<queue.size(); i++) {
    if (errors.count(i)==0) {
      polarization = queue[i].getStorage();
      break;
    }
  }

  for (i++; i<queue.size(); i++) {
    if (errors.count(i)==0)
      polarization += queue[i].getStorage();
  }
  polarization /= numSuccessful();
}

void Tracking::savePolarization()
{
  std::ofstream file;
  std::string filename = config->polFile().string();
  unsigned int w = 14;
  
  file.open(filename);
  if (!file.is_open())
    throw TrackFileError(filename);

  file << config->metadata();
  file << "# Polarization calculated as average over " << numSuccessful() << " spins" << std::endl;
  file << polarization.printHeader(w, "P") << std::endl;
  file << polarization.print(w);
  
  file.close();
  std::cout << "* Polarization written for " << polarization.size() << " steps to " << filename <<"."<< std::endl;
}

