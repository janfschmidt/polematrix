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

#ifndef __POLEMATRIX__TRACKING_HPP_
#define __POLEMATRIX__TRACKING_HPP_

#include <vector>
#include <thread>
#include <mutex>
#include <list>
#include <memory>
#include "Simulation.hpp"
#include "TrackingTask.hpp"



class Tracking : public Simulation
{
private:
  std::vector<TrackingTask> queue;
  std::vector<std::thread> threadPool;
  std::vector<TrackingTask>::iterator queueIt;
  std::list<std::vector<TrackingTask>::const_iterator> runningTasks; // to display progress
  std::mutex mutex;
  
  void processQueue();
  void printProgress() const;

  SpinMotion polarization;
  void calcPolarization();  //calculate polarization: average over all spin vectors for each time step


public:
  bool showProgressBar;

  Tracking(unsigned int nThreads=std::thread::hardware_concurrency());  // queue nP. tasks & create nT. threads
  Tracking(const Tracking& o) = delete;
  ~Tracking() {}
  
  void start();                  // start tracking (processing queued tasks)

  unsigned int numParticles() const {return config->nParticles();}
  unsigned int numThreads() const {return threadPool.size();} // number of threads (particle trackings) executed in parallel

  std::map<double,arma::colvec3> getPolarization() const {return polarization;}
  void savePolarization();
};



#endif
// __POLEMATRIX__TRACKING_HPP_
