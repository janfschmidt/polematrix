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

#ifndef __POLEMATRIX__SIMULATION_HPP_
#define __POLEMATRIX__SIMULATION_HPP_

#include <libpalattice/AccLattice.hpp>
#include <libpalattice/FunctionOfPos.hpp>
#include <thread>
#include <mutex>
#include <memory>
#include <list>
#include <vector>
#include "Configuration.hpp"
#include "Trajectory.hpp"


// abstract base class for a simulation task for a single particle
// can be used by a Simulation object, which passes its config/lattice/orbit
// to many SingleParticleSimulation

class SingleParticleSimulation {
protected:
  std::unique_ptr<Trajectory> trajectory;     // particle trajectory, implementation depends TrajectoryMode
  
public:
  const unsigned int particleId;
  const std::shared_ptr<Configuration> config; //not const Object, because SimToolInstance status can be changed
  std::shared_ptr<const pal::AccLattice> lattice;
  std::shared_ptr<const pal::FunctionOfPos<pal::AccPair>> orbit;

  SingleParticleSimulation(unsigned int id, const std::shared_ptr<Configuration> c);
  void setModel(std::shared_ptr<const pal::AccLattice> l, std::shared_ptr<const pal::FunctionOfPos<pal::AccPair>> o);
  virtual void run() =0;
};






// abstract base class for a simulation
// including a config as well as lattice and orbit, which can be set using config

template <typename T=SingleParticleSimulation>
class Simulation {
protected:
  std::shared_ptr<pal::AccLattice> lattice;
  std::shared_ptr<pal::FunctionOfPos<pal::AccPair>> orbit;

  // queue
  typedef typename std::vector<T>::iterator taskIterator;
  typedef typename std::vector<T>::const_iterator const_taskIterator;
  std::vector<T> queue;
  taskIterator queueIt;
  std::list<const_taskIterator> runningTasks; // to display progress


  // thread management
  std::vector<std::thread> threadPool;
  std::mutex mutex;
  void initThreadPool(unsigned int nThreads);
  void startThreads();
  void waitForThreads();
  void processQueue();
  
  std::map<unsigned int,std::string> errors;

public:
  const std::shared_ptr<Configuration> config;
  
  Simulation(unsigned int nThreads=std::thread::hardware_concurrency()) : queueIt(queue.begin()), config(new Configuration) {initThreadPool(nThreads);}
  Simulation(const std::shared_ptr<Configuration> c, unsigned int nThreads=std::thread::hardware_concurrency()) : queueIt(queue.begin()), config(c) {initThreadPool(nThreads);}
  Simulation(const Simulation& o) = delete;
  
  void setModel();
  
  bool modelReady() {if (lattice->size()==0 || orbit->size()==0) return false; else return true;}
  unsigned int numParticles() const {return config->nParticles();}
  unsigned int numSuccessful() const {return numParticles() - errors.size();}
    
  virtual void start() =0;
  
  void saveLattice() const {lattice->print( (config->outpath()/"lattice.dat").string() );}
  void saveOrbit() const {orbit->print( (config->outpath()/"closedorbit.dat").string() );}
};







// ---- template class implementation ----
template <typename T>
void Simulation<T>::initThreadPool(unsigned int nThreads)
{
  // use at least one thread
  if (nThreads == 0)
    nThreads = 1;
  
  // create threads
  for (auto i=0u; i<nThreads; i++) {
    threadPool.emplace_back(std::thread());
  }
}

template <typename T>
void Simulation<T>::startThreads()
{
  for (std::thread& t : threadPool) {
    t = std::thread(&Simulation::processQueue,this);
  }
}

template <typename T>
void Simulation<T>::waitForThreads()
{
  for (std::thread& t : threadPool) {
    t.join();
  }
}

template <typename T>
void Simulation<T>::setModel()
{
  auto& palattice = config->getSimToolInstance();
  lattice.reset( new pal::AccLattice(palattice) );
  orbit.reset( new pal::FunctionOfPos<pal::AccPair>(palattice) );
  config->updateSimToolSettings(*lattice);
  orbit->simToolClosedOrbit( palattice );
  config->writeRfMagnetsToLattice(*lattice);

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

template <typename T>
void Simulation<T>::processQueue()
{
  while (true) {
    mutex.lock();
    if (queueIt == queue.end()) {
      mutex.unlock();
      return;   // finish this thread
    }
    else {
      taskIterator myTask = queueIt;
      queueIt++;
      runningTasks.push_back(myTask); // to display progress
      mutex.unlock();
      try {
	myTask->setModel(lattice, orbit);
	myTask->run(); // run next queued task
      }
      //cancel thread in error case
      catch (std::exception &e) {
	std::cout << "ERROR @ particle " << myTask->particleId
	  // << " (thread_id "<< std::this_thread::get_id() << ")"
		  <<":"<< std::endl
		  << e.what() << std::endl;
	mutex.lock();
	errors.emplace(myTask->particleId, e.what());
	mutex.unlock();
      }
      mutex.lock();
      runningTasks.remove(myTask); // to display progress
      mutex.unlock();
    }//else
  }//while
}




#endif
// __POLEMATRIX__SIMULATION_HPP_
