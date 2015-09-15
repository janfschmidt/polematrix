#ifndef __POLEMATRIX__TRACKING_HPP_
#define __POLEMATRIX__TRACKING_HPP_

// Class Tracking manages the tracking of multiple spins via a thread pool.
// It contains a queue of TrackingTasks (each representing a single spin tracking)
// and configuration.

#include <vector>
#include <thread>
#include <mutex>
#include <list>
#include "Configuration.hpp"
#include "TrackingTask.hpp"

class Tracking
{
private:
  std::vector<TrackingTask> queue;
  std::vector<std::thread> threadPool;
  std::vector<TrackingTask>::iterator queueIt;
  std::list<std::vector<TrackingTask>::const_iterator> runningTasks; // to display progress
  std::mutex mutex;
  void processQueue();
  void printProgress() const;

  //pal::SimToolInstance *sim;
  const pal::AccLattice *lattice;
  const pal::FunctionOfPos<pal::AccPair> *orbit;

  SpinMotion polarization;
  void calcPolarization();  //calculate polarization: average over all spin vectors for each time step


public:
  Configuration config;
  bool showProgressBar;

  Tracking(unsigned int nThreads=std::thread::hardware_concurrency());  // queue nP. tasks & create nT. threads
  ~Tracking();

  void start();                  // start tracking (processing queued tasks)

  unsigned int numParticles() const {return config.nParticles();}
  unsigned int numThreads() const {return threadPool.size();} // number of threads (particle trackings) executed in parallel
  //pal::SimToolInstance* getSimToolInstance() const {return sim;}
  const pal::AccLattice* getLattice() const {return lattice;}
  const pal::FunctionOfPos<pal::AccPair>* getOrbit() const {return orbit;}

  void setModel() {setLattice(); setOrbit();}
  void setLattice();
  void setOrbit();

  std::map<double,arma::colvec3> getPolarization() const {return polarization;}
  void savePolarization();
};



#endif
// __POLEMATRIX__TRACKING_HPP_
