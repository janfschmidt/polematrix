#ifndef __POLEMATRIX__TRACKING_HPP_
#define __POLEMATRIX__TRACKING_HPP_

// Class Tracking manages the tracking of multiple spins via a thread pool.
// It contains a queue of TrackingTasks (each representing a single spin tracking)
// and configuration.

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
