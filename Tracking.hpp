#ifndef __POLEMATRIX__TRACKING_HPP_
#define __POLEMATRIX__TRACKING_HPP_

// Class Tracking manages the tracking of multiple spins via a thread pool.
// It contains a queue of TrackingTasks (each representing a single spin tracking)
// and configuration.

#include <vector>
#include <thread>
#include <mutex>
#include "TrackingTask.hpp"

class Tracking
{
private:
  std::vector<TrackingTask> queue;
  std::vector<std::thread> threadPool;
  std::vector<TrackingTask>::iterator queueIt;
  std::mutex mutex;
  void processQueue();

  unsigned int nParticles; // total number of particles (tasks)
  unsigned int nThreads;   // number of threads to be executed in parallel

public:
 

  Tracking(unsigned int nParticles, unsigned int nThreads=std::thread::hardware_concurrency());  // queue nP. tasks & create nT. threads
  ~Tracking() {}

  void start();                  // start tracking (processing queued tasks)

  unsigned int numParticles() const {return nParticles;}
  unsigned int numThreads() const {return nThreads;}
};

#endif
// __POLEMATRIX__TRACKING_HPP_
