#ifndef __POLEMATRIX__TRACKING_HPP_
#define __POLEMATRIX__TRACKING_HPP_

// Class Tracking manages the tracking of multiple spins via a thread pool.
// It contains a queue of TrackingTasks (each representing a single spin tracking)
// and configuration.

#include <vector>
#include <thread>
#include <mutex>
#include "Configuration.hpp"
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

  pal::SimToolInstance *sim;
  const pal::AccLattice *lattice;
  const pal::FunctionOfPos<pal::AccPair> *orbit;


public:
  Configuration config;

  Tracking(unsigned int nParticles, unsigned int nThreads=std::thread::hardware_concurrency());  // queue nP. tasks & create nT. threads
  ~Tracking() {}

  void start();                  // start tracking (processing queued tasks)

  unsigned int numParticles() const {return nParticles;}
  unsigned int numThreads() const {return nThreads;}
  pal::SimToolInstance* getSimToolInstance() const {return sim;}
  const pal::AccLattice* getLattice() const {return lattice;}
  const pal::FunctionOfPos<pal::AccPair>* getOrbit() const {return orbit;}

  void setModel(pal::SimToolInstance *s, pal::AccLattice *l, pal::FunctionOfPos<pal::AccPair> *o);
  void setSimToolInstance(pal::SimToolInstance *s) {sim = s;}
  void setLattice(pal::AccLattice *l) {lattice = l;}
  void setOrbit(pal::FunctionOfPos<pal::AccPair> *o) {orbit = o;}
};

#endif
// __POLEMATRIX__TRACKING_HPP_
