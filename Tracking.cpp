#include <iostream>
#include "Tracking.hpp"


Tracking::Tracking(unsigned int nParticles_in, unsigned int nThreads_in) : nParticles(nParticles_in), nThreads(nThreads_in), sim(NULL), lattice(NULL), orbit(NULL)
{
  // use at least one thread
  if (nThreads == 0)
    nThreads = 1;

  // fill queue
  for (unsigned int i=0; i<nParticles; i++) {
    TrackingTask toll(i,&config);
    queue.push_back(std::move(toll));
  }

  // set iterator to begin of queue
  queueIt = queue.begin();

  // create threads
  for (unsigned int i=0; i<nThreads; i++) {
    threadPool.emplace_back(std::thread());
  }
}


void Tracking::start()
{
  if (lattice==NULL || orbit==NULL)
    throw TrackError("ERROR: Tracking::start(): Cannot start tracking, if model is not specified (Lattice, Orbit).");

  auto start = std::chrono::high_resolution_clock::now();

  //start threads
  for (std::thread& t : threadPool) {
    t = std::thread(&Tracking::processQueue,this);
  }

  // wait for threads to finish
  for (std::thread& t : threadPool) {
    t.join();
  }

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout << "Tracking took " << duration.count() << " s." << std::endl;
}


void Tracking::processQueue()
{
  std::vector<TrackingTask>::iterator myTask;

  while (true) {
    mutex.lock();
    if (queueIt == queue.end()) {
      mutex.unlock();
      return;
    }
    else {
      myTask = queueIt;
      queueIt++;
      mutex.unlock();
      myTask->lattice=lattice;
      myTask->orbit=orbit;
      myTask->run(); // run next queued TrackingTask
    }
  }
}


void Tracking::setModel(pal::SimToolInstance *s, pal::AccLattice *l, pal::FunctionOfPos<pal::AccPair> *o)
{
  setSimToolInstance(s);
  setLattice(l);
  setOrbit(o);
}
