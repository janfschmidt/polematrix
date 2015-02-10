#include <iostream>
#include "Tracking.hpp"


Tracking::Tracking(unsigned int nParticles_in, unsigned int nThreads_in) : nParticles(nParticles_in), nThreads(nThreads_in)
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
      myTask->run(); // run next queued TrackingTask
    }
  }
}
