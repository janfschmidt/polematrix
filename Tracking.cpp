#include "Tracking.hpp"


Tracking::Tracking(unsigned int nThreads_in) : nThreads(nThreads_in)
{
  double x_in = 1;
  for (unsigned int i=0; i<nThreads; i++) {
    // TrackingThread t(i,x_in);
    // threads.emplace_back(std::move(t));
    threads.emplace_back(TrackingThread(i,x_in));
  }
}

// Tracking::~Tracking()
// {

// }


void Tracking::run()
{
  for (unsigned int i=0; i<nThreads; i++) {
    threads[i].start();
  }
}

