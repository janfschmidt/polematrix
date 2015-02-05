#include <iostream>
#include <cmath>
#include "TrackingThread.hpp"

TrackingThread::TrackingThread(TrackingThread &&other)
  : particleId(other.particleId), x(other.x)
{
  myThread = std::move(other.myThread);
}


void TrackingThread::threadMain()
{
  for (unsigned int i=0; i<this->particleId*10000000; i++) {
    x += std::cos(i);
  }
  std::cout << "Thread " << this->particleId << " result: " << x << std::endl;
}
