#include <iostream>
#include <cmath>
#include "TrackingTask.hpp"



void TrackingTask::run()
{
  for (unsigned int i=0; i<this->particleId*10000000; i++) {
    x += std::cos(i);
  }
  std::cout << "Thread " << this->particleId << " result: " << x << std::endl;
}
