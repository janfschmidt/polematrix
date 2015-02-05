#ifndef __POLEMATRIX__TRACKING_HPP_
#define __POLEMATRIX__TRACKING_HPP_

#include <vector>
#include "TrackingThread.hpp"

class Tracking
{
private:
  std::vector<TrackingThread> threads;

public:
  const unsigned int nThreads;

  Tracking(unsigned int nThreads_in);
  ~Tracking() {} //wird benötigt??

  void run();
};


#endif
// __POLEMATRIX__TRACKING_HPP_
