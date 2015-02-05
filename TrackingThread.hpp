#ifndef __POLEMATRIX__TRACKINGTHREAD_HPP_
#define __POLEMATRIX__TRACKINGTHREAD_HPP_

#include <thread>

class TrackingThread
{
private:
  void threadMain();   // executed in thread

public:
  std::thread myThread;
  const unsigned int particleId;
  double x;

  TrackingThread(unsigned int id, double x_in) : particleId(id), x(x_in) {}
  TrackingThread(TrackingThread &&other);
  ~TrackingThread() {if(myThread.joinable()) myThread.join();} //waits for threadMain() to finish

  void start() {myThread = std::thread(&TrackingThread::threadMain,this);}     //start execution of threadMain() in a thread
};

#endif
// __POLEMATRIX__TRACKINGTHREAD_HPP_
