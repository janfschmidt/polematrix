#include <iostream>
//#include <thread>
#include "Tracking.hpp"

static const int num_threads = 10;

//This function will be called from a thread

void call_from_thread(int tid) {
  std::cout << "Launched by thread " << tid << std::endl;
}

int main()
{
  Tracking t(6);
  std::cout << t.numThreads() << " threads, " << t.numParticles() << " tasks." << std::endl;

  t.config.setPath("./test");
  t.config.x = 1; //rotation per step in degree
  t.config.s_start = {0,1,0};
  t.config.dt_out = 1e5;

  t.start();

  return 0;
}
