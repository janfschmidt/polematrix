#ifndef __POLEMATRIX__TRACKINGTASK_HPP_
#define __POLEMATRIX__TRACKINGTASK_HPP_


class TrackingTask
{
public:
  const unsigned int particleId;
  double x;

  TrackingTask(unsigned int id, double x_in) : particleId(id), x(x_in) {}

  void run();     //run tracking task
};

#endif
// __POLEMATRIX__TRACKINGTASK_HPP_
