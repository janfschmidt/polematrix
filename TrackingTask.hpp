#ifndef __POLEMATRIX__TRACKINGTASK_HPP_
#define __POLEMATRIX__TRACKINGTASK_HPP_

#include <armadillo>
#include "Configuration.hpp"

class TrackingTask
{
public:
  const unsigned int particleId;
  const Configuration *config;

  TrackingTask(unsigned int id, const Configuration *c) : particleId(id), config(c) {}

  void run();     //run tracking task
  void matrixTracking();
  inline arma::mat33 rotxMatrix(double angle) const; //inline speeds up example by ~ factor 4!
};

#endif
// __POLEMATRIX__TRACKINGTASK_HPP_
