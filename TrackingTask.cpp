#include <iostream>
#include <cmath>
#include "TrackingTask.hpp"



void TrackingTask::run()
{
  matrixTracking();

  // double x = config->x;
  // for (unsigned int i=0; i<this->particleId*10000000; i++) {
  //   x += std::cos(i);
  // }
  // std::cout << "Thread " << this->particleId << " result: " << x << std::endl;
}


void TrackingTask::matrixTracking()
{
  // arma::mat33 A = {,2,3,4,5,6,7,8,9};
  // arma::mat33 one(arma::fill::eye);
  arma::colvec3 s = config->s_start;
  double x=config->x*M_PI/180.;

  for (unsigned int i=0; i<this->particleId*123456789; i++) {
    s = rotxMatrix(x) * s;
  }

  s.print("s = ");
}

arma::mat33 TrackingTask::rotxMatrix(double angle) const {
  double c=std::cos(angle);
  double s=std::sin(angle);

  arma::mat33 rotx = {1,0,0, 0,c,s, 0,-s,c};
  return rotx;
}
