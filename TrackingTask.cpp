#include <iostream>
#include <iomanip>
#include <cmath>
#include "TrackingTask.hpp"


TrackingTask::TrackingTask(unsigned int id, const Configuration *c)
  : particleId(id), config(c), w(14)
{
  outfile = std::unique_ptr<std::ofstream>(new std::ofstream());
}



void TrackingTask::run()
{
  outfileOpen();
  
  matrixTracking();
  
  outfileClose();
}


void TrackingTask::matrixTracking()
{
  // arma::mat33 A = {,2,3,4,5,6,7,8,9};
  // arma::mat33 one(arma::fill::eye);

  arma::colvec3 s = config->s_start;
  double x=config->x*M_PI/180.;

  for (unsigned int i=0; i<this->particleId*123456789; i++) {
    s = rotxMatrix(x) * s;

    if (i%config->dt_out == 0) //output
      storeStep(i,s);
  }
}

arma::mat33 TrackingTask::rotxMatrix(double angle) const
{
  double c=std::cos(angle);
  double s=std::sin(angle);

  arma::mat33 rotx = {1,0,0, 0,c,s, 0,-s,c};
  return rotx;
}



std::string TrackingTask::outfileName() const
{
  std::stringstream ss;
  ss << config->subfolder("spins") << "spin_" << std::setw(4)<<std::setfill('0')<<particleId << ".dat";
  return ss.str();
}



void TrackingTask::outfileOpen()
{
  outfile->open(outfileName());
  if (!outfile->is_open())
    throw TrackFileError(outfileName());
  else
    *outfile << "#"<<std::setw(w+1)<< "t / s" <<std::setw(w)<< "Sx" <<std::setw(w)<< "Sz" <<std::setw(w)<< "Ss"
	    <<std::setw(w)<< "|S|" <<std::setw(w)<< "gamma" << std::endl;
}


void TrackingTask::outfileClose()
{
  outfile->close();
  std::cout << "* Wrote " << outfileName() << "." << std::endl;
}


void TrackingTask::outfileAdd(const double &t, const arma::colvec3 &s)
{
  *outfile <<std::resetiosflags(std::ios::fixed)<<std::setiosflags(std::ios::scientific)
	   <<std::showpoint<<std::setprecision(8)<<std::setw(w+2)<< t
	   <<std::resetiosflags(std::ios::scientific)<<std::setiosflags(std::ios::fixed)<<std::setprecision(5)
	   <<std::setw(w)<< s[0] <<std::setw(w)<< s[1] <<std::setw(w)<< s[2]
	   <<std::setw(w)<< norm(s,2) <<std::setw(w)<< "gamma" << std::endl;
}



void TrackingTask::storeStep(const double &t, const arma::colvec3 &s)
{
  outfileAdd(t,s);
}
