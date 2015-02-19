#include <iostream>
#include <iomanip>
#include <cmath>
#include "TrackingTask.hpp"


TrackingTask::TrackingTask(unsigned int id, const Configuration *c)
  : particleId(id), config(c), w(14)
{
  one.eye(); // fill unit matrix
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
  arma::colvec3 s = config->s_start;
  pal::AccTriple omega;
  double pos = config->pos_start();
  double pos_stop = config->pos_stop();
  double dpos_out = config->dpos_out();

  // set current iterator and position
  pal::const_AccIterator it=lattice->nextCIt( orbit->posInTurn(pos) );
  pos = (orbit->turn(pos)-1)*lattice->circumference() + lattice->pos(it);
  //B.x=config->x*M_PI/180.;

  //  for (unsigned int i=0; i<this->particleId*123456789; i++) {
  while (pos < pos_stop) {
    omega = it->second->B_int( orbit->interp( orbit->posInTurn(pos) ) ); // field of element
    omega *= config->agamma(0.);
    s = rotMatrix(omega) * s; // spin rotation
    std::cout << rotMatrix(omega) << std::endl;

    // if (std::fmod(pos,dpos_out) == 0) //output
      storeStep(pos,s);

    // step to next element
    pos += lattice->distanceNext(it);
    it = lattice->revolve(it);
  }
}

arma::mat33 TrackingTask::rotxMatrix(double angle) const
{
  double c=std::cos(angle);
  double s=std::sin(angle);

  arma::mat33 rotx = {1,0,0, 0,c,s, 0,-s,c};
  return rotx;
}


arma::mat33 TrackingTask::rotMatrix(pal::AccTriple B_in) const
{
  arma::colvec3 B = {B_in.x,B_in.z,B_in.s};
  double angle = std::sqrt(std::pow(B(0),2) + std::pow(B(1),2) + std::pow(B(2),2)); //faster than arma::norm(B);
  if (angle < MIN_AMPLITUDE) return one;

  arma::colvec3 n = B/angle; //faster than arma::normalise(B)

  double c=std::cos(angle);
  double onemc=1-c;
  double s=std::sin(angle);

  arma::mat33 rot = {std::pow(n(0),2)*onemc+c, n(1)*n(0)*onemc+n(2)*s, n(2)*n(0)*onemc-n(1)*s,  //column 1
		     n(0)*n(1)*onemc-n(2)*s, std::pow(n(1),2)*onemc+c, n(2)*n(1)*onemc+n(0)*s,  //column 2
                     n(0)*n(2)*onemc+n(1)*s, n(1)*n(2)*onemc-n(0)*s, std::pow(n(2),2)*onemc+c}; //column 3
  return rot;
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
	   <<std::setw(w)<< arma::norm(s) <<std::setw(w)<< "gamma" << std::endl;
}



void TrackingTask::storeStep(const double &t, const arma::colvec3 &s)
{
  storage.insert(std::pair<double,arma::colvec3>(t,s));
  outfileAdd(t,s);
}
