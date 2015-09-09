#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include "TrackingTask.hpp"


void SpinMotion::operator+=(const SpinMotion &other)
{
  if (this->size() != other.size())
    throw std::runtime_error("SpinMotion::operator+= not possible for objects of different size");
  
  std::map<double,arma::colvec3>::const_iterator otherIt = other.begin();
    
  for (auto it = this->begin(); it != this->end(); ++it) {
    if (it->first == otherIt->first) {
      it->second += otherIt->second;
      otherIt++;
    }
    else
      throw std::runtime_error("SpinMotion::operator+= with incompatible tracking time steps");
  }
}

void SpinMotion::operator/=(const unsigned int &num)
{
  for (auto it = this->begin(); it!= this->end(); ++it) {
    it->second /= num;
  }
}


std::string SpinMotion::printHeader(unsigned int w, std::string name) const
{
  std::stringstream ss;
  ss << "#"<<std::setw(w+1)<< "t / s" <<std::setw(w)<< name+"x" <<std::setw(w)<< name+"z"
     <<std::setw(w)<< name+"s" <<std::setw(w)<< "|"+name+"|";
  return ss.str();
}

std::string SpinMotion::print(unsigned int w) const
{
  std::stringstream ss;
  for (auto it = this->begin(); it!= this->end(); ++it) {
    ss << std::resetiosflags(std::ios::fixed)<<std::setiosflags(std::ios::scientific)
       <<std::showpoint<<std::setprecision(8)<<std::setw(w+2)<< it->first
       <<std::resetiosflags(std::ios::scientific)<<std::setiosflags(std::ios::fixed)<<std::setprecision(5)
       <<std::setw(w)<< it->second[0] <<std::setw(w)<< it->second[2] <<std::setw(w)<< it->second[1]
       <<std::setw(w)<< arma::norm(it->second) << std::endl;
  }
    return ss.str();
}

std::string SpinMotion::printLine(unsigned int w, const double &key) const
{
  std::stringstream ss;
  ss << std::resetiosflags(std::ios::fixed)<<std::setiosflags(std::ios::scientific)
     <<std::showpoint<<std::setprecision(8)<<std::setw(w+2)<< key
     <<std::resetiosflags(std::ios::scientific)<<std::setiosflags(std::ios::fixed)<<std::setprecision(5)
     <<std::setw(w)<< this->at(key)[0] <<std::setw(w)<< this->at(key)[2] <<std::setw(w)<< this->at(key)[1]
     <<std::setw(w)<< arma::norm(this->at(key));
  return ss.str();
}

std::string SpinMotion::printAnyData(unsigned int w, const double &t, const arma::colvec3 &s) const
{
  std::stringstream ss;
  ss << std::resetiosflags(std::ios::fixed)<<std::setiosflags(std::ios::scientific)
     <<std::showpoint<<std::setprecision(8)<<std::setw(w+2)<< t
     <<std::resetiosflags(std::ios::scientific)<<std::setiosflags(std::ios::fixed)<<std::setprecision(5)
     <<std::setw(w)<< s[0] <<std::setw(w)<< s[2] <<std::setw(w)<< s[1]
     <<std::setw(w)<< arma::norm(s);
  return ss.str();
}






TrackingTask::TrackingTask(unsigned int id, const Configuration &c)
  : particleId(id), config(c), w(14), completed(false)
{
  one.eye(); // fill unit matrix
  outfile = std::unique_ptr<std::ofstream>(new std::ofstream());
}



void TrackingTask::run()
{
  outfileOpen();
  
  //std::cout << "* start tracking particle " << particleId << std::endl;

  matrixTracking();
  
  outfileClose();

  completed = true;
}


void TrackingTask::matrixTracking()
{
  arma::colvec3 s = config.s_start;
  pal::AccTriple omega;
  double pos = config.pos_start();
  double t = pos/GSL_CONST_MKSA_SPEED_OF_LIGHT;
  double pos_stop = config.pos_stop();
  double dpos_out = config.dpos_out();
  double pos_nextOut = pos + dpos_out;
  double gamma;

  // set current iterator and position
  pal::const_AccIterator it=lattice->nextCIt( orbit->posInTurn(pos) );
  pos = (orbit->turn(pos)-1)*lattice->circumference() + lattice->pos(it);
  t = pos/GSL_CONST_MKSA_SPEED_OF_LIGHT;

  //  for (unsigned int i=0; i<this->particleId*123456789; i++) {
  while (pos < pos_stop) {
    gamma = config.gamma(t);
    omega = it->second->B_int( orbit->interp( orbit->posInTurn(pos) ) ) * config.a_gyro; // field of element
    omega.x *= gamma;
    omega.z *= gamma;
    // omega.s: Precession around s is suppressed by factor gamma (->TBMT-equation)

    s = rotMatrix(omega) * s; // spin rotation
    //std::cout << rotMatrix(omega) << std::endl;

    //renormalize
    //s = s/std::sqrt(std::pow(s(0),2) + std::pow(s(1),2) + std::pow(s(2),2));

    if (pos >= pos_nextOut) {//output
      storeStep(t,s,gamma);
      pos_nextOut += dpos_out;
    }

    // step to next element
    pos += lattice->distanceNext(it);
    t = pos/GSL_CONST_MKSA_SPEED_OF_LIGHT;
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
  arma::colvec3 B = {B_in.x,B_in.s,B_in.z};
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
  //  ss << config.subfolder("spins") << "spin_" << std::setw(4)<<std::setfill('0')<<particleId << ".dat";
    ss << "spin_" << std::setw(4)<<std::setfill('0')<<particleId << ".dat";
    return ( config.spinDirectory()/ss.str() ).string();
}



void TrackingTask::outfileOpen()
{
  if ( fs::create_directory(config.spinDirectory()) )
       std::cout << "* created directory " << config.spinDirectory() << std::endl;
  outfile->open(outfileName());
  if (!outfile->is_open())
    throw TrackFileError(outfileName());
  else
    *outfile << storage.printHeader(w) <<std::setw(w)<< "gamma" << std::endl;
}


void TrackingTask::outfileClose()
{
  outfile->close();
  std::cout << "* Wrote " << storage.size() << " steps to " << outfileName()
	    <<std::setw(40)<<std::left<< "." << std::endl;
}


void TrackingTask::outfileAdd(const double &t, const arma::colvec3 &s, const double &gamma)
{
  *outfile << storage.printAnyData(w,t,s) << std::setw(w)<< gamma << std::endl;
}



void TrackingTask::storeStep(const double &t, const arma::colvec3 &s, const double &gamma)
{
  storage.insert(std::pair<double,arma::colvec3>(t,s));
  outfileAdd(t,s,gamma);
}


std::string TrackingTask::getProgressBar() const
{
  std::stringstream bar;
  double stepsTotal = 20;
  unsigned int steps = (1.0 * stepsTotal * storage.size() / config.outSteps()) + 0.5;
  unsigned int i=0;

  bar << "particle " << particleId << " [";
  
  for (; i<steps; i++)
    bar << "=";
  for (; i<stepsTotal; i++)
    bar << " ";
    
  bar << "] " << steps/stepsTotal*100 << "%";
  
  return bar.str();
}
