/* TrackingTask Class
 * It representing a single spin tracking and is managed by an object of type Tracking
 *
 * Copyright (C) 2017 Jan Felix Schmidt <janschmidt@mailbox.org>
 *   
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <gsl/gsl_spline.h>
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






TrackingTask::TrackingTask(unsigned int id, const std::shared_ptr<Configuration> c)
  : SingleParticleSimulation(id,c), w(14), completed(false),
    gammaSimTool(config->getSimToolInstance(), gsl_interp_akima),
    syliModel(config->seed()+particleId, config),
    currentElement(pal::AccLattice().begin())
    // pal::AccLattice::const_iterator currentElement is initialized with empty lattice (dirty)!
{
  one.eye(); // fill unit matrix
  outfile = std::unique_ptr<std::ofstream>(new std::ofstream());
  outfile_ps = std::unique_ptr<std::ofstream>(new std::ofstream());

  switch (config->gammaMode()) {
  case GammaMode::simtool:
    gamma = &TrackingTask::gammaFromSimTool;
    break;
  case GammaMode::simtool_plus_linear:
    gamma = &TrackingTask::gammaFromSimToolPlusConfig;
    break;
  case GammaMode::simtool_no_interpolation:
    gamma = &TrackingTask::gammaFromSimToolNoInterpolation;
    break;
  case GammaMode::radiation:
    gamma = &TrackingTask::gammaRadiation;
    break;
  case GammaMode::offset:
    gamma = &TrackingTask::gammaOffset;
    break;
      case GammaMode::oscillation:
    gamma = &TrackingTask::gammaOscillation;
    break;
  default:
    gamma = &TrackingTask::gammaFromConfig;
  }
}




void TrackingTask::run()
{
  initGamma();
  trajectory->init();
  
  outfileOpen();

  // std::cout << "* start tracking particle " << particleId << std::endl;
  matrixTracking();
  
  outfileClose();

  // clear interpolation to save memory
  gammaSimTool.clear();
  trajectory->clear();
  
  completed = true;
}

void TrackingTask::initGamma()
{
  if ( config->gammaMode()==GammaMode::simtool
       || config->gammaMode()==GammaMode::simtool_plus_linear
       || config->gammaMode()==GammaMode::simtool_no_interpolation )
    {
      // simtool: sdds import thread safe since SDDSToolKit-devel-3.3.1-2
      gammaSimTool.readSimToolParticleColumn( config->getSimToolInstance(), particleId+1, "p" );
      gammaSimToolCentral = config->getSimToolInstance().readGammaCentral();
      if (config->gammaMode() != GammaMode::simtool_no_interpolation) {
	    gammaSimTool.init();
      }
      saveGammaSimTool();
    }
  else if ( config->gammaMode()==GammaMode::radiation
	    || config->gammaMode()==GammaMode::offset
	    || config->gammaMode()==GammaMode::oscillation ) {
    syliModel.init(lattice);
  }
  //else: no init needed
}



void TrackingTask::matrixTracking()
{
  arma::colvec3 s = config->s_start();
  pal::AccTriple omega;
  double pos = config->pos_start();
  double pos_stop = config->pos_stop();
  double dpos_out = config->dpos_out();
  double pos_nextOut = pos;

  // set start lattice element and position
  currentElement = lattice->behind( orbit->posInTurn(pos), pal::Anchor::end );
  pos = (orbit->turn(pos)-1)*lattice->circumference() + currentElement.pos();

  while (pos < pos_stop) {
    currentGamma = (this->*gamma)(pos);
    auto rf = currentElement.element()->rfFactor(orbit->turn(pos));
    auto Bint = currentElement.element()->B_int( trajectory->get(pos) ) * rf;  // field of element
    // Dipole: Integral field including Bx from edge focussing (! uses vertical trajectory at "pos" for magnet entrance and exit)
    if (config->edgefoc() && currentElement.element()->type == pal::dipole) {
      Bint.x -= ( tan(currentElement.element()->e1) + tan(currentElement.element()->e2) )/(1./currentElement.element()->k0.z) * trajectory->get(pos).z;
    }
    omega = Bint * config->a_gyro;
    omega.x *= currentGamma;
    omega.z *= currentGamma;
    // omega.s: Precession around s is suppressed by factor gamma (->TBMT-equation)

    // spin rotation
    s = rotMatrix(omega) * s;

    // output
    if (pos >= pos_nextOut) {
      if (!config->outElementUsed() || currentElement.element()->name == config->outElement()) {
	storeStep(pos,s);
	pos_nextOut += dpos_out;
      }
    }
    gammaStat(currentGamma);

    //long. phase space output is in TrackingTask::gammaRadiation() -> called above via (this->*gamma)(pos)

    // step to next element
    pos += currentElement.distanceNext();
    currentElement.revolve();
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


void TrackingTask::saveGammaSimTool()
{
  if ( !config->saveGamma(particleId) )
    return;
  
  gammaSimTool.info.add("polematrix particle ID", particleId);
  std::stringstream file;
  file <<std::setw(4)<<std::setfill('0')<< particleId << ".dat";
  gammaSimTool.print( (config->outpath()/"gammaSimTool_").string() + file.str() );

  trajectory->saveSimtoolData();
}


std::string TrackingTask::outfileName() const
{
  std::stringstream ss;
  //  ss << config->subfolder("spins") << "spin_" << std::setw(4)<<std::setfill('0')<<particleId << ".dat";
    ss << "spin_" << std::setw(4)<<std::setfill('0')<<particleId << ".dat";
    return ( config->spinDirectory()/ss.str() ).string();
}
std::string TrackingTask::phasespaceOutfileName() const
{
  std::stringstream ss;
    ss << "longPhaseSpace_" << std::setw(4)<<std::setfill('0')<<particleId << ".dat";
    return ( config->outpath()/ss.str() ).string();
}



void TrackingTask::outfileOpen()
{
  if ( fs::create_directory(config->spinDirectory()) )
       std::cout << "* created directory " << config->spinDirectory() << std::endl;
  outfile->open(outfileName());
  if (!outfile->is_open())
    throw TrackFileError(outfileName());

  *outfile << config->metadata();
  *outfile << storage.printHeader(w) <<std::setw(w)<< "gamma" << std::endl;

  if (config->gammaMode()==GammaMode::radiation) {
    if (config->savePhaseSpace(particleId)) {
      outfile_ps->open(phasespaceOutfileName());
      if (!outfile_ps->is_open())
	throw TrackFileError(phasespaceOutfileName());

      *outfile_ps << config->metadata();
      *outfile_ps << "# longitudinal phase space at " << config->savePhaseSpaceElement() << ", particleId " << particleId << std::endl;
      *outfile_ps <<"# " <<std::setw(w-2)<< "t / s" <<std::setw(w)<< "dphase / rad" <<std::setw(w)<< "dgamma/gamma0" << std::endl;
    }
  }
}


void TrackingTask::outfileClose()
{
  *outfile << "# gamma statistics:" << std::endl
	   << "# mean:  " << gammaStat.mean() << std::endl
	   << "# stddev: " << gammaStat.stddev(1) << std::endl;
    
  outfile->close();
  if (config->verbose()) {
    std::cout << "* " << storage.size() << " steps written to " << outfileName()
	      <<std::setw(40)<<std::left<< "." << std::endl;
  }

  if (outfile_ps->is_open()) {
    outfile_ps->close();
    std::cout << "* " << phasespaceOutfileName() << " written" << std::endl;
  }
}


void TrackingTask::outfileAdd(const double &t, const arma::colvec3 &s)
{
  *outfile << storage.printAnyData(w,t,s) << std::setw(w)<< currentGamma;
  if (config->gammaMode() == GammaMode::radiation)
    *outfile <<std::setw(w)<< syliModel.phase();
  *outfile << std::endl;
}

void TrackingTask::outfileAdd_ps(const double &pos)
{
  double t = pos/GSL_CONST_MKSA_SPEED_OF_LIGHT;
  *outfile_ps <<std::setw(w)<< t
	      <<std::setw(w)<< syliModel.phase() - syliModel.ref_phase()
	      <<std::setw(w)<< (currentGamma - syliModel.gamma0())/syliModel.gamma0() << std::endl;
}



void TrackingTask::storeStep(const double &pos, const arma::colvec3 &s)
{
  double t = pos/GSL_CONST_MKSA_SPEED_OF_LIGHT;
  storage.insert(std::pair<double,arma::colvec3>(t,s));
  outfileAdd(t,s);
}




double TrackingTask::gammaRadiation(const double &pos)
{
  // long. phase space output
  if (config->savePhaseSpace(particleId) && currentElement.element()->name == config->savePhaseSpaceElement()) {
    outfileAdd_ps(pos);
  }
  
  syliModel.update(currentElement.element(), pos, gammaFromConfig(pos));
  return syliModel.gamma();
}

