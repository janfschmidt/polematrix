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

#ifndef __POLEMATRIX__TRACKINGTASK_HPP_
#define __POLEMATRIX__TRACKINGTASK_HPP_

#include <fstream>
#include <stdexcept>
#include <map>
#include <memory>
#include <functional>
#define ARMA_NO_DEBUG
#include <armadillo>
#include <libpalattice/AccLattice.hpp>
#include <libpalattice/FunctionOfPos.hpp>
#include <libpalattice/types.hpp>
#include "Simulation.hpp"
#include "RadiationModel.hpp"
#include "Trajectory.hpp"


// spin tracking result container (3d spin vector as function of time)
class SpinMotion : public std::map<double,arma::colvec3>
{
public:
  void operator+=(const SpinMotion &other); // operators for calculation of polarization (average over spins)
  void operator/=(const unsigned int &num);
  
  std::string printHeader(unsigned int columnWidth, std::string name="S") const; //name for vector
  std::string print(unsigned int columnWidth) const;
  std::string printLine(unsigned int columnWidth, const double &key) const;
  //print a line of external data. used for faster "online" file output during tracking
  std::string printAnyData(unsigned int columnWidth, const double &t, const arma::colvec3 &s) const;
};



class TrackingTask : public SingleParticleSimulation
{
private:
  arma::mat33 one;
  SpinMotion storage;                         // store results
  std::unique_ptr<std::ofstream> outfile;     // output file via pointer, std::ofstream not moveable in gcc 4.9
  std::unique_ptr<std::ofstream> outfile_ps;  // output file for long. phase space (gammaMode radiation only)
  unsigned int w;                             // output column width (print)
  bool completed;                             // tracking completed
  pal::FunctionOfPos<double> gammaSimTool;    // gamma(pos) from elegant
  double gammaSimToolCentral;                 // gamma central from elegant (set energy)
  LongitudinalPhaseSpaceModel syliModel;      // for gammaMode "radiation"
  
  //variables for current tracking step
  pal::AccLattice::const_iterator currentElement; // position in lattice
  double currentGamma;                            // gamma

  arma::running_stat<double> gammaStat;       // gamma statistics
  
  void outfileOpen();                         // open output file and write header
  void outfileClose();                        // write footer and close output file
  void outfileAdd(const double &t, const arma::colvec3 &s);  // append s(t) to outfile
  void storeStep(const double &pos, const arma::colvec3 &s); // append s(t) to storage and outfile
  void outfileAdd_ps(const double &pos);                     // append long. phase space(t) to outfile_ps

  void checkLongStability() const;            // check if longitudinal motion is stable (gammaMode "radiation")

  
public:
  TrackingTask(unsigned int id, const std::shared_ptr<Configuration> c);
  TrackingTask(const TrackingTask& other) = delete;
  TrackingTask(TrackingTask&& other) = default;
  ~TrackingTask() {}

  void run();                                 //run tracking task
  void matrixTracking();

  // particle energy gamma(pos), implementation depends GammaMode
  double (TrackingTask::*gamma)(const double&);
  void initGamma();
  void saveGammaSimTool();
  
  // gamma(pos) implementations:
  double gammaFromConfig(const double &pos) {return config->gamma(pos/GSL_CONST_MKSA_SPEED_OF_LIGHT);}
  double gammaFromSimTool(const double &pos) {return gammaSimTool.interpPeriodic(pos-config->pos_start());}
  double gammaFromSimToolPlusConfig(const double &pos) {return gammaFromSimTool(pos) - gammaSimToolCentral + gammaFromConfig(pos); }
  double gammaFromSimToolNoInterpolation(const double &pos) {return gammaSimTool.infrontof(pos-config->pos_start());}
  double gammaRadiation(const double &pos);
  double gammaOffset(const double &pos) {return gammaFromConfig(pos) + syliModel.gammaMinusGamma0();}
  double gammaOscillation(const double &pos) {return gammaFromConfig(pos) + syliModel.gammaMinusGamma0()*cos(2*M_PI*syliModel.synchrotronFreq_current()*pos/GSL_CONST_MKSA_SPEED_OF_LIGHT + particleId);} // uses particleId for individual start phases

  inline arma::mat33 rotxMatrix(double angle) const;
  inline arma::mat33 rotMatrix(pal::AccTriple B) const;
  
  std::string outfileName() const;            // output file name
  std::string phasespaceOutfileName() const; // phase space output file name

  SpinMotion getStorage() const {return storage;}
  double getProgress() const {return (double)storage.size() / config->outSteps();}
  bool isCompleted() const {return completed;}
  
};



// exceptions
class TrackError : public std::runtime_error {
public:
  TrackError(std::string msg) : std::runtime_error(msg) {}
};

class TrackFileError : public TrackError {
public:
  TrackFileError(std::string file) : TrackError("Cannot open "+file) {}
};



#endif
// __POLEMATRIX__TRACKINGTASK_HPP_
