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
#include "Configuration.hpp"


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



class TrackingTask
{
private:
  arma::mat33 one;
  SpinMotion storage;                         // store results
  std::unique_ptr<std::ofstream> outfile;     // output file via pointer, std::ofstream not moveable in gcc 4.9
  unsigned int w;                             // output column width (print)
  bool completed;                             // tracking completed
  pal::FunctionOfPos<double> gammaSimTool;    // gamma(pos) from elegant
  double gammaCentralSimTool;                // gamma central from elegant (set energy)
  
  void outfileOpen();                         // open output file and write header
  void outfileClose();                        // write footer and close output file
  void outfileAdd(const double &t, const arma::colvec3 &s, const double &agamma); // append s(t) to outfile
  void storeStep(const double &t, const arma::colvec3 &s, const double &agamma);  // append s(t) to storage and outfile

  
public:
  const unsigned int particleId;
  Configuration &config;
  const pal::AccLattice *lattice;
  const pal::FunctionOfPos<pal::AccPair> *orbit;
  
  TrackingTask(unsigned int id, Configuration &c);
  TrackingTask(const TrackingTask& other) = delete;
  TrackingTask(TrackingTask&& other) = default;
  ~TrackingTask() {}
  
  void run();                                 //run tracking task
  void matrixTracking();
  
  double (TrackingTask::*gamma)(double) const;
  //gamma modes:
  double gammaFromConfig(double pos) const {return config.gamma(pos/GSL_CONST_MKSA_SPEED_OF_LIGHT);}
  double gammaFromSimTool(double pos) const {return gammaSimTool.interpPeriodic(pos-config.pos_start());}
  double gammaFromSimToolPlusConfig(double pos) const {return gammaFromSimTool(pos) - gammaCentralSimTool + gammaFromConfig(pos); }
  
  void saveGammaSimTool();

  inline arma::mat33 rotxMatrix(double angle) const;
  inline arma::mat33 rotMatrix(pal::AccTriple B) const;
  
  std::string outfileName() const;            // output file name

  SpinMotion getStorage() const {return storage;}
  std::string getProgressBar(unsigned int barWidth=20) const;
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
