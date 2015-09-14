#ifndef __POLEMATRIX__CONFIGURATION_HPP_
#define __POLEMATRIX__CONFIGURATION_HPP_

#include <string>
#include <armadillo>
#include <gsl/gsl_const_mksa.h>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <libpalattice/SimTools.hpp>

namespace pt = boost::property_tree;
namespace fs = boost::filesystem;


enum GammaMode{linear, simtool};

class Configuration
{
private:
  //palattice
  pal::SimToolInstance *palattice;

  pal::SimTool toolFromTree(pt::ptree tree, std::string key) const;
  pal::SimToolMode modeFromTree(pt::ptree tree, std::string key) const;
  void setSimToolInstance(pt::ptree &tree);
  void setGammaMode(pt::ptree &tree);

  
public:
  fs::path outpath;

  //spintracking
  arma::colvec3 s_start;
  double t_start;          // time / s
  double t_stop;
  double dt_out;           // output(!) step width / s
  double E0;               // energy at t=0 / GeV
  double dE;               // dE/dt / GeV/s
  unsigned int nParticles; // number of tracked particles
  GammaMode gammaMode;

  //constants / internal configuration (constructor)
  const double E_rest;               // electron rest energy / GeV
  const double a_gyro;               // electron gyromagnetic anomaly a = (g-2)/2
  const unsigned int default_steps;  // number of output steps if not specified by dt_out
  const std::string spinDirName;     // directory name for tracking output files (outpath/spinDirName)
  const std::string polFileName;     // file name for polarization output file (Tracking::savePolarization())
  const std::string confOutFileName; // file name for config output file (written by Tracking::start())

  Configuration(std::string path=".");
  ~Configuration();

  pal::SimToolInstance& getSimToolInstance() const {return *palattice;}  
  fs::path subDirectory(std::string folder) const {return outpath/folder;}
  fs::path spinDirectory() const {return outpath/spinDirName;}
  fs::path polFile() const {return outpath/polFileName;}
  fs::path confOutFile() const {return outpath/confOutFileName;}
  double pos_start() const {return GSL_CONST_MKSA_SPEED_OF_LIGHT * t_start;}
  double pos_stop() const {return GSL_CONST_MKSA_SPEED_OF_LIGHT * t_stop;}
  double dpos_out() const {return GSL_CONST_MKSA_SPEED_OF_LIGHT * dt_out;}
  unsigned int outSteps() const {return (t_stop-t_start)/dt_out;}

  void printSummary() const;

  double gamma(double t) const;
  double agamma(double t) const;

  // save & load configuration to/from file
  void save(const std::string &filename) const;
  void load(const std::string &filename);
  
};




#endif
// __POLEMATRIX__CONFIGURATION_HPP_
