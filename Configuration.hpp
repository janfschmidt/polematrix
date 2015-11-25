#ifndef __POLEMATRIX__CONFIGURATION_HPP_
#define __POLEMATRIX__CONFIGURATION_HPP_

#include <string>
#include <armadillo>
#include <gsl/gsl_const_mksa.h>
#include <fstream>
#include <vector>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>
#include <libpalattice/SimTools.hpp>

namespace pt = boost::property_tree;
namespace fs = boost::filesystem;


enum GammaMode{linear, simtool, simtool_plus_linear, radiation};

class Configuration
{
private:
  //palattice
  std::shared_ptr<pal::SimToolInstance> palattice;
  std::vector<bool> _saveGamma;

  pal::SimTool toolFromTree(pt::ptree tree, std::string key) const;
  pal::SimToolMode modeFromTree(pt::ptree tree, std::string key) const;
  void setSimToolInstance(pt::ptree &tree);
  void setGammaMode(pt::ptree &tree);

  fs::path _outpath;

  //spintracking
  arma::colvec3 _s_start;
  double _t_start;          // time / s
  double _t_stop;
  double _dt_out;           // output(!) step width / s
  double _E0;               // energy at t=0 / GeV
  double _dE;               // dE/dt / GeV/s
  unsigned int _nParticles; // number of tracked particles
  GammaMode _gammaMode;
  int _seed;                 // random number seed for gammaMode "radiation"

  
public:
  //constants / internal configuration (constructor)
  const double E_rest_GeV;           // electron rest energy / GeV
  const double E_rest_keV;           // electron rest energy / keV
  const double a_gyro;               // electron gyromagnetic anomaly a = (g-2)/2
  const unsigned int default_steps;  // number of output steps if not specified by dt_out
  const std::string spinDirName;     // directory name for tracking output files (outpath/spinDirName)
  const std::string polFileName;     // file name for polarization output file (Tracking::savePolarization())
  const std::string confOutFileName; // file name for config output file (written by Tracking::start())

  Configuration(std::string path=".");

  //getter
  fs::path outpath() const {return _outpath;}
  arma::colvec3 s_start() const {return _s_start;}
  double t_start() const {return _t_start;}
  double t_stop() const {return _t_stop;}
  double dt_out() const {return _dt_out;}
  double E0() const {return _E0;}
  double dE() const {return _dE;}
  unsigned int nParticles() const {return _nParticles;}
  GammaMode gammaMode() const {return _gammaMode;}
  int seed() const {return _seed;}
  pal::SimToolInstance& getSimToolInstance() {return *palattice;}
  bool saveGamma(unsigned int particleId) const {return _saveGamma.at(particleId);}
  
  //setter
  void set_outpath(fs::path p) {_outpath=p;}
  void set_s_start(arma::colvec3 s) {_s_start=arma::normalise(s);}
  void set_t_start(double t) {_t_start=t;}
  void set_t_stop(double t) {_t_stop=t;}
  void set_dt_out(double dt) {_dt_out=dt;}
  void set_E0(double E) {_E0=E;}
  void set_dE(double dEin) {_dE=dEin;}
  void set_nParticles(unsigned int n);
  void set_gammaMode(GammaMode g) {_gammaMode=g;}
  void set_saveGamma(std::string particleList);
  void set_seed(int s) {_seed=s;}

  double duration() const {return t_stop() - t_start();}
  fs::path subDirectory(std::string folder) const {return outpath()/folder;}
  fs::path spinDirectory() const {return outpath()/spinDirName;}
  fs::path polFile() const {return outpath()/polFileName;}
  fs::path confOutFile() const {return outpath()/confOutFileName;}
  double pos_start() const {return GSL_CONST_MKSA_SPEED_OF_LIGHT * t_start();}
  double pos_stop() const {return GSL_CONST_MKSA_SPEED_OF_LIGHT * t_stop();}
  double dpos_out() const {return GSL_CONST_MKSA_SPEED_OF_LIGHT * dt_out();}
  double gamma_start() const {return gamma(t_start());}
  double gamma_stop() const {return gamma(t_stop());}
  double agamma_start() const {return agamma(t_start());}
  double agamma_stop() const {return agamma(t_stop());}
  unsigned int outSteps() const {return (t_stop()-t_start())/dt_out();}
  std::string saveGammaList() const;

  void printSummary() const;

  double gamma(double t) const;
  double agamma(double t) const;

  // save & load configuration to/from file
  void save(const std::string &filename) const;
  void load(const std::string &filename);

protected:
  int randomSeed() const;
};




#endif
// __POLEMATRIX__CONFIGURATION_HPP_
