/* Configuration Class
 * stores program configuration and can be loaded from & saved to xml file
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
#include <libpalattice/AccLattice.hpp>
#include <libpalattice/Metadata.hpp>
#include "version.hpp"

namespace pt = boost::property_tree;
namespace fs = boost::filesystem;


enum class GammaMode{linear, offset, oscillation, radiation, simtool, simtool_plus_linear, simtool_no_interpolation};
enum class TrajectoryMode{closed_orbit, simtool};



// parse configuration of multiple rf magnets from comma separated lists
// in a polematrix config file (boost property tree)
// and write them to a pal::AccLattice
class RfMagnetConfig {
protected:
  std::vector<std::string> elements; // element names
  std::vector<double> Q1;            //  tunes for turn 1 (AccElement::Qrf1)
  std::vector<double> dQ;            //  tune changes per turn (AccElement::dQrf)
  std::vector<unsigned int> period;  //  sweep lengths in turns (AccElement::rfPeriod)

  template<typename T> std::string getStringList(const std::vector<T> v) const;
  template<typename T> void setFromStringList(std::vector<T>& v, std::string list);
  
public:
  RfMagnetConfig() {};
  ~RfMagnetConfig() {};

  void set(const pt::ptree &tree);
  void writeToLattice(pal::AccLattice& lattice) const;
  void writeToConfig(pt::ptree &tree) const;

  // comma separated string getters to write config
  std::string getElements() const {return getStringList<std::string>(elements);}
  std::string getQ1() const {return getStringList<double>(Q1);}
  std::string getDQ() const {return getStringList<double>(dQ);}
  std::string getPeriod() const {return getStringList<unsigned int>(period);}
};







class Configuration
{
private:
  //palattice
  std::shared_ptr<pal::SimToolInstance> palattice;
  std::vector<bool> _saveGamma;
  double _circumference;
  bool _simToolRamp;
  unsigned int _simToolRampSteps;

  pal::SimTool toolFromTree(pt::ptree tree, std::string key) const;
  pal::SimToolMode modeFromTree(pt::ptree tree, std::string key) const;
  void setSimToolInstance(pt::ptree &tree);
  void setGammaMode(pt::ptree &tree);
  void setTrajectoryMode(pt::ptree &tree);

  //not in config file (cmdline options)
  fs::path _outpath;
  bool _verbose;

  //spintracking
  arma::colvec3 _s_start;
  double _t_start;          // time / s
  double _t_stop;
  double _dt_out;           // output(!) step width / s
  double _E0;               // energy at t=0 / GeV
  double _dE;               // dE/dt / GeV/s
  double _Emax;             // ramp end energy / GeV
  unsigned int _nParticles; // number of tracked particles
  GammaMode _gammaMode;
  TrajectoryMode _trajectoryMode;
  bool _edgefoc;            // edge focussing field (Bx) of Dipoles included ?
  std::string _outElement;  // output at the lattice element with this name only
                            // (wait for next occurrence after dt_out)
  bool _outElementUsed;

  //rf magnets
  RfMagnetConfig rf;
  
  //radiation (used with gammaMode radiation only)
  int _seed;                // random number seed
  double _q;                 // over voltage factor
  double _alphac;           // momentum compaction factor
  double _alphac2;          // 2nd order momentum compaction factor
  unsigned int _h;           // harmonic number (number of buckets, h = f_rf/f_rev)
  double _R;                // (mean) dipole bending radius / m
  double _Js;               // longitudinal damping partition number
  std::vector<bool> _savePhaseSpace;  // configure particles to write long. phase space to file
  std::string _savePhaseSpaceElement; // Name of Lattice Element at which long. phase space is recorded
  double _sigmaPhaseFactor; // start value for sigma_phase in units of equilibrium value
  double _sigmaGammaFactor; // start value for sigma_gamma in units of equilibrium value

  //Resonance Strengths
  double _agammaMin;    // output spin tune range
  double _agammaMax;
  double _dagamma;      // spin tune output step width
  unsigned int _nTurns; // number of turns for res.strengths calculation

  pal::Metadata info;
  
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
  bool verbose() const {return _verbose;}
  arma::colvec3 s_start() const {return _s_start;}
  double t_start() const {return _t_start;}
  double t_stop() const {return _t_stop;}
  double dt_out() const {return _dt_out;}
  std::string outElement() const {return _outElement;}
  bool outElementUsed() const {return _outElementUsed;}
  double E0() const {return _E0;}
  double dE() const {return _dE;}
  double Emax() const {return _Emax;}
  unsigned int nParticles() const {return _nParticles;}
  GammaMode gammaMode() const {return _gammaMode;}
  std::string gammaModeString() const;
  std::string trajectoryModeString() const;
  TrajectoryMode trajectoryMode() const {return _trajectoryMode;}
  bool edgefoc() const {return _edgefoc;}
  int seed() const {return _seed;}
  double q() const {return _q;}
  double alphac() const {return _alphac;}
  double alphac2() const {return _alphac2;}
  unsigned int h() const {return _h;}
  double R() const {return _R;}
  double Js() const {return _Js;}
  pal::SimToolInstance& getSimToolInstance() {return *palattice;}
  bool simToolRamp() const {return _simToolRamp;}
  unsigned int simToolRampSteps() const {return _simToolRampSteps;}
  bool saveGamma(unsigned int particleId) const {return _saveGamma.at(particleId);}
  bool savePhaseSpace(unsigned int particleId) const {return _savePhaseSpace.at(particleId);}
  std::string savePhaseSpaceElement() const {return _savePhaseSpaceElement;}
  double sigmaPhaseFactor() const {return _sigmaPhaseFactor;}
  double sigmaGammaFactor() const {return _sigmaGammaFactor;}
  double agammaMin() const {return _agammaMin;}
  double agammaMax() const {return _agammaMax;}
  std::string metadata() const {return info.out("#");}
  // #turns for resonance strengths calc are calculated from tracking duration() if not set
  unsigned int numTurns() const;
  // spin tune output step width for resonance strengths is calculated from #turns if not set
  double dagamma() const;
  
  //setter
protected:
  void set_saveList(const std::string &particleList, std::vector<bool> &list, const std::string &optionName);
  void set_metadata(const std::string &configfile);
public:
  void set_outpath(fs::path p) {_outpath=p;}
  void set_verbose(bool v=true) {_verbose=v;}
  void set_s_start(arma::colvec3 s) {_s_start=arma::normalise(s);}
  void set_t_start(double t) {_t_start=t;}
  void set_t_stop(double t) {_t_stop=t;}
  void set_dt_out(double dt) {_dt_out=dt;}
  void set_outElement(std::string name) {_outElement=name; _outElementUsed=true;}
  void set_E0(double E) {_E0=E;}
  void set_dE(double dEin) {_dE=dEin;}
  void set_Emax(double E) {_Emax=E;}
  void set_nParticles(unsigned int n);
  void set_gammaMode(GammaMode g) {_gammaMode=g;}
  void set_trajectoryMode(TrajectoryMode t) {_trajectoryMode=t;}
  void set_edgefoc(bool e) {_edgefoc = e;}
  void set_saveGamma(std::string particleList) {set_saveList(particleList,_saveGamma,"saveGamma");}
  void set_seed(int s) {_seed=s;}
  void set_q(double q) {_q=q;}
  void set_alphac(double ac) {_alphac = ac;}
  void set_alphac2(double ac2) {_alphac2 = ac2;}
  void set_h(unsigned int h) {_h=h;}
  void set_R(double R) {_R=R;}
  void set_Js(double Js) {_Js=Js;}
  void set_simToolRamp(bool r) {_simToolRamp=r;}
  void set_simToolRampSteps(unsigned int n) {_simToolRampSteps=n;}
  void set_savePhaseSpaceElement(std::string name) {_savePhaseSpaceElement=name;}
  void set_savePhaseSpace(std::string particleList) {set_saveList(particleList,_savePhaseSpace,"savePhaseSpace");}
  void set_sigmaPhaseFactor(double sP) {_sigmaPhaseFactor = sP;}
  void set_sigmaGammaFactor(double sG) {_sigmaGammaFactor = sG;}
  void set_agammaMin(double a) {_agammaMin = a;}
  void set_agammaMax(double a) {_agammaMax = a;}
  void set_dagamma(double a) {_dagamma = a;}
  void set_nTurns(unsigned int n) {_nTurns = n;}
  // set parameters, which are currently unset, from SimToolInstance and given lattice
  void autocomplete(const pal::AccLattice& lattice);

protected:
  std::string getSaveList(std::vector<bool> list) const;
  int randomSeed() const;
  
public:
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
  std::string saveGammaList() const {return getSaveList(_saveGamma);}
  std::string savePhaseSpaceList() const {return getSaveList(_savePhaseSpace);}

  void printSummary() const;

  double gamma(double t) const;
  double agamma(double t) const;

  // save & load configuration to/from file
  void save(const std::string &filename) const;
  void load(const std::string &filename);

  // write energy and, if needed, number of turns to SimToolInstance
  void updateSimToolSettings(const pal::AccLattice& lattice);

  void writeRfMagnetsToLattice(pal::AccLattice& lattice) const {rf.writeToLattice(lattice);}

};







template<typename T>
std::string RfMagnetConfig::getStringList(const std::vector<T> v) const
{
  std::stringstream ss;
  for (auto i=0u; i<v.size(); i++) {
    ss << v[i];
    if (i != v.size()-1)
      ss << ",";
  }
  return ss.str();
}


template<typename T>
void RfMagnetConfig::setFromStringList(std::vector<T>& v, std::string list)
{
  std::istringstream ss(list);
  T tmp;
   while ( (ss >> tmp) ) {
     v.push_back(tmp);
     if (ss.peek() == ',' || ss.peek() == ' ')
       ss.ignore();
   }
}

template<> void RfMagnetConfig::setFromStringList<std::string>(std::vector<std::string>& v, std::string list);


#endif
// __POLEMATRIX__CONFIGURATION_HPP_
