#ifndef __POLEMATRIX__CONFIGURATION_HPP_
#define __POLEMATRIX__CONFIGURATION_HPP_

#include <string>
#include <armadillo>
#include <gsl/gsl_const_mksa.h>
#include <fstream>
// #include <boost/serialization/base_object.hpp>
// #include <boost/archive/xml_oarchive.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace pt = boost::property_tree;

class Configuration
{
private:
  std::string p; //path
  const unsigned int default_steps;

  
public:
  arma::colvec3 s_start;
  double t_start;         // time / s
  double t_stop;
  double dt_out;
  double E_start;         // energy / GeV
  double dE;              // dE/dt / GeV/s

  //constants
  const double E_rest;    // electron rest energy / GeV
  const double a_gyro;    // electron gyromagnetic anomaly a = (g-2)/2

  Configuration(std::string path="~");
  ~Configuration() {}

  std::string path() const {return p;}
  void setPath(std::string path);
  std::string subfolder(std::string folder) const {return p + folder + "/";}

  double pos_start() const {return GSL_CONST_MKSA_SPEED_OF_LIGHT * t_start;}
  double pos_stop() const {return GSL_CONST_MKSA_SPEED_OF_LIGHT * t_stop;}
  double dpos_out() const {return GSL_CONST_MKSA_SPEED_OF_LIGHT * dt_out;}

  double gamma(double t) const;
  double agamma(double t) const;

  // save & load configuration to/from file
  void save(const std::string &filename) const;
  void load(const std::string &filename);
  
};


#endif
// __POLEMATRIX__CONFIGURATION_HPP_
