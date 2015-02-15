#ifndef __POLEMATRIX__CONFIGURATION_HPP_
#define __POLEMATRIX__CONFIGURATION_HPP_

#include <string>
#include <armadillo>

class Configuration
{
private:
  std::string p; //path

public:
  double x;
  arma::colvec3 s_start;
  unsigned int dt_out;

  Configuration(std::string path="~");
  ~Configuration() {}

  std::string path() const {return p;}
  void setPath(std::string path);
  std::string subfolder(std::string folder) const {return p + folder + "/";}
};


#endif
// __POLEMATRIX__CONFIGURATION_HPP_
