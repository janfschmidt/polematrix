#include "Configuration.hpp"

Configuration::Configuration(std::string pathIn) : E_rest(0.000511), a_gyro(0.001159652)//E_rest(GSL_CONST_MKSA_MASS_ELECTRON*pow(GSL_CONST_MKSA_SPEED_OF_LIGHT,2)/GSL_CONST_MKSA_ELECTRON_CHARGE/1e9)
{
  setPath(pathIn);
}


void Configuration::setPath(std::string path)
{
  if (path.back() != '/')
    path += "/";
  this->p = path;
}

double Configuration::gamma(double t) const
{
  return (E_start + dE * t) / E_rest;
}

double Configuration::agamma(double t) const
{
  return a_gyro * gamma(t);
}
