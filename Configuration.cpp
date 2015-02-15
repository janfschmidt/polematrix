#include "Configuration.hpp"

Configuration::Configuration(std::string pathIn)
{
  setPath(pathIn);
}


void Configuration::setPath(std::string path)
{
  if (path.back() != '/')
    path += "/";
  this->p = path;
}
