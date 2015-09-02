#include "Configuration.hpp"

Configuration::Configuration(std::string pathIn) : default_steps(200), E_rest(0.000511), a_gyro(0.001159652)//E_rest(GSL_CONST_MKSA_MASS_ELECTRON*pow(GSL_CONST_MKSA_SPEED_OF_LIGHT,2)/GSL_CONST_MKSA_ELECTRON_CHARGE/1e9)
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


void Configuration::save(const std::string &filename) const
{
  pt::ptree tree;
  tree.put("config.path", path());
  tree.put("config.t_start", t_start);
  tree.put("config.t_stop", t_stop);
  tree.put("config.dt_out", dt_out);
  tree.put("config.E_start", E_start);
  tree.put("config.dE", dE);
  tree.put("config.s_start.x", s_start[0]);
  tree.put("config.s_start.z", s_start[2]);
  tree.put("config.s_start.s", s_start[1]);
  pt::xml_writer_settings<char> settings(' ', 4); //indentation
  pt::write_xml(filename, tree, std::locale(), settings);
}

void Configuration::load(const std::string &filename)
{
  pt::ptree tree;
  pt::read_xml(filename, tree);

  //obligatory config
  try {
    setPath( tree.get<std::string>("config.path") );
    t_start = tree.get<double>("config.t_start");
    t_stop = tree.get<double>("config.t_stop"); 
    E_start = tree.get<double>("config.E_start");
    dE = tree.get<double>("config.dE", dE);
  } catch (pt::ptree_error &e) {
    std::cout << "Error loading configuration file:" << std::endl
	      << e.what() << std::endl;
    exit(1);
  }
    
  // optional config with default values
  dt_out = tree.get("config.dt_out", (t_stop-t_start)/default_steps);
  s_start[0] = tree.get("config.s_start.x", 0.0);
  s_start[2] = tree.get("config.s_start.z", 1.0);
  s_start[1] = tree.get("config.s_start.s", 0.0);
}

