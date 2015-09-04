#include "Configuration.hpp"

Configuration::Configuration(std::string pathIn) : E_rest(0.000511), a_gyro(0.001159652), default_steps(200), spinDirName("spins"), polFileName("polarization.dat")//E_rest(GSL_CONST_MKSA_MASS_ELECTRON*pow(GSL_CONST_MKSA_SPEED_OF_LIGHT,2)/GSL_CONST_MKSA_ELECTRON_CHARGE/1e9)
{
  //setPath(pathIn);
  outpath = pathIn;
  nParticles = 1;
}



double Configuration::gamma(double t) const
{
  return (E0 + dE * t) / E_rest;
}

double Configuration::agamma(double t) const
{
  return a_gyro * gamma(t);
}


void Configuration::save(const std::string &filename) const
{
  pt::ptree tree;
  tree.put("config.numParticles", nParticles);
  tree.put("config.t_start", t_start);
  tree.put("config.t_stop", t_stop);
  tree.put("config.dt_out", dt_out);
  tree.put("config.E0", E0);
  tree.put("config.dE", dE);
  tree.put("config.s_start.x", s_start[0]);
  tree.put("config.s_start.z", s_start[2]);
  tree.put("config.s_start.s", s_start[1]);
  tree.put("config.simTool.file", simFile);
  pt::xml_writer_settings<char> settings(' ', 4); //indentation
  pt::write_xml(filename, tree, std::locale(), settings);

  std::cout << "* current configuration saved to " << filename << std::endl;
  return;
}

void Configuration::load(const std::string &filename)
{
  pt::ptree tree;
  pt::read_xml(filename, tree);

  //obligatory config
  try {
    t_stop = tree.get<double>("config.t_stop"); 
    E0 = tree.get<double>("config.E0");
    dE = tree.get<double>("config.dE");
    s_start[0] = tree.get<double>("config.s_start.x");
    s_start[2] = tree.get<double>("config.s_start.z");
    s_start[1] = tree.get<double>("config.s_start.s");
    simFile = tree.get<std::string>("config.simTool.file");
  } catch (pt::ptree_error &e) {
    std::cout << "Error loading configuration file:" << std::endl
	      << e.what() << std::endl;
    exit(1);
  }
    
  // optional config with default values
  nParticles = tree.get("config.numParticles", 1);
  t_start = tree.get("config.t_start", 0.0);
  dt_out = tree.get("config.dt_out", (t_stop-t_start)/default_steps);

  std::cout << "* configuration loaded from " << filename << std::endl;
  return;
}

