#include "Configuration.hpp"

Configuration::Configuration(std::string pathIn)
  : E_rest(0.000511), a_gyro(0.001159652), default_steps(200), spinDirName("spins"), polFileName("polarization.dat"), confOutFileName("currentconfig.pole")
{
  //setPath(pathIn);
  _outpath = pathIn;
  _nParticles = 1;

  _t_start = _t_stop = 0.;
  _dt_out = 1e-5;
  _E0 = 1;
  _dE = 0;
  _s_start.zeros();
  _s_start[2] = 1;
  _gammaMode = linear;

  palattice.reset(new pal::SimToolInstance(pal::madx, pal::offline, ""));  
}


double Configuration::gamma(double t) const
{
  return (E0() + dE() * t) / E_rest;
}

double Configuration::agamma(double t) const
{
  return a_gyro * gamma(t);
}



void Configuration::save(const std::string &filename) const
{
  pt::ptree tree;
  tree.put("spintracking.numParticles", _nParticles);
  tree.put("spintracking.t_start", _t_start);
  tree.put("spintracking.t_stop", _t_stop);
  tree.put("spintracking.dt_out", _dt_out);
  tree.put("spintracking.E0", _E0);
  tree.put("spintracking.dE", _dE);
  tree.put("spintracking.s_start.x", _s_start[0]);
  tree.put("spintracking.s_start.z", _s_start[2]);
  tree.put("spintracking.s_start.s", _s_start[1]);
  tree.put("palattice.simTool", palattice->tool_string());
  tree.put("palattice.mode", palattice->mode_string());
  tree.put("palattice.file", palattice->inFile());

  if (_gammaMode==linear) tree.put("spintracking.gammaMode", "linear");
  else if (_gammaMode==simtool) tree.put("spintracking.gammaMode", "simtool");

  pt::xml_writer_settings<std::string> settings(' ', 2); //indentation
  pt::write_xml(filename, tree, std::locale(), settings);

  std::cout << "* current configuration saved in " << filename << std::endl;
  return;
}




void Configuration::load(const std::string &filename)
{
  pt::ptree tree;
  pt::read_xml(filename, tree);

  //obligatory config
  try {
    _t_stop = tree.get<double>("spintracking.t_stop"); 
    _E0 = tree.get<double>("spintracking.E0");
    _dE = tree.get<double>("spintracking.dE");
    _s_start[0] = tree.get<double>("spintracking.s_start.x");
    _s_start[2] = tree.get<double>("spintracking.s_start.z");
    _s_start[1] = tree.get<double>("spintracking.s_start.s");
    
    setSimToolInstance(tree);
    setGammaMode(tree); //optional, but throws if invalid value
  }
  catch (pt::ptree_error &e) {
    std::cout << "Error loading configuration file:" << std::endl
	      << e.what() << std::endl;
    exit(1);
  }

  
  // optional config with default values
  _nParticles = tree.get("spintracking.numParticles", 1);
  try {
    palattice->setNumParticles(_nParticles);
  }
  catch (pal::palatticeError &e) {
    std::cout << "ignoring numParticles for madx tracking: " <<std::endl << e.what() << std::endl;
  }
  
  _t_start = tree.get("spintracking.t_start", 0.0);
  _dt_out = tree.get("spintracking.dt_out", duration()/default_steps);
  
  std::cout << "* configuration loaded from " << filename << std::endl;
  return;
}




void Configuration::printSummary() const
{
  std::stringstream s;

  s << "-----------------------------------------------------------------" << std::endl;
  s << "Tracking " << _nParticles << " Spins" << std::endl
    << "time      " <<  _t_start << " s   -------------------->   " << _t_stop << " s" << std::endl
    << "energy    " << gamma(_t_start)*E_rest << " GeV   ----- " << _dE << " GeV/s ----->   " << gamma(_t_stop)*E_rest << " GeV" << std::endl
    << "spin tune " << agamma(_t_start) << "   -------------------->   " << agamma(_t_stop) << std::endl;
  s << "start polarization: Px = " << _s_start[0] << ", Ps = " << _s_start[1] << ", Pz = " << _s_start[2] << std::endl;
  s << "-----------------------------------------------------------------" << std::endl;

  std::cout << s.str();
}



pal::SimTool Configuration::toolFromTree(pt::ptree tree, std::string key) const
{
  std::string tmp = tree.get<std::string>(key);
  if (tmp == "madx")
    return pal::madx;
  else if (tmp == "elegant")
    return pal::elegant;
  else {
    throw pt::ptree_error("Invalid pal::SimTool "+tmp);
  }
}

pal::SimToolMode Configuration::modeFromTree(pt::ptree tree, std::string key) const
{
  std::string tmp = tree.get<std::string>(key);
  if (tmp == "online")
    return pal::online;
  else if (tmp == "offline")
    return pal::offline;
  else
    throw pt::ptree_error("Invalid pal::SimToolMode "+tmp);
}

void Configuration::setSimToolInstance(pt::ptree &tree)
{
  pal::SimTool tool = toolFromTree(tree, "palattice.simTool");
  pal::SimToolMode mode = modeFromTree(tree, "palattice.mode");
  std::string file = tree.get<std::string>("palattice.file");

  //do not change SimToolInstance if nothing has to be changed
  if (tool==palattice->tool && mode==palattice->mode && file==palattice->inFile()) {
    std::cout << "DEBUG: Configuration::setSimToolInstance(): no changes" << std::endl;
    return;
  }

  // delete palattice;
  // palattice = new pal::SimToolInstance(tool, mode, file);
  palattice.reset(new pal::SimToolInstance(tool, mode, file));
 

}

void Configuration::setGammaMode(pt::ptree &tree)
{
  std::string s = tree.get<std::string>("spintracking.gammaMode", "linear");
  
  if (s == "linear")
    _gammaMode = linear;
  else if (s == "simtool")
    _gammaMode = simtool;
  else
    throw pt::ptree_error("Invalid gammaMode "+s);
}
