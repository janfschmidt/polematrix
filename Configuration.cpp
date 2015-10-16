#include "Configuration.hpp"

Configuration::Configuration(std::string pathIn)
  : E_rest(0.000511), a_gyro(0.001159652), default_steps(200), spinDirName("spins"), polFileName("polarization.dat"), confOutFileName("currentconfig.pole")
{
  _outpath = pathIn;
  _nParticles = 1;
  _saveGamma.assign(1, false);

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
  tree.put("palattice.saveGamma", saveGammaList());

  if (_gammaMode==linear) tree.put("spintracking.gammaMode", "linear");
  else if (_gammaMode==simtool) tree.put("spintracking.gammaMode", "simtool");
  else if (_gammaMode==simtool_plus_linear) tree.put("spintracking.gammaMode", "simtool+linear");

  pt::xml_writer_settings<char> settings(' ', 2); //indentation
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
    _s_start = arma::normalise(_s_start);
    
    setSimToolInstance(tree);
    setGammaMode(tree); //optional, but fails if invalid value
  }
  catch (pt::ptree_error &e) {
    std::cout << "Error loading configuration file:" << std::endl
	      << e.what() << std::endl;
    exit(1);
  }

  
  // optional config with default values
  set_nParticles( tree.get("spintracking.numParticles", 1) );
  _t_start = tree.get("spintracking.t_start", 0.0);
  _dt_out = tree.get("spintracking.dt_out", duration()/default_steps);
  set_saveGamma( tree.get<std::string>("palattice.saveGamma", "") );
  
  std::cout << "* configuration loaded from " << filename << std::endl;
  return;
}


void Configuration::set_nParticles(unsigned int n)
{
  _nParticles=n;
  _saveGamma.assign(_nParticles, false);
  try {
    palattice->setNumParticles(_nParticles);
  }
  catch (pal::palatticeError &e) {
    std::cout << "ignoring numParticles for madx tracking: " <<std::endl << e.what() << std::endl;
  }
}




void Configuration::printSummary() const
{
  std::stringstream s;

  s << "-----------------------------------------------------------------" << std::endl;
  s << "Tracking " << _nParticles << " Spins" << std::endl
    << "time      " <<  _t_start << " s   -------------------->   " << _t_stop << " s" << std::endl
    << "energy    " << gamma(_t_start)*E_rest << " GeV   ----- " << _dE << " GeV/s ----->   " << gamma(_t_stop)*E_rest << " GeV" << std::endl
    << "spin tune " << agamma(_t_start) << "   -------------------->   " << agamma(_t_stop) << std::endl;
  s << "start spin direction: Sx = " << _s_start[0] << ", Ss = " << _s_start[1] << ", Sz = " << _s_start[2] << std::endl;
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

  palattice.reset(new pal::SimToolInstance(tool, mode, file));
}

void Configuration::setGammaMode(pt::ptree &tree)
{
  std::string s = tree.get<std::string>("spintracking.gammaMode", "linear");
  
  if (s == "linear")
    _gammaMode = linear;
  else if (s == "simtool")
    _gammaMode = simtool;
  else if (s == "simtool+linear")
    _gammaMode = simtool_plus_linear;
  
  else
    throw pt::ptree_error("Invalid gammaMode "+s);
}


void Configuration::set_saveGamma(std::string particleList)
{
 std::stringstream ss(particleList);
  unsigned int tmp;
  while ( (ss >> tmp) ) {
    try {
    _saveGamma.at(tmp) = true;
    }
    catch (std::out_of_range) {
      std::cout << "* ignore saveGamma config option of particle ID " << tmp
		<< ", which is out of range" << std::endl;
    }
      
    if (ss.peek() == ',')
      ss.ignore();
  }
}

std::string Configuration::saveGammaList() const
{
  std::stringstream s;
  for (auto i=0u; i<_saveGamma.size(); i++) {
    if (_saveGamma[i]) s << i << ",";
  }
  return s.str();
}
