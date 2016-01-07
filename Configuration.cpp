#include <chrono>
#include "Configuration.hpp"

Configuration::Configuration(std::string pathIn)
  : E_rest_GeV(0.000511), E_rest_keV(E_rest_GeV*1e6), a_gyro(0.001159652), default_steps(1000), spinDirName("spins"), polFileName("polarization.dat"), confOutFileName("currentconfig.pole")
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
  
  _seed = randomSeed();
  _q = 0.;
  _alphac = 0.;
  _h = 0.;

  palattice.reset(new pal::SimToolInstance(pal::madx, pal::offline, ""));
}


double Configuration::gamma(double t) const
{
  return (E0() + dE() * t) / E_rest_GeV;
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
  tree.put("radiation.seed", seed());
  tree.put("radiation.overvoltage_factor", q());
  tree.put("radiation.momentum_compaction_factor", alphac());
  tree.put("radiation.harmonic_number", h());
  
  if (_gammaMode==linear) tree.put("spintracking.gammaMode", "linear");
  else if (_gammaMode==simtool) tree.put("spintracking.gammaMode", "simtool");
  else if (_gammaMode==simtool_plus_linear) tree.put("spintracking.gammaMode", "simtool+linear");
  else if (_gammaMode==radiation) tree.put("spintracking.gammaMode", "radiation");


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
    set_t_stop( tree.get<double>("spintracking.t_stop") );
    set_E0( tree.get<double>("spintracking.E0") );
    set_dE( tree.get<double>("spintracking.dE") );
    arma::colvec3 tmp;
    tmp[0] = tree.get<double>("spintracking.s_start.x");
    tmp[2] = tree.get<double>("spintracking.s_start.z");
    tmp[1] = tree.get<double>("spintracking.s_start.s");
    set_s_start(tmp);
    
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
  set_t_start( tree.get("spintracking.t_start", 0.0) );
  set_dt_out( tree.get("spintracking.dt_out", duration()/default_steps) );
  set_saveGamma( tree.get<std::string>("palattice.saveGamma", "") );
  set_seed( tree.get<int>("radiation.seed", randomSeed()) );
  set_alphac( tree.get("radiation.momentum_compaction_factor", getSimToolInstance().readAlphaC()) );
  set_q( tree.get("radiation.overvoltage_factor", 10) );
  set_h( tree.get("radiation.harmonic_number", 274) );
  
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
  unsigned int w=8;

  s << "-----------------------------------------------------------------" << std::endl;
  s << "Tracking " << _nParticles << " Spins" << std::endl
    << "time      " <<std::setw(w-2)<<  _t_start << " s   -------------------->   " <<std::setw(w-2)<< _t_stop << " s" << std::endl;
  if(_gammaMode == simtool)
    s << "energy from " << palattice->tool_string() << std::endl;
  else
    s << "energy    " <<std::setw(w-4)<< gamma_start()*E_rest_GeV << " GeV   ----- " <<std::setw(3)<< _dE << " GeV/s ---->   " <<std::setw(w-4)<< gamma_stop()*E_rest_GeV << " GeV" << std::endl
      << "spin tune " <<std::setw(w)<< agamma_start() << "   -------------------->   " <<std::setw(w)<< agamma_stop() << std::endl;
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
  else if (s == "radiation")
    _gammaMode = radiation;
  
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

int Configuration::randomSeed() const
{
  auto now = std::chrono::system_clock::now();
  return std::move( now.time_since_epoch().count() );
}
