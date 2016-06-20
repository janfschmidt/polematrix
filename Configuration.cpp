#include <chrono>
#include <boost/version.hpp>
#include "Configuration.hpp"

Configuration::Configuration(std::string pathIn)
  : E_rest_GeV(0.0005109989), E_rest_keV(E_rest_GeV*1e6), a_gyro(0.001159652), default_steps(1000), spinDirName("spins"), polFileName("polarization.dat"), confOutFileName("currentconfig.pole")
{
  _outpath = pathIn;
  _nParticles = 1;
  _saveGamma.assign(1, false);
  _savePhaseSpace.assign(1, false);

  _t_start = _t_stop = 0.;
  _dt_out = 1e-5;
  _E0 = 1;
  _dE = 0;
  _s_start.zeros();
  _s_start[2] = 1;
  _gammaMode = GammaMode::linear;
  _trajectoryMode = TrajectoryMode::closed_orbit;
  
  _seed = randomSeed();
  _q = 0.;
  _alphac = _alphac2 = 0.;
  _h = 0;
  _R = 0.;
  _Js = 0.;
  _savePhaseSpaceElement = "";

  palattice.reset(new pal::SimToolInstance(pal::elegant, pal::online, ""));
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
  tree.put("radiation.momentum_compaction_factor_2", alphac2());
  tree.put("radiation.harmonic_number", h());
  tree.put("radiation.bending_radius", R());
  tree.put("radiation.longitudinal_damping_partition_number", Js());
  tree.put("radiation.savePhaseSpace", savePhaseSpaceList());
  tree.put("radiation.savePhaseSpaceElement", savePhaseSpaceElement());
  
  if (_gammaMode==GammaMode::linear) tree.put("spintracking.gammaMode", "linear");
  else if (_gammaMode==GammaMode::simtool) tree.put("spintracking.gammaMode", "simtool");
  else if (_gammaMode==GammaMode::simtool_plus_linear) tree.put("spintracking.gammaMode", "simtool+linear");
    else if (_gammaMode==GammaMode::simtool_no_interpolation) tree.put("spintracking.gammaMode", "simtool_no_interpolation");
  else if (_gammaMode==GammaMode::radiation) tree.put("spintracking.gammaMode", "radiation");

  if (_trajectoryMode==TrajectoryMode::closed_orbit) tree.put("spintracking.trajectoryMode", "closed orbit");
  else if (_trajectoryMode==TrajectoryMode::simtool) tree.put("spintracking.trajectoryMode", "simtool");

  #if BOOST_VERSION < 105600
  pt::xml_writer_settings<char> settings(' ', 2); //indentation
  #else
  pt::xml_writer_settings<std::string> settings(' ', 2); //indentation
  #endif
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
    
    //optional, but throws if invalid value
    setGammaMode(tree);
    setTrajectoryMode(tree);
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
  set_alphac( tree.get("radiation.momentum_compaction_factor", 0.0) );
  set_alphac2( tree.get("radiation.momentum_compaction_factor_2", 0.0) );
  set_q( tree.get("radiation.overvoltage_factor", 0.0) );
  set_h( tree.get("radiation.harmonic_number", 0) );
  set_R( tree.get("radiation.bending_radius", 0.0) );
  set_Js( tree.get("radiation.longitudinal_damping_partition_number", 0.0) );
  set_savePhaseSpace( tree.get<std::string>("radiation.savePhaseSpace", "") );
  set_savePhaseSpaceElement( tree.get<std::string>("radiation.savePhaseSpaceElement", "") );
  
  std::cout << "* configuration loaded from " << filename << std::endl;
  return;
}


void Configuration::set_nParticles(unsigned int n)
{
  _nParticles=n;
  _saveGamma.assign(_nParticles, false);
  _savePhaseSpace.assign(_nParticles, false);
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
  if(_gammaMode == GammaMode::simtool) {
    s << "energy ";
    if(palattice->tool==pal::SimTool::elegant)
      s << E0()*1000. <<" MeV";
    s << " from " << palattice->tool_string() << std::endl;
  }
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
    _gammaMode = GammaMode::linear;
  else if (s == "simtool")
    _gammaMode = GammaMode::simtool;
  else if (s == "simtool+linear")
    _gammaMode = GammaMode::simtool_plus_linear;
    else if (s == "simtool_no_interpolation")
    _gammaMode = GammaMode::simtool_no_interpolation;
  else if (s == "radiation")
    _gammaMode = GammaMode::radiation;
  
  else
    throw pt::ptree_error("Invalid gammaMode "+s);
}

void Configuration::setTrajectoryMode(pt::ptree &tree)
{
  std::string s = tree.get<std::string>("spintracking.trajectoryMode", "closed orbit");
  
  if (s == "closed orbit")
    _trajectoryMode = TrajectoryMode::closed_orbit;
  else if (s == "simtool")
    _trajectoryMode = TrajectoryMode::simtool;
  
  else
    throw pt::ptree_error("Invalid trajectoryMode "+s);
}


// parse particleIds from comma separated string
// also ranges (e.g. 0-99) can be parsed
void Configuration::set_saveList(const std::string &particleList, std::vector<bool> &list, const std::string &optionName)
{
 std::stringstream ss(particleList);
  unsigned int tmp;
  while ( (ss >> tmp) ) {
    try {
      list.at(tmp) = true;
    }
    catch (std::out_of_range) {
      std::cout << "* ignore "<<optionName<<" config option of particle ID " << tmp
		<< ", which is out of range" << std::endl;
    }
    if (ss.peek() == ',')
      ss.ignore();

    // parse range
    if (ss.peek() == '-') {
      ss.ignore();
      unsigned int max;
      ss >> max;
      if (ss.peek() == ',')
	ss.ignore();
      for (; tmp<=max; tmp++) {
	try {
	  list.at(tmp) = true;
	}
	catch (std::out_of_range) {
	  std::cout << "* ignore "<<optionName<<" config option of particle IDs " << tmp <<"-"<<max
		    << ", which are out of range" << std::endl;
	  break;
	}
      }
    }
    
  }
}


std::string Configuration::getSaveList(std::vector<bool> list) const
{
  std::stringstream s;
  for (auto i=0u; i<list.size(); i++) {
    if (list[i]) {
      if (i==0u || !list[i-1] || !list[i+1]) {
	s << i;
	if (i<list.size()-1) {
	  if (list[i+1])
	    s << "-";
	  else
	    s << ",";
	}
      }
    }
  }
  return s.str();
}

int Configuration::randomSeed() const
{
  auto now = std::chrono::system_clock::now();
  return std::move( now.time_since_epoch().count() );
}
