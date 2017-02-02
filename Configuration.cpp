/* Configuration Class
 * stores program configuration and can be loaded from & saved to xml file
 *
 * Copyright (C) 2017 Jan Felix Schmidt <janschmidt@mailbox.org>
 *   
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <chrono>
#include <boost/version.hpp>
#include "Configuration.hpp"
#include "debug.hpp"

Configuration::Configuration(std::string pathIn)
  : E_rest_GeV(0.0005109989), E_rest_keV(E_rest_GeV*1e6), a_gyro(0.001159652), default_steps(1000), spinDirName("spins"), polFileName("polarization.dat"), confOutFileName("currentconfig.pole")
{
  _outpath = pathIn;
  _verbose = false;
  _nParticles = 1;
  _saveGamma.assign(1, false);
  _savePhaseSpace.assign(1, false);

  _t_start = _t_stop = 0.;
  _dt_out = 1e-5;
  _E0 = 1;
  _dE = 0;
  _s_start.zeros();
  _s_start[2] = 1;
  _gammaMode = GammaMode::radiation;
  _trajectoryMode = TrajectoryMode::closed_orbit;
  _edgefoc = false;
  
  _seed = randomSeed();
  _q = 0.;
  _alphac = _alphac2 = 0.;
  _h = 0;
  _R = 0.;
  _Js = 0.;
  _savePhaseSpaceElement = "";
  _sigmaPhaseFactor = 1.;
  _sigmaGammaFactor = 1.;
  _agammaMin = 0.;
  _agammaMax = 10.;
  _nTurns = 0;
  _dagamma = 0.;

  palattice.reset(new pal::SimToolInstance(pal::elegant, pal::online, ""));

  // metadata
  info.add("and polematrix version", polemversion());
}

std::string Configuration::gammaModeString() const
{
  if (_gammaMode==GammaMode::linear) return "linear";
  else if (_gammaMode==GammaMode::offset) return "offset";
  else if (_gammaMode==GammaMode::oscillation) return "oscillation";
  else if (_gammaMode==GammaMode::simtool) return "simtool";
  else if (_gammaMode==GammaMode::simtool_plus_linear) return "simtool+linear";
  else if (_gammaMode==GammaMode::simtool_no_interpolation) return "simtool_no_interpolation";
  else if (_gammaMode==GammaMode::radiation) return "radiation";
  else
    return "Please implement this GammaMode in Configuration::gammaModeString()!";
}

std::string Configuration::trajectoryModeString() const
{
  if (_trajectoryMode==TrajectoryMode::closed_orbit) return "closed orbit";
  else if (_trajectoryMode==TrajectoryMode::simtool) return "simtool";
  else
    return "Please implement this TrajectoryMode in Configuration::trajectoryModeString()!";
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
  tree.put("spintracking.edgeFocussing", _edgefoc);
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
  tree.put("radiation.startDistribution.sigmaPhaseFactor", sigmaPhaseFactor());
  tree.put("radiation.startDistribution.sigmaGammaFactor", sigmaGammaFactor());
  tree.put("resonancestrengths.minSpintune", agammaMin());
  tree.put("resonancestrengths.maxSpintune", agammaMax());
  tree.put("resonancestrengths.spintuneStep", _dagamma);
  tree.put("resonancestrengths.turns", _nTurns);
  
  tree.put("spintracking.gammaMode", gammaModeString());
  tree.put("spintracking.trajectoryMode", trajectoryModeString());

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
  set_edgefoc( tree.get<bool>("spintracking.edgeFocussing", false) );
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
  set_sigmaPhaseFactor( tree.get<double>("radiation.startDistribution.sigmaPhaseFactor", 1.0) );
  set_sigmaGammaFactor( tree.get<double>("radiation.startDistribution.sigmaGammaFactor", 1.0) );
  set_agammaMin( tree.get<double>("resonancestrengths.minSpintune", 0.) );
  set_agammaMax( tree.get<double>("resonancestrengths.maxSpintune", 10.) );
  set_nTurns( tree.get<unsigned int>("resonancestrengths.turns", 0) );
  set_dagamma( tree.get<double>("resonancestrengths.spintuneStep", 0.) );

  set_metadata(filename);
  
  std::cout << "* configuration loaded from " << filename << std::endl;
  return;
}

void Configuration::set_metadata(const std::string& configfile)
{
  info.add("configuration file", configfile);
  info.add("gammaMode", gammaModeString());
  info.add("trajectoryMode", trajectoryModeString());
  info.add("simtool", palattice->tool_string());
  info.add("simtool file", palattice->inFile());
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
  s << "longitudinal phase space model (GammaMode): \"" << gammaModeString() << "\"" << std::endl;
  if (edgefoc())
    s << "horizontal dipole edge focussing field used" << std::endl;
  s << "output for each spin vector to " << spinDirectory().string() <<"/"<< std::endl;
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
    polematrix::debug(__PRETTY_FUNCTION__, "no changes");
    return;
  }

  palattice.reset(new pal::SimToolInstance(tool, mode, file));
}

void Configuration::setGammaMode(pt::ptree &tree)
{
  std::string s = tree.get<std::string>("spintracking.gammaMode", "radiation");
  
  if (s == "linear")
    _gammaMode = GammaMode::linear;
  else if (s == "offset")
    _gammaMode = GammaMode::offset;
    else if (s == "oscillation")
    _gammaMode = GammaMode::oscillation;
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



// set parameters, which are currently unset, from SimToolInstance and given lattice
void Configuration::autocomplete(const pal::AccLattice& lattice)
{
  if (q()==0.) {
    set_q(lattice.overvoltageFactor(gamma_start()));
    std::cout << "* set overvoltage factor from lattice"
	      << ": q=" << q() << std::endl;
  }
  if (h()==0) {
    set_h(lattice.harmonicNumber());
    std::cout << "* set harmonic number from lattice"
	      << ": h=" << h() << std::endl;
  }
  if (R()==0.) {
    set_R(lattice.integralDipoleRadius());
    std::cout << "* set dipole bending radius from lattice"
	      << ": R=" << R() << std::endl;
  }
  if (alphac()==0.) {
    set_alphac(palattice->readAlphaC());
    std::cout << "* set momentum compaction factor from " << palattice->tool_string()
	      << ": alphac=" << alphac() << std::endl;
  }
  if (alphac2()==0.) {
    set_alphac2(palattice->readAlphaC2());
    std::cout << "* set 2nd order momentum compaction factor from " << palattice->tool_string()
	      << ": alphac2=" << alphac2() << std::endl;
  }
  if (Js()==0.) {
    set_Js(palattice->readDampingPartitionNumber_syli().s);
    std::cout << "* set long. damping partition number from " << palattice->tool_string()
	      << ": Js=" << Js() << std::endl;
  }

  _circumference = lattice.circumference();
}


// write energy and, if needed, number of turns to SimToolInstance
void Configuration::updateSimToolSettings(const pal::AccLattice& lattice)
{
  //set energy to E0 (ramp not considered!)
  double p_MeV = E0()*1000.;
  palattice->setMomentum_MeV(p_MeV);

  
  //set number of turns based on tracking time & circumference
  if (gammaMode() == GammaMode::simtool
      || gammaMode()==GammaMode::simtool_plus_linear
      || gammaMode()==GammaMode::simtool_no_interpolation
      || trajectoryMode() == TrajectoryMode::simtool) {
    unsigned int turns = (duration()*GSL_CONST_MKSA_SPEED_OF_LIGHT / lattice.circumference()) + 1;
    palattice->verbose = true;
    palattice->setTurns(turns);
    std::cout << "* " << palattice->tool_string() <<" tracking " << turns <<" turns to get single particle trajectories" << std::endl;
  }

}

// #turns for resonance strengths calc are calculated from tracking duration() if not set
unsigned int Configuration::numTurns() const
{
  if (_nTurns != 0)
    return _nTurns;
  else {
    if (std::fabs(_circumference) < COMPARE_DOUBLE_EQUAL)
      throw std::runtime_error("Configuration::numTurns() cannot guess number of turns, because it's neither set in config file nor loaded from model (Configuration::autocomlete()");
    return (duration()*GSL_CONST_MKSA_SPEED_OF_LIGHT / _circumference) + 1;
  }
}

// spin tune output step width for resonance strengths is calculated from #turns if not set
double Configuration::dagamma() const
{
  auto dag = std::fabs(_dagamma);
  if ( dag > COMPARE_DOUBLE_EQUAL )
    return dag;
  else
    return 1./double(numTurns());
}
