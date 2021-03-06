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
#include <algorithm>
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
  _simToolRamp = true;
  _simToolRampSteps = 200;

  _t_start = _t_stop = 0.;
  _dt_out = 1e-4;
  _E0 = 1;
  _dE = 0;
  _Emax = 1e10;
  _s_start.zeros();
  _s_start[2] = 1;
  _gammaMode = GammaMode::radiation;
  _trajectoryMode = TrajectoryMode::closed_orbit;
  _edgefoc = false;
  _outElementUsed = false;
  
  _seed = randomSeed();
  _q = 0.;
  _alphac = _alphac2 = 0.;
  _h = 0;
  _R = 0.;
  _Js = 0.;
  _savePhaseSpaceElement = "";
  _sigmaPhaseFactor = 1.;
  _sigmaGammaFactor = 1.;
  _checkStability = true;
  _agammaMin = 0.;
  _agammaMax = 10.;
  _nTurns = 0;
  _dagamma = 1.;

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
    return "Please implement this GammaModel in Configuration::gammaModeString()!";
}

std::string Configuration::trajectoryModeString() const
{
  if (_trajectoryMode==TrajectoryMode::closed_orbit) return "closed orbit";
  else if (_trajectoryMode==TrajectoryMode::simtool) return "simtool";
  else if (_trajectoryMode==TrajectoryMode::oscillation) return "oscillation";
  else
    return "Please implement this TrajectoryModel in Configuration::trajectoryModeString()!";
}

double Configuration::gamma(double t) const
{
  double E = (E0() + dE() * t);
  if (E > Emax())
    E = Emax();
  
  return E / E_rest_GeV;
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
  tree.put("palattice.simToolRamp.set", simToolRamp());
  tree.put("palattice.simToolRamp.steps", simToolRampSteps());
  tree.put("radiation.seed", seed());
  tree.put("radiation.savePhaseSpace.list", savePhaseSpaceList());
  tree.put("radiation.savePhaseSpace.elementName", savePhaseSpaceElement());
  tree.put("radiation.startDistribution.sigmaPhaseFactor", sigmaPhaseFactor());
  tree.put("radiation.startDistribution.sigmaGammaFactor", sigmaGammaFactor());
  tree.put("radiation.checkStability", checkStability());
  tree.put("resonancestrengths.spintune.min", agammaMin());
  tree.put("resonancestrengths.spintune.max", agammaMax());
  tree.put("resonancestrengths.spintune.step", _dagamma);
  tree.put("resonancestrengths.turns", _nTurns);
  tree.put("oscillation.emittance.x", emittance().x);
  tree.put("oscillation.emittance.z", emittance().z);

  tree.put("spintracking.gammaModel", gammaModeString());
  tree.put("spintracking.trajectoryModel", trajectoryModeString());
  rf.writeToConfig(tree);

  // options, which are only saved if not default value
  if (Emax() < 1e10) {
    tree.put("spintracking.Emax", Emax());
  }
  if (outElementUsed()) {
    tree.put("spintracking.outElement", outElement());
  }
  if (tune().x != 0. || tune().z != 0.) {
    tree.put("oscillation.tune.x", tune().x);
    tree.put("oscillation.tune.z", tune().z);
  }
  if (q() != 0. || h() != 0.) {
    tree.put("radiation.overvoltage_factor", q());
    tree.put("radiation.harmonic_number", h());
  }
  if (alphac() != 0. || alphac2() != 0.) {
    tree.put("radiation.momentum_compaction_factor", alphac());
    tree.put("radiation.momentum_compaction_factor_2", alphac2());
  }
  if (R() != 0.) {
    tree.put("radiation.bending_radius", R());
  }
  if (Js() != 0.) {
    tree.put("radiation.longitudinal_damping_partition_number", Js());
  }


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
  try {
    pt::read_xml(filename, tree);
  }
  catch (pt::xml_parser_error &e){
    std::cout << "Error loading configuration file:" << std::endl
	      << e.what() << std::endl;
    exit(1);
  }

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


  //obligatory for trajectoryMode oscillation
  if (trajectoryMode() == TrajectoryMode::oscillation) {
    try {
      set_emittance_x( tree.get<double>("oscillation.emittance.x") );
      set_emittance_z( tree.get<double>("oscillation.emittance.z") );
    }
    catch (pt::ptree_error &e) {
      std::cout << "Error loading configuration file:" << std::endl
		<< "emittance has to be set for trajectoryModel \"oscillation\"" << std::endl
		<< " --> " << e.what() << std::endl;
      exit(1);
    }
    //do not accept emittance zero
    if (emittance().x==0.0 || emittance().z==0.0) {
      std::cout << "Error loading configuration file:" << std::endl
		<< "emittance > 0 has to be set for trajectoryModel \"oscillation\"" << std::endl;
      exit(1);
    }
  }


  
  // optional config with default values
  set_nParticles( tree.get("spintracking.numParticles", 1) );
  set_t_start( tree.get("spintracking.t_start", 0.0) );
  set_dt_out( tree.get("spintracking.dt_out", duration()/default_steps) );
  set_Emax( tree.get("spintracking.Emax", 1e10) );
  set_edgefoc( tree.get<bool>("spintracking.edgeFocussing", false) );
  set_saveGamma( tree.get<std::string>("palattice.saveGamma", "") );
  set_simToolRamp( tree.get<bool>("palattice.simToolRamp.set", true) );
  set_simToolRampSteps( tree.get<unsigned int>("palattice.simToolRamp.steps", 200) );
  set_seed( tree.get<int>("radiation.seed", randomSeed()) );
  set_alphac( tree.get("radiation.momentum_compaction_factor", 0.0) );
  set_alphac2( tree.get("radiation.momentum_compaction_factor_2", 0.0) );
  set_q( tree.get("radiation.overvoltage_factor", 0.0) );
  set_h( tree.get("radiation.harmonic_number", 0) );
  set_R( tree.get("radiation.bending_radius", 0.0) );
  set_Js( tree.get("radiation.longitudinal_damping_partition_number", 0.0) );
  set_savePhaseSpace( tree.get<std::string>("radiation.savePhaseSpace.list", "") );
  set_savePhaseSpaceElement( tree.get<std::string>("radiation.savePhaseSpace.elementName", "") );
  set_sigmaPhaseFactor( tree.get<double>("radiation.startDistribution.sigmaPhaseFactor", 1.0) );
  set_sigmaGammaFactor( tree.get<double>("radiation.startDistribution.sigmaGammaFactor", 1.0) );
  set_checkStability( tree.get<bool>("radiation.checkStability", true) );
  set_tune_x( tree.get<double>("oscillation.tune.x", 0.0) );
  set_tune_z( tree.get<double>("oscillation.tune.z", 0.0) );
  set_agammaMin( tree.get<double>("resonancestrengths.spintune.min", 0.) );
  set_agammaMax( tree.get<double>("resonancestrengths.spintune.max", 10.) );
  set_dagamma( tree.get<double>("resonancestrengths.spintune.step", 1.) );
  set_nTurns( tree.get<unsigned int>("resonancestrengths.turns", 0) );
  rf.set(tree);

  try {
    set_outElement( tree.get<std::string>("spintracking.outElement") );
  }
  catch (pt::ptree_error &e) {
    _outElementUsed = false;
  }
  
  set_metadata(filename);
  
  std::cout << "* configuration loaded from " << filename << std::endl;
  return;
}

void Configuration::set_metadata(const std::string& configfile)
{
  info.add("configuration file", configfile);
  info.add("gammaModel", gammaModeString());
  info.add("trajectoryModel", trajectoryModeString());
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
  s << "Tracking " << _nParticles << " Spins" << std::endl;
  s << "time      " <<std::setw(w-2)<<  _t_start << " s   -------------------->   " <<std::setw(w-2)<< _t_stop << " s" << std::endl;
  if(_gammaMode == GammaMode::simtool
     && (palattice->tool==pal::SimTool::madx || palattice->mode==pal::offline)) {
    s << "energy from " << palattice->tool_string() << std::endl;
  }
  else {
    s << "energy    " <<std::setw(w-4)<< E_GeV(t_start()) << " GeV   ----- " <<std::setw(3)<< _dE << " GeV/s ---->   " <<std::setw(w-4)<< E_GeV(t_stop()) << " GeV" << std::endl
      << "spin tune " <<std::setw(w)<< agamma_start() << "   -------------------->   " <<std::setw(w)<< agamma_stop() << std::endl;
  }
  if (palattice->mode == pal::online)
    s << "lattice: " << palattice->inFile() << std::endl;
  s << "start spin direction: Sx = " << _s_start[0] << ", Ss = " << _s_start[1] << ", Sz = " << _s_start[2] << std::endl;
  s << "longitudinal phase space model (GammaModel): \"" << gammaModeString() << "\"" << std::endl;
  s << "transversal phase space model (TrajectoryModel): \"" << trajectoryModeString() << "\"" << std::endl;
  if (edgefoc())
    s << "horizontal dipole edge focussing field used" << std::endl;
  s << "output for each spin vector to " << spinDirectory().string() <<"/"<< std::endl;
  if (outElementUsed())
    s << "output at lattice element " << outElement() << " only "<< std::endl;
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
  std::string s;
  // rename Mode->Model, accept "Mode" if "Model" found for compatibility
  try {
    s = tree.get<std::string>("spintracking.gammaModel");
  }
  catch (pt::ptree_error &e) {
    s = tree.get<std::string>("spintracking.gammaMode", "radiation");
    std::cout << "WARNING: option gammaMode is deprecated. It has been renamed as gammaModel." << std::endl;
  }
  
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
    throw pt::ptree_error("Invalid gammaModel "+s);
}

void Configuration::setTrajectoryMode(pt::ptree &tree)
{
  std::string s;
  // rename Mode->Model, accept "Mode" if "Model" found for compatibility
  try {
    s = tree.get<std::string>("spintracking.trajectoryModel");
  }
  catch (pt::ptree_error &e) {
    s = tree.get<std::string>("spintracking.trajectoryMode", "closed orbit");
    std::cout << "WARNING: option trajectoryMode is deprecated. It has been renamed as trajectoryModel." << std::endl;
  }
  
  if (s == "closed orbit")
    _trajectoryMode = TrajectoryMode::closed_orbit;
  else if (s == "simtool")
    _trajectoryMode = TrajectoryMode::simtool;
  else if (s == "oscillation")
      _trajectoryMode = TrajectoryMode::oscillation;
  else
    throw pt::ptree_error("Invalid trajectoryModel "+s);
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
    set_R(lattice.avgDipoleRadius());
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
  if (trajectoryMode() == TrajectoryMode::oscillation) {
    if (tune().x==0.) {
      set_tune_x(palattice->readTune().x);
      std::cout << "* set horizontal tune from " << palattice->tool_string()
		<< ": Qx=" << tune().x << std::endl;
    }
    if (tune().z==0.) {
      set_tune_z(palattice->readTune().z);
      std::cout << "* set vertical tune from " << palattice->tool_string()
		<< ": Qz=" << tune().z << std::endl;
    }
  }

  _circumference = lattice.circumference();
}


// write energy and, if needed, number of turns to SimToolInstance
void Configuration::updateSimToolSettings(const pal::AccLattice& lattice)
{
  //set energy to E0
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

  //set energy ramp
  if (gammaMode() == GammaMode::simtool
      || gammaMode() == GammaMode::simtool_no_interpolation
      || trajectoryMode() == TrajectoryMode::simtool) {
    if (simToolRamp()) {
      if (palattice->tool==pal::SimTool::elegant) {
	palattice->elegantEnergyRamp.tStart = t_start();
	palattice->elegantEnergyRamp.tStop = t_stop();
	palattice->elegantEnergyRamp.nSteps = simToolRampSteps();
	palattice->elegantEnergyRamp.set([&](double t) {return gamma(t)*E_rest_GeV/E0();} );
	std::cout << "* " << palattice->tool_string() << " energy ramp set" << std::endl;
      }
      else
	std::cout << "WARNING: Setting SimTool energy ramp is not implemented for " << palattice->tool_string() << std::endl
		  << "         Option <simToolRamp> is ignored." << std::endl;
    }
  }
}

// #turns for resonance strengths calc are calculated from stepwidth dagamma() if not set
unsigned int Configuration::numTurns() const
{
  if (_nTurns != 0)
    return _nTurns;
  else {
    return (unsigned int) (1./std::fabs(dagamma()) + 0.5);
  }
}







// class RfMagnetConfig
// parse configuration of multiple rf magnets from comma separated lists
// in a polematrix config file (boost property tree)
// and write them to a pal::AccLattice

void RfMagnetConfig::set(const pt::ptree &tree)
{
  setFromStringList<std::string>(elements, tree.get<std::string>("palattice.rfMagnets.elements", ""));
  setFromStringList<double>(Q1, tree.get<std::string>("palattice.rfMagnets.Q1", ""));
  setFromStringList<double>(dQ, tree.get<std::string>("palattice.rfMagnets.dQ", ""));
  setFromStringList<unsigned int>(period, tree.get<std::string>("palattice.rfMagnets.period", ""));

  auto n = elements.size();
  if (Q1.size()!=n || dQ.size()!=n || period.size()!=n)
    throw std::runtime_error("Cannot set up RF magnets from config file! Unequal number of entries.");
}


void RfMagnetConfig::writeToLattice(pal::AccLattice& lattice) const
{
  for (auto i=0u; i<elements.size(); i++) {
    lattice[elements[i]].element()->Qrf1 = Q1[i];
    lattice[elements[i]].element()->dQrf = dQ[i];
    lattice[elements[i]].element()->rfPeriod = period[i];
    std::cout << "* set up " << elements[i] << " as RF magnet (Q1=" << Q1[i]
	      << ", dQ=" << dQ[i] << ", period=" << period[i] << ")" << std::endl;
  }
}


void RfMagnetConfig::writeToConfig(pt::ptree &tree) const
{
  if (elements.size() == 0)
    return;

  tree.put("palattice.rfMagnets.elements", getElements());
  tree.put("palattice.rfMagnets.Q1", getQ1());
  tree.put("palattice.rfMagnets.dQ", getDQ());
  tree.put("palattice.rfMagnets.period", getPeriod());
}


template<>
void RfMagnetConfig::setFromStringList<std::string>(std::vector<std::string>& v, std::string list)
{
  list.erase(std::remove_if(list.begin(), list.end(), ::isspace), list.end());
  std::istringstream ss(list);
  std::string tmp;
  while ( std::getline(ss,tmp,',') ) {
     v.push_back(tmp);
   }
}
