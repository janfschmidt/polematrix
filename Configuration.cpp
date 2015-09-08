#include "Configuration.hpp"

Configuration::Configuration(std::string pathIn)
  : E_rest(0.000511), a_gyro(0.001159652), default_steps(200), spinDirName("spins"), polFileName("polarization.dat"), confOutFileName("currentconfig.pole")
{
  //setPath(pathIn);
  outpath = pathIn;
  nParticles = 1;

  t_start = t_stop = 0.;
  dt_out = 1e-5;
  E0 = 1;
  dE = 0;
  s_start.zeros();
  s_start[2] = 1;

  palattice = new pal::SimToolInstance(pal::madx, pal::offline, "");
}

Configuration::~Configuration()
{
  delete palattice;
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
  tree.put("spintracking.numParticles", nParticles);
  tree.put("spintracking.t_start", t_start);
  tree.put("spintracking.t_stop", t_stop);
  tree.put("spintracking.dt_out", dt_out);
  tree.put("spintracking.E0", E0);
  tree.put("spintracking.dE", dE);
  tree.put("spintracking.s_start.x", s_start[0]);
  tree.put("spintracking.s_start.z", s_start[2]);
  tree.put("spintracking.s_start.s", s_start[1]);
  tree.put("palattice.simTool", palattice->tool_string());
  tree.put("palattice.mode", palattice->mode_string());
  tree.put("palattice.file", palattice->inFile());

  pt::xml_writer_settings<char> settings(' ', 4); //indentation
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
    t_stop = tree.get<double>("spintracking.t_stop"); 
    E0 = tree.get<double>("spintracking.E0");
    dE = tree.get<double>("spintracking.dE");
    s_start[0] = tree.get<double>("spintracking.s_start.x");
    s_start[2] = tree.get<double>("spintracking.s_start.z");
    s_start[1] = tree.get<double>("spintracking.s_start.s");
    
    setSimToolInstance(tree);    
  }
  catch (pt::ptree_error &e) {
    std::cout << "Error loading configuration file:" << std::endl
	      << e.what() << std::endl;
    exit(1);
  }
    
  // optional config with default values
  nParticles = tree.get("spintracking.numParticles", 1);
  t_start = tree.get("spintracking.t_start", 0.0);
  dt_out = tree.get("spintracking.dt_out", (t_stop-t_start)/default_steps);

  std::cout << "* configuration loaded from " << filename << std::endl;
  return;
}

void Configuration::printSummary() const
{
  std::stringstream s;

  s << "-----------------------------------------------------------------" << std::endl;
  s << "Tracking " << nParticles << " Spins" << std::endl
    << "time      " <<  t_start << " s   -------------------->   " << t_stop << " s" << std::endl
    << "energy    " << gamma(t_start)*E_rest << " GeV   ----- " << dE << " GeV/s ----->   " << gamma(t_stop)*E_rest << " GeV" << std::endl
    << "spin tune " << agamma(t_start) << "   -------------------->   " << agamma(t_stop) << std::endl;
  s << "start polarization: Px = " << s_start[0] << ", Ps = " << s_start[1] << ", Pz = " << s_start[2] << std::endl;
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

  delete palattice;
  palattice = new pal::SimToolInstance(tool, mode, file);

}
