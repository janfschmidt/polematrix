#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <libpalattice/AccLattice.hpp>
#include <libpalattice/FunctionOfPos.hpp>
#include "Tracking.hpp"

namespace po = boost::program_options;


int main(int argc, char *argv[])
{
  unsigned int nThreads;
  std::string toll;
  
  // Declare the supported options
  po::options_description desc("Allowed options");
  desc.add_options()
    ("config,c", po::value<std::string>(&toll), "the configuration file for the tracking")
    ("help,h", "display this help message")
    ("threads,t", po::value<unsigned int>(&nThreads)->default_value(std::thread::hardware_concurrency()), "number of threads used for tracking")
    ("template,T", "create config file template (template.pole) and quit")
    ;
  po::positional_options_description pd;
  pd.add("config", 1);

  po::variables_map args;
  //  po::store(po::parse_command_line(argc, argv, desc), args);
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), args);
  po::notify(args);
  
  if (args.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }
  if (args.count("template")) {
    Configuration empty;
    empty.save("template.pole");
    return 0;
  }

  Tracking t(nThreads);
  std::cout << t.numThreads() << " threads" << std::endl;
    std::cout << toll << std::endl;
  
  t.config.setPath("/home/schmidt/pole/MyProjects/depolVStracking");

  // t.config.s_start = {0,0.7,0.714};
  // t.config.t_start = 0.;
  // t.config.t_stop = 548e-8;
  // t.config.E0 = 1.32194;
  // t.config.dE = 0.;
  // t.config.dt_out = 1e-9;

  t.config.nParticles = 2;
  t.config.s_start = {0.,0.,1.};
  t.config.t_start = 0.085;
  t.config.t_stop =  0.104;
  t.config.E0 = 1.2;
  t.config.dE = 6.0;
  t.config.dt_out = 2e-5;
  // t.config.save("test.conf");
  t.config.load( args["config"].as<std::string>() );

  std::cout << "dpos: " <<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<< t.config.dpos_out() << std::endl;

  pal::SimToolInstance sim(pal::madx, pal::offline, t.config.subfolder("madx")+ "madx.twiss");
  pal::AccLattice lattice("polematrix",sim);
  pal::FunctionOfPos<pal::AccPair> orbit(sim);
  orbit.simToolClosedOrbit(sim);

  t.setModel(&sim,&lattice,&orbit);

  try{
    t.start();
  }
  catch (TrackError &e) {
    std::cout << e.what() << std::endl << "Quit.";
    return 1;
  }

  return 0;
}
