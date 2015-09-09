#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <libpalattice/AccLattice.hpp>
#include <libpalattice/FunctionOfPos.hpp>
#include "Tracking.hpp"

namespace po = boost::program_options;



void usage(po::options_description &desc)
{
  std::cout << std::endl
	    << "polematrix [CONFIGURATION FILE] [options]" <<std::endl<<std::endl
	    << "[CONFIGURATION FILE] holds the tracking parameters." <<std::endl
	    << "A template config file can be generated with option --template" <<std::endl<<std::endl<<std::endl
	    << "Allowed options:" <<std::endl;
  std::cout << desc << std::endl;
  return;
}



int main(int argc, char *argv[])
{
  unsigned int nThreads;
  std::string configfile;
    std::string outpath;
  
  // Declare the supported options
  po::options_description modes("Program modes");
  modes.add_options()
    ("help,h", "display this help message")
    ("template,T", "create config file template (template.pole) and quit")
    ;

  po::options_description confs("Configuration options");
  confs.add_options()  
    ("threads,t", po::value<unsigned int>(&nThreads)->default_value(std::thread::hardware_concurrency()), "number of threads used for tracking")
    ("output-path,o", po::value<std::string>(&outpath)->default_value("."), "path for output files")
    ("verbose,v", "show progress bar during tracking")
    ;
  
  po::options_description hidden("Hidden options");
  hidden.add_options()  
    ("config", po::value<std::string>(&configfile), "the configuration file for the tracking")
    ;
	
  po::positional_options_description pd;
  pd.add("config", 1);

  po::options_description all;
  all.add(modes).add(confs).add(hidden);
  po::options_description visible;
  visible.add(modes).add(confs);

  po::variables_map args;
  //  po::store(po::parse_command_line(argc, argv, all), args);
  po::store(po::command_line_parser(argc, argv).options(all).positional(pd).run(), args);
  po::notify(args);


  // special program modes and abortion due to input errors
  if (args.count("help")) {
    usage(visible);
    return 0;
  }
  
  if (args.count("template")) {
    Configuration empty;
    empty.save("template.pole");
    return 0;
  }

  if (!args.count("config")) {
    std::cout << "ERROR: No configuration file given. Use -h for help." << std::endl;
    return 1;
  }


  // initialize tracking & config
  Tracking t(nThreads);
  std::cout << "* " << t.numThreads() << " threads used." << std::endl;
  
  t.config.load( configfile );
  t.config.outpath = outpath; // NOT in config file
  
  t.config.printSummary();

  if (args.count("verbose"))
    t.showProgressBar = true;

  //tool/mode/file muss in constructor gesetzt werden
  //muss weiterexistieren für alle initialisierungen, NICHT für ganze lebensdauer von lattice/orbit etc
  //lattice&orbit müssen in Tracking const sein, möglichst nicht kopieren
  
  //trajectory müsste beim thread-starten .simToolTrajectory(sim) aufrufen
  // => sim muss weiterleben (in tracking oder config oder task?)
  // => trajectory muss in task leben (nicht pointer) da individuell verschieden

  t.setModel();
  
  try{
    t.start();
  }
  catch (TrackError &e) {
    std::cout << e.what() << std::endl << "Quit.";
    return 2;
  }
  
  t.savePolarization();

  return 0;
}




  // t.config.s_start = {0,0.7,0.714};
  // t.config.t_start = 0.;
  // t.config.t_stop = 548e-8;
  // t.config.E0 = 1.32194;
  // t.config.dE = 0.;
  // t.config.dt_out = 1e-9;

  // t.config.nParticles = 2;
  // t.config.s_start = {0.,0.,1.};
  // t.config.t_start = 0.085;
  // t.config.t_stop =  0.104;
  // t.config.E0 = 1.2;
  // t.config.dE = 6.0;
  // t.config.dt_out = 2e-5;
