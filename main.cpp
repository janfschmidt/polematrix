/* polematrix - matrix spin tracking code
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

#include <iostream>
#include <iomanip>
#include <boost/program_options.hpp>
#include <libpalattice/AccLattice.hpp>
#include <libpalattice/FunctionOfPos.hpp>
#include "Tracking.hpp"
#include "ResStrengths.hpp"
#include "version.hpp"

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
  
  po::options_description modes("Program modes");
  modes.add_options()
    ("help,h", "display this help message")
    ("version,v", "display version")
    ("template,T", "create config file template (template.pole) and quit")
    ("resonance-strengths,R", "estimate strengths of depolarizing resonances")
    ;

  po::options_description confs("Configuration options");
  confs.add_options()  
    ("threads,t", po::value<unsigned int>(&nThreads)->default_value(std::thread::hardware_concurrency()), "number of threads used for tracking")
    ("output-path,o", po::value<std::string>(&outpath)->default_value("."), "path for output files")
    ("verbose,V", "more output, e.g. each written spin file")
    ("no-progressbar,n", "do not show progress bar during tracking")
    ("all,a", "write all output (e.g. lattice and orbit)")
    ("spintune,s", po::value<double>(), "in resonance-strengths mode: calculate for given spin tune only")
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
  try {
    po::store(po::command_line_parser(argc, argv).options(all).positional(pd).run(), args);
    po::notify(args);
  }
  catch(po::error_with_option_name e){
    std::cout << "ERROR: " << e.what() << std::endl;
    std::cout << "use -h for help." << std::endl;
    return 1;
  }


  // special program modes and abortion due to input errors
  if (args.count("help")) {
    usage(visible);
    return 0;
  }

  if (args.count("version")) {
    std::cout << "polematrix " << polemversion() << std::endl;
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
  
  t.config->load( configfile );
  // alternatively config can be set hardcoded- e.g.:
  // t.config->nParticles = 2;
  // t.config->s_start = {0.,0.,1.};
  // t.config->t_start = 0.085;

  
  // options from command line, NOT stored in config file
  t.config->set_outpath(outpath);
  if (args.count("verbose"))
    t.config->set_verbose(true);
  if (args.count("no-progressbar"))
    t.showProgressBar = false;


  
  // resonance strengths mode (no tracking)
  if (args.count("resonance-strengths")) {
    ResStrengths r(t.config);
    // initialize model from simtool
    try {
      r.setModel();
    }
    catch (pal::palatticeError &e) {
      std::cout << e.what() << std::endl << "Quit." << std::endl;
      return 3;
    }
    if (args.count("all")) {
      r.saveLattice();
      r.saveOrbit();
    }
    try{
      if (args.count("spintune")) {
	std::cout << r.getSingle( args["spintune"].as<double>() ) << std::endl;
      }
      else {
	r.start();
	r.save();
      }
    }
    catch (std::exception &e) {
      std::cout << e.what() << std::endl << "Quit." << std::endl;
      return 2;
    }
    return 0;
  }


  
  t.config->printSummary();
  
  // initialize model from simtool
  try {
    t.setModel();
  }
  catch (pal::palatticeError &e) {
    std::cout << e.what() << std::endl << "Quit." << std::endl;
    return 3;
  }

  if (args.count("all")) {
    t.saveLattice();
    t.saveOrbit();
  }


  // start tracking
  try{
    t.start();
  }
  catch (TrackError &e) {
    std::cout << e.what() << std::endl << "Quit." << std::endl;
    return 2;
  }

  t.savePolarization();

  return 0;
}
