#include <iostream>
#include "Tracking.hpp"


Tracking::Tracking(unsigned int nThreads) : error(false), lattice(0., pal::Anchor::end), orbit(0., gsl_interp_akima_periodic), gammaCentral(0.), showProgressBar(true)
{
  // use at least one thread
  if (nThreads == 0)
    nThreads = 1;

  // set iterator to begin of queue
  queueIt = queue.begin();

  // create threads
  for (unsigned int i=0; i<nThreads; i++) {
    threadPool.emplace_back(std::thread());
  }
}



void Tracking::start()
{ 
  if (lattice.size()==0 || orbit.size()==0)
    throw TrackError("Cannot start tracking, if model is not specified (Lattice, Orbit).");

  if (config.t_stop() <= config.t_start()) {
    std::stringstream msg;
    msg << "Cannot track backwards: t_stop=" << config.t_stop() << " < t_start=" << config.t_start();
    throw TrackError(msg.str());
  }

  // fill queue
  for (unsigned int i=0; i<config.nParticles(); i++) {
    queue.emplace_back( TrackingTask(i,config) );
  }
  // set iterator to begin of queue
  queueIt = queue.begin();

  // write current config to file
  config.save( config.confOutFile().string() );

  auto start = std::chrono::high_resolution_clock::now();

  //start threads
  for (std::thread& t : threadPool) {
    t = std::thread(&Tracking::processQueue,this);
  }

  //start extra thread for progress bars
  if (showProgressBar) {
  sleep(1);
  std::thread progress(&Tracking::printProgress,this);
  progress.join();
  }
  else {
    std::cout << "* start tracking..." << std::endl;
  }

  // wait for threads to finish
  for (std::thread& t : threadPool) {
    t.join();
  }

  auto stop = std::chrono::high_resolution_clock::now();
  auto secs = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  auto mins = std::chrono::duration_cast<std::chrono::minutes>(stop-start);
  std::cout << "-----------------------------------------------------------------" << std::endl;
  if (error)
    std::cout << "An error occured during tracking!\nAt least one Spin was not tracked successfully.\nStopped after ";
  else
    std::cout << "Tracking "<<numParticles()<< " Spins done. Tracking took ";
  std::cout << secs.count() << " s = "<< mins.count() << " min." << std::endl;
  std::cout << "-----------------------------------------------------------------" << std::endl;
  
  if (!error)
    calcPolarization();
}


void Tracking::processQueue()
{
  std::vector<TrackingTask>::iterator myTask;

  while (true) {
    mutex.lock();
    if (queueIt == queue.end()) {
      mutex.unlock();
      return;   // finish this thread
    }
    else {
      myTask = queueIt;
      queueIt++;
      runningTasks.push_back(myTask); // to display progress
      mutex.unlock();
      myTask->lattice=&lattice;
      myTask->orbit=&orbit;
      try {
	myTask->initGamma(gammaCentral);  // simtool: sdds import thread safe since SDDSToolKit-devel-3.3.1-2
	myTask->initTrajectory();
	myTask->run(); // run next queued TrackingTask
      }
      //cancel thread in error case
      catch (std::exception &e) {
	std::cout << "ERROR @ particle " << myTask->particleId
		  << " (thread_id "<< std::this_thread::get_id() <<"):"<< std::endl
		  << e.what() << std::endl;
	mutex.lock();
	error = true;
	runningTasks.remove(myTask); // to display progress
	mutex.unlock();
	return;
      }
      mutex.lock();
      runningTasks.remove(myTask); // to display progress
      mutex.unlock();
    }//else
  }//while
}


void Tracking::printProgress() const
{
  std::list<std::vector<TrackingTask>::const_iterator> tmp;
  unsigned int barWidth;
  if (runningTasks.size() < 5)
    barWidth = 20;
  else
    barWidth = 15;
  //shorter looks ugly
  
  while (runningTasks.size() > 0) {
    tmp = runningTasks;
    unsigned int n=0;
    for (std::vector<TrackingTask>::const_iterator task : tmp) {
      if (n<2)
	std::cout << task->getProgressBar(barWidth) << "  ";
      else
	std::cout << task->getProgressBar(0) << " ";
      n++;
    }
    std::cout <<"\r"<< std::flush;
    sleep(1);
  }
}

void Tracking::setModel()
{
  setLattice();
  setOrbit();

  //set energy in SimTool to E0 (ramp not considered!)
  double p_MeV = config.E0()*1000.;
  config.getSimToolInstance().setMomentum_MeV(p_MeV);

  
  //set number of turns for SimTool based on tracking time
  if (config.gammaMode() == GammaMode::simtool || config.gammaMode()==GammaMode::simtool_plus_linear || config.trajectoryMode() == TrajectoryMode::simtool) {
    unsigned int turns = (config.duration()*GSL_CONST_MKSA_SPEED_OF_LIGHT / lattice.circumference()) + 1;
    config.getSimToolInstance().verbose = true;
    config.getSimToolInstance().setTurns(turns);
    std::cout << "* Elegant tracking " << turns <<" turns to get single particle trajectories" << std::endl;
  }
  

  //set physical quantities from SimTool if not set by config
  if (config.gammaMode() == GammaMode::simtool || config.gammaMode() == GammaMode::simtool_plus_linear) {
    gammaCentral = config.getSimToolInstance().readGammaCentral();
  }
  if (config.gammaMode() == GammaMode::radiation) {
    if (config.q()==0.) {
      config.set_q(lattice.overvoltageFactor(config.gamma_start()));
      std::cout << "* set overvoltage factor from lattice"
		<< ": q=" << config.q() << std::endl;
	}
    if (config.h()==0) {
      config.set_h(lattice.harmonicNumber());
      std::cout << "* set harmonic number from lattice"
		<< ": h=" << config.h() << std::endl;
    }
    if (config.R()==0.) {
      config.set_R(lattice.integralDipoleRadius());
      std::cout << "* set dipole bending radius from lattice"
		<< ": R=" << config.R() << std::endl;
    }
    if (config.alphac()==0.) {
      config.set_alphac(config.getSimToolInstance().readAlphaC());
      std::cout << "* set momentum compaction factor from " << config.getSimToolInstance().tool_string()
		<< ": alphac=" << config.alphac() << std::endl;
    }
    if (config.Js()==0.) {
      config.set_Js(config.getSimToolInstance().readDampingPartitionNumber_syli().s);
      std::cout << "* set long. damping partition number from " << config.getSimToolInstance().tool_string()
		<< ": Js=" << config.Js() << std::endl;
    }
  }
}

void Tracking::setLattice()
{
  lattice = pal::AccLattice(config.getSimToolInstance());
}

void Tracking::setOrbit()
{
  orbit = pal::FunctionOfPos<pal::AccPair>( config.getSimToolInstance() );
  orbit.simToolClosedOrbit( config.getSimToolInstance() );
}


//calculate polarization: average over all spin vectors for each time step
void Tracking::calcPolarization()
{
  polarization = queue[0].getStorage();
  
  for (unsigned int i=1; i<queue.size(); i++) {
    polarization += queue[i].getStorage();
  }
  polarization /= numParticles();
}

void Tracking::savePolarization()
{
  std::ofstream file;
  std::string filename = config.polFile().string();
  unsigned int w = 14;
  
  file.open(filename);
  if (!file.is_open())
    throw TrackFileError(filename);

  file << polarization.printHeader(w, "P") << std::endl;
  file << polarization.print(w);
  
  file.close();
  std::cout << "* Wrote polarization for " << polarization.size() << " steps to " << filename <<"."<< std::endl;
}
