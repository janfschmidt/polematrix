#include <iostream>
#include "Tracking.hpp"
#include "version.hpp"


Tracking::Tracking(unsigned int nThreads) : gammaCentral(0.), config(new Configuration), showProgressBar(true)
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
  if (lattice->size()==0 || orbit->size()==0)
    throw TrackError("Cannot start tracking, if model is not specified (Lattice, Orbit).");

  if (config->t_stop() <= config->t_start()) {
    std::stringstream msg;
    msg << "Cannot track backwards: t_stop=" << config->t_stop() << " < t_start=" << config->t_start();
    throw TrackError(msg.str());
  }

  // fill queue
  for (unsigned int i=0; i<config->nParticles(); i++) {
    queue.emplace_back( TrackingTask(i,config) );
  }
  // set iterator to begin of queue
  queueIt = queue.begin();

  // write current config to file
  config->save( config->confOutFile().string() );

  std::cout << "Start tracking "<<config->nParticles()<<" Spins..." << std::endl;
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

  // wait for threads to finish
  for (std::thread& t : threadPool) {
    t.join();
  }

  auto stop = std::chrono::high_resolution_clock::now();
  auto secs = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout << std::endl
	    << "-----------------------------------------------------------------" << std::endl;
  if (errors.size()>0) {
    std::cout << "ERRORS occured during tracking!" << std::endl;
  }
  std::cout << "Tracking "<<numSuccessful()<< " Spins done. Tracking took ";
  std::cout << secs.count() << " s = "<< int(secs.count()/60.+0.5) << " min." << std::endl;
  std::cout << "Thanks for using polematrix " << polemversion() << std::endl;
  for (auto& it : errors)
    std::cout << "ERROR @ particle " << it.first << ": " << it.second << std::endl;
  std::cout << "-----------------------------------------------------------------" << std::endl;
  
  if (numSuccessful() > 0)
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
      try {
	myTask->setModel(lattice,orbit);
	myTask->initGamma(gammaCentral);  // simtool: sdds import thread safe since SDDSToolKit-devel-3.3.1-2
	myTask->initTrajectory();
	myTask->run(); // run next queued TrackingTask
      }
      //cancel thread in error case
      catch (std::exception &e) {
	std::cout << "ERROR @ particle " << myTask->particleId
	  // << " (thread_id "<< std::this_thread::get_id() << ")"
		  <<":"<< std::endl
		  << e.what() << std::endl;
	mutex.lock();
	errors.emplace(myTask->particleId, e.what());
	mutex.unlock();
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
  auto numTasks = runningTasks.size();
  if ( numTasks < 5)
    barWidth = 20;
  else
    barWidth = 15;
  //shorter looks ugly
  
  while (runningTasks.size() > 0) {
    tmp = runningTasks;
    unsigned int n=0;
    for (std::vector<TrackingTask>::const_iterator task : tmp) {
      if (n<2) // first 2 with progress bar
	std::cout << task->getProgressBar(barWidth) << "  ";
      else     // others percentage only
	std::cout << task->getProgressBar(0) << " ";
      n++;
    }
    while (n<numTasks) { // clear finished tasks
      std::cout << "       ";
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
  double p_MeV = config->E0()*1000.;
  config->getSimToolInstance().setMomentum_MeV(p_MeV);

  
  //set number of turns for SimTool based on tracking time
  if (config->gammaMode() == GammaMode::simtool
      || config->gammaMode()==GammaMode::simtool_plus_linear
      || config->gammaMode()==GammaMode::simtool_no_interpolation
      || config->trajectoryMode() == TrajectoryMode::simtool) {
    unsigned int turns = (config->duration()*GSL_CONST_MKSA_SPEED_OF_LIGHT / lattice->circumference()) + 1;
    config->getSimToolInstance().verbose = true;
    config->getSimToolInstance().setTurns(turns);
    std::cout << "* Elegant tracking " << turns <<" turns to get single particle trajectories" << std::endl;
  }
  

  //set physical quantities from SimTool if not set by config
  if (config->gammaMode() == GammaMode::simtool
      || config->gammaMode() == GammaMode::simtool_plus_linear) {
    gammaCentral = config->getSimToolInstance().readGammaCentral();
  }
  else if (config->gammaMode() == GammaMode::linear)
    {
      // no model setup needed
    }
  else {
    if (config->q()==0.) {
      config->set_q(lattice->overvoltageFactor(config->gamma_start()));
      std::cout << "* set overvoltage factor from lattice"
		<< ": q=" << config->q() << std::endl;
	}
    if (config->h()==0) {
      config->set_h(lattice->harmonicNumber());
      std::cout << "* set harmonic number from lattice"
		<< ": h=" << config->h() << std::endl;
    }
    if (config->R()==0.) {
      config->set_R(lattice->integralDipoleRadius());
      std::cout << "* set dipole bending radius from lattice"
		<< ": R=" << config->R() << std::endl;
    }
    if (config->alphac()==0.) {
      config->set_alphac(config->getSimToolInstance().readAlphaC());
      std::cout << "* set momentum compaction factor from " << config->getSimToolInstance().tool_string()
		<< ": alphac=" << config->alphac() << std::endl;
    }
    if (config->alphac2()==0.) {
      config->set_alphac2(config->getSimToolInstance().readAlphaC2());
      std::cout << "* set 2nd order momentum compaction factor from " << config->getSimToolInstance().tool_string()
		<< ": alphac2=" << config->alphac2() << std::endl;
    }
    if (config->Js()==0.) {
      config->set_Js(config->getSimToolInstance().readDampingPartitionNumber_syli().s);
      std::cout << "* set long. damping partition number from " << config->getSimToolInstance().tool_string()
		<< ": Js=" << config->Js() << std::endl;
    }
  }
}

void Tracking::setLattice()
{
  lattice.reset( new pal::AccLattice(config->getSimToolInstance()) );
}

void Tracking::setOrbit()
{
  orbit.reset( new pal::FunctionOfPos<pal::AccPair>(config->getSimToolInstance()) );
  orbit->simToolClosedOrbit( config->getSimToolInstance() );
}


//calculate polarization: average over all successfully tracked spin vectors for each time step
void Tracking::calcPolarization()
{
  unsigned int i=0;
  for (; i<queue.size(); i++) {
    if (errors.count(i)==0) {
      polarization = queue[i].getStorage();
      break;
    }
  }

  for (i++; i<queue.size(); i++) {
    if (errors.count(i)==0)
      polarization += queue[i].getStorage();
  }
  polarization /= numSuccessful();
}

void Tracking::savePolarization()
{
  std::ofstream file;
  std::string filename = config->polFile().string();
  unsigned int w = 14;
  
  file.open(filename);
  if (!file.is_open())
    throw TrackFileError(filename);

  file << "# Polarization calculated as average over " << numSuccessful() << " spins" << std::endl;
  file << polarization.printHeader(w, "P") << std::endl;
  file << polarization.print(w);
  
  file.close();
  std::cout << "* Polarization written for " << polarization.size() << " steps to " << filename <<"."<< std::endl;
}
