#include <iostream>
#include "Tracking.hpp"


Tracking::Tracking(unsigned int nThreads) : error(false), lattice("tut", 0, pal::end), orbit(0., gsl_interp_akima_periodic), showProgressBar(true)
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
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout << "-----------------------------------------------------------------" << std::endl;
  if (error)
    std::cout << "An error occured during tracking! Stopped after " << duration.count() << " s." << std::endl;
  else
    std::cout << "Tracking "<<numParticles()<< " Spins done. Tracking took " << duration.count() << " s." << std::endl;
  std::cout << "-----------------------------------------------------------------" << std::endl;

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
      myTask->lattice=&lattice;
      myTask->orbit=&orbit;
      myTask->initGammaSimTool(); // has to be mutexed, because all particles are read from the same sdds files
      mutex.unlock();
      try {
	myTask->run(); // run next queued TrackingTask
      }
      //cancel thread in error case
      catch (std::exception &e) {
	std::cout << "ERROR:"<< std::endl
		  << e.what() << " (thread_id " << std::this_thread::get_id()
		  << ", particle " << myTask->particleId << ")"<< std::endl;
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
  
  while (runningTasks.size() > 0) {
    tmp = runningTasks;
    unsigned int n=0;
    for (std::vector<TrackingTask>::const_iterator task : tmp) {
      n++;
      if (n<5)
	std::cout << task->getProgressBar(barWidth) << "   ";
      else
	std::cout << task->getProgressBar(0) << "   ";
    }
    std::cout <<"\r"<< std::flush;
    sleep(1);
  }
}

void Tracking::setModel(bool resetTurns)
{
  if (resetTurns)
    config.getSimToolInstance().setTurns(0);
  
  setLattice();
  setOrbit();

  //set number of turns for SimTool based on tracking time
  if (config.gammaMode() == simtool) {
    unsigned int turns = (config.duration()*GSL_CONST_MKSA_SPEED_OF_LIGHT / lattice.circumference()) + 1;
    std::cout << "* tracking " << turns <<" turns" << std::endl;
    config.getSimToolInstance().setTurns(turns);
  }
  else {
    config.getSimToolInstance().setTurns(0);
  }

}

void Tracking::setLattice()
{
  lattice = pal::AccLattice("polematrix", config.getSimToolInstance());
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
