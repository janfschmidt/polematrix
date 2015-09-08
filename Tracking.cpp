#include <iostream>
#include "Tracking.hpp"


Tracking::Tracking(unsigned int nThreads) : lattice(NULL), orbit(NULL)
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

Tracking::~Tracking()
{
  delete lattice;
  delete orbit;
}


void Tracking::start()
{ 
  if (lattice==NULL || orbit==NULL)
    throw TrackError("ERROR: Cannot start tracking, if model is not specified (Lattice, Orbit).");

  if (config.t_stop <= config.t_start) {
    std::stringstream msg;
    msg << "ERROR: Cannot track backwards: t_stop=" << config.t_stop << " < t_start=" << config.t_start;
    throw TrackError(msg.str());
  }

  // fill queue
  for (unsigned int i=0; i<config.nParticles; i++) {
    TrackingTask toll(i,config);
    queue.push_back(std::move(toll));
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

  // wait for threads to finish
  for (std::thread& t : threadPool) {
    t.join();
  }

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
  std::cout << "Tracking took " << duration.count() << " s." << std::endl;

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
      mutex.unlock();
      myTask->lattice=lattice;
      myTask->orbit=orbit;
      try {
	myTask->run(); // run next queued TrackingTask
      }
      catch (TrackError &e) {
	std::cout << "ERROR: " << e.what() << " (thread_id " << std::this_thread::get_id()
		  << ", particle " << myTask->particleId << ")"<< std::endl;
	return;
      }
    }//else
  }//while
}


void Tracking::setLattice()
{
  lattice = new pal::AccLattice("polematrix", config.getSimToolInstance());
}

void Tracking::setOrbit()
{
  pal::FunctionOfPos<pal::AccPair> *tmp = new pal::FunctionOfPos<pal::AccPair>(config.getSimToolInstance());
  tmp->simToolClosedOrbit(config.getSimToolInstance());
  orbit = tmp;
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
