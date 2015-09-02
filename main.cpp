#include <iostream>
#include <iomanip>
#include <libpalattice/AccLattice.hpp>
#include <libpalattice/FunctionOfPos.hpp>
#include "Tracking.hpp"


int main()
{
  Tracking t(2,2);
  std::cout << t.numThreads() << " threads, " << t.numParticles() << " tasks." << std::endl;

  t.config.setPath("/home/schmidt/pole/MyProjects/depolVStracking");

  // t.config.s_start = {0,0.7,0.714};
  // t.config.t_start = 0.;
  // t.config.t_stop = 548e-8;
  // t.config.E_start = 1.32194;
  // t.config.dE = 0.;
  // t.config.dt_out = 1e-9;

  t.config.s_start = {0.,0.,1.};
  t.config.t_start = 0.085;
  t.config.t_stop =  0.104;
  t.config.E_start = 1.2;
  t.config.dE = 6.0;
  t.config.dt_out = 2e-5;
  t.config.save("test.conf");
  //t.config.load("test.conf");

  std::cout << "dpos: " <<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<< t.config.dpos_out() << std::endl;

  pal::SimToolInstance sim(pal::madx, pal::offline, t.config.subfolder("madx")+ "madx.twiss");
  pal::AccLattice lattice("polematrix",sim);
  pal::FunctionOfPos<pal::AccPair> orbit(sim);
  orbit.simToolClosedOrbit(sim);

  t.setModel(&sim,&lattice,&orbit);

  t.start();

  return 0;
}
