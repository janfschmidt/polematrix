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
  t.config.t_start = 0.01;
  t.config.t_stop =  0.03;
  t.config.E_start = 1.2;
  t.config.dE = 6.0;
  t.config.dt_out = 2e-5;

  std::cout << "dpos: " <<std::setiosflags(std::ios::scientific)<<std::setprecision(8)<< t.config.dpos_out() << std::endl;

  pal::SimToolInstance sim(pal::madx, pal::online, t.config.subfolder("madx")+ "elsa.madx");
  pal::AccLattice lattice("polematrix",sim);
  pal::FunctionOfPos<pal::AccPair> orbit(sim);
  orbit.simToolClosedOrbit(sim);

  t.setModel(&sim,&lattice,&orbit);

  t.start();

  return 0;
}
