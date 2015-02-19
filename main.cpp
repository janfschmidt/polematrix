#include <iostream>
#include <libpal/AccLattice.hpp>
#include <libpal/FunctionOfPos.hpp>
#include "Tracking.hpp"


int main()
{
  Tracking t(2,1);
  std::cout << t.numThreads() << " threads, " << t.numParticles() << " tasks." << std::endl;

  t.config.setPath("/home/schmidt/pole/MyProjects/libpalfield");
  t.config.x = 1; //rotation per step in degree
  t.config.s_start = {0,0.7,0.714};
  t.config.t_start = 0.;
  t.config.t_stop = 548e-9;
  t.config.dt_out = 1e-9;
  t.config.E_start = 1.32194;
  t.config.dE = 0.;

  pal::SimToolInstance sim(pal::madx, pal::online, t.config.subfolder("madx")+ "elsa_harmcorr.seq");
  pal::AccLattice lattice("polematrix",sim);
  pal::FunctionOfPos<pal::AccPair> orbit(sim);
  orbit.simToolClosedOrbit(sim);

  t.setModel(&sim,&lattice,&orbit);

  t.start();

  return 0;
}
