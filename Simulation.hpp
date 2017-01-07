#ifndef __POLEMATRIX__SIMULATION_HPP_
#define __POLEMATRIX__SIMULATION_HPP_

#include <libpalattice/AccLattice.hpp>
#include <libpalattice/FunctionOfPos.hpp>
#include "Configuration.hpp"

class Simulation {
protected:
  std::shared_ptr<pal::AccLattice> lattice;
  std::shared_ptr<pal::FunctionOfPos<pal::AccPair>> orbit;
  
  std::map<unsigned int,std::string> errors;

public:
  const std::shared_ptr<Configuration> config;
  
  Simulation() : config(new Configuration) {}
  void setModel();
  
  bool modelReady() {if (lattice->size()==0 || orbit->size()==0) return false; else return true;}
  unsigned int numParticles() const {return config->nParticles();}
  unsigned int numSuccessful() const {return numParticles() - errors.size();}
    
  virtual void start() =0;
  
  void saveLattice() const {lattice->print( (config->outpath()/"lattice.dat").string() );}
  void saveOrbit() const {orbit->print( (config->outpath()/"closedorbit.dat").string() );}
};

#endif
// __POLEMATRIX__SIMULATION_HPP_
