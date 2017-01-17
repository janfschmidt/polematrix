/* Simulation base classes
 * used by Tracking/TrackingTask and ResStrengths/ParticleResStrengths
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

#ifndef __POLEMATRIX__SIMULATION_HPP_
#define __POLEMATRIX__SIMULATION_HPP_

#include <libpalattice/AccLattice.hpp>
#include <libpalattice/FunctionOfPos.hpp>
#include "Configuration.hpp"
#include "Trajectory.hpp"


// abstract base class for a simulation
// including a config as well as lattice and orbit, which can be set using config

class Simulation {
protected:
  std::shared_ptr<pal::AccLattice> lattice;
  std::shared_ptr<pal::FunctionOfPos<pal::AccPair>> orbit;
  
  std::map<unsigned int,std::string> errors;

public:
  const std::shared_ptr<Configuration> config;
  
  Simulation() : config(new Configuration) {}
  Simulation(const std::shared_ptr<Configuration> c) : config(c) {}
  void setModel();
  
  bool modelReady() {if (lattice->size()==0 || orbit->size()==0) return false; else return true;}
  unsigned int numParticles() const {return config->nParticles();}
  unsigned int numSuccessful() const {return numParticles() - errors.size();}
    
  virtual void start() =0;
  
  void saveLattice() const {lattice->print( (config->outpath()/"lattice.dat").string() );}
  void saveOrbit() const {orbit->print( (config->outpath()/"closedorbit.dat").string() );}
};



// abstract base class for a simulation task for a single particle
// can be used by a Simulation object, which passes its config/lattice/orbit
// to many SingleParticleSimulation

class SingleParticleSimulation {
protected:
  std::unique_ptr<Trajectory> trajectory;     // particle trajectory, implementation depends TrajectoryMode
  
public:
  const unsigned int particleId;
  const std::shared_ptr<Configuration> config; //not const Object, because SimToolInstance status can be changed
  std::shared_ptr<const pal::AccLattice> lattice;
  std::shared_ptr<const pal::FunctionOfPos<pal::AccPair>> orbit;

  SingleParticleSimulation(unsigned int id, const std::shared_ptr<Configuration> c);
  void setModel(std::shared_ptr<const pal::AccLattice> l, std::shared_ptr<const pal::FunctionOfPos<pal::AccPair>> o);
  virtual void run() =0;
};

#endif
// __POLEMATRIX__SIMULATION_HPP_
