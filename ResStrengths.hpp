/* ResStrengths Class
 * estimate the complex strengths of depolarizing resonances from a lattice & orbit
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
 *
 *
 * Based on formulas by Courant & Ruth from:
 * Courant, E. D., and Ronald D. Ruth.
 * "The acceleration of polarized protons in circular accelerators."
 * BNL–51270 and UC–28 and ISA–80–5 (1980).
 *
 */

#ifndef __POLEMATRIX__RESSTRENGTHS_HPP_
#define __POLEMATRIX__RESSTRENGTHS_HPP_

#include <map>
#include <complex>
#include <memory>
#include <libpalattice/AccLattice.hpp>
#include <libpalattice/FunctionOfPos.hpp>
#include <libpalattice/Metadata.hpp>
#include "Simulation.hpp"
#include "Trajectory.hpp"


class ResStrengthsData {
protected:
  std::map<double,std::complex<double> > cache;
  const std::complex<double> im;                  // imaginary unit i

  void cacheIt(double agamma, std::complex<double>& epsilon);

public:
  ResStrengthsData() : im(std::complex<double> (0,1)) {}
  std::complex<double> operator[](double agamma);        // get res. strength (from cache)
};




// Resonance Strengths for a single particle

class ParticleResStrengths : public SingleParticleSimulation, public ResStrengthsData {
protected:
  std::complex<double> calculate(double agamma);    // calculate res. strength freq. omega=agamma
  
public:
  ParticleResStrengths(unsigned int id, const std::shared_ptr<Configuration> c,std::shared_ptr<const pal::AccLattice> l, std::shared_ptr<const pal::FunctionOfPos<pal::AccPair>> o);
  ParticleResStrengths(const ParticleResStrengths& other) = delete;
  ParticleResStrengths(ParticleResStrengths&& other) = default;
  void run();

  unsigned int nTurns() const {return config->numTurns(lattice->circumference());}
  double spintuneStep() const {return 1./double(nTurns());}

};




class ResStrengths : public Simulation, public ResStrengthsData {
protected:
  std::vector<ParticleResStrengths> particles;

  std::complex<double> calculate(double agamma);
  
public:
  Metadata info;
  
  ResStrengths();
  ResStrengths(const std::shared_ptr<Configuration> c) : Simulation(c) {}
  void start();                                    // calculate & print all res. strengths according to config

  //print all res. strengths in cache
  void print(string filename="");
  void save() {print( (config->outpath()/"resonance-strengths.dat").string() );}

  unsigned int nTurns() const {return config->numTurns(lattice->circumference());}
  double spintuneStep() const {return 1./double(nTurns());}
};


#endif
/*__POLEMATRIX__RESSTRENGTHS_HPP_*/
