/* ResStrengths Classes
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
#include "version.hpp"


class ResStrengthsData {
protected:
  std::map<double,std::complex<double> > cache;
  const std::complex<double> im;                  // imaginary unit i

  std::string header(unsigned int w=16) const;
  void cacheIt(double agamma, const std::complex<double>& epsilon);
  virtual std::complex<double> calculate(double agamma) =0;
  std::complex<double> get(double agamma);               // get res. strength from cache or calculate it

public:
  ResStrengthsData() : im(std::complex<double> (0,1)) {}
  std::complex<double> operator[](double agamma);        // get res. strength from cache

  std::string printSingle(double agamma, std::complex<double> epsilon) const; // formated output of entry
  std::string printSingle(const std::pair<double,std::complex<double> >& it) const {return printSingle(it.first,it.second);}
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
  void runSingle();
};



class ResStrengths : public Simulation<ParticleResStrengths>, public ResStrengthsData {
protected:
  std::complex<double> calculate(double agamma);
  void init();
  
public:
  Metadata info;
  
  ResStrengths(unsigned int nThreads=std::thread::hardware_concurrency()) : Simulation(nThreads) {}
  ResStrengths(const std::shared_ptr<Configuration> c, unsigned int nThreads=std::thread::hardware_concurrency()) : Simulation(c,nThreads) {}

  void start();                                  // calculate & print all res. strengths according to config
  std::string getSingle(double agamma);          // calculate & print single resonance strength

  void print(string filename="");                //print all res. strengths in cache
  void save() {print( (config->outpath()/"resonance-strengths.dat").string() );}
};


#endif
/*__POLEMATRIX__RESSTRENGTHS_HPP_*/
