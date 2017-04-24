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

#include "ResStrengths.hpp"
#include "debug.hpp"


// get res. strength from cache
std::complex<double> ResStrengthsData::operator[](double agamma)
{
  std::map<double,std::complex<double> >::const_iterator it = cache.find(agamma);
  if (it != cache.end())
    return it->second;
  else
    throw std::runtime_error("Resonance Strength not known for requested spin tune");
}


// get res. strength from cache or calculate it
std::complex<double> ResStrengthsData::get(double agamma)
{
  auto it = cache.find(agamma);
  if (it != cache.end())
    return it->second;
  else
    return calculate(agamma);
}


void ResStrengthsData::cacheIt(double agamma, const std::complex<double>& epsilon)
{
  cache.insert( std::pair<double,std::complex<double> >(agamma,epsilon) );
}


std::string ResStrengthsData::header(unsigned int w) const
{
  std::stringstream s;
  s << "#"<<std::setw(w)<<"agamma"<<std::setw(w)<< "real(epsilon)" <<std::setw(w)<< "imag(epsilon)" <<std::setw(w)<< "abs(epsilon)";
  return s.str();
}

std::string ResStrengthsData::printSingle(double agamma, std::complex<double> epsilon) const
{
  const unsigned int w=16;
  std::stringstream s;
  s <<resetiosflags(ios::scientific)<<setiosflags(ios::fixed)<<setprecision(4);
  s <<std::setw(1+w)<< agamma;
  s <<resetiosflags(ios::fixed)<<setiosflags(ios::scientific)<<showpoint<<setprecision(5);
  s <<std::setw(w)<< epsilon.real() <<std::setw(w)<< epsilon.imag() <<std::setw(w)<< std::abs(epsilon);
  return s.str();
}




void ResStrengths::init()
{
  if (queue.size() > 0)
    return;
  
  std::cout << "Estimate Resonance Strengths using "
	    << config->numTurns() << " turns for "
    	    << config->nParticles() << " particles:" << std::endl;
  for (unsigned int i=0; i<config->nParticles(); i++) {
    queue.emplace_back( ParticleResStrengths(i,config,lattice,orbit) );
  }
}

std::complex<double> ResStrengths::calculate(double agamma)
{
  std::stringstream msg;
  msg << "average over particles for gamma*a=" << agamma;
  polematrix::debug(__PRETTY_FUNCTION__, msg.str());
  
  std::complex<double> epsilon (0,0);
  for (auto &p : queue) {
    epsilon += p[agamma];
  }
  epsilon /= numParticles();
  cacheIt(agamma,epsilon);
  return epsilon;
}


void ResStrengths::start()
{
  // fill particle queue
  init();

  // set iterator to begin of queue
  queueIt = queue.begin();

  // write current config to file
  config->save( config->confOutFile().string() );

  auto start = std::chrono::high_resolution_clock::now();
  
  // start threads
  startThreads();
  
  waitForThreads();

  // finished: average ResStrengths
  std::cout << printErrors();
  for (double agamma=config->agammaMin(); agamma<=config->agammaMax(); agamma+=config->dagamma()) {
    calculate(agamma);
  }
  if (numSuccessful() > 0) {
    auto stop = std::chrono::high_resolution_clock::now();
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(stop-start);
    std::cout << std::endl
	      << "-----------------------------------------------------------------" << std::endl;
    std::cout << "Resonance Strengths estimated via "<<numSuccessful()<< " particles in ";
    std::cout << secs.count() << " s = "<< int(secs.count()/60.+0.5) << " min." << std::endl;
    std::cout << "Thanks for using polematrix " << polemversion() << std::endl;
    std::cout << "-----------------------------------------------------------------" << std::endl;
  }
  else
    std::cout << "Aborted due to ERRORs." << std::endl;
}


std::string ResStrengths::getSingle(double agamma)
{
  config->set_agammaMin(agamma);
  config->set_agammaMax(agamma);
  start();
  
  std::stringstream s;
  s << header() << std::endl << printSingle(agamma, operator[](agamma));
  return s.str();
}




// print all res. strength in cache
void ResStrengths::print(string filename)
{
  std::fstream file;
  std::stringstream s;
  const unsigned int w = 16;
  std::complex<double> epsilon;

 //metadata
  info.add("and polematrix version", polemversion());
  info.add("Description", "strengths of depolarizing resonances (complex numbers)");
  info.add("turns used for res. strength calc.", config->numTurns());
  info += lattice->info;
  s << info.out("#");

  // header
  s << header(w) << std::endl;

  // data
  for (auto& it : cache) {
    s << printSingle(it) << std::endl;
  }

 // output of s
 if (filename == "")
   cout << s.str();
 else {
   file.open(filename.c_str(), ios::out);
   if (!file.is_open()) {
     throw palatticeFileError(filename);
   }
   file << s.str();
   file.close();
   cout << "* Wrote " << filename  << std::endl;
 }

}





ParticleResStrengths::ParticleResStrengths(unsigned int id, const std::shared_ptr<Configuration> c,std::shared_ptr<const pal::AccLattice> l, std::shared_ptr<const pal::FunctionOfPos<pal::AccPair>> o)
  : SingleParticleSimulation(id,c)
{
  setModel(l,o);
}



// calculate res. strength freq. omega=agamma
// based on Courant-Ruth formalism using the magnetic fields B(orbit)
// Fields are NOT expressed by linear approx. of particle motion as by Courant-Ruth and DEPOL code,
// but the magnetic fields are used directly.
// !!! the orbit/field inside a magnet is assumed to be constant.
// !!! edge focusing is not included
std::complex<double> ParticleResStrengths::calculate(double agamma)
{
  std::stringstream msg;
  msg << "calculate gamma*a=" << agamma;
  polematrix::debug(__PRETTY_FUNCTION__, msg.str());
  
  std::complex<double> epsilon (0,0);

  for (unsigned int turn=0; turn<config->numTurns(); turn++) {
    for (AccLattice::const_iterator it=lattice->begin(); it!=lattice->end(); ++it) {
      double pos = it.pos() + turn*lattice->circumference();
      //field from Thomas-BMT equation:
      // omega = (1+agamma) * B_x - i * (1+a) * B_s
      // assume particle velocity parallel to s-axis, B is already normalized to rigidity (BR)_0 = p_0/e
      std::complex<double> omega = (1+agamma)*it.element()->B(trajectory->get(pos)).x - im * (1+config->a_gyro)*it.element()->B(trajectory->get(pos)).s;

      // dipole
      if (it.element()->type == dipole) {
	double R = ((Dipole*)it.element())->R(); // bending radius
	// calculate for dipole: epsilon = 1/2pi * omega * R/(i*agamma) * (e^{i*agamma*theta2}-e^{i*agamma*theta1})
	epsilon += 1/(2*M_PI) * omega * R/(im*agamma) * (std::exp(im*agamma*(lattice->theta(it.end())+turn*2*M_PI)) - std::exp(im*agamma*(lattice->theta(it.begin())+turn*2*M_PI)));
      }
      else {
	// calculate for all others: epsilon = 1/2pi * e^{i*agamma*theta} *  omega * l
	epsilon += 1/(2*M_PI) * std::exp(im*agamma*(lattice->theta(it.pos())+turn*2*M_PI)) * omega * it.element()->length;
      }
    }//lattice
  }//turn

  epsilon /= double(config->numTurns());
  cacheIt(agamma,epsilon);
  return epsilon;
}



void ParticleResStrengths::run()
{
  trajectory->init();

  for (double agamma=config->agammaMin(); agamma<=config->agammaMax(); agamma+=config->dagamma()) {
    calculate(agamma);
  }
  
  trajectory->clear();
}


