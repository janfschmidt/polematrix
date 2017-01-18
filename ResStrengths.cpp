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


void ResStrengthsData::cacheIt(double agamma, std::complex<double>& epsilon)
{
  cache.insert( std::pair<double,std::complex<double> >(agamma,epsilon) );
}






std::complex<double> ResStrengths::calculate(double agamma)
{
  std::cout << "average over particles for gamma*a=" << agamma << std::endl;
  std::complex<double> epsilon (0,0);
  for (auto &p : particles) {
    epsilon += p[agamma];
  }
  epsilon /= numParticles();
  cacheIt(agamma,epsilon);
  return epsilon;
}


void ResStrengths::start()
{
  
  std::cout << "Estimate Resonance Strengths using "
    	    << config->nParticles() << " particles and "
	    << nTurns() << " turns:" << std::endl;
  for (unsigned int i=0; i<config->nParticles(); i++) {
    particles.emplace_back( ParticleResStrengths(i,config,lattice,orbit) );
  }
  for (auto& p : particles) {
    p.run();
  }
  for (double agamma=config->agammaMin(); agamma<=config->agammaMax(); agamma+=spintuneStep()) {
    calculate(agamma);
  }
  std::cout << "Done." << std::endl;
}




// print all res. strength in cache
void ResStrengths::print(string filename)
{
  std::fstream file;
  std::stringstream s;
  const int w = 16;
  std::complex<double> epsilon;

 //metadata
  info.add("and polematrix version", polemversion());
  info.add("Description", "strengths of depolarizing resonances (complex numbers)");
  info.add("turns used for res. strength calc.", nTurns());
  info += lattice->info;
  s << info.out("#");

 // header
  s <<"#"<<std::setw(w)<<"agamma"<<std::setw(w)<< "real(epsilon)" <<std::setw(w)<< "imag(epsilon)" <<std::setw(w)<< "abs(epsilon)" << std::endl;
 
  for (auto it : cache) {
    auto epsilon = it.second;
    s <<resetiosflags(ios::scientific)<<setiosflags(ios::fixed)<<setprecision(4);
    s <<std::setw(1+w)<< it.first;
    s <<resetiosflags(ios::fixed)<<setiosflags(ios::scientific)<<showpoint<<setprecision(5);
    s <<std::setw(w)<< epsilon.real() <<std::setw(w)<< epsilon.imag() <<std::setw(w)<< std::abs(epsilon)
      << std::endl;
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
// !!! At the moment the orbit/field inside a magnet is assumed to be constant.
// !!! edge focusing of dipoles is not included
std::complex<double> ParticleResStrengths::calculate(double agamma)
{
  std::cout << "calculate gamma*a=" << agamma << std::endl;
  std::complex<double> epsilon (0,0);

  for (unsigned int turn=0; turn<nTurns(); turn++) {
    for (AccLattice::const_iterator it=lattice->begin(); it!=lattice->end(); ++it) {
      double pos = it.pos() + turn*lattice->circumference();
      //field from Thomas-BMT equation:
      // omega = (1+agamma) * B_x - i * (1+a) * B_s
      // assume particle velocity parallel to s-axis, B is already normalized to rigidity (BR)_0 = p_0/e
      std::complex<double> omega = (1+agamma)*it.element()->B(trajectory->get(pos)).x - im * (1+config->a_gyro)*it.element()->B(trajectory->get(pos)).s;

      // dipole
      if (it.element()->type == dipole) {
	double R = ((Dipole*)it.element())->R(); // bending radius
	// dipole: epsilon = 1/2pi * omega * R/(i*agamma) * (e^{i*agamma*theta2}-e^{i*agamma*theta1})
	epsilon += 1/(2*M_PI) * omega * R/(im*agamma) * (std::exp(im*agamma*(lattice->theta(it.end())+turn*2*M_PI)) - std::exp(im*agamma*(lattice->theta(it.begin())+turn*2*M_PI)));
      }
      // all others: epsilon = 1/2pi * e^{i*agamma*theta} *  omega * l
      else {
	epsilon += 1/(2*M_PI) * std::exp(im*agamma*(lattice->theta(it.pos())+turn*2*M_PI)) * omega * it.element()->length;
      }
    }//lattice
  }//turn

  epsilon /= double(nTurns());
  cacheIt(agamma,epsilon);
  return epsilon;
}



void ParticleResStrengths::run()
{
  trajectory->init();
  polematrix::debug(__PRETTY_FUNCTION__, "init done");

  for (double agamma=config->agammaMin(); agamma<=config->agammaMax(); agamma+=spintuneStep()) {
    calculate(agamma);
  }
  
  trajectory->clear();
}
