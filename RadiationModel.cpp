/* Classes to model longitudinal phase space
 * implementation of gamma(t) with stochastic radiation of photons
 * used by gammaMode 'radiation'
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

#include "RadiationModel.hpp"
#include "gsl/gsl_sf_synchrotron.h"

// photon spectrum. used for probabilities of photon energies
double SynchrotronRadiationModel::nPhoton(double u_per_uc) const
{
  // integrated modified bessel function K_5/3 - Implementation from GSL
  return gsl_sf_synchrotron_1(u_per_uc)/u_per_uc;

  // --------------
  // This is an implementation of the integrated modified bessel function K_5/3
  // from V. O. Kostroun "Simple numerical Evaluation of modified bessel functions of fractional order [...]"
  // double term;
  // double h=0.4;
  // double result = std::exp(-u_per_uc)/2.0;
  // unsigned int r = 1;
  // do {
  //   term = std::exp(-u_per_uc*std::cosh(r*h)) * std::cosh(r*h*5./3.) / std::cosh(r*h);
  //   result += term;
  //   r++;
  // } while (term/result > 1e-5);
  // return result;
  // --------------
}


// initialization of photon energy distribution
SynchrotronRadiationModel::SynchrotronRadiationModel(int _seed) : seed(_seed), rng(seed)
{ 
  std::vector<double> intervals; // energies u, normalized to critical energy
  std::vector<double> weights;   // number of photons emitted at these energies
  for(double u=1e-20; u<=100.; u*=1.1) {
    intervals.push_back(u);
    weights.push_back(nPhoton(u));
  }
  //std::cout << intervals.size() <<" energy spectrum sampling points"<< std::endl;

  // photon energy distribution
  photonEnergy = boost::random::piecewise_linear_distribution<>(intervals.begin(), intervals.end(), weights.begin());
}


// energy radiated by particle passing given element IN UNITS OF GAMMA. entering with energy given by gammaIn
double SynchrotronRadiationModel::radiatedEnergy(const pal::AccElement* element, const double& gamma0, const double& gammaIn)
{
  double g = gammaIn; // current particle energy in units of gamma
  double grad = 0;    // radiated energy in units of gamma

  // poisson distribution for number of emitted photons
  boost::random::poisson_distribution<unsigned int> photonsPerDipole(element->syli_meanPhotons(g));
  unsigned int n = photonsPerDipole(rng);

  // for each photon: get radiated energy from photon energy distribution (normalized to critical energy)
  // critical energy has to be corrected by gamma0/gamma, because 1/R decreases with gamma (R prop. to gamma),
  // since the magnetic field is set for gamma0.
  for(auto photon=0u; photon<n; photon++) {
    double dg = photonEnergy(rng) * element->syli_Ecrit_gamma(g) * gamma0/g;
    grad += dg;
    g -= dg;
  }

  // return total radiated energy in units of gamma
  return std::move( grad );
}






void LongitudinalPhaseSpaceModel::init(std::shared_ptr<const pal::AccLattice> l)
{  
  //gamma0 & pos start from config
  lastPos = config->pos_start();
  lattice = l;
  nCavities = lattice->size(pal::cavity);
  set_gamma0(config->gamma_start());

  // if (config->verbose()) {
  //   std::cout <<"gamma0=" << gamma0() << ", sigma_gamma=" << sigma_gamma()
  // 	      << "\nref_phase=" << ref_phase() << ", sigma_phase=" << sigma_phase()
  // 	      << "\nq=" << config->q()<< ", syncr_freq=" << synchrotronFreq() << std::endl;
  // }

  // init statistical distributions:
  boost::random::normal_distribution<> phaseDistribution(ref_phase(), sigma_phase());
  boost::random::normal_distribution<> gammaDistribution(gamma0(), sigma_gamma());
  
  //initial phase space coordinate for this particle
  boost::random::mt11213b initrng(seed);
  _gamma =  gammaDistribution(initrng);
  _phase = phaseDistribution(initrng);
  
  // if (config->verbose()) {
  //   std::cout << "longitudinal phase space initialized with "
  // 	      << "phase="<< phase() <<", gamma="<< gamma() << std::endl;
  // }
}


void LongitudinalPhaseSpaceModel::update(const pal::AccElement* element, const double& pos, const double& newGamma0)
{
  if(element->type == pal::dipole) {
    // phase change from momentum compaction (1st + 2nd order!)
    // calculate for whole turn and use percentage of bent length
    _phase += 2*M_PI * config->h() * (config->alphac() + config->alphac2()*delta()) * delta()  * (element->length/lattice->bentLength());
    // energy loss in dipole: radiate
    _gamma -= radModel.radiatedEnergy(element, gamma0(), gamma());
  }
  else if(element->type == pal::cavity) {
    set_gamma0(newGamma0); // update reference energy (energy ramp)
    // energy gain in cavity
    double tmp =  gammaU0()/nCavities * std::sin(phase());
    _gamma += tmp;
  }
    
  lastPos = pos;
}

// bunch length, calculated as time, converted to rf-phase (factor 2pi cancels)
double LongitudinalPhaseSpaceModel::sigma_phase() const
{
  return config->sigmaPhaseFactor() * config->alphac()/synchrotronFreq() * sigma_gamma()/gamma0()
    * config->h()*GSL_CONST_MKSA_SPEED_OF_LIGHT/lattice->circumference();
}

// energy spread in units of gamma
double LongitudinalPhaseSpaceModel::sigma_gamma() const
{
  return config->sigmaGammaFactor() * std::pow(gamma0(),2) * std::sqrt( 3.84e-13/(config->Js()*config->R()) );
}

// synchrotron frequency in Hz
double LongitudinalPhaseSpaceModel::synchrotronFreq_formula(const double& gammaIn) const
{
  return GSL_CONST_MKSA_SPEED_OF_LIGHT/lattice->circumference()
    * std::sqrt( -U0_keV()*config->h() / (2*M_PI*gammaIn*config->E_rest_keV)
		 * std::cos(ref_phase())*config->alphac() );
}
