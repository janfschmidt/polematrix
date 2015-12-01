#include "RadiationModel.hpp"
#include "gsl/gsl_sf_synchrotron.h"

  // photon spectrum. used for probabilities of photon energies
double SynchrotronRadiationModel::nPhoton(double u_per_uc) const
{
  // This is an implementation of the integrated modified bessel function K_5/3
  // from V. O. Kostroun "Simple numerical Evaluation of modified bessel functions of fractional order [...]"
  // -----
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

  // Implementation from GSL
  return gsl_sf_synchrotron_1(u_per_uc)/u_per_uc;
}

SynchrotronRadiationModel::SynchrotronRadiationModel(int _seed) : seed(_seed), rng(seed)
{
  std::vector<double> intervals;
  std::vector<double> weights;
  for(double u=1e-10; u<=10.; u*=2) {
    intervals.push_back(u);
    weights.push_back(nPhoton(u));
    std::cout << u <<"\t"<< nPhoton(u) << std::endl;
  }
  // std::cout << intervals.size() <<" energy spectrum sampling points"<< std::endl;
  photonEnergy = boost::random::piecewise_linear_distribution<>(intervals.begin(), intervals.end(), weights.begin());
}

double SynchrotronRadiationModel::radiatedEnergy(const pal::AccElement* element, const double& gamma)
{
  double Erad = 0;
  boost::random::poisson_distribution<> photonsPerDipole(element->meanPhotons_syli(gamma));
  unsigned int n = photonsPerDipole(rng);
  for(auto photon=0u; photon<n; photon++) {
    Erad += photonEnergy(rng);
  }
  return std::move( Erad *  element->Ecrit_keV_syli(gamma) );
}



void LongitudinalPhaseSpaceModel::init(const pal::AccLattice* l)
{  
  //gamma & pos start from config
  lastPos = config.pos_start();
  lattice = l;
  nCavities = lattice->size(pal::cavity);
  set_gamma0(config.gamma_start());

  //sigma_phase -> bunch length
  boost::random::normal_distribution<> phaseDistribution(M_PI-std::asin(1/config.q()), 0.6);     //electron beam: stable phase on falling slope of sine
  //J_s (& wieder R aus lattice)
  boost::random::normal_distribution<> gammaDistribution(gamma0(), std::pow(gamma0(),2)*std::sqrt(3.84e-13/(1.994*11.)));
    //initial phase space coordinate
  boost::random::mt11213b initrng(seed);
  double tmp = gammaDistribution(initrng);
  updateCavityVoltage();
  _gamma =  tmp;
  phase = phaseDistribution(initrng);
  //std::cout << phase <<"\t"<< gamma() << std::endl;
}


void LongitudinalPhaseSpaceModel::update(const pal::AccElement* element, const double& pos)
{
  phase += 2*M_PI*(stepDistance(pos)/lattice->circumference()) * config.h() * (config.alphac()-std::pow(gamma(),-2)) * delta();
  
  if(element->type == pal::dipole) {
    double tmp = radModel.radiatedEnergy(element, gamma()) / E_rest_keV;
    //std::cout << tmp << std::endl;
    _gamma -= tmp;
  }

  else if(element->type == pal::cavity) {
    double tmp =  gammaU0() * std::sin(phase) / nCavities;
    // std::cout << "cavity: "<< tmp  <<"\t"<< phase<< std::endl;
    // std::cout << nCavities <<" cavities, U0="<< q()*dGamma_ref*E_rest_keV << " keV" << std::endl;
    _gamma += tmp;
  }
    
  lastPos = pos;
}
