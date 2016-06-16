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
  for(double u=1e-20; u<=100.; u*=2) {  // ??? what is a reasonable range ???
    intervals.push_back(u);
    weights.push_back(nPhoton(u));
  }
  // std::cout << intervals.size() <<" energy spectrum sampling points"<< std::endl;

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
  // critical energy is calculated for reference energy gamma0, corrected here by (gamma/gamma0)^2:
  // ()^3 due to gamma^3 and ()^-1 due to 1/R also depending on energy
  for(auto photon=0u; photon<n; photon++) {
    double dg = photonEnergy(rng) * element->syli_Ecrit_gamma(gamma0) * std::pow(g/gamma0, 2);
    grad += dg;
    g -= dg;
  }

  // return total radiated energy in units of gamma
  return std::move( grad );
}






void LongitudinalPhaseSpaceModel::init(const pal::AccLattice* l)
{  
  //gamma0 & pos start from config
  lastPos = config.pos_start();
  lattice = l;
  nCavities = lattice->size(pal::cavity);
  set_gamma0(config.gamma_start());

  // init statistical distributions:
  boost::random::normal_distribution<> phaseDistribution(ref_phase(), sigma_phase());
  boost::random::normal_distribution<> gammaDistribution(gamma0(), sigma_gamma());
  // std::cout << "ref_phase=" << ref_phase() << "\nsigma_phase=" << sigma_phase()
  // 	    << "\nsigma_gamma=" << sigma_gamma() << "\nfreq=" << synchrotronFreq() << std::endl;
  
  //initial phase space coordinate for this particle
  boost::random::mt11213b initrng(seed);
  _gamma =  gammaDistribution(initrng);
  _phase = phaseDistribution(initrng);
  //std::cout << phase() <<"\t"<< gamma() << std::endl;
}


void LongitudinalPhaseSpaceModel::update(const pal::AccElement* element, const double& pos)
{
  // phase change from momentum compaction (1st + 2nd order!)
  // calculate for whole turn and use percentage (by distance s)
  _phase += 2*M_PI*(stepDistance(pos)/lattice->circumference()) * config.h() * (config.alphac() + config.alphac2()*delta()) * delta();


  // energy loss in dipole: radiate
  if(element->type == pal::dipole) {
    _gamma -= radModel.radiatedEnergy(element, gamma0(), gamma());
  }
  // energy gain in cavity
  else if(element->type == pal::cavity) {
    double tmp =  gammaU0()/nCavities * std::sin(phase());
    // std::cout << "cavity: "<< tmp  <<"\t"<< phase()<< std::endl;
    // std::cout << nCavities <<" cavities, U0="<< q()*dGamma_ref*E_rest_keV << " keV" << std::endl;
    _gamma += tmp;
  }
    
  lastPos = pos;
}

// bunch length, calculated as time, converted to rf-phase (factor 2pi cancels)
double LongitudinalPhaseSpaceModel::sigma_phase() const
{
  return config.alphac()/synchrotronFreq() * sigma_gamma()/gamma0()
    * config.h()*GSL_CONST_MKSA_SPEED_OF_LIGHT/lattice->circumference();
}

// energy spread in units of gamma
double LongitudinalPhaseSpaceModel::sigma_gamma() const
{
  return std::pow(gamma0(),2) * std::sqrt( 3.84e-13/(config.Js()*config.R()) );
}

// synchrotron frequency in Hz
double LongitudinalPhaseSpaceModel::synchrotronFreq() const
{
  return GSL_CONST_MKSA_SPEED_OF_LIGHT/lattice->circumference()
    * std::sqrt( -U0_keV()*config.h() / (2*M_PI*gamma0()*config.E_rest_keV)
		 * std::cos(ref_phase())*config.alphac() );
}
