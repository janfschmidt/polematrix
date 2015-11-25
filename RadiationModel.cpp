#include "RadiationModel.hpp"

double SynchrotronRadiationModel::nPhoton(double u_per_uc) const
{
  double term;
  double h=0.4;
  double result = std::exp(-u_per_uc)/2.0;
  unsigned int r = 1;
  do {
    term = std::exp(-u_per_uc*std::cosh(r*h)) * std::cosh(r*h*5./3.) / std::cosh(r*h);
    result += term;
    r++;
  } while (term/result > 1e-5);
  return result;
}

SynchrotronRadiationModel::SynchrotronRadiationModel(int _seed) : seed(_seed), rng(seed)
{
  std::vector<double> intervals;
  std::vector<double> weights;
  for(double u=1e-10; u<=10.; u*=2) {
    intervals.push_back(u);
    weights.push_back(nPhoton(u));
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

