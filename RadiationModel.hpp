#include <vector>
#include <cmath>
#include <memory>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/piecewise_linear_distribution.hpp>
#include "libpalattice/AccLattice.hpp"
#include "Configuration.hpp"


class SynchrotronRadiationModel {
protected:
  int seed;
  boost::random::mt11213b rng;
  boost::random::piecewise_linear_distribution<> photonEnergy;

  // photon spectrum. used for probabilities of photon energies
  double nPhoton(double u_per_uc) const;

public:
  SynchrotronRadiationModel(int _seed=1);
  int getSeed() const {return seed;}
  // photon energy radiated within this element (e.g. Dipole) by an electron entering with energy gammaIn
  // at a reference beam energy given by gamma0.
  // returns energy in units of gamma
  double radiatedEnergy(const pal::AccElement* element, const double& gamma0, const double& gammaIn);
};


class LongitudinalPhaseSpaceModel {
protected:
  int seed;
  SynchrotronRadiationModel radModel;        // stochastical model for radiation
  std::shared_ptr<const pal::AccLattice> lattice;
  const std::shared_ptr<const Configuration> config;
  unsigned int nCavities;
  double _gamma0;  // reference energy in units of gamma
  double _gammaU0; // cavity amplitude U0 in units of gamma
  
  double _phase;   //current synchrotron phase
  double _gamma;   //current energy
  double lastPos;  //total distance currently traveled in m (to calc distance since last step)

  void updateCavityVoltage() {_gammaU0 = U0_keV() / config->E_rest_keV;}
  double synchrotronFreq_formula(const double& gammaIn) const;


public:
  LongitudinalPhaseSpaceModel(int _seed, std::shared_ptr<const Configuration> c)
    : seed(_seed),radModel(seed), config(c) {lastPos=_phase=_gamma=_gamma0=_gammaU0=0;}
  double gammaU0() const {return _gammaU0;}
  double gamma0() const {return _gamma0;}
  double phase() const {return _phase;}
  double gamma() const {return _gamma;}
  
  void set_gamma0(double x) {_gamma0=x; updateCavityVoltage();}
  
  double stepDistance(const double& pos) const {return pos - lastPos;}
  double delta() const {return (gamma()-gamma0())/gamma0();}
  double gammaMinusGamma0() const {return gamma()-gamma0();}

  void init(std::shared_ptr<const pal::AccLattice> l);
  void update(const pal::AccElement* element, const double& pos, const double& newGamma0);

  //cavity voltage in keV
  double U0_keV() const {return config->q() * lattice->Erev_keV_syli(gamma0());}

  //reference phase -> stable phase for given overvoltage q
  // (! electron beam: stable phase on falling slope of sine !)
  double ref_phase() const {return M_PI - std::asin( 1/config->q() );}
  
  double sigma_phase() const;     //bunch length as phase in units of radian
  double sigma_gamma() const;     //energy spread in units of gamma
  
  //synchrotron frequency in Hz (for reference energy):
  double synchrotronFreq() const {return synchrotronFreq_formula(gamma0());}
  //synchrotron frequency  in Hz for current energy of this particle:
  double synchrotronFreq_current() const {return synchrotronFreq_formula(gamma());}

};
