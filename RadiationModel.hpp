#include <vector>
#include <cmath>
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
  // photon energy radiated within this element (e.g. Dipole) in keV
  // at electron beam energy given by gamma
  double radiatedEnergy(const pal::AccElement* element, const double& gamma);
};


class LongitudinalPhaseSpaceModel {
protected:
  int seed;
  SynchrotronRadiationModel radModel;        // stochastical model for radiation
  const pal::AccLattice* lattice;
  Configuration &config;
  unsigned int nCavities;
  double _gamma0;  // reference energy in units of gamma
  double _gammaU0; // cavity amplitude U0 in units of gamma
  
  double phase;    //current synchrotron phase
  double _gamma;   //current energy
  double lastPos;  //total distance currently traveled in m (to calc distance since last step)

  void updateCavityVoltage()
  {
    _gammaU0 = config.q() * lattice->Erev_keV_syli(gamma0()) / E_rest_keV;
  }


public:
  const double E_rest_keV;
  
  LongitudinalPhaseSpaceModel(int _seed, Configuration &c)
    : seed(_seed),radModel(seed),lattice(NULL), config(c), E_rest_keV(511) {lastPos=phase=_gamma=_gamma0=_gammaU0=0;}
  double gammaU0() const {return _gammaU0;}
  double gamma() const {return _gamma;}
  double gamma0() const {return _gamma0;}
  
  
  void set_gamma0(double x) {_gamma0=x; updateCavityVoltage();}
  
  double stepDistance(const double& pos) const {return pos - lastPos;}
  double delta() const {return (gamma()-gamma0())/gamma0();}

  
  void init(const pal::AccLattice* l);
  void update(const pal::AccElement* element, const double& pos);
};
