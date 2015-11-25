#include <vector>
#include <cmath>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/piecewise_linear_distribution.hpp>
#include "libpalattice/AccLattice.hpp"
#include "Configuration.hpp"


class SynchrotronRadiationModel {
protected:
  int seed;
  boost::random::mt11213b rng;
  boost::random::piecewise_linear_distribution<> photonEnergy;
  //boost::random::poisson_distribution photonsPerDipole;

  double nPhoton(double u_per_uc) const;

public:
  SynchrotronRadiationModel(int _seed=1);
  int getSeed() const {return seed;}
  // photon energy radiated within this element (e.g. Dipole) in keV
  // at electron beam energy given by gamma
  double radiatedEnergy(const pal::AccElement* element, const double& gamma);
};


class longitudinalPhaseSpaceModel {
protected:
  SynchrotronRadiationModel radModel;        // stochastical model for radiation
  const pal::AccLattice* lattice;
  unsigned int nCavities;
  double _q;       //overvoltage factor
  double _ac;      //momentum compaction factor
  unsigned int _h; //harmonic number
  double _gamma_0; //center energy
  double dGamma_ref; //energy gain per turn for reference phase (= radiation energy loss per turn), in unit of gamma
  
  double phase;    //current synchrotron phase
  double _gamma;   //current energy
  double lastPos;  //total distance currently traveled in m (to calc distance since last step)

public:
  const double E_rest_keV;
  
  longitudinalPhaseSpaceModel(int seed)
    : radModel(seed),lattice(NULL), E_rest_keV(511) {lastPos=_q=_ac=_h=phase=_gamma=dGamma_ref=0;}
  double q() const {return _q;}
  double ac() const {return _ac;}
  unsigned int h() const {return _h;}
  double gamma() const {return _gamma;}
  double gammaCentral() const {return _gamma_0;}
  
  void set_q(double x) {_q=x;}
  void set_ac(double x) {_ac=x;}
  void set_h(unsigned int x) {_h=x;}

  double stepDistance(const double& pos) const {return pos - lastPos;}
  double delta() const {return (gamma()-gammaCentral())/gammaCentral();}

  void init(const pal::AccLattice* l, const Configuration& config)
  {
    //gamma & pos start from config
    lastPos = config.pos_start();
    _gamma = config.gamma_start();
    lattice = l;
    nCavities = lattice->size(pal::cavity);
    set_gammaCentral(gamma());
    //load from config: q,ac,h (oder SimTool? oder config entscheidet woher?)
    _q=10; _ac=0.0601; _h=274;
    phase = M_PI - std::asin(1/q()); //electron beam: stable phase on falling slope of sine
  }
  
  void set_gammaCentral(const double& gamma)
  {
    _gamma_0 = gamma;
    dGamma_ref = lattice->Erev_keV_syli(gamma) / E_rest_keV;
  }

  
  void update(const pal::AccElement* element, const double& pos)
  {
    phase += (stepDistance(pos)/lattice->circumference()) * h() * (ac()-std::pow(gamma(),-2)) * delta();

    if(element->type == pal::dipole) {
      double tmp = radModel.radiatedEnergy(element, gamma()) / E_rest_keV;
      //std::cout << "dipole: "<< tmp <<"\t"<< phase << std::endl;
      _gamma -= tmp;
    }

    else if(element->type == pal::cavity) {
      double tmp = q()*dGamma_ref * std::sin(phase) / nCavities;
      // std::cout << "cavity: "<< tmp  <<"\t"<< phase<< std::endl;
      // std::cout << nCavities <<" cavities, U0="<< q()*dGamma_ref*E_rest_keV << " keV" << std::endl;
      _gamma += tmp;
    }
    
    lastPos = pos;
  }
};
