#include "gtest/gtest.h"
#include "RadiationModel.hpp"

#include <sstream>
#include <fstream>
#include <map>
#include <random>
#include <gsl/gsl_spline.h>


class PhotonEnergy : public ::testing::Test {
private:
  gsl_interp_accel *acc;
  gsl_spline *spline;
  unsigned int n_ele;
  
public:
  SynchrotronRadiationModel m;
  std::vector<double> u_ele, wcum_ele;
  const unsigned int Nstat;

  PhotonEnergy() : n_ele(200), m(47891), Nstat(1e7)
  {
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_akima, n_ele);
  }

  ~PhotonEnergy()
  {
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }

  void init_ele()
  {
    std::ifstream ele;
    ele.open("elegant_photonEnergy.dat");
    if (!ele.is_open())
      throw std::runtime_error("Error opening elegant_photonEnergy.dat");
    
    while (ele.good()) {
      double tmp=0.;
      for (auto v : {&u_ele, &wcum_ele}) {
	ele >> tmp;
	if (ele.eof()) break;
	v->push_back(tmp);
      }
    }
    ele.close();
    if (u_ele.size()!=n_ele || wcum_ele.size()!=n_ele)
      throw std::runtime_error("Error reading elegant_photonEnergy.dat");

    gsl_spline_init (spline, &wcum_ele[0], &u_ele[0], n_ele);
  }

  std::map<double,double> getBoostDist()
  {
    std::map<double,double> dist;
    
    for (auto& u : u_ele)
      dist.emplace(u,0);
    
    for (auto i=0u; i<Nstat; i++) {
      double tmp = m.getPhotonEnergy();
      auto it = dist.lower_bound(tmp);
      it->second+=1.;
    }
    return std::move(dist);
  }

  // get random photon energy from uniform random number [0,1] by interpolation of u_ele(wcum_ele)
  double interp_ele(double rn) const
  {
    return gsl_spline_eval (spline, rn, acc);
  }

};





TEST(PhotonNumber, Mean) {
  auto l = 2.875;
  auto R = 10.98;
  pal::AccTriple k0; k0.z = 1/R;
  pal::Dipole d("M2", l, k0);
  auto gamma = 4599;
  
  boost::random::poisson_distribution<unsigned int> photonsPerDipole(d.syli_meanPhotons(gamma));
  boost::random::mt11213b rng(47891);
  
  unsigned int n = 100000;
  double p = 0.;
  for (auto i=0u; i<n; i++) {
    p += photonsPerDipole(rng);
  }
  p/=n;

  ASSERT_NEAR(5/(2*std::sqrt(3))*gamma/137.*l/R, p, 0.001);
}





TEST_F(PhotonEnergy, Mean) {
  double u = 0.;
  for (auto i=0u; i<Nstat; i++) {
    u += m.getPhotonEnergy();
  }
  u/=Nstat;
  std::cout <<"avg(u)="<< u << std::endl;

  // average photon energy - from M.Sands "physics of electron storage rings" eq. (5.14)
  ASSERT_NEAR(8/(15*std::sqrt(3)), u, 0.005); 
}



// calculate cumulative probabilities from Bessel function (m.nPhoton())
// and compare with elegant
TEST_F(PhotonEnergy, ElegantVsBessel) {
  init_ele();
  
  std::ofstream f;
  f.open("photonEnergyBessel.dat");
  std::vector<double> wcum_polem;
  wcum_polem.push_back(0.);
  for (auto i=1u; i<u_ele.size(); i++) {
    wcum_polem.push_back( wcum_polem.back() + m.nPhoton((u_ele[i]+u_ele[i-1])/2.) * (u_ele[i]-u_ele[i-1]) );
  }

  for (auto i=1u; i<u_ele.size(); i++) {
    wcum_polem[i] /= wcum_polem.back();
    f <<std::scientific<<std::setprecision(15)<< u_ele[i] <<"\t"<< wcum_ele[i] <<"\t"<< wcum_polem[i] << std::endl; 
    EXPECT_NEAR(wcum_ele[i], wcum_polem[i], 0.001);
  }
  f.close();
}


// calculate cumulative probabilities from boost::random::piecewise_linear_distribution (m.getPhotonEnergy())
// and compare with elegant
TEST_F(PhotonEnergy, ElegantVsDistCum) {
  init_ele();
  
  // get energy probabilities from random number distribution
  auto dist = getBoostDist();
  
  std::vector<double> wcum_polem;
  wcum_polem.push_back(0.);
  for (auto i=1u; i<u_ele.size(); i++) {
    wcum_polem.push_back( wcum_polem.back() + dist.at(u_ele[i]) );
  }

  std::ofstream f;
  f.open("photonEnergyCum.dat");
  for (auto i=1u; i<u_ele.size(); i++) {
    wcum_polem[i] /= wcum_polem.back();
    f <<std::scientific<<std::setprecision(15)<< u_ele[i] <<"\t"<< wcum_ele[i] <<"\t"<< wcum_polem[i] << std::endl; 
    EXPECT_NEAR(wcum_ele[i], wcum_polem[i], 0.001);
  }
  f.close();
}


// calculate probability distribution from elegant
// and compare with boost::random::piecewise_linear_distribution
TEST_F(PhotonEnergy, ElegantVsDistNoncum) {

  // get elegant probabilities from cumulative data
  init_ele();
  std::map<double,double> dist_ele;
  
  for (auto& u : u_ele)
    dist_ele.emplace(u,0);
  
  std::default_random_engine rng;
  std::uniform_real_distribution<double> uniform(0.0,1.0);
  for (auto i=0u; i<Nstat; i++) {
    double tmp = interp_ele(uniform(rng));
    
    auto it = dist_ele.lower_bound(tmp);
    it->second+=1./Nstat;
  }
  
  // get energy probabilities from random number distribution
  auto dist = getBoostDist();
  for (auto &it : dist) it.second /= Nstat;
  
   std::ofstream f;
  f.open("photonEnergyDist.dat");
  for (auto i=1u; i<u_ele.size(); i++) {
    f <<std::scientific<<std::setprecision(15)<< u_ele[i] <<"\t"<< dist_ele.at(u_ele[i]) <<"\t"<< dist.at(u_ele[i]) << std::endl; 
    EXPECT_NEAR(dist.at(u_ele[i]),dist_ele.at(u_ele[i]), 0.001);
  }
  f.close();
}




int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
