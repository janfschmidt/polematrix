#include "gtest/gtest.h"
#include "RadiationModel.hpp"

#include <sstream>

TEST(StatisticsTests, Spectrum) {
  SynchrotronRadiationModel m(47891);

  unsigned int n = 100000;
  double u = 0.;
  for (auto i=0u; i<n; i++) {
    u += m.getPhotonEnergy();
  }
  u/=n;

  // average photon energy - from M.Sands "physics of electron storage rings" eq. (5.14)
  ASSERT_NEAR(8/(15*std::sqrt(3)),u,0.005); 
 
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
