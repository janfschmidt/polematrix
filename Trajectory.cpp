#include "Trajectory.hpp"

SimtoolTrajectory::SimtoolTrajectory(unsigned int id, const std::shared_ptr<Configuration> c)
  : Trajectory(id,c), simtoolTrajectory(config->getSimToolInstance(), gsl_interp_akima) {}


void SimtoolTrajectory::init()
{
  // simtool: sdds import thread safe since SDDSToolKit-devel-3.3.1-2
  simtoolTrajectory.simToolTrajectory( config->getSimToolInstance(), particleId+1 );
  //trajectorySimTool.init();  --> not used to improve performance, trajectory is accessed by infrontof()
}


void SimtoolTrajectory::saveSimtoolData()
{
  if ( !config->saveGamma(particleId) )
    return;
  
  simtoolTrajectory.info.add("polematrix particle ID", particleId);
  std::stringstream file;
  file <<std::setw(4)<<std::setfill('0')<< particleId << ".dat";
  simtoolTrajectory.print( (config->outpath()/"trajectorySimtool_").string() + file.str() );
}
