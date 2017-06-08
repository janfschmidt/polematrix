/* Trajectory Classes
 * It represents the trajectory of a single particle, which can be
 * an Elegant/MadX trajectory or just a closed orbit
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

#ifndef __POLEMATRIX__TRAJECTORY_HPP_
#define __POLEMATRIX__TRAJECTORY_HPP_

#include <memory>
#include <libpalattice/FunctionOfPos.hpp>
#include "Configuration.hpp"


class Trajectory
{
public:
  const unsigned int particleId;
  std::shared_ptr<const pal::FunctionOfPos<pal::AccPair>> orbit;

protected:
  const std::shared_ptr<Configuration> config;
  bool initDone;

public:
  Trajectory(unsigned int id, const std::shared_ptr<Configuration> c) : particleId(id), config(c), initDone(false) {}
  virtual ~Trajectory() {}

  void setOrbit(std::shared_ptr<const pal::FunctionOfPos<pal::AccPair>> o) {orbit = o;}
  
  virtual pal::AccPair get(const double& pos) =0;
  
  void init();
  virtual void saveSimtoolData() {}
  virtual void clear() {} // clear trajectory data

protected:
  virtual void initImplementation() =0;
};



class Orbit : public Trajectory
{
public:
  Orbit(unsigned int id, const std::shared_ptr<Configuration> c) : Trajectory(id,c) {}
    virtual ~Orbit() {}
  virtual pal::AccPair get(const double& pos) {return orbit->interp( orbit->posInTurn(pos) );}

protected:
  virtual void initImplementation() {}
};



class SimtoolTrajectory : public Trajectory
{
private:
  pal::FunctionOfPos<pal::AccPair> simtoolTrajectory;
  
public:
  SimtoolTrajectory(unsigned int id, const std::shared_ptr<Configuration> c);
  virtual ~SimtoolTrajectory() {}
  virtual pal::AccPair get(const double& pos) {return simtoolTrajectory.interpPeriodic(pos-config->pos_start());}
  virtual void clear() {simtoolTrajectory.clear();}
  virtual void saveSimtoolData();

  protected:
  virtual void initImplementation();
};


class Oscillation : public Trajectory
{
private:
  pal::FunctionOfPos<pal::AccPair> beta;  // beta function (twiss)
  AccPair emittance;                      // single particle emittance
  AccPair phase0;                         // phase(0) start value
  AccPair freq;                           // oscillation frequency / 1/m
  
public:
  Oscillation(unsigned int id, const std::shared_ptr<Configuration> c);
  virtual ~Oscillation() {}
  virtual pal::AccPair get(const double& pos);

protected:
  virtual void initImplementation();
};


#endif
// __POLEMATRIX__TRAJECTORY_HPP_
