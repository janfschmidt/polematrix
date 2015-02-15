#ifndef __POLEMATRIX__TRACKINGTASK_HPP_
#define __POLEMATRIX__TRACKINGTASK_HPP_

#include <fstream>
#include <stdexcept>
#include <map>
#include <memory>
#define ARMA_NO_DEBUG
#include <armadillo>
#include <libpal/types.hpp>
#include "Configuration.hpp"

class TrackingTask
{
public:
  const unsigned int particleId;
  const Configuration *config;
  
  TrackingTask(unsigned int id, const Configuration *c);
  TrackingTask(const TrackingTask& other) = delete;
  TrackingTask(TrackingTask&& other) = default;
  ~TrackingTask() {}
  
  void run();                                 //run tracking task
  void matrixTracking();
  inline arma::mat33 rotxMatrix(double angle) const;
  inline arma::mat33 rotMatrix(pal::AccTriple B) const;
  std::string outfileName() const;            // output file name
  
private:
  std::map<double,arma::colvec3> storage;     // store results
  std::unique_ptr<std::ofstream> outfile;     // output file via pointer, std::ofstream not moveable in gcc 4.9
  unsigned int w;                             // output file column width
  void outfileOpen();                         // open output file and write header
  void outfileClose();                        // write footer and close output file
  void outfileAdd(const double &t, const arma::colvec3 &s); // append s(t) to outfile
  void storeStep(const double &t, const arma::colvec3 &s);  // append s(t) to storage and outfile
};


// exceptions
class TrackError : public std::runtime_error {
public:
  TrackError(std::string msg) : std::runtime_error(msg) {}
};

class TrackFileError : public TrackError {
public:
  TrackFileError(std::string file) : TrackError("Cannot open "+file) {}
};



#endif
// __POLEMATRIX__TRACKINGTASK_HPP_
