#ifndef LOOPINFO_H
#define LOOPINFO_H

#include <Chrono.h>
#include <string>

class LoopInfo
{
public:
  /// @brief Constructor. This also sets the initial_time
  LoopInfo(std::string name);
  /// @brief Constructor for convenience, this calls the other constructor
  LoopInfo(const char* name);
  
  /// @brief use this object again in another iteration. This resets the 
  /// initial_time as well.
  /// @details You can call this even if the loop has never started.
  void restart(std::string name, double initial_residual);
  
  /// @brief write out some information during the loop
  void print(unsigned int loop_index, double current_residual);
  /// @brief write out some information after the loop has finished
  void finish(unsigned int loop_index, double current_residual);
  /// @brief print the time for each iteration within LoopInfo::print.
  bool print_time_every_step = false;
  /// @brief print the reduction factors for each (except the first) iteration
  /// within LoopInfo::print
  bool print_reduction_rates = true;
  /// @brief the verbosity determines if anything is printed at all during 
  /// LoopInfo::print. Small verbosity_threshold means much output
  unsigned int verbosity_threshold = 3;
  /// @brief return the initial residual. This can be used to compute the 
  /// reduction of the residual during the loop.
  double get_initial_residual() const;
  /// @brief return the residual from the previous iteration. This can be used
  /// to check for divergence.
  double get_previous_residual() const;
  
private:
  /// @brief string to be preceded during output
  std::string name;
  /// @brief the residual before the iteration starts
  double initial_residual;
  /// @brief the residual of the previous iteration
  double old_residual;
  /// @brief the time before the iteration starts
  Chrono initial_time;
  /// @brief the time at the end of the last iteration
  Chrono old_time;
};

#endif // LOOPINFO_H