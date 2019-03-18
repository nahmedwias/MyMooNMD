#include <LoopInfo.h>
#include <MooNMD_Io.h>
#include <limits>
#include <cmath>

#ifdef _MPI
#include <mpi.h>
#endif

LoopInfo::LoopInfo(const std::string& name)
 : name(name), initial_residual(std::numeric_limits<double>::max()),
   old_residual(std::numeric_limits<double>::max()), timer(),
   n_previous_iterations(0)
{

}

LoopInfo::LoopInfo(const char* name)
 : LoopInfo(std::string(name))
{
  
}

LoopInfo::LoopInfo(const std::string& name,
                   bool print_time_every_step_in,
                   bool print_reduction_rates_in,
                   size_t verbosity_threshold_in)
  : LoopInfo(name)
{
  print_time_every_step = print_time_every_step_in;
  print_reduction_rates = print_reduction_rates_in;
  verbosity_threshold   = verbosity_threshold_in;
}

void LoopInfo::restart(const std::string& name, double initial_residual)
{
  this->name = name;
  this->initial_residual = initial_residual;
  this->old_residual = std::numeric_limits<double>::max();
  this->timer.reset();
}

void LoopInfo::print(unsigned int loop_index, double current_residual)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif

  using namespace std;
  std::stringstream s;
  
  s << this->name << " iteration: " << setw(3) << loop_index
    << "  residual: " << left << setprecision(10) << setw(15)
    << current_residual;
  if(loop_index > 0)
  {
    if(this->print_reduction_rates)
    {
      s << "  red/step: " << setprecision(10) << setw(15) 
        << current_residual/old_residual
        << "  red: " << setprecision(10) << setw(15)
        << current_residual/initial_residual;
    }
    if(this->print_time_every_step)
    {
      s << "  t[s]/step: " << setprecision(4) << setw(9)
        << timer.time_since_last_start()
        << "  t[s]: " << setprecision(4) << setw(9) 
        << timer.elapsed_time();
    }
  }
  else
  {
    // should we reset the initial_time here or rather not?
    this->restart(this->name, current_residual);
  }
  old_residual = current_residual;
  timer.stop();
  timer.start();
  // print only if verbosity is high enough.
  if(this->verbosity_threshold <= Output::getVerbosity() && my_rank == 0)
    Output::print<1>(s.str());
}


void LoopInfo::finish(unsigned int loop_index, double current_residual)
{
  this->n_previous_iterations += loop_index;
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif
  
  timer.stop();
  using namespace std;
  std::stringstream s;
  
  s << this->name << " iteration finished: " << setw(3) << loop_index
    << " residual: " << left << setprecision(10) << setw(15)
    << current_residual;
  if(loop_index > 0)
  {
    s << "  red: " << setprecision(10) << setw(15) 
      << current_residual/this->initial_residual 
      << "  red/step: " << setprecision(10) << setw(15)
      << std::pow(current_residual/this->initial_residual, 1./loop_index);
    s << "  t[s]: " << setprecision(4) << setw(9) 
      << timer.elapsed_time()
      << "  t[s]/step: " << setprecision(4) << setw(9)
      << this->timer.elapsed_time()/loop_index;
  }
  if(loop_index != this->n_previous_iterations)
  {
    s << "  total iterations: " << setw(3) <<this-> n_previous_iterations;
  }
  if(my_rank == 0)
    Output::print<1>(s.str());
}

double LoopInfo::get_initial_residual() const
{
  return this->initial_residual;
}

double LoopInfo::get_previous_residual() const
{
  return this->old_residual;
}

unsigned int LoopInfo::get_n_previous_iterations() const
{
  return n_previous_iterations;
}
