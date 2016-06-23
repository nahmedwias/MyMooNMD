#include <LoopInfo.h>
#include <MooNMD_Io.h>
#include <limits>
#include <cmath>

LoopInfo::LoopInfo(std::string name)
 : name(name), initial_residual(std::numeric_limits<double>::max()),
   old_residual(std::numeric_limits<double>::max()), initial_time(), old_time()
{

}

LoopInfo::LoopInfo(const char* name)
 : LoopInfo(std::string(name))
{
  
}

void LoopInfo::restart(std::string name, double initial_residual)
{
  this->name = name;
  this->initial_residual = initial_residual;
  this->old_residual = std::numeric_limits<double>::max();
  this->initial_time.reset();
  this->old_time.reset();
}

void LoopInfo::print(unsigned int loop_index, double current_residual)
{
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
        << old_time.elapsed_time()
        << "  t[s]: " << setprecision(4) << setw(9) 
        << initial_time.elapsed_time();
    }
  }
  else
  {
    // should we reset the initial_time here or rather not?
    this->restart(this->name, current_residual);
  }
  old_residual = current_residual;
  old_time.reset();
  // print only if verbosity is high enough.
  if(this->verbosity_threshold <= Output::getVerbosity())
    Output::print<1>(s.str());
}


void LoopInfo::finish(unsigned int loop_index, double current_residual)
{
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
      << initial_time.elapsed_time()
      << "  t[s]/step: " << setprecision(4) << setw(9)
      << this->initial_time.elapsed_time()/loop_index;
  }
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
