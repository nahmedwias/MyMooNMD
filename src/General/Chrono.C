/**
 * @file Chrono.C
 * Implementation of class Chrono declared in Chrono.h.
 *
 * @date 2016/04/18
 */
#include <Chrono.h>
#include <MooNMD_Io.h>


Chrono::Chrono()
{
  reset();
}

void Chrono::print_time(std::string program_part)
{
#ifdef _MPI
  Output::print("WARNING -- No MPI time printing enabled yet.");
#elif _OMP
  Output::print("WARNING -- No OpenMP time printing enabled yet.")
#else //sequential execution case
  double time = get_exec_time_s() - start_time;
  Output::print("--- time for ", program_part,": ", time, " s"  );
#endif
}

double Chrono::get_exec_time_s()
{
  rusage usage;

  if(getrusage(RUSAGE_SELF, &usage) == -1)
    ErrThrow("Error in GetTime!");

  // user mode time and system mode time in s
  double user_mode_time = usage.ru_utime.tv_sec;
  double kernel_mode_time =  usage.ru_stime.tv_sec;

  //add microseconds
  user_mode_time += ((double) usage.ru_utime.tv_usec)/1000000;
  kernel_mode_time += ((double) usage.ru_stime.tv_usec)/1000000;

  //return sum of both times
  return user_mode_time + kernel_mode_time;
}

void Chrono::reset()
{
  struct rusage usage;

  if(getrusage(RUSAGE_SELF, &usage) == -1)
    ErrThrow("Error in GetTime!");

  //add up system and user time since program start and set them as start_time
  double tv_sec = usage.ru_stime.tv_sec + usage.ru_utime.tv_sec;
  double tv_usec = (usage.ru_stime.tv_usec + usage.ru_utime.tv_usec)/(double)1000000;

  start_time = tv_sec + tv_usec;
}


