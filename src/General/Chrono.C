/**
 * @file Chrono.C
 * Implementation of class Chrono declared in Chrono.h.
 *
 * @date 2016/04/18
 */
#include <Chrono.h>
#include <MooNMD_Io.h>
#include <time.h>

#ifdef _MPI
#include <mpi.h>
#endif

#ifdef _OMP
#include <omp.h>
#endif


Chrono::Chrono()
{
  reset();
}

void Chrono::print_time(const std::string& program_part)
{
#ifdef _MPI
  print_time_mpi(program_part);
  return;
#endif
#ifdef _OMP
  print_time_omp(program_part);
  return;
#endif

  print_time_seq(program_part);
}

void Chrono::print_time_seq(const std::string& program_part) const
{
  //rusage time
  double time_rusage = get_rusage_time() - start_time_rusage;
  Output::print("--- time for ", program_part,": ", time_rusage, " s (sequential rusage time)");

  //wall clock time
  double time_wall = get_wall_time() - start_time_wall;
  Output::print("--- time for ", program_part,": ", time_wall, " s (sequential wall time)");
}
#ifdef _MPI
void Chrono::print_time_mpi(const std::string& program_part) const
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int my_rank;
  MPI_Comm_rank(comm, &my_rank);

  //rusage time in root
  double time_rusage = get_rusage_time() - start_time_rusage;
  if(my_rank==0)
    Output::print("--- time for ", program_part,": ", time_rusage, " s (mpi root rusage time)");

  //wall clock time in root
  double time_wall = get_wall_time() - start_time_wall;
  if(my_rank==0)
    Output::print("--- time for ", program_part,": ", time_wall, " s (mpi root wall time)");

  //mpi_wtime
  double time_spent = MPI_Wtime() - start_time_mpi_wtime;
  double time_max;
  double time_min;
  MPI_Reduce(&time_spent, &time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
  MPI_Reduce(&time_spent, &time_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
  if(my_rank == 0)
  {
   Output::print("--- time for ", program_part,": ", time_max, " s (max wtime in a process)");
   Output::print("--- time for ", program_part,": ", time_min, " s (min wtime in a process)");
  }
  return;
}
#endif

#ifdef _OMP
void Chrono::print_time_omp(const std::string& program_part) const
{
  //rusage time
  double time_rusage = get_rusage_time() - start_time_rusage;
  Output::print("--- time for ", program_part,": ", time_rusage, " s (openmp rusage time)");

  //wall clock time
  double time_wall = get_wall_time() - start_time_wall;
  Output::print("--- time for ", program_part,": ", time_wall, " s (openmp wall time)");
}
#endif


double Chrono::get_rusage_time()
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

double Chrono::get_wall_time()
{
#ifdef _OMP
  return omp_get_wtime();
#endif
  //wall time
  struct timeval wall_time;
  if (gettimeofday(&wall_time,NULL))
    ErrThrow("Error in gettimeofday!");
  return wall_time.tv_sec + wall_time.tv_usec/(double)1000000;


}

void Chrono::reset()
{
  //rusage time (system + cpu)
  struct rusage usage;
  if(getrusage(RUSAGE_SELF, &usage) == -1)
    ErrThrow("Error in getrusage!");
  //add up system and user time since program start and set them as start_time_rusage
  double tv_sec = usage.ru_stime.tv_sec + usage.ru_utime.tv_sec;
  double tv_usec = (usage.ru_stime.tv_usec + usage.ru_utime.tv_usec)/(double)1000000;
  start_time_rusage = tv_sec + tv_usec;

  //wall time
  struct timeval wall_time;
  if (gettimeofday(&wall_time,NULL))
    ErrThrow("Error in gettimeofday!");
  start_time_wall = wall_time.tv_sec + wall_time.tv_usec/(double)1000000;
#ifdef _OMP
  start_time_wall = omp_get_wtime();
#endif

  //mpi wtime
#ifdef _MPI
  start_time_mpi_wtime = MPI_Wtime();
#endif
}


