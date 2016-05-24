/**
 * @file Chrono.h
 *
 * This class is intended for time measuring. Interface is pretty easy. Use
 * constructor to set up a timer, which immediately starts time measuring.
 * To start a new measuring with the same timer, call reset.
 * To print out measured time to console and outputfile call print_time.
 *
 * Since we are still a bit undecisive which time to measure, this class
 * uses more than one method (per parallel type) and prints out all. Those are
 *
 *    in SEQUENTIAL case:
 *      - cpu time (user plus system time) with "rusage"
 *      - wall clock time with "gettimeofday"
 *    in OPENMP case:
 *      - cpu time (user plus system time) with "rusage", summed up over all threads
 *      - wall clock time with "omp_get_wtime"
 *    in MPI case:
 *      - cpu time (user plus system time) with "rusage", measured in root only
 *      - wall clock time with "gettimeofday", measured in root only
 *      - mpi_wtime maximum over all processes
 *      - mpi_wtime minimum over all processes
 *
 * @date 2016/04/18
 *
 * @ruleof0
 */

#ifndef INCLUDE_GENERAL_CHRONO_H_
#define INCLUDE_GENERAL_CHRONO_H_

#include <sys/time.h>
#include <sys/resource.h>
#include <string>

class Chrono
{
  public:
    /// Create a Chrono object with start times set to the moment of its creation.
    Chrono();

    /// Reset the start times.
    void reset();

    /**
     * Print out the measured times. The output as well as the means of time
     * measuring depend on the parallel type.
     *
     * @param program_part will be put into the output string
     */
    void print_time(const std::string& program_part);

  private:
    //// The starting time for time measurement with means of rusage
    double start_time_rusage;

    /// The wall clock starting time
    double start_time_wall;

#ifdef _MPI
    /// The starting time as gained with MPI_Wtime
    double start_time_mpi_wtime;
#endif

    void print_time_seq(const std::string& program_part) const;

    void print_time_mpi(const std::string& program_part) const;

    void print_time_omp(const std::string& program_part) const;


    /**
     * Evaluates system time spent in user mode and in kernel mode so far.
     * Returns sum of both.
     *
     * @return  added up time spent in user and kernel mode since start of the program
     */
    static double get_rusage_time();

    /// @return Current wall time.
    static double get_wall_time();

};



#endif /* INCLUDE_GENERAL_CHRONO_H_ */
