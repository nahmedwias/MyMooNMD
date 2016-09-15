/**
 * @file Chrono.h
 *
 * This class is intended for time measuring. Interface is pretty easy. Use
 * constructor to set up a timer, which immediately starts time measuring.
 * To start a new measuring with the same timer, call reset.
 * To print out measured time to console and outputfile call print_time.
 * 
 * You can stop and start the timer whenever you like. This class adds all the 
 * durations spent between successive calls to start() and stop(). Between such
 * calls the bool `running` is set to true, otherwise false.
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
    /// @brief Create a Chrono object with start times set to the moment of its 
    /// creation.
    Chrono();

    /// Reset the start times. 
    void reset();
    
    /// start the clock
    void start();
    
    /// @brief stop the clock
    /// 
    /// The time since the last call to start() is added to the total time. If
    /// this Chrono object was not running, nothing is done. 
    /// @return the wall time since last start() or zero if running was false
    double stop();

    /**
     * @brief Print out the measured times since construction or last call to 
     * reset().
     * 
     * The CPU time and the wall time are printed. In MPI mode the maximum, 
     * minimum and average over all processes of these times are printed.
     *
     * @param program_part will be put into the output string
     */
    void print_time(const std::string& program_part) const;
    void print_time(const char* program_part) const
    { print_time(std::string(program_part)); }
    
    /**
     * @brief Print out the measured times since last call to start().
     * 
     * This invokes stop() and start().
     * 
     * The CPU time and the wall time are printed. In MPI mode the maximum, 
     * minimum and average over all processes of these times are printed.
     *
     * @param program_part will be put into the output string
     */
    void print_time_since_last_start(const std::string& program_part);
    void print_time_since_last_start(const char* program_part)
    { print_time_since_last_start(std::string(program_part)); }
    
    /// @brief return the CPU time accumulated over all start/stop cycles
    ///
    /// This includes the elapsed time since the last call to start(), even if
    /// stop() has not yet been called.
    double elapsed_time() const;
    /// @brief return the wall time accumulated over all start/stop cycles
    ///
    /// This includes the elapsed time since the last call to start(), even if
    /// stop() has not yet been called.
    double elapsed_wall_time() const;
    
    /// @brief return the CPU time elapsed since last start()
    ///
    /// Returns zero if not running.
    double time_since_last_start() const;
    
    /// @brief return the wall time elapsed since last start()
    ///
    /// Returns zero if not running.
    double wall_time_since_last_start() const;

  private:
    /// The starting time for time measurement with means of rusage
    double start_time_rusage;
    
    /// accumulated time over multiple start/stop cycles (using rusage)
    double cumulative_time_rusage;

    /// The wall clock starting time
    double start_time_wall;
    
    /// accumulated time over multiple start/stop cycles (using wall time)
    double cumulative_time_wall;
    
    /// true after calling start(), false after calling stop()
    bool running;
};



#endif /* INCLUDE_GENERAL_CHRONO_H_ */
