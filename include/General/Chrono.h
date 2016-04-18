/**
 * @file StopWatch.h
 *
 *  A classe whci takes care of perormance time measuring.
 *  Kind of rivals the Method GetTime() which is implemented in MainUtilities
 *  but I would encourage using it.
 *
 *  TODO So far it only operates properly in sequential.
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
    /**
     * Create a Chrono object with start time set to the user plus system
     * time at the moment of its creation.
     */
    Chrono();

    /** Reset the start time. */
    void reset();

    /** Print out the added up user time and kernel time since last reset.
     * @param program_part will be put into the output string
     */
    void print_time(std::string program_part);

  private:
    /** The starting time for time measurement as used in sequential (s). */
    double start_time;

    /**
     * Evaluates system time spent in user mode and in kernel mode so far.
     * Returns sum of both.
     *
     * TODO Should not be used in MPI - what about OpenMP?
     * @return  added up time spent in user and kernel mode since start of the program
     */
    static double get_exec_time_s();
};



#endif /* INCLUDE_GENERAL_CHRONO_H_ */
