/** ************************************************************************ 
*
* @name       SnapshotsCollector
* @brief      write snapshots (finite elements coefficients of the solution) into a file
*
* @author     Swetlana Giere & Alfonso Caiazzo
* @date       08.03.2017. Restarted on 15.1.2019
*
****************************************************************************/


#ifndef __SNAPSHOTSCOLLECTOR__
#define __SNAPSHOTSCOLLECTOR__


#include <MooNMD_Io.h>
#include <ParameterDatabase.h>
#include <BlockVector.h>


#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

class SnapshotsCollector{
  
  public:
    
    /** 
    * @brief constructor
    * 
    * The filename for storing snapshots is constructed as
    * filename = db["snaps_directory"] + "/" + db["snaps_basename"] + "snap".
    * Then, the stream is created and stored as a member of the class.
    *
    */
    SnapshotsCollector( const ParameterDatabase& param_db );
    
    /** 
    * @brief destructor closes the datafile
    * 
    */
    ~SnapshotsCollector();

    /** 
    * @brief Write s snapshot into file
    * 
    * Write finite elemenent coefficients of the solution into the file
    * represented by the class member 'datastream'. Moreover, the snapshots
    * will be written only if the database parameter db["write_snaps"] is set
    * to true and the condition time_step % db["steps_per_snap"] == 0 is
    * satisfied.
    * Note: For 3D problems a binary format could be more appropriate to avoid
    * huge storage volume requirements (to be implemented!)
    * 
    * @param solution   vector containing the FE coeffs of the solution (snapshot)
    * @param time_step  count of the time step from the time loop (main program)
    *
    */
    void write_data( const BlockVector &solution, size_t time_step=0);
    
  private:

    /* stream for writing snapshots */
    fstream datastream;
    ParameterDatabase db;
    int snap_count;
    string snapshot_filename;
};

#endif
