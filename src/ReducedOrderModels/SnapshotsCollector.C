/** **************************************************************************** 
*
* @name       SnapshotsCollector
* @brief      write/read snapshot data into/from a file
*
* @author     Swetlana Giere & Alfonso Caiazzo
* @date       08.05.2012 (start of implementaion). Restarted on 15.1.2019
*
*******************************************************************************/

#include <SnapshotsCollector.h>
#include <sys/stat.h>

/** ***************************************************************************/
ParameterDatabase SnapshotsCollector::default_snapshots_database()
{
  ParameterDatabase db("Default ParMooN parameter database for storing "
                       "snapshots");

  db.add("snaps_directory", ".",
         "This directory is where the snapshots will be written. This "
         "directory will be created, if it does not exist already. Files in "
         "this directory will be overwritten without any warning.");

  db.add("snaps_basename", "parmoon_snapshots",
         "Basename for file where the snapshots are stored.");

  db.add("steps_per_snap", (size_t) 5,
         "This integer specifies how many time steps are performed "
         "before a snapshot has to be written into file.");

  return db;
}

/** ***************************************************************************/
SnapshotsCollector::SnapshotsCollector(const ParameterDatabase& param_db):
                      db(SnapshotsCollector::default_snapshots_database()),
                      snap_count(0)
{
  this->db.merge(param_db, true);
  snapshot_filename =  this->db["snaps_directory"].get<std::string>() + "/";
  snapshot_filename += this->db["snaps_basename"].get<std::string>();
//  snapshot_filename += ".snap";

  if (db["write_snaps"])
  {
    // create directory db["snaps_directory"]
    std::string directory_name = this->db["snaps_directory"].get<std::string>();
    mkdir(directory_name.c_str(), 0777);

    Output::print<1>("Filename for storing snapshots: ", snapshot_filename);

    this->datastream.open(snapshot_filename.c_str(),
                          ios::out | ios::trunc | ios::in);

    if( ! this->datastream.good() )
    {
      ErrThrow("Error: File ", snapshot_filename,
               " could not be opened in SnapshotsCollector.\n"
               "(No read access to file or file already open)");
    }

    this->datastream << setprecision(12);
  }
}

/** ***************************************************************************/
SnapshotsCollector::~SnapshotsCollector(){
  if( this->datastream.is_open() ) this->datastream.close();
}

/** ***************************************************************************/
void SnapshotsCollector::write_data(const BlockVector& solution,
                                    size_t             time_step)
{
  if( (time_step % (int)db["steps_per_snap"]) == 0 )
  {
    Output::print<1>("Time step ", time_step, ", writing snapshot n. ",
                     snap_count+1, " on ", snapshot_filename);
    snap_count++;

    if( ! this->datastream.is_open() )
    {
      ErrThrow("Error: Snapfile is not open.");
    }

    for( size_t i = 0; i < solution.length(); ++i )
    {
      this->datastream << solution.get_entries()[i] << " ";
    }

    this->datastream << endl;
    this->datastream.flush();
  }
}
