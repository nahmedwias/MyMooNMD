/** ************************************************************************ 
*
* @name       SNAPS
* @brief      write/read snapshot data into/from a file
*
* @author     Swetlana Schyschlowa
* @date       08.05.2012 (start of implementaion)
*
****************************************************************************/

#include <SNAPS.h>

ParameterDatabase get_default_snaps_parameters()
{
  Output::print<5>("Creating a default parameter database for managing snapshots...");
  ParameterDatabase db("Default ParMooN parameter database for storing snapshots");

  db.add("snaps_directory", ".",
    	         "This directory is where the snapshots are written. This "
    	         "directory will be created, if it does not exist already. Files in "
    	         "this directory will be overwritten without any warning.");

  db.add("snaps_basename", "parmoon_snapshots",
  		  "Basename for file where the snapshots are stored. The basename should end with a dot.");

  db.add("write_snaps", false,
             "This is the flag whether the snapshots should be written to a file. ",
  		   {true,false});

  db.add("steps_per_snap", (size_t) 5,
    	     "This integer specifies how many time steps are performed "
    		 "before a snapshot has to be written into file. ");
  return db;
}

/** ***********************************************************************/
SNAPS::SNAPS( const ParameterDatabase& param_db ) : db(get_default_snaps_parameters())
{
  this->db.merge(param_db, false);
  string filename =  this->db["snaps_directory"].get<std::string>() + "/";
  filename += this->db["snaps_basename"].get<std::string>() + "snap";
  
  Output::print<1>( "Filename for storing snapshots: ", filename);

  this->datastream.open( filename.c_str(), ios::out | ios::trunc | ios::in);
  
  if( ! this->datastream.good() )
  {
    ErrThrow("Error: File ", filename, " could not be opened in SNAPS.\n"
    		"(No read access to file or file already open)");
  }
  
  this->datastream << setprecision( 12 );
}

/** ***********************************************************************/
SNAPS::~SNAPS(){
  if( this->datastream.is_open() ) this->datastream.close();
}

/** ***********************************************************************/
void SNAPS::write_data(const BlockVector &solution, size_t time_step )
{
  if(db["write_snaps"] && (time_step % (int)db["steps_per_snap"] == 0))
  {
    Output::print<1>("Writing snapshot into file...");
  
    if( ! this->datastream.is_open() )
      ErrThrow( "Error: Snapfile is not open." );
  
    for(size_t i = 0; i < solution.length(); ++i)
      this->datastream << solution.get_entries()[ i ] << " ";
    this->datastream << endl;
    this->datastream.flush();
  }
}
