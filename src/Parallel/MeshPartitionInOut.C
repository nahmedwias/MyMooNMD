/*
 * Read_write_metis.C
 *
 *  Created on: Mar 13, 2017
 *      Author: rominger
 */
#ifdef _MPI

#include <Domain.h>
#include <MeshPartitionInOut.h>
#include <MooNMD_Io.h>

#include "mpi.h"



void MeshPartitionInOut::read_file (const TDomain& Domain,int size,
                                    int N_Cells, idx_t *Cell_Rank)
{
    std::string file_name = Domain.get_database()["read_metis_file"];

    Output::info("PARTITIONING", "Reading mesh partition information "
         "from file ",file_name);

	  std::ifstream ifs; //Deklarieren von ifs zum Oeffnen der Textdatei
	  ifs.open(file_name); //Textdatei wird geoeffnet
	  if(ifs.is_open())
	      {
		      std::string line;
		      int l =0; //a line counter
		      while (std::getline(ifs, line)) //Reading line by line
		      {
		    	  if(l==0)//we expect this line to hold "ParMooN Mesh Partition File"
		    	  {
		    		  if (line!="ParMooN Mesh Partition File")
		    		  {
		    		    throw std::runtime_error("The first line of the "
		    		        "written-file does not match.");
		    		  }
		    		  ++l;
		    		  continue;
		    	  }
		    	  else if(l==1)
		    		  // Reads the total Number of Cells and the Number of Partitions
		    		  //  of the written file and compares them to current configuration
		    		  // throws error if they don't match
		    	  {
		    		  int position = line.find_first_of(',');
		    		  std::string number_partition = line.substr(position+1);
		    		  std::string number_cells = line.substr(0, position);
		    		  number_partition.erase (0,14);
		    		  number_cells.erase (0,9);
		    		  int total_partition_number;
		    		  try
		    		  {
		    			  total_partition_number=std::stoi(number_partition);
		    			  if (total_partition_number != size)
		    			  {
		    				  throw std::runtime_error("The total number of partitions of "
		    				      "the written file does not match the total number of "
		    				      "partitions that is used");
		    			  }
		    		  }
		    		  catch (std::exception& e)
		    		  {
		    			  Output::print(e.what());
		    		  }
		    		  int total_number_of_cells;
		    		  try
		    		  {
		    			  total_number_of_cells=std::stoi(number_cells);
		    			  if (total_number_of_cells!=N_Cells)
		    			  {
		    				  throw std::runtime_error("The total number of cells of the "
		    				      "written file does not match the total number of "
		    				      "cells that are used.");
		    			  }
		    		  }
		    		  catch (std::exception& e)
		    		  {
		    			  Output::print(e.what());
		    		  }
		    		  ++l;
		    		  continue;
		    	  }
		    	  else
		    	  {
		    		  int Comment = line.find_first_of('#'); //Looks for lines that start with #
		    		  if (Comment == 0){cout<<line<<endl;continue;} // If a Line starts with # it prints the comment
		    		  else
		    			  //Here we write out the Cellnumber and the matching Cell_Rank.
		    			  //If the partition number or number of cells is too big, than an error is thrown
		    		  {
		    			  int pos = line.find_first_of(',');
		    		  	  std::string Rank = line.substr(pos+1);
		    		  	  std::string Cell = line.substr(0, pos);
		    		  	  Rank.erase (0,5);
		    		  	  Cell.erase (0,5);
		    		  	  int Rank_Number;
		    		  	  try
		    		  	  {
		    		  		  Rank_Number=std::stoi(Rank);
		    		  		  if (Rank_Number >= size)
		    		  		  {
		    		  			  ErrThrow("****************************************"
		    		  			      "The input file contains a rank greater than mpi_size.");
		    		  		  }
		    		  	  }
		    		  	  catch (std::exception& e)
		    		  	  {
		    		  		  Output::print(e.what());
		    		  	  }

		    		  	  int Cellnumber;
		    		  	  try
		    		  	  {
		    		  		  Cellnumber=std::stoi(Cell);
		    		  		  if (Cellnumber>=N_Cells)
		    		  		  {
		    		  		    ErrThrow("****************************************"
		    		  			    "The input file contains a cell number greater than the "
		    		  			    "maximal cell number. Additionally the cell number does "
		    		  			    "not match to the used line number l-2.");
		    		  		  }
		    		  	  }
		    		  	  catch (std::exception& e)
		    		  	  {
		    		  		Output::print(e.what());
		    		  	  }
		    		  	  Cell_Rank[l-2]=Rank_Number;
		    			 //cout<<"The "<<l-2<<" entry of the CellRank is: "<< Cell_Rank[l-2]<<endl;

		    			  if (l>N_Cells+1)
		    			  {
		    			    ErrThrow("****************************************"
		    			      "The file has more lines than it should.");
		    			  }
		    			  ++l;
		    		  }
		    	  }

		      }
		      if (l-2!=N_Cells)
		      {
		        ErrThrow("****************************************"
		            "The file has less lines than it should.");
		      }
	      }
	  else
	  {
	     ErrThrow("**********************************The file was not read. ");
	  }

}
void MeshPartitionInOut::write_file(const TDomain& Domain, int size,
                                    int N_Cells, idx_t *Cell_Rank)
{
	  std::string file_name = Domain.get_database()["write_metis_file"];

    Output::info("MeshPartitionInOut", "Writing mesh partition information "
        "to file ",file_name);

	  std::ofstream ofs (file_name); //Create write-txt-file
	  ofs << "ParMooN Mesh Partition File\n" << "N_Cells: "<< N_Cells << " ,N_Processors: " << size <<"\n";

	  for(int i =0; i < N_Cells ; i++) // Loop cell entries of cell rank
	  {
	    //Write entries of each Cell of Cell_Rank in seperate line in txt-file
	    ofs <<"Cell "<< i<< " ,Rank "<< Cell_Rank[i] <<endl;
	  }
	  ofs.close(); //closes write in file
	}

#endif
