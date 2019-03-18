#ifndef INCLUDE_PARALLEL_MESH_PARTITION_IN_OUT_INCLUDED
#define INCLUDE_PARALLEL_MESH_PARTITION_IN_OUT_INCLUDED

#ifdef _MPI
#include "metis.h" //for the id_x typedef

//forward declaration
class TDomain;

/**
 * A namespace holding to methods. They are used for reading/writing the
 * information, which cell of a domain belongs to which processor in a parallel
 * run of the program to/from a text file. The respective filename must be known
 * to the database of the domain.
 * Both methods need be called by root only, this is a sequential part of the
 * code (still).
 *
 *
 * TODO: Describe the used file format.
 */
namespace MeshPartitionInOut
{
  /// Read domain partitioning information from a file.
	void read_file(const TDomain& Domain, int size,
	               int N_Cells, idx_t *Cell_Rank);

	/// Write domain partitioning information to a file.
	void write_file(const TDomain& Domain, int size,
	                int N_Cells, idx_t *Cell_Rank);
}
#endif

#endif //INCLUDE_PARALLEL_MESH_PARTITION_IN_OUT_INCLUDED
