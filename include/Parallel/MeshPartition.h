// =======================================================================
// @(#)MeshPartition.h
// 
// Purpose:     partition the domain into "npart" parts for parallel computing
// 
// Author:      Sashikumaar Ganesan
// History:      start of implementation  07/09/09 (Sashikumaar Ganesan)
// =======================================================================
#ifndef __MESH_PARTITION__
#define __MESH_PARTITION__

class TDomain;
class TVertex;

namespace MeshPartitionAuxFunctions
{
void Sort(TVertex **Array, int length);

int GetIndex(TVertex **Array, int Length, TVertex *Element);
}

#ifdef  _MPI

/** @brief Not implemented.*/
void Partition_Mesh(MPI_Comm comm, TDomain *Domain, int &MaxCpV);

/**@brief Calls one of two Metis methods to distribute the mesh among
 * participating processes. Which one is determined by the control parameter
 * Par_P2. Setting that to 0 means call of METIS_PartMeshNodal, 1 a call of METIS_PartMeshDual.
 * The method further sets the datastructure of cells, domain
 * and collection to store the information on own and halo cells.
 * Works on finest refinement level only.
 * Can so far deal with homogeneous quad or hexahedra meshes only.
 *
 * @note One should keep in mind that after calling this function the Domain
 * is reduced to the process' subdomain - no process knows the whole domain anymore.
 *
 * @param[in] comm The used MPI Communicator.
 * @param[in,out] Domain The domain whose finest collection is to be distributed.
 * @param[out] MaxCpV The globally unique maximum number of cells per vertex.
 *
 * @return 0 if partitioning was successful, 1 if unsuccessful.
 *           At the moment no strong guarantee can be given, because
 *           we don't know enough about the method yet.
 */
int Partition_Mesh3D(MPI_Comm comm, TDomain *Domain, int &MaxCpV);

/**
 * @brief Removes unneeded halo cells which get introduced by refining a process' subdomain
 * (so after the initial mesh has been partitioned).
 *
 * This method changes the Cell Tree of the Domain it works on in such a way,
 * that all remaining cells are set as root cells. This makes the GetCollection method
 * relatively useless, for after the call only one cell level (0) is known to the Domains tree.
 * \todo Rework Domain_Crop as to modify the Domain's cell tree less brutally!
 *
 * @note The method does not physicallly delete any cells, edges or nodes at the moment.
 * That code was commented out for some reason.
 *
 *  @param[in] comm The MPI Communicator. In ParMooN this is almost always MPI_COMM_WORLD.
 *  @param[in,out] Domain The domain on whose finest level superfluous halo cells have to be removed.
 */
void Domain_Crop(MPI_Comm comm, TDomain *Domain);

#endif
#endif

