/** ***************************************************************************
 *
 * @name       TParFEMapper3D
 * @brief      Class containing communication routines between processes which all
 * 			   administer their own subdomain and corresponding degrees of freedom of
 * 			   an FE system.
 *
 *             Objects of this type are set up after the domain is split, so
 *             there exists a particular one in each process. The data members
 *             refer to the data of the processes' Domain and FESpace.
 *
 *			   All the information the class needs is known to an object of class
 *			   TParFEMapper3D.
 *
 * @author     Sashikumaar Ganesan
 * @date       2015/04/24
 *
 * @ruleof0
 *
 ******************************************************************************/

#ifdef _MPI

#ifndef __PARFECOMMUNICATOR3D__
#define __PARFECOMMUNICATOR3D__

#include "mpi.h"

#include <ParFEMapper3D.h>
#include <string>

class TParFECommunicator3D
{
  public:
    /**
     * Set up a ParFECommunicator belonging to a certain FESpace.
     *
     * TODO This should be replaced by a ctor taking an FESpace - ParFEMapper
     * should only work as the underlying data structure and not be visible from
     * the outside.
     */
    TParFECommunicator3D(TParFEMapper3D *mapper);
    
    /// Default constructor. Creates an almost empty object, which should not be used.
    TParFECommunicator3D();
    
    /// Gather information about this communicator in root and print it
    /// to console and output file.
    void print_info() const;

    /**
     * Restores the required level of consistency of a vector.
     * For the notion of consistency levels refer to our ParMooN paper:
     * TODO cite the ParMooN paper correctly...
     *
     * Not that this method is not range checked - you must make sure that
     * "vector" whcih is supposed to be updated is as long as N_Dof.
     */
    void consistency_update(double* vector, size_t level) const;
    
    //TODO comment
    // Verdacht: die Methode addiert die Werte von Interface-Slaves und Interface-
    // Mastern auf und stellt mit diesen addierten Werten Level-1-Konsistenz her.
    // Dies scheint nach Multigrid-Gridtransfers noetig zu sein ( zumindest wird
    // die Methode dort aufgerufen).
    void CommUpdateReduce(double *rhs) const;

    //TODO Comment the usage of this!
    const int *GetMaster() const
    {return Master;}
    
    //TODO Comment the usage of this!
    int GetNDof() const
    {return N_Dof;}
    
    /// Return number of dimensions the communicator is used for.
    int get_n_dim() const
    {
      return Mapper->get_n_dim();
    }

    //TODO Comment the usage of this!
    const int* Get_Local2Global() const
    { return Mapper->Get_Local2Global();}
    
    /// Returns the number of master d.o.f. on this process.
    int GetN_Master() const
    { return Mapper->GetN_Master();}
    
    /// Return a pointer to the fe space this communicator belongs to.
    const TFESpace3D* get_fe_space() const
    {
      return Mapper->get_fe_space();
    }

    //Special member functions. Rule of zero applied (class does not manage ressources).

    //Declaration of special member functions - rule of zero

    //! Default copy constructor.
    TParFECommunicator3D(const TParFECommunicator3D&) = default;

    //! Default move constructor.
    TParFECommunicator3D(TParFECommunicator3D&&) = default;

    //! Default copy assignment operator. Performs deep copy.
    TParFECommunicator3D& operator=(const TParFECommunicator3D&) = default;

    //! Default move assignment operator
    TParFECommunicator3D& operator=(TParFECommunicator3D&&) = default;

    //! Default destructor.
    ~TParFECommunicator3D() = default;

  private:
    TParFEMapper3D *Mapper;

    int N_Dim;

    int N_Dof;

    int N_SendDof, N_SendDofMS, N_SendDofH1, N_SendDofH2;

    double *Send_Info, *Send_InfoMS, *Send_InfoH1, *Send_InfoH2;

    double *Recv_Info, *Recv_InfoMS, *Recv_InfoH1, *Recv_InfoH2;

    int *N_DofSend, *N_DofSendMS, *N_DofSendH1, *N_DofSendH2;

    int *N_DofRecv, *N_DofRecvMS, *N_DofRecvH1, *N_DofRecvH2;

    int *sdispl, *sdisplMS, *sdisplH1, *sdisplH2;

    int *rdispl, *rdisplMS, *rdisplH1, *rdisplH2;

    int *DofSend, *DofSendMS, *DofSendH1, *DofSendH2;

    int *DofRecv, *DofRecvMS, *DofRecvH1, *DofRecvH2;

    int N_Slave, N_InterfaceS, N_Halo1, N_Halo2;

    MPI_Comm Comm;

    int *Master;

    /**
     * Sets all interface slave dof values in fe_vec to the value of their
     * interface master dof on other processors.
     *
     * @param[in] fe_vec The fe_vector to be updated. Its length must be equal
     * to the number of dofs stored in this communicator, else: garbage-in-garbage-out.
     */
    void CommUpdateMS(double* fe_vec) const;

    /**
     * Sets all halo 1 dof values in fe_vec to the value of their
     * master dof on other processors.
     *
     * @param[in] fe_vec The fe_vector to be updated. Its length must be equal
     * to the number of dofs stored in this communicator, else: garbage-in-garbage-out.
     */
    void CommUpdateH1(double* fe_vec) const;

    /**
     * Sets all halo 2 dof values in fe_vec to the value of their
     * master dof on other processors.
     *
     * @param[in] fe_vec The fe_vector to be updated. Its length must be equal
     * to the number of dofs stored in this communicator, else: garbage-in-garbage-out.
     */
    void CommUpdateH2(double* fe_vec) const;

};

#endif
#endif






















