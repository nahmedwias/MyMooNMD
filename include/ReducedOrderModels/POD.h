/** ************************************************************************ 
*
* @name       POD
* @brief      Computation/manipulation of POD basis independently of problem dimension
*
* @author     Swetlana Giere & Alfonso Caiazzo
* @date       08.03.2017 (start of implementation), 15.1.2019 (restart)
*
****************************************************************************/

#ifndef __POD__
#define __POD__

// boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>

#include <MooNMD_Io.h>
#include <Matrix.h>
#include <ParameterDatabase.h>
#include <fstream>
#include <memory>


// lapack
// #include "clapack.h"

using namespace std;
namespace ublas = boost::numeric::ublas;

class POD
{
    public:
    /** 
    * @brief constructor
    *
    */ 
    POD(const ParameterDatabase& param_db);
    
    /** 
    * @brief default destructor
    *
    */
    ~POD();
    
    /** 
    * @brief set gramian matrix (matrix describing inner product for POD)
    * 
    * Set the gramian matrix describing the inner product with respect to which
    * POD basis will be computed.
    * NOTE: If no gramian matrix is set, then it is automatically set to identity,
    * i.e. POD basis will be computed with respect to the euclidean inner product.
    *
    * @param mat ParMooN gramian matrix
    */
    void set_gramian( std::shared_ptr<TMatrix> mat );

    /**
    * @brief Convert parMooN matrix to ublas::compressed_matrix<double, row_major>
    *
    *
    * @param mat ParMooN matrix
    * @param res_mat Resulting ublas CSR matrix
    */
    void convert_to_ublas( std::shared_ptr<TMatrix> mat,
    		               ublas::compressed_matrix<double> & res );

    
    /**
    * @brief compute POD basis
    *
    * Compute POD basis from snapshots stored in
    * rom_db["snaps_directory"] + '/' + rom_db["snaps_basename"] + "snap".
    * If this->pod_fluct (rom_db["pod_fluct"]) is set to true then POD
    * basis will be computed from the fluctuating part of the snapshots and the
    * time average of the snapshots will be stored in the class member
    * variable this->snaps_mean. The rank of the basis is given by the class
    * member variable this->rank (or rom_db["pod_rank"]). If this->rank<=0, then
    * all possible POD basis functions will be computed and stored.
    */
    void compute_basis();
    
    /** 
    * @brief read basis from file
    * 
    * Read POD basis from the file with the name
    * this->rom_db["pod_directory"] + "/" + this->rom_db["pod_basename"] + ".pod"
    * If this->rom_db["pod_fluct"] is set to true, then the average
    * of the snapshots will be automatically read from the file 
    * this->rom_db["pod_directory"] + "/" + this->rom_db["pod_directory"] + ".mean"
    * and stored in the member variable snaps_mean.
    */
    void read_basis();
    
    /** 
    * @brief Write POD data into file
    * 
    * Write POD basis functions(row-wise) into
    * this->rom_db["pod_directory"] + "/" + this->rom_db["pod_basename"] + this->rom_db["pod_inner_product"] + ".pod",
    * time average of snapshots (if this->rom_db["pod_fluct"]==true) into
    * this->rom_db["pod_directory"] + "/" + this->rom_db["pod_basename"] + this->rom_db["pod_inner_product"] + ".mean",
    * [time, POD eigenvalues, missing energy ratio] into
    * this->rom_db["pod_directory"] + "/" + this->rom_db["pod_basename"] + this->rom_db["pod_inner_product"] + ".eigs".
    *
    * @param basename = this->rom_db["pod_basename"] + this->rom_db["pod_inner_product"] + "."
    * TODO at the momnet basename is set in class POD_TCDR2D. It is not necessary. DO it in this class.
    */
    void write_pod( std::string basename );
    
    /** 
    * @brief Write matrix into file
    * 
    * Write the matrix _mat into the file _filename
    *
    * @param _mat matrix which will be written into file
    * @param _filename full name of the file
    */
    void write_data( ublas::matrix<double> &_mat, string _filename );
    
    /* getta functions */
    int get_rank() const {
      return rank;
    }  
    int get_length() const {
      return length;
    }
    const ParameterDatabase & get_db() const
           { return rom_db; }

    const ublas::matrix<double> get_basis() const {
          return pod_basis;
    }

    const ublas::vector<double> get_snaps_avr() const {
              return snaps_mean;
    }
  
    protected:
    // memeber variables

    /* parameter database */
    ParameterDatabase rom_db;
    /* rank of pod basis */
    int rank;
    /* dof of pod basis functions */
    int length;
    /* threshold value for eigenvalues of autocorrelation matrix */
    double eigen_threshold;
    /* number of eigenvalues greater than eigen_threshold */
    int valid_eigs;
    /* dof of snapshots */
    int length_snaps;
    /* number of snapshots */
    int number_snaps;
    /* array for pod eigenvalues */
    double* eigs;
    /* inner product matrix for pod computation */ 
    ublas::compressed_matrix<double> gramian_mat;
    
    /* matrix with snapshots values (DO WE NEED IT? TOO MUCH OF STORAGE?) */
    ublas::matrix<double> snaps_mat;
    /* matrix for pod basis */
    ublas::matrix<double> pod_basis;
    /* time_mean of the snapshots */
    ublas::vector<double> snaps_mean;

  private:

    /* functions */
    
    /** 
    * @brief read snapshots
    *
    * Read snapshots from file
    * this->rom_db["snaps_directory"] + "/" + this->rom_db["snaps_basename"] + "snap".
    * This function is called by the constructor.
    * The snapshots will be stored into the class member this->snaps_mat colomn-wise.
    *  
    */
    void read_snapshots();
    
    /**
    * @brief Decompose snapshots by substracting their mean values
    *
    * Decompose snapshots into the time average of the snapshots and the fluctuating
    * part of the snapshots stored in this->snaps_mat.
    */
    void decompose_snaps();
    
    /**
    * @brief Compute auto-correlation matrix for eigenvalue problem
    *  
    * Compute autocorrelation matrix for POD computation which is computed as 
    * U^T*S*U, where U is the snapshot matrix (this->snaps_mat) and S is the inner
    * product matrix (gramian_mat). If S is not specified (i.e. set_gramian() was
    * not called), then the POD basis will be computed with respect to the 
    * euclidian inner product with S = Id.
    *
    * @param corr_mat resulting auto-correlation matrix
    */
    void compute_autocorr_mat( ublas::matrix<double> &corr_mat );
    
    /** 
    * @brief Write time average of snapshots into file
    * 
    * See documentation for write_pod(string) for more info.
    */
    void write_averages( std::string basename );
    
    /** 
    * @brief Write POD eigenvalue data into file
    * 
    * See documentation for write_pod(string) for more info.
    */
    void write_eigenvalues( std::string basename );
    
    /** 
    * @brief Read time average of snapshots from file
    *
    * Read snapshots mean from the file
    * this->rom_db["pod_directory"] + "/" + this->rom_db["pod_basename"] + ".mean".
    * It is automatically called from read_basis() if rom_db["pod_fluct"]==true.
    */
    void read_averages();
    
    /**
    * @brief Read data from file
    *
    * Read data from the file _filename and save it into the matrix _mat
    *
    * @param _filename full name of the file
    * @param _mat matrix  
    */
    void read_data( std::string _filename, vector < vector<double> > &_mat );
    
    /**
    * @brief read data from file
    *
    * Read data from the file _filename and save it into the vector 'data'
    *
    * @param _filename full name of the file
    * @param data vector  
    */
    void read_data( std::string _filename, vector<double> &data);

    /**
    * @brief matrix-matrix product
    *
    * C = A * B
    */
    void mat_prod(ublas::matrix<double>& C, ublas::matrix<double>& A, ublas::matrix<double>& B);
};

#endif
