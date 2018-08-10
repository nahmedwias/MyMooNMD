/** ************************************************************************ 
*
* @name       POD
* @brief      all about computation and manipulation of POD basis
*
* @author     Swetlana Schyschlowa
* @date       08.05.2012 (start of implementation)
*
****************************************************************************/

#ifndef __POD_2D__
#define __POD_2D__

// boost
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <MooNMD_Io.h>
#include <Matrix.h>
#include <POD.h>
#include <PostProcessing2D.h>

#include <fstream>

// lapack
#include "clapack.h"

using namespace std;
namespace ublas = boost::numeric::ublas;


class POD_2D : public POD
{
  protected:
	/** @brief Finite Element space */
	TFESpace2D fe_space;
	/** @brief Finite Element function */
	TFEFunction2D fe_function;
	/** @brief the system matrix */
	BlockFEMatrix gramian_matrix;
	/** @brief solution vector with one component. */
	BlockVector pod_mode;

	/** @brief Definition of the used example */
	const Example_TimeCD2D example;

	/** @brief a local parameter database which constrols this class
	* The database given to the constructor will be merged into this one. Only
	* parameters which are of interest to this class are stored (and the
    * defualt ParMooN parameters).
	*/
	ParameterDatabase db;

	/** @brief class for output handling */
	PostProcessing2D outputWriter;

	/** @brief set parameters in database
	*
	* This functions checks if the parameters in the database are meaningful
	* and resets them otherwise. The hope is that after calling this function
	* this class is fully functional.
	*
	* If some parameters are set to unsupported values, an error occurs and
	* throws an exception.
	*/
	void set_parameters();

  public:
      /** @brief constructor
       * This constructor calls the other constructor creating an Example_TimeCD2D
       * object.
       */
      POD_2D(const TDomain& domain, const ParameterDatabase& param_db);

      /** @brief constructor
       *
       * The domain must have been refined a couple of times already. On the
       * finest level the finite element spaces and functions as well as
       * gramian matrix, pod mode and right hand side vectors are initialized.
       */
      POD_2D(const TDomain& domain, const ParameterDatabase& param_db,
      		const Example_TimeCD2D& ex);

      /**
      * @brief Compute POD basis
      *
      * Compute POD basis from snapshots with the compute_basis() routine
      * from the class POD. Within the function, appropriate gramian matrix
      * according to the pod_inner_product parameter is assembled and involved
      * into the computation of the POD basis.
      */
      void compute_basis();

      /**
      * @brief Read POD basis from file
      *
      * Read POD basis from file. The functionality could be interesting if POD
      * basis is already computed and atored into file and one
      * only wants to write vtk files of the POD basis functions.
      */
      void read_basis();

      /// @brief write POD basis to file and write vtk files
      void output();

       // getters and setters
      const Example_TimeCD2D& get_example() const
      { return example; }
      const TFEFunction2D & get_function() const
      { return this->systems.front().fe_function; }
      TFEFunction2D & get_function()
      { return this->systems.front().fe_function; }
      const BlockFEMatrix & get_stiff_matrix() const
      { return this->systems.front().stiff_matrix; }
      const BlockVector & get_rhs() const
      { return this->systems.front().rhs; }
      BlockVector & get_rhs()
      { return this->systems.front().rhs; }
      const BlockVector & get_solution() const
      { return this->systems.front().solution; }
      const TFESpace2D & get_space() const
      { return this->systems.front().fe_space; }
      const ParameterDatabase & get_db() const
      { return db; }

      /**
      * @brief return the computed errors at each discre time point
      *
      */
      std::array<double, int(3)> get_errors() const;
  
    public:
      
    /** 
    * @brief constructor
    * 
    * This constructor is suitable if one wishes to compute a POD basis.
    * The snapshots will be read from the specified snapshots file.
    * 
    * @param _snap_filename full name of the file where the snapshots 
    *                       are stored row-wise.
    * @param _rank desired rank of pod basis. if _rank<=0 is set, then 
    *              all available pod modes will be processed.
    * @param _pod_fluct a flag denoting whether POD basis has to be 
    *                   computed from the fluctuating part of the snapshots
    *
    */ 
    //POD(string _snap_filename, int _rank, bool _pod_fluct);
    POD(string _snap_filename, int _rank, bool _pod_fluct);
      
    /** 
    * @brief constructor
    * 
    * This constructor is suitable if one doesn't wish to compute POD basis
    * but rather read POD basis from a file in order to e.g. plot it. Use 
    * read_basis(..) to read the POD basis.
    * 
    * @param _rank desired rank of pod basis. if _rank<=0 is set, then 
    *              all available pod modes will be used.
    * @param _pod_fluct a flag denoting whether POD basis is
    *                   computed from the fluctuating part of the snapshots
    *                   
    */   
    POD(int _rank, bool _pod_fluct);
    
    /** 
    * @brief default destructor
    */
    ~POD();
    
    /** 
    * @brief set Gramian (inner product) matrix 
    * 
    * Set a Gramian matrix describing the inner product with respect to which 
    * POD basis will be computed.
    * NOTE: If no Gramian matrix is set, then it is set to identity, i.e. POD
    * basis will be computed with respect to Euclidean inner product.
    *
    * @param mat Gramian matrix  
    */
    void set_gramian( TMatrix &mat );
    
    /**
    * @brief change snapshot matrix
    *
    *  Overwrite member snapshots matrix snaps_mat  by the new one. 
    *  The members of the class length_snaps" and number_snaps 
    *  will be correspondingly changed. This could be necessary if 
    *  the snapshots have to be projected on a different space.
    *
    * @param new_snaps_mat new snapshots matrix
    */
    void change_snaps(ublas::matrix<double> &new_snaps_mat);
    
    /**
    * @brief compute POD basis
    *
    * Compute POD basis from snapshots stored in snaps_mat.
    * If the member pod_fluct is set to true then POD basis will 
    * be computed from the fluctuating part of the snapshots and the 
    * time average of the snapshots will be stored in the member 
    * variable average_snaps. The rank of the basis is given by the
    * member variable rank_basis.
    *
    */
    void compute_basis();
    
    /** 
    * @brief read basis from file
    * 
    * Read POD basis basis from the file tieh the name basename+pod. 
    * If POD basis doesn't have to be computed then use constructor 
    * POD(int, bool). If pod_fluct is set to true, then the average 
    * of the snapshots will be automatically read from the file 
    * basename+avr and stored in the member variable average_snaps. 
    *
    * @param basename name of the file without extension
    */
    void read_basis( string basename );
    
    /** 
    * @brief write POD data into file
    * 
    * Write POD basis (row-wise), time average of the snapshots
    * (if pod_fluct==true), POD eigenvalues into files with the 
    * names basename+pod, basename+avr, basename+svd, respectively.
    *
    * @param basename name of the file without extension
    */
    void write_pod( string basename );
    
    /** 
    * @brief write matrix into file _filename
    * 
    * Write the matrix _mat into the file 
    * TODO: delete this function and put it somewhere else, e.g.
    * main program (no connection to POD)
    *
    * @param _mat matrix which will be written into file
    * @param _filename full name of the file
    */
    void write_data( ublas::matrix<double> &_mat, string _filename );
    
    
    /* variables */
  
    /* matrix with snapshots values (DO WE NEED IT? TOO MUCH OF STORAGE?) */
    ublas::matrix<double> snaps_mat;
    /* matrix for pod basis */
    ublas::matrix<double> pod_basis;
    /* time_mean of the snapshots */
    ublas::vector<double> average_snaps;
    /* flag whether pod basis is computed from the fluctuating part of snapshots */
    bool pod_fluct;
  
  private:
    
    /* functions */
    
    /** 
    * @brief read snapshots
    *
    * Read snapshots from the member file snap_filename. This function is 
    * called by the constructor POD(string, int, bool)
    * The snapshots will be stored into the member matrix snaps_mat colomn-wise   
    *  
    */
    void read_snapshots();
    
    /**
    * @brief decompose snapshots 
    * Decompose the snapshots from the member file snap_filename
    * into the time average of the snapshots stored into the member variable
    * average_snaps and the fluctuating part of the snapshots stored into 
    * the member variable snaps_mat
    */
    void decompose_snaps();
    
    /**
    * @brief compute auto-correlation matrix
    *  
    * Compute autocorrelation matrix for POD computation which is computed as 
    * U^T*S*U, where U is the snapshot matrix (snaps_mat) and S is the inner 
    * product matrix (Gram_mat). If S is not specified (i.e. set_gramian() was 
    * not called), then the POD basis will be computed with respect to the 
    * Euclidian inner product with S = Id.  
    *
    * @param corr_mat resulting auto-correlation matrix
    */
    void compute_autocorr_mat( ublas::matrix<double> &corr_mat );
    
    /** 
    * @brief write time averages of snapshots into file
    * 
    * Write time average of the snapshots into the file
    * basename+avr. See also write_pod( string ).
    *
    * @param basename name of the file without extension
    */
    void write_averages( string basename );
    
    /** 
    * @brief write eigenvalue data into file
    * 
    * Write the POD eigenvalue data into the file
    * basename+svd. The data consists of three colomns with the numbers of the 
    * mode, corresponding eigenvalues, and the energy content, respectively.
    * See also write_pod( string ).
    *
    * @param basename name of the file without extension
    */
    void write_eigenvalues( string basename );
    
    /** 
    * @brief read time averages of snapshots
    *
    * Read the time averages of the snapshots from the file 
    * basename+avr. It will be automatically called in read_basis( string )
    * if pod_fluct==true.
    *
    * @param basename name of the file without extension 
    *  
    */
    void read_averages( string basename );
    
    /**
    * @brief read data from file
    *
    * Read data from the file _filename and save it into the matrix _mat
    *
    * @param _filename full name of the file
    * @param _mat matrix  
    */
    void read_data( string _filename, vector < vector<double> > &_mat );
    
    /**
    * @brief read data from file
    *
    * Read data from the file _filename and save it into the vector 'data'
    *
    * @param _filename full name of the file
    * @param data vector  
    */
    void read_data( string _filename, vector<double> &data);
};

#endif
