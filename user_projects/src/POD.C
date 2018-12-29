/** ************************************************************************ 
*
* @name       POD
* @brief      Computation/manipulation of POD basis independently of problem dimension
*
* @author     Swetlana Giere
* @date       08.03.2017 (start of implementation)
*
****************************************************************************/

#include <POD.h>

extern "C"
{
    void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda,
        double* w, double* work, int* lwork, int* info );
}

/** ***********************************************************************/
ParameterDatabase get_default_rom_parameters()
{
  Output::print<5>("Creating a default POD-ROM parameter database...");
  ParameterDatabase db("Default ParMooN parameter database for POD-based ROM problems");

  db.add("pod_directory", ".",
      	         "This directory is where the POD basis and Co. are written. This "
      	         "directory will be created, if it does not exist already. Files in "
      	         "this directory will be overwritten without any warning.");

  db.add("pod_basename", "parmoon_pod",
  		  "Basename for file where the pod basis and Co. are stored.");

  db.add("pod_rank", (size_t) 0,
      	         "This integer specifies the dimension of the POD space to be computed. "
      		 "If pod_rank <= 0, then all possible POD modes will be computed. ");

  db.add("pod_fluct", true,
               "This is the flag whether the POD basis should be computed from the "
  		  "part of the snapshots (also also central trajectory method).",
    		   {true,false});

  db.add("pod_inner_product", "eucl",
               "Specification of the inner product which is used to compute POD basis."
  		  "Besides default value, only 'l2' is possible at the moment.",
  		  {"eucl", "l2"});

  db.add("rom_init_regularized", false,
  		 "This is the flag whether the the ROM initial condition should be regularized.",
  		 {true,false});

  db.add("differential_filter_width", 1.0,
           "Filter width for the differential filter (Helmoltz equation) for the computation "
           "of of the regularized ROM initial condition.",
             0., 10.);

  // Merge with other databases
  db.merge(ParameterDatabase::default_output_database(), true);

  return db;
}

/** ***********************************************************************/
POD::POD(const ParameterDatabase& param_db) :
		rom_db(get_default_rom_parameters()),
		length(0),
		rank(0),
		length_snaps(0),
		number_snaps(0),
		eigen_threshold(1.0e-10),
		valid_eigs(0),
		eigs(NULL)
{
  this->rom_db.merge(param_db, true);
  rank = rom_db["pod_rank"];
  //TODO set gramian_mat, snaps_mat, pod_basis, snaps_mean to zero
}

/** ***********************************************************************/
POD::~POD() {
  if( eigs != NULL )
  {
    delete[] eigs;
  }
}

/** ***********************************************************************/
void POD::read_snapshots() {

  string snap_filename = this->rom_db["snaps_directory"].get<std::string>();
  snap_filename += "/";
  snap_filename += this->rom_db["snaps_basename"].get<std::string>();
  snap_filename += "snap";

  Output::print<1>("Reading snapshots from file: ", snap_filename);

  std::vector < std::vector<double> > tmp_snaps;
  read_data( snap_filename, tmp_snaps);
  
  this->length_snaps = tmp_snaps[0].size(); //total dof
  this->number_snaps = tmp_snaps.size();
  
  this->snaps_mat.resize(this->length_snaps,this->number_snaps);
  
  // store snapshots into member matrix snaps_mat 
  for(int i=0; i<this->number_snaps;i++)
  {
    for(int j=0; j<this->length_snaps;j++)
    	this->snaps_mat(j,i)=tmp_snaps[i][j];
  }
  
  Output::print<1>("Length of snapshots : ", this->length_snaps);
  Output::print<1>("Number of snapshots : ", this->number_snaps);
}

/** ***********************************************************************/
void POD::compute_basis() {
  Output::print<1>("Computing POD basis by method of snapshots...");

  read_snapshots();

  this->snaps_mean.resize(this->length_snaps);
  this->snaps_mean.clear();
  
  if( this->rom_db["pod_fluct"] ) {
    Output::print<1>("POD will be computed from fluctuating part of snapshots.");
    decompose_snaps();
  }
  else            
    Output::print<1>("POD will be computed from full/raw snapshots.");

  ublas::matrix<double> autocorr_mat;
  compute_autocorr_mat(autocorr_mat);

  ublas::vector<double> vec_eig, vec_Seig;
  int max_rank;
  
  /* parameters for LAPACK dsyev */
  double** snaps;
  double* work;
  double wkopt;
  char arg1 = 'V'; /* compute eigenvalues and eigenvectors */
  char arg2 = 'U'; /* for upper triangular matrix */
  int info = 0, lwork;
  int n = (int) this->number_snaps;
  if(this-> eigs != NULL ) delete[] this->eigs;
  /* array for eigenvalues */
  this->eigs = new double[ this->number_snaps ];
  
  snaps = new double*[ this->number_snaps ];
  snaps[ 0 ] = new double[ this->number_snaps*number_snaps ];

  for( int i = 0; i< this->number_snaps; i++ ) {
    snaps[ i ] = snaps[ 0 ]+ i*this->number_snaps;
  }

  for( int i = 0 ; i < this->number_snaps; i++ )
  {
    for( int j = 0; j < this->number_snaps; j++ )
    {
      snaps[ j ][ i ] = autocorr_mat( i, j );
    }
  } 
  
  Output::print<1>("Calling lapack... ");
  lwork = -1;
  dsyev_( &arg1, &arg2, &n, *snaps, &n, this->eigs, &wkopt, &lwork, &info );
  lwork = (int)wkopt;
  work = (double*)malloc( lwork*sizeof(double) );
  /* Solve eigenproblem */
  dsyev_( &arg1, &arg2, &n, *snaps, &n, this->eigs, work, &lwork, &info );
  /* Check for convergence */
  if( info > 0 )
  {
    ErrThrow("The algorithm failed to compute eigenvalues.");
  }
  
  for( int i=this->number_snaps-1; i >=0 ; i-- ) {
    if( this->eigs[ i ] > this->eigen_threshold ) {
    	this->valid_eigs ++;
    }
  }
  
  /* check the rank of the pod basis */
  
  if(this->rank<=0)
	this->rank=valid_eigs; // take all
  
  if(this->rank > this->valid_eigs)
  {
    ErrThrow("ERROR: rank of POD basis can't be greater than ", this->valid_eigs, "! "
    		"Check parameter 'pod_rank'.");
  }
  
  Output::print<1>("Number of valid POD eigenvalues: ", this->valid_eigs);
  Output::print<1>("Dimension of POD basis         : ", this->rank);

  this->pod_basis.resize( this->length_snaps, this->rank );
  vec_eig.resize( this->number_snaps );
  vec_Seig.resize( this->length_snaps );

  for( int i = this->number_snaps-1; i >= this->number_snaps - this->rank && i >=0 ; i-- )
  { 
    int count = this->number_snaps - 1 - i;
    for( int j = 0; j < this->number_snaps ; j++ ) {
      vec_eig( j ) = 1.0 / sqrt( this->eigs[ i ]) * snaps[ i ][ j ];
    }
    axpy_prod( this->snaps_mat, vec_eig, vec_Seig, true );
    for( int j = 0; j < this->length_snaps ; j++ ) {
    	this->pod_basis( j, count ) = vec_Seig( j );
    }
  } 
  for( int i= 0 ; i < this->number_snaps/2 ; i++ ){
    double help = this->eigs[ i ];
    this->eigs[ i ] = this->eigs[ this->number_snaps-1-i ] ;
    this->eigs[ this->number_snaps-1-i ] = help;
  }
  
  this->length = this->length_snaps;
  delete [] *snaps;
  delete [] snaps;
  delete [] work;

}

/** ***********************************************************************/
void POD::set_gramian( std::shared_ptr<TMatrix> mat)
{
	convert_to_ublas(mat, this->gramian_mat);
}

/** ***********************************************************************/
void POD::convert_to_ublas( std::shared_ptr<TMatrix> mat,
    		                ublas::compressed_matrix<double> & res_mat )
{
	int nrow = mat->GetN_Rows();
	int ncol = mat->GetN_Columns();
	int nvals = mat->GetN_Entries();

	const int* col_idx = mat->GetKCol();
	const int* row_ptr = mat->GetRowPtr();
	const double* vals = mat->GetEntries();

	int begin, end, pos = 0;

	res_mat.resize(nrow, ncol, false);
	for (int i=0; i<nrow; ++i)
	{
		begin = row_ptr[i];
		end   = row_ptr[i+1];
		for (int j=begin; j<end; ++j)
		{
			res_mat(i,col_idx[pos]) = vals[pos];
			++pos;
		}
	}
}

/** ***********************************************************************/
void POD::write_pod( std::string basename ) {
  
  /* writing into basis into file */
  std::string pod_filename = this->rom_db["pod_directory"].get<std::string>();
  pod_filename += "/";
  pod_filename += basename;
  pod_filename += "pod";
  std::ofstream ofile;
  ofile.open( pod_filename.c_str() , ios::out | ios::trunc );
  if( !ofile.is_open() )
  {
    ErrThrow( "Error: File for POD basis ", pod_filename," could not be created." );
  }
  ofile << setprecision( 12 );
  Output::print<1>( "Writing POD basis into file: ", pod_filename );
  if(rom_db["pod_fluct"]){
    ofile << "POD basis: fluctuating field" << "\n";
    /* write averages of snapshots into file */
    write_averages( basename );
  }
  else{
    ofile << "POD basis: full field" << "\n";
  }

  for( int i = 0; i< pod_basis.size2() ; i++ )
  {
    for( int j = 0; j < pod_basis.size1(); j++ )
      ofile << pod_basis( j, i ) << " ";
    ofile << "\n";
  }
  ofile.close();
 
  /*write pod eigenvalues into file */
  write_eigenvalues( basename );
}

/** ***********************************************************************/
void POD::read_basis() {
  std::string filename = this->rom_db["pod_directory"].get<std::string>();
  filename += "/" + this->rom_db["pod_basename"].get<std::string>();
  filename += "pod";
  vector < vector<double> > tmp_basis;
  ifstream podfile(filename.c_str());
  string line;
  int a_length;
  if( !podfile )
  {
    Output::print<1>("Error: POD file ", filename ," could not be opened.");
    exit(4711);
  }
  getline(podfile,line); // get the first line of file filename
  if(!line.compare("POD basis: fluctuating field"))
  {
    Output::print<1>("POD basis is computed out of fluctuating part of snapshots!");
    read_averages();
  }
  else if(!line.compare("POD basis: full field"))
  {
	Output::print<1>( "POD basis is computed out of entire snapshots!");
  }
  else{
	ErrThrow("No information about type of building up of POD basis. "
			 "Unsupported format.\n First row of file should contain "
			 "'POD basis: full field' or 'POD basis: fluctuating field'.");
  }
  Output::print<1>("Reading POD basis from file ", filename);
  while( getline(podfile, line) ){
    vector<double> data;
    double value;
    istringstream iss(line);
      while (iss >> value){  
      data.push_back(value);
      }
      tmp_basis.push_back(data);
  }
  podfile.close();

  if( this->rank <= 0 )
	  this->rank=tmp_basis.size();
  
  if( this->rank > tmp_basis.size() )
  {
	Output::warn<1>("Rank of POD basis can't be greater than ", tmp_basis.size(), "! "
			        "Parameter 'pod_rank' changed to ", tmp_basis.size());
	this->rank = tmp_basis.size();
  }
  
  this->length = tmp_basis[0].size(); //total dof
  this->pod_basis.resize(this->length,this->rank);
  
  Output::print<1>("Length of POD basis     : ", this->length);
  Output::print<1>("Dimension of POD basis  : ", this->rank);
  
  //pod basis (for all components)
  for(int i=0; i<this->rank;i++)
    for(int j=0; j<this->length;j++)
      this->pod_basis(j,i)=tmp_basis[i][j];
}

/** ***********************************************************************/
void POD::write_data( ublas::matrix<double> &mat, std::string filename) {
  std::ofstream ofile;
  ofile.open( filename.c_str() , ios::out | ios::trunc );
  if( !ofile.is_open() )
  {
    ErrThrow("Error: File ", filename, " could not be created!");
  }
  /** write elements */
  for(int i=0;i< mat.size2();i++)
  {
    for(int j=0;j<mat.size1();j++)
      ofile << mat(j,i) << " ";
    ofile << "\n";
  }
  ofile.close();
}

/* used only internally by the class */
void POD::write_averages( string basename ) {
  std::string avr_filename = this->rom_db["pod_directory"].get<std::string>();
  avr_filename += "/";
  avr_filename += basename;
  avr_filename += "mean";
  ofstream ofile;
  ofile.open( avr_filename.c_str() , ios::out | ios::trunc );
  if( !ofile.is_open() ) {
    ErrThrow("Error: File for averages of snapshots ", avr_filename," could not be created.");
  }
  ofile << setprecision( 12 );
  
  if(this->snaps_mean.size()!=this->length)
  {
    ErrThrow( "Error: Vector for averages of snapshots has the wrong length!" );
  }
  // write elements
  for(int i=0; i<this->snaps_mean.size() ; ++i)
  {
      ofile << this->snaps_mean(i) << " ";
  }
  ofile.close();
}

/** ***********************************************************************/
void POD::write_eigenvalues( string basename ) {
  std::string eig_filename = this->rom_db["pod_directory"].get<std::string>();
  eig_filename += "/";
  eig_filename += basename;
  eig_filename += "eigs";
  int valid_eigen = 0;
  ofstream ofile;
  double sum_eigen = 0.0, cumulative=0.0;
  int no_eigen = snaps_mat.size2();
  if( this->eigs == NULL ) {
    ErrThrow("Error: Eigenvalues not found. POD not created?");
  }
  for( int i = 0 ; i< no_eigen ; i ++ ) {
    if( this->eigs[ i ] >= this->eigen_threshold ) {
      valid_eigen ++;
    }
  }
  ofile.open( eig_filename.c_str() , ios::out | ios::trunc );
  ofile << setprecision( 12 );
  for( int i = 0 ; i< no_eigen ; i++ ) {
    sum_eigen += this->eigs[ i ];
  }
  for( int i = 0 ; i < min(no_eigen,valid_eigen); i++ ) {
    cumulative += this->eigs[ i ] / sum_eigen;
    ofile << i+1 << " "<< this->eigs[ i ] << " \t" << cumulative << endl;
  }
  ofile.close();
}

/** ***********************************************************************/
void POD::read_averages() {
  std::string filename = this->rom_db["pod_directory"].get<std::string>();
  filename += "/";
  filename += this->rom_db["pod_basename"].get<std::string>();
  filename += "mean";
  Output::print<1>("Reading averages of snapshots from file ", filename);

  vector<double> data;

  read_data( filename , data );
  this->snaps_mean.resize(data.size());

  for(int i=0 ; i<this->snaps_mean.size() ; ++i)
	  this->snaps_mean( i ) = data[ i ];
}

/** ***********************************************************************/
void POD::decompose_snaps() {
  if( this->length_snaps == 0 ) {
    ErrThrow("Error: Snapshots are not available. Before computing POD basis "
    		 "read_snapshots() has to be called!");
  }
  this->snaps_mean.resize(this->length_snaps);
  this->snaps_mean.clear();

  /* compute snaps' mean */
  for(int i = 0; i < length_snaps; i++)
  {
	//get ith row of snapshot matrix
    ublas::matrix_row<ublas::matrix<double>> snap_row(this->snaps_mat,i);
    this->snaps_mean(i) = (1./this->number_snaps) * sum(snap_row);//compute the time-mean in ith node
  }
  
  /* substract the mean from snapshots */
  for(int j = 0; j < this->number_snaps; j++) {
    ublas::matrix_column<ublas::matrix<double>> snap_column(this->snaps_mat,j);
    snap_column -= this->snaps_mean;
  }     
}

/** ***********************************************************************/
void POD::compute_autocorr_mat( ublas::matrix<double> &corr_mat ) {
  
  Output::print<1>("Computing autocorrelation matrix...");
  
  corr_mat.resize( this->snaps_mat.size2(),  this->snaps_mat.size2());

  if( this->length_snaps == 0 )
  {
    ErrThrow("Error: Snapshots are not available! Before computing POD "
    		 "basis read_snapshots() must be called!");
  }
  
  if( this->rom_db["pod_inner_product"].get<std::string>() == "eucl" ){
    noalias(corr_mat) = prod( trans( this->snaps_mat ), this->snaps_mat );
  }
  else{
    if( this->gramian_mat.size1() != this->length_snaps )
    {
      ErrThrow("ERROR: Dimension of inner product matrix and length of "
               "snapshots must coincide!\nCurrently: dim of matrix : ",
               this->gramian_mat.size1(), " length of snapshots : ",
               this->length_snaps);
    }
    ublas::matrix<double> tmp_mat (this->gramian_mat.size1(), this->snaps_mat.size2());
    noalias(tmp_mat)  = ublas::prod( this->gramian_mat, this->snaps_mat );
    noalias(corr_mat) = ublas::prod(trans(this->snaps_mat), tmp_mat);
  } // end else
  Output::print<1>("Computation of autocorrelation matrix completed...");
}

/** ***********************************************************************/
void POD::read_data( string _filename, vector < vector<double> > &_mat ) {
  
  ifstream file(_filename.c_str());
  if( !file ){
    ErrThrow("Error: File ", _filename, " could not be opened");
  }
  string line;
  while( getline(file, line) )
  { 
    vector<double> data;
    double value;
    istringstream iss(line);
    while (iss >> value)
    {
      data.push_back(value);
    }
    _mat.push_back(data);
  }
  file.close();
}

/** ***********************************************************************/
void POD::read_data( string _filename, vector<double> &data) {
  
  string line;
  ifstream file(_filename.c_str());
  if( !file )
  {
    ErrThrow("Error: File ", _filename, " could not be opened.");
  }

  while( getline(file, line) ){
    double value;
    istringstream iss(line);
    while (iss >> value){  
    data.push_back(value);
    }
  }
  file.close();
}

