#include "Variational_multiscale.h"
#include "MooNMD_Io.h"
#include "list"

void matrices_a_tilde(std::vector<std::shared_ptr<FEMatrix>> &a, 
                        std::shared_ptr<const FEMatrix> g, std::shared_ptr<const FEMatrix> gt, 
                        std::shared_ptr<const FEMatrix> m, std::string flag)
{
  int nActive=0, nDof=0;
#ifdef __2D__
   nActive = a[0]->GetAnsatzSpace2D()->GetN_ActiveDegrees();
   nDof = a[0]->GetAnsatzSpace2D()->GetN_DegreesOfFreedom();
#else
   nActive = a[0]->GetAnsatzSpace3D()->GetN_ActiveDegrees();
   nDof = a[0]->GetAnsatzSpace3D()->GetN_DegreesOfFreedom();
#endif
  
  std::vector<double> value(nDof,0);
  int begin, end;
  double* entries_a11=nullptr;
  double* entries_a22=nullptr;
  double* entries_a33=nullptr;
  double* entries_a=nullptr;
  // matrix entries
  if(a.size() == 3)
  {
    entries_a11 = a[0]->GetEntries(); // entries of matrix a11
    entries_a22 = a[1]->GetEntries(); // entries of matrix a22
    entries_a33 = a[2]->GetEntries();// entries of matrix a33
  }
  else
  {
    entries_a = a[0]->GetEntries(); // off diagonal matrix a
  }
  // 
  const double* entries_gt = gt->GetEntries(); // entries of matrix G tilde
  const double* entries_g=g->GetEntries(); // entries of matrix G 
  const double* entries_m_vms = m->GetEntries(); // mass matrix from vms 
  for(int i=0; i<nActive; ++i)
  {
    // multiplication of matrices 
    {
      // i-th row of  Matrix tildeG
      begin = gt->GetRowPtr()[i];
      end = gt->GetRowPtr()[i+1];
      for(int j=begin; j<end; ++j)
      {
        // row of the matrix G11
        int index = gt->GetKCol()[j];
        int begin_1 = g->GetRowPtr()[index];
        int end_1 = g->GetRowPtr()[index+1];
        double val = entries_gt[j]/entries_m_vms[index];
        for(int k=begin_1; k<end_1; ++k)
        {
          int index1 = g->GetKCol()[k];
          value.at(index1) += val * entries_g[k];
        }
      }
    }
    // add matrix entries to A
    begin=a[0]->GetRowPtr()[i];
    end = a[0]->GetRowPtr()[i+1];
    for(int j=begin; j<end; ++j)
    {
      int index=a[0]->GetKCol()[j];
      if(flag.compare("a11_tilde") == 0)
      {
        entries_a11[j] -= value.at(index);
        entries_a22[j] -= value.at(index)/2.;
        entries_a33[j] -= value.at(index)/2.;
      }
      else if(flag.compare("a22_tilde") == 0)
      {
        entries_a11[j] -= value.at(index)/2.;
        entries_a22[j] -= value.at(index);
        entries_a33[j] -= value.at(index)/2.;
      }
      else if(flag.compare("a33_tilde") == 0)
      {
        entries_a11[j] -= value.at(index)/2.;
        entries_a22[j] -= value.at(index)/2.;
        entries_a33[j] -= value.at(index);
      }
      else
      {
        entries_a[j] -= value.at(index)/2.;
      }
      value.at(index) = 0.;
    }
  }
}
//====================================================================================

template<int d>
void LumpMassMatrixToDiagonalMatrix(std::shared_ptr<FEMatrix> & matrix)
{
  if(d==2)
    ErrThrow("Variational MultiScale method is not supported in 2D yet");
 double* entries = matrix->GetEntries();
  int* rowptr = matrix->GetRowPtr();
  int* kcol = matrix->GetKCol();
  
  int rows = matrix->GetN_Rows();
  for(int i=0; i<rows; ++i)
  {
    int begin = rowptr[i];
    int end   = rowptr[i+1];
    for(int j=begin; j<end; ++j)
    {
      // diagonal entry 
      if(kcol[j] == i)
      {
        rowptr[i] = i;
        entries[i] = entries[j];
        kcol[i] = i;
        break;
      }
    }
    if(fabs(entries[i]) < 1e-10)
      entries[i] = 1.;
  }
  rowptr[rows] = rows; 
}

template<int d>
void VMS_ProjectionUpdateMatrices(std::vector< std::shared_ptr< FEMatrix > >& blocks, 
                     std::vector< std::shared_ptr< FEMatrix > > matrices_vms)
{
  if(d==2)
    ErrThrow("Variational MultiScale method is not supported in 2D yet");
  std::vector<std::shared_ptr<FEMatrix>> A(3);
  A.at(0) = blocks.at(0); // A11 
  A.at(1) = blocks.at(5); // A22 
  A.at(2) = blocks.at(10);// A33
  
  // std::shared_ptr<FEMatrix> M    = matrices_vms.at(6);
  auto M = std::next(matrices_vms.begin(), 6); 
  //GT = matrices_vms.at(0);
  auto GT = std::next(matrices_vms.begin(), 0); 
  // std::shared_ptr<FEMatrix> G    = matrices_vms.at(3);  
  auto G = std::next(matrices_vms.begin(), 3); 
  // diagonal A tilde matrices 
  // (tilde_A11, tilde_A22, tilde_A33)
  matrices_a_tilde(A, *G, *GT, *M, "a11_tilde");
  matrices_a_tilde(A, *G, *GT, *M, "a22_tilde");
  matrices_a_tilde(A, *G, *GT, *M, "a33_tilde");
  
// off diagonal A tilde matrices
  // tilde_A12
  A.resize(1); A.at(0) = blocks.at(1); // A12 
  GT = std::next(matrices_vms.begin(), 1); 
  G = std::next(matrices_vms.begin(), 3); 
  matrices_a_tilde(A, *G, *GT, *M, "a12_tilde");
  // tilde_A1213
  A.resize(1); A.at(0)=blocks.at(2);
  GT = std::next(matrices_vms.begin(), 2); 
  G = std::next(matrices_vms.begin(), 3); 
  matrices_a_tilde(A, *G, *GT, *M, "a13_tilde");
  //tilde_A21
  A.resize(1); A.at(0)=blocks.at(4);
  GT = std::next(matrices_vms.begin(), 0); 
  G = std::next(matrices_vms.begin(), 4); 
  matrices_a_tilde(A, *G, *GT, *M, "a21_tilde");
  //tilde_A23
  A.resize(1); A.at(0)=blocks.at(6);
  GT = std::next(matrices_vms.begin(), 2); 
  G = std::next(matrices_vms.begin(), 4); 
  matrices_a_tilde(A, *G, *GT, *M, "a23_tilde");
  // tilde_A31
  A.resize(1); A.at(0)=blocks.at(8);
  GT = std::next(matrices_vms.begin(), 0); 
  G = std::next(matrices_vms.begin(), 5); 
  matrices_a_tilde(A, *G, *GT, *M, "a31_tilde");
  // tilde_A32
  A.resize(1); A.at(0)=blocks.at(9);
  GT = std::next(matrices_vms.begin(), 1); 
  G = std::next(matrices_vms.begin(), 5); 
  matrices_a_tilde(A, *G, *GT, *M, "a32_tilde"); 
}

#ifdef __2D__
template void LumpMassMatrixToDiagonalMatrix<2>(std::shared_ptr<FEMatrix> & matrix);
template void VMS_ProjectionUpdateMatrices<2>(std::vector< std::shared_ptr< FEMatrix > >& blocks,
                     std::vector< std::shared_ptr< FEMatrix > > matrices_vms);
#endif
#ifdef __3D__
template void LumpMassMatrixToDiagonalMatrix<3>(std::shared_ptr<FEMatrix> & matrix);
template void VMS_ProjectionUpdateMatrices<3>(std::vector< std::shared_ptr< FEMatrix > >& blocks,
                     std::vector< std::shared_ptr< FEMatrix > > matrices_vms);
#endif
