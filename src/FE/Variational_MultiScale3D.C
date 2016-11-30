#include <Variational_MultiScale3D.h>
#include <LinAlg.h>

void matrices_a_tilde(std::vector<std::shared_ptr<FEMatrix>> &a, 
                        std::shared_ptr<FEMatrix> g, std::shared_ptr<FEMatrix> gt, 
                        std::shared_ptr<FEMatrix> m, std::string flag)
{
  int nActive = a[0]->GetAnsatzSpace3D()->GetN_ActiveDegrees();
  int nDof = a[0]->GetAnsatzSpace3D()->GetN_DegreesOfFreedom();
  
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
  double* entries_gt = gt->GetEntries(); // entries of matrix G tilde
  double* entries_g=g->GetEntries(); // entries of matrix G 
  double* entries_m_vms = m->GetEntries(); // mass matrix from vms 
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
void LumpMassMatrixToDiagonalMatrix3D(std::shared_ptr< FEMatrix >& matrix)
{
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

//====================================================================================
void VMS_ProjectionUpdateMatrices3D(std::vector<std::shared_ptr<FEMatrix>> &blocks,
  std::array< std::shared_ptr< FEMatrix >, int(7) > matrices_vms)
{
  std::vector<std::shared_ptr<FEMatrix>> A(3);
  A.at(0) = blocks.at(0); // A11 
  A.at(1) = blocks.at(5); // A22 
  A.at(2) = blocks.at(10);// A33
  
  std::shared_ptr<FEMatrix> M    = matrices_vms.at(6);
  
  std::shared_ptr<FEMatrix> GT = matrices_vms.at(0);
  std::shared_ptr<FEMatrix> G    = matrices_vms.at(3);  
  // diagonal A tilde matrices 
  // (tilde_A11, tilde_A22, tilde_A33)
  matrices_a_tilde(A, G, GT, M, "a11_tilde");
  GT = matrices_vms.at(1);
  G    = matrices_vms.at(4);
  matrices_a_tilde(A, G, GT, M, "a22_tilde");
  GT = matrices_vms.at(1);
  G    = matrices_vms.at(4);
  matrices_a_tilde(A, G, GT, M, "a33_tilde");
  
  // off diagonal A tilde matrices
  // tilde_A12
  A.resize(1); A.at(0) = blocks.at(1); // A12 
  GT = matrices_vms.at(1); // GT22 
  G = matrices_vms.at(3); // G11 
  matrices_a_tilde(A, G, GT, M, "a12_tilde");
  // tilde_A1213
  A.resize(1); A.at(0)=blocks.at(2);
  GT=matrices_vms.at(2); // GT33
  G=matrices_vms.at(3);  // G11
  matrices_a_tilde(A, G, GT, M, "a13_tilde");
  //tilde_A21
  A.resize(1); A.at(0)=blocks.at(4);
  GT=matrices_vms.at(0);// GT11
  G=matrices_vms.at(4); // G22
  matrices_a_tilde(A, G, GT, M, "a21_tilde");
  //tilde_A23
  A.resize(1); A.at(0)=blocks.at(6);
  GT=matrices_vms.at(2); // GT33
  G=matrices_vms.at(4);  // G22
  matrices_a_tilde(A, G, GT, M, "a23_tilde");
  // tilde_A31
  A.resize(1); A.at(0)=blocks.at(8);
  GT=matrices_vms.at(0); // GT11
  G=matrices_vms.at(5);  // G22
  matrices_a_tilde(A, G, GT, M, "a31_tilde");
  // tilde_A32
  A.resize(1); A.at(0)=blocks.at(9);
  GT=matrices_vms.at(1); // GT22
  G=matrices_vms.at(5);  // G33
  matrices_a_tilde(A, G, GT, M, "a32_tilde"); 
}

void ComputeVMSProjection(
  const std::array< std::shared_ptr< FEMatrix >, int(4) > matrixG, 
  const TFEVectFunct3D& velocity, TFEVectFunct3D& projection)
{
  double* proj = projection.GetValues();
  double* u1=velocity.GetComponent(0)->GetValues();
  double* u2=velocity.GetComponent(1)->GetValues();
  double* u3=velocity.GetComponent(2)->GetValues();
  
  int n=matrixG.at(3)->GetN_Rows();
  std::vector<double> value(n, 0.);
  // G11 * u1
  matrixG.at(0)->multiply(u1, proj,     1.0);
  // G22 * u2
  matrixG.at(1)->multiply(u2, proj+3*n, 1.0);
  // G33 * u3 
  matrixG.at(2)->multiply(u3, proj+5*n, 1.0);
  
  // off diagonals
  matrixG.at(0)->multiply(u2, proj, 1.0);
  matrixG.at(1)->multiply(u1, &value[0], 1.0);
  Daxpy(n, 1., &value[0], proj+n);
  
  matrixG.at(0)->multiply(u3, proj+2*n, 1.0);
  matrixG.at(2)->multiply(u1, &value[0], 1.0);
  Daxpy(n, 1., &value[0], proj + 2*n);
  
  matrixG.at(1)->multiply(u3, proj+4*n, 1.0);
  matrixG.at(2)->multiply(u2, &value[0], 1.0);
  Daxpy(n, 1., &value[0], proj+4*n);
  
  if(matrixG.at(3)->GetRowPtr()[n] == n)
  {
    for(int i=0; i<6; ++i)
    {
      for(int j=0; j<n; ++j)
        proj[j+i*n] /= -matrixG.at(3)->GetEntries()[j];
    }
  }
  else
  {
    ErrThrow("not a diagonal matrix " );
  }
}