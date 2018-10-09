#include "LPS_scott_zhang.h"
#include "Enumerations.h"
#include "FEDatabase2D.h"
#include "FEMatrix.h"
#include "FEFunction2D.h"

double ComputePSPGParameterLPS(double hK, bool velocity, double nu,
                               size_t lps_coeff_type, double delta0,
                               double delta1)
{
  double val;
  switch(lps_coeff_type)
  {
    // [DGJN17], Section 3
    case 1:
      val = delta0 * hK * hK / nu;
      break;
    // [DGJN17], Section 4, implicit
    case 2:
    // [DGJN17], Section 4, explicit
    case 3:
      if ( velocity )
        val = delta1;
      else
        val = delta0 * hK * hK / nu;
      break;
    // [DGJN17], Section 6, implicit
    case 4:
    // [DGJN17], Section 6, explicit
    case 5:
      if ( velocity )
        val = delta1 * hK / nu;
      else
        val = delta0 * hK / nu;
      break;
    default:
      val = delta0 * hK * hK / nu;
      break;
  }
  return val;
}

// ======================================================================
// implementation like Badia, CMAME 2012
// ======================================================================
std::shared_ptr<FEMatrix> LPS_for_pressure_Scott_Zhang(
  std::shared_ptr<const FEMatrix> C, bool velocity, double nu,
  LPS_parameter_set lps_ps)
{
  int end, i, index, index_j, index_k, j, jj, k, l, N_, N_Entries, N_P, N_Cells, start;
  int N_Entries_new, index_jj, N_neigh_j, neigh_no_j, index_ansatz, index_test;
  int N_neigh_k, neigh_no_k, kk, index_kk, N_Points, N_Edges, N_Active;
  const int *DOF, *GlobalNumbers, *BeginIndex;
  int *N_BaseFunct;
  const int *DOF_neigh_j, *DOF_neigh_k;
  double *EntriesC_new;
  double *weights, *xi, *eta, *Orig, **OrigFEValues;
  double X[MaxN_QuadPoints_2D], Y[MaxN_QuadPoints_2D], x[4], y[4];
  double AbsDetjk[MaxN_QuadPoints_2D], value, value_grad_grad, value_fct_gradx, value_fct_grady, h, delta;
  double val[6];
  TCollection *coll;
  FE2D CurrentElement, CurrentElement_neigh_j, CurrentElement_neigh_k;
  TBaseCell *cell, *neigh_j, *neigh_k;
  BaseFunct2D *BaseFuncts;
  BaseFunct2D BaseFunct;
  MultiIndex2D der[3] = {D00, D10, D01};

  int MaxN_EntriesPerRow = 200;
  int N_MaxBaseFuncts = 25;
  int N_MaxQuadPoints = 100;
  double values_in_quad_points[N_MaxBaseFuncts][3][N_MaxQuadPoints];
  double integral_fct_fct[N_MaxBaseFuncts][N_MaxBaseFuncts];
  double integral_grad_grad[N_MaxBaseFuncts][N_MaxBaseFuncts];
  double integral_fct_gradx[N_MaxBaseFuncts][N_MaxBaseFuncts];
  double integral_fct_grady[N_MaxBaseFuncts][N_MaxBaseFuncts];
  double X_dof[N_MaxBaseFuncts], Y_dof[N_MaxBaseFuncts];

  if ( !velocity )
  {
    OutPut ( "LPS for pressure (Scott/Zhang)" << endl );
  }
  else
  {
    OutPut ( "LPS for velocity (Scott/Zhang)" << endl );
  }
  const TFESpace2D *pressure_space = C->GetFESpace2D();
  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFunct = TFEDatabase2D::GetN_BaseFunctFromFE2D();
  // get arrays with the numbering of the dof
  GlobalNumbers = pressure_space->GetGlobalNumbers();
  BeginIndex = pressure_space->GetBeginIndex();

  // allocate array for node-to-cell map
  N_P = pressure_space->GetN_DegreesOfFreedom();
  N_Active = pressure_space->GetActiveBound();
  std::vector<int> node_to_cell(N_P, -1);

  // get collection and number of cells
  coll = pressure_space->GetCollection();
  N_Cells = coll->GetN_Cells();

//   OutPut ( "set node-to-cell map" << endl );
  // loop over all cells
  for ( i = 0; i < N_Cells; i++ )
  {
    cell = coll->GetCell ( i );
    DOF = GlobalNumbers + BeginIndex[i];
    CurrentElement = pressure_space->GetFE2D ( i, cell );
    // # basis functions
    N_ = N_BaseFunct[CurrentElement];
    if ( N_ > N_MaxBaseFuncts )
    {
      OutPut ( "LPS_for_pressure_Scott_Zhang: N_MaxBaseFuncts too small !! " << N_ << endl );
      exit ( 4711 );
    }
    for ( j = 0; j < N_; j++ )
    {
      index = DOF[j];
      if ( node_to_cell[index] < 0 )
        // assign the cell number
        node_to_cell[index] = i;
    }
  }
// for (i=0;i<N_P;i++)
//   OutPut("ntc " << i << " " << node_to_cell[i] << endl);

  std::vector<int> RowPtr = C->get_row_array();
  N_Entries = C->GetN_Entries();
//   OutPut ( "original matrix " << N_Entries << endl )

  // compute new sparsity pattern
  // allocate array for storage
  std::vector<int> EntriesPerRow(N_P * MaxN_EntriesPerRow, -1);
  N_Entries_new = 0;

  for ( i = 0; i < N_Cells; i++ )
  {
    cell = coll->GetCell ( i );
    DOF = GlobalNumbers + BeginIndex[i];
    CurrentElement = pressure_space->GetFE2D ( i, cell );
    // # basis functions
    N_ = N_BaseFunct[CurrentElement];
    // loop over the dofs
    for ( j = 0; j < N_; j++ )
    {
      index_j = DOF[j];
      // the corresponding row
      start = index_j * MaxN_EntriesPerRow;
      // loop over all dofs of same cell
      for ( k = 0; k < N_; k++ )
      {
        index_k = DOF[k];
        // check whether the index index_j*MaxN_EntriesPerRow;is already in the list
        for ( l = start; l < start + MaxN_EntriesPerRow; l++ )
        {
          if ( EntriesPerRow[l] == -1 )
          {
            EntriesPerRow[l] = index_k;
            N_Entries_new++;
            break;
          }
          if ( EntriesPerRow[l] == index_k )
            break;
        }
        if ( l == start + MaxN_EntriesPerRow )
        {
          OutPut ( "LPS_for_pressure_Scott_Zhang: MaxN_EntriesPerRow too small !!!" << endl );
          exit ( 4711 );
        }
      }

      // take the cell that is assigned to index_j
      neigh_no_j = node_to_cell[index_j];
      neigh_j = coll->GetCell ( neigh_no_j );
      DOF_neigh_j = GlobalNumbers + BeginIndex[neigh_no_j];
      CurrentElement_neigh_j = pressure_space->GetFE2D ( neigh_no_j, neigh_j );
      N_neigh_j = N_BaseFunct[CurrentElement_neigh_j];
      // loop over the neighbor dofs
      // this loop is just for safety, it should not give new entries
      for ( jj = 0; jj < N_neigh_j; jj++ )
      {
        index_jj = DOF_neigh_j[jj];
        // the indices of this cell might be coupled to index_k
        for ( k = 0; k < N_; k++ )
        {
          index_k = DOF[k];
          start = index_k * MaxN_EntriesPerRow;
          // check whether the index is already in the list
          for ( l = start; l < start + MaxN_EntriesPerRow; l++ )
          {
            if ( EntriesPerRow[l] == -1 )
            {
              EntriesPerRow[l] = index_jj;
              N_Entries_new++;
              break;
            }
            if ( EntriesPerRow[l] == index_jj )
              break;
          }
          if ( l == start + MaxN_EntriesPerRow )
          {
            OutPut ( "LPS_for_pressure_Scott_Zhang: MaxN_EntriesPerRow too small !!!" << endl );
            exit ( 4711 );
          }

          // now check the transposed entry
          start = index_jj * MaxN_EntriesPerRow;
          // check whether the index is already in the list
          for ( l = start; l < start + MaxN_EntriesPerRow; l++ )
          {
            if ( EntriesPerRow[l] == -1 )
            {
              EntriesPerRow[l] = index_k;
              N_Entries_new++;
              break;
            }
            if ( EntriesPerRow[l] == index_k )
              break;
          }
          if ( l == start + MaxN_EntriesPerRow )
          {
            OutPut ( "LPS_for_pressure_Scott_Zhang: MaxN_EntriesPerRow too small !!!" << endl );
            exit ( 4711 );
          }
        }
      } // end jj

//       start =8*MaxN_EntriesPerRow;
//       for (l=start;l<start+MaxN_EntriesPerRow;l++)
//  OutPut(EntriesPerRow[l] << " " );
//       OutPut(endl);

      // check the far connections
      for ( k = 0; k < N_; k++ )
      {
        // take the cell that is assigned to index_k
        index_k = DOF[k];
        neigh_no_k = node_to_cell[index_k];
        neigh_k = coll->GetCell ( neigh_no_k );
        DOF_neigh_k = GlobalNumbers + BeginIndex[neigh_no_k];
        CurrentElement_neigh_k = pressure_space->GetFE2D ( neigh_no_k, neigh_k );
        N_neigh_k = N_BaseFunct[CurrentElement_neigh_k];
        // loop over the dofs assigned to index_j and index_k
        for ( jj = 0; jj < N_neigh_j; jj++ )
        {
          index_jj = DOF_neigh_j[jj];
          start = index_jj * MaxN_EntriesPerRow;
          for ( kk = 0; kk < N_neigh_k; kk++ )
          {
            index_kk = DOF_neigh_k[kk];
            // check whether the index is already in the list
            for ( l = start; l < start + MaxN_EntriesPerRow; l++ )
            {
              if ( EntriesPerRow[l] == -1 )
              {
                EntriesPerRow[l] = index_kk;
                N_Entries_new++;
                break;
              }
              if ( EntriesPerRow[l] == index_kk )
                break;
            }
            if ( l == start + MaxN_EntriesPerRow )
            {
              OutPut ( "LPS_for_pressure_Scott_Zhang: MaxN_EntriesPerRow too small !!!" << endl );
              exit ( 4711 );
            }
          } // end kk
        } // end jj
      } // end k
    } // end j
  } // end i
//   OutPut("modified matrix " << N_Entries_new << endl)

  ////////////////////////////////////////////
  // new column pointer and entry pointer
  std::vector<int> ColInd_new(N_Entries_new);
  //EntriesC_new = new double[N_Entries_new];
  //memset(EntriesC_new, 0.0, N_Entries_new*SizeOfDouble);
  // reset and fill the pointers
  // loop over all dofs
  l = 0;
  RowPtr[0] = 0;
  for ( i = 0; i < N_P; i++ )
  {
    start = i * MaxN_EntriesPerRow;
    k = 0;
    for ( j = start; j < start + MaxN_EntriesPerRow; j++ )
    {
      if ( EntriesPerRow[j] > -1 )
      {
        ColInd_new[l] = EntriesPerRow[j];
        l++;
        k++;
      }
      else
        break;
    }
    RowPtr[i + 1] = RowPtr[i] + k;
  }

// for (i=0;i<=N_P;i++)
//    OutPut(RowPtr[i] << " ");
//  OutPut(endl);
//  OutPut(N_Entries_new << endl );
//   for (i=0;i<=N_Entries_new;i++)
//     OutPut(ColInd_new[i] << " ");
//   OutPut(endl);
// OutPut(l << endl);
// define new matrix

  auto sqstructureC = std::make_shared<TStructure>(N_P, N_Entries_new,
                                                   ColInd_new.data(),
                                                   RowPtr.data());
  ColInd_new = sqstructureC->get_columns();
  //sqstructureC = new TSquareStructure2D(N_P, N_Entries_new, ColInd_new, RowPtr, NULL);
  auto C_new = std::make_shared<FEMatrix>(pressure_space, sqstructureC);
  EntriesC_new = C_new->GetEntries();

  char PString[] = "pressure";
  std::vector<double> tmp_pres(N_P, 0.);
  TFEFunction2D tmp_pressure(pressure_space, PString,  PString, tmp_pres.data(),
                             N_P);
  /////////////////////////////////////////////////
  // compute matrix entries
  std::vector<double> FEFunctValues(MaxN_EntriesPerRow);

  for ( i = 0; i < N_Cells; i++ )
  {
    cell = coll->GetCell ( i );
    h = cell->GetDiameter();
    delta = ComputePSPGParameterLPS(h, velocity, 1, lps_ps.lps_coeff_type,
                                    lps_ps.delta0, lps_ps.delta1);
    //delta = compute_PSPG_delta(0.1, h, 0.001);
    DOF = GlobalNumbers + BeginIndex[i];
    CurrentElement = pressure_space->GetFE2D ( i, cell );
    // # basis functions
    N_ = N_BaseFunct[CurrentElement];
    // basis functions
    BaseFunct = BaseFuncts[CurrentElement];

    bool SecondDer[1];
    SecondDer[0] = false;
    // get reference transformation, quadrature points
    TFEDatabase2D::GetOrig(1, &CurrentElement, coll, cell, SecondDer, N_Points,
                           xi, eta, weights, X, Y, AbsDetjk);

    if ( N_Points > N_MaxQuadPoints )
    {
      OutPut ( "LPS_for_pressure_Scott_Zhang: N_MaxQuadPoints too small !! " << N_Points << endl );
      exit ( 4711 );
    }

    // assign geometric positions to the dofs
    // compute vertices
    N_Edges = cell->GetN_Edges();
    for ( j = 0; j < N_Edges; j++ )
    {
      x[j] = cell->GetVertex ( j )->GetX();
      y[j] = cell->GetVertex ( j )->GetY();
    }

    switch ( CurrentElement )
    {
      case C_P1_2D_T_A:
        X_dof[0] = x[0];
        Y_dof[0] = y[0];
        X_dof[1] = x[1];
        Y_dof[1] = y[1];
        X_dof[2] = x[2];
        Y_dof[2] = y[2];
        break;
      case C_P2_2D_T_A:
        X_dof[0] = x[0];
        Y_dof[0] = y[0];
        X_dof[1] = ( x[0] + x[1] ) / 2.0;
        Y_dof[1] = ( y[0] + y[1] ) / 2.0;
        X_dof[2] = x[1];
        Y_dof[2] = y[1];
        X_dof[3] = ( x[0] + x[2] ) / 2.0;
        Y_dof[3] = ( y[0] + y[2] ) / 2.0;
        X_dof[4] = ( x[1] + x[2] ) / 2.0;
        Y_dof[4] = ( y[1] + y[2] ) / 2.0;
        X_dof[5] = x[2];
        Y_dof[5] = y[2];
        break;
      case C_P3_2D_T_A:
        X_dof[0] = x[0];
        Y_dof[0] = y[0];
        X_dof[1] = 2.0 * x[0] / 3.0 + x[1] / 3.0;
        Y_dof[1] = 2.0 * y[0] / 3.0 + y[1] / 3.0;
        X_dof[2] = x[0] / 3.0 + 2.0 * x[1] / 3.0;
        Y_dof[2] = y[0] / 3.0 + 2.0 * y[1] / 3.0;
        X_dof[3] = x[1];
        Y_dof[3] = y[1];
        X_dof[4] = 2.0 * x[0] / 3.0 + x[2] / 3.0;
        Y_dof[4] = 2.0 * y[0] / 3.0 + y[2] / 3.0;
        X_dof[5] = ( x[0] + x[1] + x[2] ) / 3.0;
        Y_dof[5] = ( y[0] + y[1] + y[2] ) / 3.0;
        X_dof[6] = 2.0 * x[1] / 3.0 + x[2] / 3.0;
        Y_dof[6] = 2.0 * y[1] / 3.0 + y[2] / 3.0;
        X_dof[7] = x[0] / 3.0 + 2.0 * x[2] / 3.0;
        Y_dof[7] = y[0] / 3.0 + 2.0 * y[2] / 3.0;
        X_dof[8] = x[1] / 3.0 + 2.0 * x[2] / 3.0;
        Y_dof[8] = y[1] / 3.0 + 2.0 * y[2] / 3.0;
        X_dof[9] = x[2];
        Y_dof[9] = y[2];
        break;
      default:
        ErrThrow("LPS_for_pressure_Scott_Zhang: element not implemented!");
    }

    for ( l = 0; l < N_; l++ )
      for ( j = 0; j < N_; j++ )
        integral_grad_grad[l][j] = integral_fct_fct[l][j] = integral_fct_gradx[l][j] = integral_fct_grady[l][j] = 0.0;

    for ( l = 0; l < N_; l++ )
      for ( j = 0; j < 3; j++ )
        for ( k = 0; k < N_Points; k++ )
          values_in_quad_points[l][j][k] = 0.0;

    // loop over each basis function
    for ( j = 0; j < N_; j++ )
    {
      // reset the function values
      for ( k = 0; k < N_; k++ )
        FEFunctValues[k] = 0;
      FEFunctValues[j] = 1.0;

      // compute values for all derivatives
      // in all quadrature points
      // in original mesh cell
      for ( k = 0; k < 3; k++ )         // for all derivatives
      {
        // get values in original cell
        OrigFEValues = TFEDatabase2D::GetOrigElementValues ( BaseFunct, der[k] );
        for ( jj = 0; jj < N_Points; jj++ )            // for all quadrature points
        {
          Orig = OrigFEValues[jj];                   // value in original cell
          value = 0;
          for ( l = 0; l < N_; l++ )                // for all basis functions
            value += FEFunctValues[l] * Orig[l];    // accumulate value of derivative in point j
          values_in_quad_points[j][k][jj] = value;
          //OutPut(i << " " << j <<" " << k << " " << value<<" :: ");
        }  // endfor jj
        //OutPut(endl);
      }                                             // endfor k
    }
    // values in the quad points are stored
    // compute integrals

    // loop over each basis function
    for ( j = 0; j < N_; j++ ) // test fct
    {
      for ( k = 0; k < N_; k++ )
      {
        value = value_grad_grad = value_fct_gradx = value_fct_grady = 0.0;
        // loop over quad points
        for ( jj = 0; jj < N_Points; jj++ )
        {
          value += values_in_quad_points[j][0][jj] * values_in_quad_points[k][0][jj] * weights[jj] * AbsDetjk[jj];
          value_grad_grad += ( values_in_quad_points[j][1][jj] * values_in_quad_points[k][1][jj]
                               + values_in_quad_points[j][2][jj] * values_in_quad_points[k][2][jj] ) * weights[jj] * AbsDetjk[jj];
          value_fct_gradx += values_in_quad_points[j][0][jj] * values_in_quad_points[k][1][jj] * weights[jj] * AbsDetjk[jj];
          value_fct_grady += values_in_quad_points[j][0][jj] * values_in_quad_points[k][2][jj] * weights[jj] * AbsDetjk[jj];
        }
        integral_grad_grad[j][k] += value_grad_grad;
        integral_fct_fct[j][k] += value;
        integral_fct_gradx[j][k] += value_fct_gradx;
        integral_fct_grady[j][k] += value_fct_grady;
      }
    }

//    for (j=0;j<N_;j++)
//     {
//       for (k=0;k<N_;k++)
//  OutPut(j << " " << k << " " << integral_grad_grad[j][k] << " :: ")
//     }
//     OutPut(endl);

    // fill the matrix with standard Brezzi-Pitkaeranta term
    for ( j = 0; j < N_; j++ ) // test fct.
    {
      index_test = DOF[j];
//      if (index_test >= N_Active)
//  continue;
      start = RowPtr[index_test];
      end = RowPtr[index_test + 1];
      for ( k = 0; k < N_; k++ )
      {
        index_ansatz = DOF[k];
//  if (index_ansatz >=  N_Active)
//    continue;
        for ( l = start; l < end; l++ )
        {
          if ( index_ansatz == ColInd_new[l] )
          {
            EntriesC_new[l] -= delta * integral_grad_grad[j][k];
            break;
          }
        }
      }
    } //  standard Brezzi-Pitkaeranta term done

    // medium far entries
    for ( j = 0; j < N_; j++ ) // this is c in Badia (2012)
    {
      index_j = DOF[j];
//      if (index_j >= N_Active)
//  continue;
      // take the cell that is assigned to index_j
      neigh_no_j = node_to_cell[index_j];
      neigh_j = coll->GetCell ( neigh_no_j );
      DOF_neigh_j = GlobalNumbers + BeginIndex[neigh_no_j];
      CurrentElement_neigh_j = pressure_space->GetFE2D ( neigh_no_j, neigh_j );
      N_neigh_j = N_BaseFunct[CurrentElement_neigh_j];
      // loop over the neighbor dofs
      for ( jj = 0; jj < N_neigh_j; jj++ ) // this is a in Badia (2012)
      {
        index_jj = DOF_neigh_j[jj];
//        if (index_jj >= N_Active)
//    continue;
        // memset(tmp_pres,0.0,N_P*SizeOfDouble);
        // set the basis fct. in jj
        tmp_pres[index_jj] = 1.0;
        tmp_pressure.FindGradientLocal ( neigh_j, neigh_no_j, X_dof[j], Y_dof[j], val );
        tmp_pres[index_jj] = 0.0;

        // the indices of this cell might be coupled to index_k
        for ( k = 0; k < N_; k++ ) // this is d in Badia (2012)
        {
          index_k = DOF[k];
//         if (index_k >= N_Active)
//      continue;
          value = delta * ( val[1] * integral_fct_gradx[j][k] + val[2] * integral_fct_grady[j][k] );

          start = RowPtr[index_k];
          end = RowPtr[index_k + 1];
          // check whether the index is already in the list
          for ( l = start; l < end; l++ )
          {
            if ( index_jj == ColInd_new[l] )
            {
              EntriesC_new[l] += value;
              break;
            }
          }
          // now do the transposed entry
          start = RowPtr[index_jj];
          end = RowPtr[index_jj + 1];
          // check whether the index is already in the list
          for ( l = start; l < end; l++ )
          {
            if ( index_k == ColInd_new[l] )
            {
              EntriesC_new[l] += value;
              break;
            }
          }
        }
      } // end jj
    } // end j

    // far entries
    for ( j = 0; j < N_; j++ ) // this is c in Badia (2012)
    {
      index_j = DOF[j];
//    if (index_j >= N_Active)
//      continue;
      // take the cell that is assigned to index_j
      neigh_no_j = node_to_cell[index_j];
      neigh_j = coll->GetCell ( neigh_no_j );
      DOF_neigh_j = GlobalNumbers + BeginIndex[neigh_no_j];
      CurrentElement_neigh_j = pressure_space->GetFE2D ( neigh_no_j, neigh_j );
      N_neigh_j = N_BaseFunct[CurrentElement_neigh_j];
      // loop over the neighbor dofs
      for ( jj = 0; jj < N_neigh_j; jj++ ) // this is a in Badia (2012)
      {
        index_jj = DOF_neigh_j[jj];
//        if (index_jj >= N_Active)
//    continue;
        //memset(tmp_pres,0.0,N_P*SizeOfDouble);
        // set the basis fct. in jj
        tmp_pres[index_jj] = 1.0;
        tmp_pressure.FindGradientLocal ( neigh_j, neigh_no_j, X_dof[j], Y_dof[j], val );
        tmp_pres[index_jj] = 0.0;

        // the indices of this cell might be coupled to index_k
        for ( k = 0; k < N_; k++ ) // this is d in Badia (2012)
        {
          // take the cell that is assigned to index_k
          index_k = DOF[k];
//          if (index_k >= N_Active)
//          continue;
          neigh_no_k = node_to_cell[index_k];
          neigh_k = coll->GetCell ( neigh_no_k );
          DOF_neigh_k = GlobalNumbers + BeginIndex[neigh_no_k];
          CurrentElement_neigh_k = pressure_space->GetFE2D ( neigh_no_k, neigh_k );
          N_neigh_k = N_BaseFunct[CurrentElement_neigh_k];
          // loop over the neighbor dofs
          for ( kk = 0; kk < N_neigh_k; kk++ ) // the b in Badia (2012)
          {
            index_kk = DOF_neigh_k[kk];
//            if (index_kk >= N_Active)
//              continue;
            //memset(tmp_pres,0.0,N_P*SizeOfDouble);
            //tmp_pres[index_set_in_tmp_pres] = 0.0;
            // set the basis fct. in kk
            tmp_pres[index_kk] = 1.0;
            tmp_pressure.FindGradientLocal ( neigh_k, neigh_no_k, X_dof[k], Y_dof[k], val + 3 );
            tmp_pres[index_kk] = 0.0;
            // the value to add
            value = delta * ( val[1] * val[4] + val[2] * val[5] ) * integral_fct_fct[j][k];
            //if ((index_j==0) && (index_k==1))
            //{OutPut(neigh_no_j<< " " << X_dof[j] << " "  << Y_dof[j]<<  " " << neigh_no_k << " " <<  X_dof[k] << " "  << Y_dof[k] << endl);
            //OutPut(val[1] << " " << val[2] <<  " " << val[4] <<  " " << val[5] <<  " " << integral_fct_fct[j][k] << " " << value << endl);
            //}
            start = RowPtr[index_jj];
            end = RowPtr[index_jj + 1];
            // check whether the index is already in the list
            for ( l = start; l < end; l++ )
            {
              if ( index_kk == ColInd_new[l] )
              {
                EntriesC_new[l] -= value;
                break;
              }
            }
          }
        }
      } // end jj
    } // end j
  } // end i

  // correct Diriclet boundary conditions
  if ( ( velocity ) && ( 0 ) )
  {
    OutPut ( "dof " << N_P << " active " << N_Active << endl );
    int *row_ptr = C_new->GetRowPtr();
    /*for (i=0;i<N_Active;i++)
    {
       start = row_ptr[i];
       end = row_ptr[i+1];
       for(j=start;j<end;j++)
       {
    if (ColInd_new[j] >= N_Active)
    {
     EntriesC_new[j] = 0.0;
    }
       }

    }*/
    for ( i = N_Active; i < N_P; i++ )
    {
      start = row_ptr[i];
      end = row_ptr[i + 1];
      for ( j = start; j < end; j++ )
        EntriesC_new[j] = 0.0;
    }
  }
//   OutPut ( "modified matrix " << N_Entries_new << endl )
//   int *row_ptr = C_new->GetRowPtr();
//   for (i=0;i<C->GetN_Rows();i++)
//   //for (i=0;i<10;i++)
//     {
//        start = row_ptr[i];
//        end = row_ptr[i+1];
//        for(j=start;j<end;j++)
//          Output::print(i, " ", ColInd_new[j], " ", EntriesC_new[j], " :: ");
//     }
//   OutPut(endl);

  return C_new;
}
