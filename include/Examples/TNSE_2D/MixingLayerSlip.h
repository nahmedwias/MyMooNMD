#ifndef MIXINGLAYERSLIP_H
#define MIXINGLAYERSLIP_H

double DIMENSIONLESS_VISCOSITY;

#define U_INFTY 1

#include <vector>
#include <Time_NSE2D_Merged.h>
#include <algorithm>

namespace Mixing_layer
{
  std::vector<double> layers;
  size_t nlayers;
  std::vector<double>xDofs, yDofs;
  size_t nDofs;
  
  std::vector< double > temporal_mean_u1; 
  std::vector< double > temporal_mean_u2; 
  std::vector< double > R11; 
  std::vector< double > R22;
  
  void fill_arrays(const Time_NSE2D_Merged& tnse2d)
  {
     const TFESpace2D& space = tnse2d.get_velocity_space();
     nDofs=space.GetN_DegreesOfFreedom();
     
     xDofs.resize(nDofs); yDofs.resize(nDofs); 
    
    for(size_t i=0; i<nDofs; ++i)
      space.GetDOFPosition(i, xDofs.at(i), yDofs.at(i));
    
    
    nlayers =0;
    layers.resize(10000);
    
    for(size_t i=0; i<nDofs; i++)
    {
      if(fabs(xDofs.at(i)) < 1e-6 )
      {
        layers.at(nlayers) = yDofs.at(i);
        nlayers++;
      }
    }
    layers.shrink_to_fit();
    
    std::sort(layers.begin(), layers.begin()+nlayers);
  
    temporal_mean_u1.resize(nlayers, 0.);
    temporal_mean_u2.resize(nlayers, 0.);
    R11.resize(nlayers, 0.);
    R22.resize(nlayers, 0.);
  }
  
  void compute_mean_velocity(const Time_NSE2D_Merged& tnse2d);
  
};

void Mixing_layer::compute_mean_velocity(const Time_NSE2D_Merged& tnse2d)
{
  
  
  const TFEVectFunct2D& U = tnse2d.get_velocity();
  size_t nuDofs = U.GetComponent(0)->GetLength();
  const std::vector<double> u1(U.GetComponent(0)->GetValues(), 
                           U.GetComponent(0)->GetValues()+nuDofs);
  const std::vector<double> u2(U.GetComponent(1)->GetValues(), 
                           U.GetComponent(1)->GetValues()+nuDofs);
 
  std::vector<double> u1_spt_mean(nlayers), u2_spt_mean(nlayers);
  std::vector<int> counts_ndofs_per_layer(nlayers);
  for(size_t i=0; i<u1.size(); ++i)
  {
    for(size_t j=0; j<nlayers; ++j)
    {
      if(fabs(yDofs[i]-layers[j]) < 1e-6)
      {
        u1_spt_mean.at(j) = u1_spt_mean.at(j) + u1.at(i);
        u2_spt_mean.at(j) = u2_spt_mean.at(j) + u2.at(i);
        
        counts_ndofs_per_layer.at(j)++;
        break;
      }
    }
  }
  // compute the average
  for(size_t i=0; i<nlayers; ++i)
  {
    u1_spt_mean.at(i) /= counts_ndofs_per_layer.at(i);
    u2_spt_mean.at(i) /= counts_ndofs_per_layer.at(i);
  }
  
  // temporal mean velocity profile
  double step_length=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double currentTime = TDatabase::TimeDB->CURRENTTIME;
  double tstart = TDatabase::TimeDB->T0;
  
  double factor;
  if(currentTime >= tstart)
    factor= step_length/(currentTime - tstart + step_length);
  else
    factor = step_length/(currentTime+step_length);
  
  for(size_t i=0; i<u1_spt_mean.size(); i++)
  {
    temporal_mean_u1.at(i) = temporal_mean_u1.at(i)
        + factor*(u1_spt_mean.at(i) - temporal_mean_u1.at(i));
    temporal_mean_u2.at(i) = temporal_mean_u2.at(i)
        + factor*(u2_spt_mean.at(i) - temporal_mean_u2.at(i));
  }
  
  // computation of rms velocity profile
  // 1. add the product of velocities 
  std::vector<double> u1u1_spt_mean(nlayers);
  std::vector<double> u2u2_spt_mean(nlayers);
  for(size_t i=0; i<u1.size(); ++i)
  {
    for(size_t j=0; j<nlayers; ++j)
    {
      if(fabs(yDofs[i]-layers[j]) < 1e-6)
      {
        u1u1_spt_mean[j] = u1u1_spt_mean[j] + u1[i]*u1[i];
        u2u2_spt_mean[j] = u2u2_spt_mean[j] + u2[i]*u2[i];
        break;
      }
    }
  }
  
  //2. compute the average
  for(size_t i=0; i<nlayers; ++i)
  {
    u1u1_spt_mean.at(i) /= counts_ndofs_per_layer.at(i);
    u2u2_spt_mean.at(i) /= counts_ndofs_per_layer.at(i);
  }
  //3. compute the Reynold stress 
  for(size_t i=0; i<u1u1_spt_mean.size(); ++i)
  {
    R11.at(i) = R11.at(i) + factor*(u1u1_spt_mean.at(i) - R11.at(i));
    R22.at(i) = R22.at(i) + factor*(u2u2_spt_mean.at(i) - R22.at(i));
  }
  
  std::vector<double> rms_u1(nlayers), rms_u2(nlayers);
  // now compute the rms velocities
  for(size_t i=0; i<nlayers; ++i)
  {
    rms_u1.at(i) = R11[i] - pow(temporal_mean_u1.at(i), 2);
    rms_u2.at(i) = R22[i] - pow(temporal_mean_u2.at(i), 2);
  }
  
  // print out the 
  double rms_u, rms_v;
  double t = TDatabase::TimeDB->CURRENTTIME; 
  for(size_t i=0; i<nlayers; ++i)
  {
    rms_u = sqrt(fabs(rms_u1.at(i)));
    rms_v = sqrt(fabs(rms_u2.at(i)));
    Output::print("t ", std::scientific, t," ", setw(8), layers.at(i), 
                  " mu ", setw(8), temporal_mean_u1.at(i), " mv ", setw(8),
                  temporal_mean_u2.at(9), " rms_u ", setw(8), rms_u, " rms_v ", 
                  setw(8), rms_v);
  }
}

void ExampleFile()
{
  Output::print("Example: MixingLayerSlip.h (correct vorticity thickness,");
  Output::print(" scale old one by 4), U_INFTY ", U_INFTY);
}

// ========================================================================
// initial solution
// ========================================================================
void InitialU1(double x, double y, double *values)
{
  double z, sigma;

  sigma = 1/TDatabase::ParamDB->P8;
  z = 2*y/sigma;
  if (z>=0)
    values[0] = U_INFTY * (1-exp(-2*z))/(1+exp(-2*z));
  else
    values[0] = U_INFTY * (exp(2*z)-1)/(exp(2*z)+1);
  values[0] -= 0.001* U_INFTY *exp(-z*z)*cos(8*Pi*x)*8*y/(sigma*sigma);
  values[0] -= 0.001* U_INFTY *exp(-z*z)*cos(20*Pi*x)*8*y/(sigma*sigma);
}

void InitialU2(double x, double y, double *values)
{
  double z, sigma;

  sigma = 1/TDatabase::ParamDB->P8;
  z = 2*y/sigma;
  values[0] = 0.001*U_INFTY*exp(-z*z)*sin(8*Pi*x)*8*Pi;
  values[0] += 0.001*U_INFTY*exp(-z*z)*sin(20*Pi*x)*20*Pi;
}

void InitialP(double x, double y, double *values)
{
  values[0] = 0;
}


// ========================================================================
// exact solution
// ========================================================================
void ExactU1(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactU2(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

void ExactP(double x, double y, double *values)
{
  values[0] = 0;
  values[1] = 0;
  values[2] = 0;
  values[3] = 0;
}

// ========================================================================
// boundary conditions
// ========================================================================
void BoundCondition(int i, double t, BoundCond &cond)
{
  cond = SLIP_FRICTION_PENETRATION_RESISTANCE;
  TDatabase::ParamDB->INTERNAL_SLIP_WITH_FRICTION = 1;
  TDatabase::ParamDB->INTERNAL_PROJECT_PRESSURE=1;
}

void U1BoundValue(int BdComp, double Param, double &value)
{
  if(BdComp>3)
  {
    ErrThrow( "ERROR in file " , __FILE__ , ", line: ",  __LINE__ ,
              ": wrong boundary part number: " , BdComp);
  }
  value =0.;
}

void U2BoundValue(int BdComp, double Param, double &value)
{
  value = 0;
}

// ========================================================================
// coefficients for Stokes form: A, B1, B2, f1, f2
// ========================================================================
void LinCoeffs(int n_points, double *X, double *Y,
               double **parameters, double **coeffs)
{
  for(int i=0;i<n_points;i++)
  {
    coeffs[i][0] = DIMENSIONLESS_VISCOSITY;
    coeffs[i][1] = 0;
    coeffs[i][2] = 0;
  }
}

void ComputeVorticiyThickness(const TFEFunction2D *Vorticity, double &thickness)
{
  int i,j,k,l,found,max_lines,N_Cells,N_Edges,N_Vort;
  TCollection *Coll;
  TBaseCell *cell;
  double x0,x1,x,y0,y1,val[3],max;
  double *global_y_coord, *aver_vort;
  const TFESpace2D *vorticity_space;


  vorticity_space=Vorticity->GetFESpace2D();
  N_Vort = vorticity_space->GetN_DegreesOfFreedom();
  global_y_coord = new double[N_Vort];
  aver_vort =  new double[N_Vort];

  // get pointer to set of mesh cells which define the fe space
  Coll = vorticity_space->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  max_lines = 0;
  // loop over all mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);
    // get number of edges
    N_Edges = cell->GetN_Edges();
    // for all vertices
    for (j=0;j<N_Edges;j++)
    {
      // compute coordinate of vertex
      cell->GetVertex(j)->GetCoords(x0, y0);
      // no vertices on lower and upper boundary
      if ((fabs(y0-1)< 1e-6)||(fabs(y0+1)< 1e-6))
        continue;
      // for all edges which has a larger number
      for (k=j+1; k < N_Edges;k++)
      {
        // compute coordinate of vertex
        cell->GetVertex(k)->GetCoords(x1, y1);
        // vertices do not posses the same y coordinate, continue
        if (fabs(y0-y1)> 1e-6)
          continue;
        // compute integral by edge midpoint rule
        x = (x0+x1)/2;
        // compute value of fe function in midpoint of edge
        Vorticity->FindGradientLocal(cell,i,x,y0,val);
        //OutPut(x << " " << y0 <<  " " << val[0] << endl);
        // add value to entry in array
        found = 0;
        for (l=0;l<max_lines;l++)
        {
          // not the correct entry
          if (fabs(y0-global_y_coord[l])> 1e-6)
            continue;
          // correct entry found
          found = 1;
          aver_vort[l] += fabs(x0-x1)*val[0];
        }
        // new entry has to be created
        if (!found)
        {
          global_y_coord[max_lines] = y0;
          aver_vort[max_lines] = fabs(x0-x1)*val[0];
          max_lines++;
          if (max_lines > N_Vort)
          {
            OutPut("Increase N_Vort in ComputeVorticiyThickness !!!"<<endl);
            exit(4711);
          }
        }
      }
    }
  }

  // compute averages and vorticity thickness
  max = -1;
  for (i=0;i<max_lines;i++)
  {
    // each integral is computed twice
    aver_vort[i]/=2.0;
    // OutPut(i << " " << aver_vort[i] << " " << global_y_coord[i] << endl);
    if (fabs(aver_vort[i]) > max)
      max = fabs(aver_vort[i]);
  }

  thickness =  4*U_INFTY/max; // (2*U_INFTY)/(max/2) (divide max by length of
                               // the integration domain)
  delete global_y_coord;
  delete aver_vort;
  return;
}

void EvaluateSolution(const Time_NSE2D_Merged &tnse2d, double & zero_vort)
{
  const TFEFunction2D& vorticity(tnse2d.get_vorticity_funct());
  double thickness;
  ComputeVorticiyThickness(&vorticity, thickness);

  if(zero_vort < 0)
    zero_vort = thickness;
  double t=TDatabase::TimeDB->CURRENTTIME;
  Output::print( t, " ", "vorticity thickness: ", thickness,
                 " ", thickness/zero_vort);
  
  /// computation of the mean velocity profile
  
}
#endif // MIXINGLAYERSLIP_H
