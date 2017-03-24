#include <MeanVelocity.h>
#include <math.h>
#include <algorithm>
#include <Database.h>
#include <MooNMD_Io.h>

#ifdef __2D__

void MeanVelocity::fill_arrays(const Time_NSE2D_Merged& tnse2d)
{
  const TFESpace2D& space = tnse2d.get_velocity_space();
  nDofs=space.GetN_DegreesOfFreedom();
  
  xDofs.resize(nDofs); yDofs.resize(nDofs); 
  
  for(size_t i=0; i<nDofs; ++i)
    space.GetDOFPosition(i, xDofs.at(i), yDofs.at(i));
  
  n_xlayers =0;
  n_ylayers = 0;
  xlayers.resize(nDofs);
  ylayers.resize(nDofs);
  
  for(size_t i=0; i<nDofs; i++)
  {
    if(fabs(xDofs.at(i)) < 1e-6 )
    {
      ylayers.at(n_ylayers) = yDofs.at(i);
      n_ylayers++;
    }
    if(fabs(yDofs.at(i)) < 1e-6 )
    {
      xlayers.at(n_xlayers) = xDofs.at(i);
      n_xlayers++;
    }
  }
  
  xlayers.shrink_to_fit();
  ylayers.shrink_to_fit();
  
  std::sort(xlayers.begin(), xlayers.begin()+n_xlayers);
  std::sort(ylayers.begin(), ylayers.begin()+n_ylayers);
  
  temporal_mean_u1.resize(n_ylayers, 0.);
  temporal_mean_u2.resize(n_xlayers, 0.);
  R11.resize(n_ylayers, 0.);
  R22.resize(n_xlayers, 0.);
}
void MeanVelocity::compute_mean_velocity(const Time_NSE2D_Merged& tnse2d)
{
  const TFEVectFunct2D& U = tnse2d.get_velocity();
  size_t nuDofs = U.GetComponent(0)->GetLength();
  const std::vector<double> u1(U.GetComponent(0)->GetValues(), 
                           U.GetComponent(0)->GetValues()+nuDofs);
  const std::vector<double> u2(U.GetComponent(1)->GetValues(), 
                           U.GetComponent(1)->GetValues()+nuDofs);
  
  std::vector<double> u1_spt_mean(n_ylayers), u2_spt_mean(n_xlayers);

  std::vector<int> counts_ndofs_per_xlayer(n_xlayers);
  std::vector<int> counts_ndofs_per_ylayer(n_ylayers);
  
  for(size_t i=0; i<u1.size(); ++i)
  {
    for(size_t j=0; j<n_ylayers; ++j)
    {
      if(fabs(yDofs[i]-ylayers[j]) < 1e-6)
      {
        u1_spt_mean.at(j) = u1_spt_mean.at(j) + u1.at(i);
        
        counts_ndofs_per_ylayer.at(j)++;
        break;
      }
    }
  }
  
  for(size_t i=0; i<u2.size(); ++i)
  {
    for(size_t j=0; j<n_xlayers; ++j)
    {
      if(fabs(xDofs[i]-xlayers[j]) < 1e-6)
      {
        u2_spt_mean.at(j) = u2_spt_mean.at(j) + u2.at(i);
        
        counts_ndofs_per_xlayer.at(j)++;
        break;
      }
    }
  }

  // compute the average
  for(size_t i=0; i<n_ylayers; ++i)
    u1_spt_mean.at(i) /= counts_ndofs_per_ylayer.at(i);
  
  for(size_t i=0; i<n_xlayers; ++i)
    u2_spt_mean.at(i) /= counts_ndofs_per_xlayer.at(i);
  
  double tau=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double ct = TDatabase::TimeDB->CURRENTTIME;
  double tstart = TDatabase::TimeDB->T0;
  
  double factor;
  if(ct >= tstart)
    factor= tau/(ct - tstart + tau);
  else
    factor = tau/(ct+tau);
  
  for(size_t i=0; i<u1_spt_mean.size(); i++)
  {
    temporal_mean_u1.at(i) = temporal_mean_u1.at(i)
        + factor*(u1_spt_mean.at(i) - temporal_mean_u1.at(i));    
  }
  for(size_t i=0; i<u2_spt_mean.size(); ++i)
  {
    temporal_mean_u2.at(i) = temporal_mean_u2.at(i)
        + factor*(u2_spt_mean.at(i) - temporal_mean_u2.at(i));
  }
  
  // computation of rms velocity profile
  // 1. add the product of velocities 
  std::vector<double> u1u1_spt_mean(n_ylayers,0.);
  for(size_t i=0; i<u1.size(); ++i)
  {
    for(size_t j=0; j<n_ylayers; ++j)
    {
      if(fabs(yDofs[i]-ylayers[j]) < 1e-6)
      {
        u1u1_spt_mean[j] = u1u1_spt_mean[j] + u1[i]*u1[i];
        break;
      }
    }
  }
  std::vector<double> u2u2_spt_mean(n_xlayers,0.);
  for(size_t i=0; i<u1.size(); ++i)
  {
    for(size_t j=0; j<n_xlayers; ++j)
    {
      if(fabs(xDofs[i]-xlayers[j]) < 1e-6)
      {
        u2u2_spt_mean[j] = u2u2_spt_mean[j] + u2[i]*u2[i];
        break;
      }
    }
  }
  //2. compute the average
  for(size_t i=0; i<n_ylayers; ++i)
    u1u1_spt_mean.at(i) /= counts_ndofs_per_ylayer.at(i);
  for(size_t i=0; i<n_xlayers; ++i)
    u2u2_spt_mean.at(i) /= counts_ndofs_per_xlayer.at(i);
  //3. compute the Reynold stress 
  for(size_t i=0; i<u1u1_spt_mean.size(); ++i)
    R11.at(i) = R11.at(i) + factor*(u1u1_spt_mean.at(i) - R11.at(i));
  //3. compute the Reynold stress 
  for(size_t i=0; i<u2u2_spt_mean.size(); ++i)
    R22.at(i) = R22.at(i) + factor*(u2u2_spt_mean.at(i) - R22.at(i));
  
  std::vector<double> rms_u1(n_ylayers), rms_u2(n_xlayers,0.);
  // now compute the rms velocities
  for(size_t i=0; i<n_ylayers; ++i)
    rms_u1.at(i) = R11[i] - pow(temporal_mean_u1.at(i), 2);
  
  for(size_t i=0; i<n_xlayers; ++i)
    rms_u2.at(i) = R22[i] - pow(temporal_mean_u2.at(i), 2);
  
  // print out the 
  double rms_u, rms_v;
  for(size_t i=0; i<n_ylayers; ++i)
  {
    rms_u = sqrt(fabs(rms_u1.at(i)));
    
    Output::print("t ", std::scientific, ct," y ", setw(8), ylayers.at(i), 
                  " mu ", setw(8), temporal_mean_u1.at(i), 
                  " rms_u ", setw(8), rms_u, 
                  " R11 ", setw(8), R11.at(i));
  }
  
  for(size_t i=0; i<n_xlayers; ++i)
  {
    rms_v = sqrt(fabs(rms_u2.at(i)));
    Output::print("t ", std::scientific, ct," x ", setw(8), xlayers.at(i), 
                  " mv ", setw(8), temporal_mean_u2.at(i), 
                  " rms_v ", setw(8), rms_v, 
                  " R22 ", setw(8), R22.at(i) );
  }
}

#endif