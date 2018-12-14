#include <MeanVelocity.h>
#include <math.h>
#include <algorithm>
#include <Database.h>
#include <MooNMD_Io.h>

#ifdef __2D__

// void MeanVelocity::fill_arrays(const Time_NSE2D& tnse2d)
// {
//   const TFESpace2D& space = tnse2d.get_velocity_space();
//   nDofs=space.GetN_DegreesOfFreedom();
//   
//   xDofs.resize(nDofs); yDofs.resize(nDofs); 
//   
//   for(size_t i=0; i<nDofs; ++i)
//     space.GetDOFPosition(i, xDofs.at(i), yDofs.at(i));
//   
//   n_xlayers =0;
//   n_ylayers = 0;
//   xlayers.resize(nDofs);
//   ylayers.resize(nDofs);
//   
//   for(size_t i=0; i<nDofs; i++)
//   {
//     if(fabs(xDofs.at(i)) < 1e-6 )
//     {
//       ylayers.at(n_ylayers) = yDofs.at(i);
//       n_ylayers++;
//     }
//     if(fabs(yDofs.at(i)) < 1e-6 )
//     {
//       xlayers.at(n_xlayers) = xDofs.at(i);
//       n_xlayers++;
//     }
//   }
//   
//   xlayers.shrink_to_fit();
//   ylayers.shrink_to_fit();
//   
//   std::sort(xlayers.begin(), xlayers.begin()+n_xlayers);
//   std::sort(ylayers.begin(), ylayers.begin()+n_ylayers);
//   
//   temporal_mean_u1.resize(n_ylayers, 0.);
//   temporal_mean_u2.resize(n_xlayers, 0.);
//   R11.resize(n_ylayers, 0.);
//   R22.resize(n_xlayers, 0.);
// }
// void MeanVelocity::compute_mean_velocity(const Time_NSE2D& tnse2d)
// {
//   const TFEVectFunct2D& U = tnse2d.get_velocity();
//   size_t nuDofs = U.GetComponent(0)->GetLength();
//   const std::vector<double> u1(U.GetComponent(0)->GetValues(), 
//                            U.GetComponent(0)->GetValues()+nuDofs);
//   const std::vector<double> u2(U.GetComponent(1)->GetValues(), 
//                            U.GetComponent(1)->GetValues()+nuDofs);
//   
//   std::vector<double> u1_spt_mean(n_ylayers), u2_spt_mean(n_xlayers);
// 
//   std::vector<int> counts_ndofs_per_xlayer(n_xlayers);
//   std::vector<int> counts_ndofs_per_ylayer(n_ylayers);
//   
//   for(size_t i=0; i<u1.size(); ++i)
//   {
//     for(size_t j=0; j<n_ylayers; ++j)
//     {
//       if(fabs(yDofs[i]-ylayers[j]) < 1e-6)
//       {
//         u1_spt_mean.at(j) = u1_spt_mean.at(j) + u1.at(i);
//         
//         counts_ndofs_per_ylayer.at(j)++;
//         break;
//       }
//     }
//   }
//   
//   for(size_t i=0; i<u2.size(); ++i)
//   {
//     for(size_t j=0; j<n_xlayers; ++j)
//     {
//       if(fabs(xDofs[i]-xlayers[j]) < 1e-6)
//       {
//         u2_spt_mean.at(j) = u2_spt_mean.at(j) + u2.at(i);
//         
//         counts_ndofs_per_xlayer.at(j)++;
//         break;
//       }
//     }
//   }
// 
//   // compute the average
//   for(size_t i=0; i<n_ylayers; ++i)
//     u1_spt_mean.at(i) /= counts_ndofs_per_ylayer.at(i);
//   
//   for(size_t i=0; i<n_xlayers; ++i)
//     u2_spt_mean.at(i) /= counts_ndofs_per_xlayer.at(i);
//   
//   double tau=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
//   double ct = TDatabase::TimeDB->CURRENTTIME;
//   double tstart = TDatabase::TimeDB->T0;
//   
//   double factor;
//   if(ct >= tstart)
//     factor= tau/(ct - tstart + tau);
//   else
//     factor = tau/(ct+tau);
//   
//   for(size_t i=0; i<u1_spt_mean.size(); i++)
//   {
//     temporal_mean_u1.at(i) = temporal_mean_u1.at(i)
//         + factor*(u1_spt_mean.at(i) - temporal_mean_u1.at(i));    
//   }
//   for(size_t i=0; i<u2_spt_mean.size(); ++i)
//   {
//     temporal_mean_u2.at(i) = temporal_mean_u2.at(i)
//         + factor*(u2_spt_mean.at(i) - temporal_mean_u2.at(i));
//   }
//   
//   // computation of rms velocity profile
//   // 1. add the product of velocities 
//   std::vector<double> u1u1_spt_mean(n_ylayers,0.);
//   for(size_t i=0; i<u1.size(); ++i)
//   {
//     for(size_t j=0; j<n_ylayers; ++j)
//     {
//       if(fabs(yDofs[i]-ylayers[j]) < 1e-6)
//       {
//         u1u1_spt_mean[j] = u1u1_spt_mean[j] + u1[i]*u1[i];
//         break;
//       }
//     }
//   }
//   std::vector<double> u2u2_spt_mean(n_xlayers,0.);
//   for(size_t i=0; i<u1.size(); ++i)
//   {
//     for(size_t j=0; j<n_xlayers; ++j)
//     {
//       if(fabs(xDofs[i]-xlayers[j]) < 1e-6)
//       {
//         u2u2_spt_mean[j] = u2u2_spt_mean[j] + u2[i]*u2[i];
//         break;
//       }
//     }
//   }
//   //2. compute the average
//   for(size_t i=0; i<n_ylayers; ++i)
//     u1u1_spt_mean.at(i) /= counts_ndofs_per_ylayer.at(i);
//   for(size_t i=0; i<n_xlayers; ++i)
//     u2u2_spt_mean.at(i) /= counts_ndofs_per_xlayer.at(i);
//   //3. compute the Reynold stress 
//   for(size_t i=0; i<u1u1_spt_mean.size(); ++i)
//     R11.at(i) = R11.at(i) + factor*(u1u1_spt_mean.at(i) - R11.at(i));
//   //3. compute the Reynold stress 
//   for(size_t i=0; i<u2u2_spt_mean.size(); ++i)
//     R22.at(i) = R22.at(i) + factor*(u2u2_spt_mean.at(i) - R22.at(i));
//   
//   std::vector<double> rms_u1(n_ylayers), rms_u2(n_xlayers,0.);
//   // now compute the rms velocities
//   for(size_t i=0; i<n_ylayers; ++i)
//     rms_u1.at(i) = R11[i] - pow(temporal_mean_u1.at(i), 2);
//   
//   for(size_t i=0; i<n_xlayers; ++i)
//     rms_u2.at(i) = R22[i] - pow(temporal_mean_u2.at(i), 2);
//   
//   // print out the 
//   double rms_u, rms_v;
//   for(size_t i=0; i<n_ylayers; ++i)
//   {
//     rms_u = sqrt(fabs(rms_u1.at(i)));
//     
//     Output::print("t ", std::scientific, ct," y ", setw(8), ylayers.at(i), 
//                   " mu ", setw(8), temporal_mean_u1.at(i), 
//                   " rms_u ", setw(8), rms_u, 
//                   " R11 ", setw(8), R11.at(i));
//   }
//   
//   for(size_t i=0; i<n_xlayers; ++i)
//   {
//     rms_v = sqrt(fabs(rms_u2.at(i)));
//     Output::print("t ", std::scientific, ct," x ", setw(8), xlayers.at(i), 
//                   " mv ", setw(8), temporal_mean_u2.at(i), 
//                   " rms_v ", setw(8), rms_v, 
//                   " R22 ", setw(8), R22.at(i) );
//   }
// }


void MeanVelocity::compute_mean_velocity_on_points(const TimeNavierStokes<2>& tnse2d,
  const std::vector<double>& vec, const std::vector<TBaseCell *> cells)
{
  const TFEFunction2D *u1 = tnse2d.get_velocity().GetComponent(0);
  const TFEFunction2D *u2 = tnse2d.get_velocity().GetComponent(1);
  
  std::vector<double> meanu(vec.size()/2,0);
  std::vector<double> meanv(vec.size()/2,0);
  double temp[3];
  
  TCollection *coll = tnse2d.get_velocity_space().GetCollection();
  for(size_t i=0; i<cells.size(); ++i)
  {
    TBaseCell *c = cells.at(i);
    double x = vec.at(2*i);
    double y = vec.at(2*i+1);
    u1->FindGradientLocal(c, coll->GetIndex(c), x, y, temp);
    meanu.at(i) = temp[0];
    
    u2->FindGradientLocal(c, coll->GetIndex(c), x, y, temp);
    meanv.at(i) = temp[0];   
  }
  
  // computation of temporal mean 
  double tau=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double ct = TDatabase::TimeDB->CURRENTTIME;
  double tstart = TDatabase::TimeDB->T0;
  
  
  if(temporal_mean_u1.empty())
  {
    temporal_mean_u1.resize(meanu.size(), 0.);
    temporal_mean_u2.resize(meanv.size(), 0.);
    R11.resize(meanu.size(), 0.);
    R22.resize(meanv.size(), 0.);
  }
  
  
  double factor;
  if(ct >= tstart)
    factor= tau/(ct - tstart + tau);
  else
    factor = tau/(ct+tau);
  
  if(ct >= tstart)
  {
    for(size_t i=0; i<meanu.size(); i++)
    {
      temporal_mean_u1.at(i) = temporal_mean_u1.at(i)
        + factor*(meanu.at(i) - temporal_mean_u1.at(i));    
    }
    for(size_t i=0; i<meanv.size(); ++i)
    {
      temporal_mean_u2.at(i) = temporal_mean_u2.at(i)
        + factor*(meanv.at(i) - temporal_mean_u2.at(i));
    }
  }
  else
  {
    for(size_t i=0; i<meanu.size(); i++)
    {
      temporal_mean_u1.at(i) = temporal_mean_u1.at(i)
        + factor*(meanu.at(i) - temporal_mean_u1.at(i));    
    }
  }
  
  std::vector<double> u1u1(meanu.size(),0.);
  
  for(size_t i=0; i<meanu.size(); ++i)
    u1u1[i] += meanu[i]*meanu[i];
  
  std::vector<double> u2u2(meanv.size(),0.);
  for(size_t i=0; i<meanu.size(); ++i)
    u2u2[i] += meanv[i]*meanv[i];
  
  //3. compute the Reynold stress 
  for(size_t i=0; i<u1u1.size(); ++i)
    R11.at(i) = R11.at(i) + factor*(u1u1.at(i) - R11.at(i));
  //3. compute the Reynold stress 
  for(size_t i=0; i<u2u2.size(); ++i)
    R22.at(i) = R22.at(i) + factor*(u2u2.at(i) - R22.at(i));
  
  std::vector<double> rms_u1(R11.size(),0), rms_u2(R22.size(),0.);
  // now compute the rms velocities
  for(size_t i=0; i<R11.size(); ++i)
    rms_u1.at(i) = R11[i] - pow(temporal_mean_u1.at(i), 2);
  
  for(size_t i=0; i<R22.size(); ++i)
    rms_u2.at(i) = R22[i] - pow(temporal_mean_u2.at(i), 2);
  
  if(ct>=tstart)
  {
    double rms_u, rms_v;
    for(size_t i=0; i<meanu.size(); ++i)
    {
      rms_u = sqrt(fabs(rms_u1.at(i)));
      
      Output::print("t ", std::scientific, ct," y ", setw(8), vec.at(2*i+1), 
                    " mu ", setw(8), temporal_mean_u1.at(i), 
                    " rms_u ", setw(8), rms_u, 
                    " R11 ", setw(8), R11.at(i));
    }
    Output::print("");
    for(size_t i=0; i<meanv.size(); ++i)
    {
      rms_v = sqrt(fabs(rms_u2.at(i)));
      Output::print("t ", std::scientific, ct," x ", setw(8), vec.at(2*i), 
                    " mv ", setw(8), temporal_mean_u2.at(i), 
                    " rms_v ", setw(8), rms_v, 
                    " R22 ", setw(8), R22.at(i) );
    }
  }  
}
void MeanVelocity::compute_velocity_on_points(const TimeNavierStokes<2>& tnse2d,
                                              const std::vector<double>& vec_x_y,
                                              const std::vector<TBaseCell *> cells,
                                              std::string basename)
{
  const TFEFunction2D *u1 = tnse2d.get_velocity().GetComponent(0);
  const TFEFunction2D *u2 = tnse2d.get_velocity().GetComponent(1);
  
  std::vector<double> u1_at_xy(vec_x_y.size()/2,0);
  std::vector<double> u2_at_xy(vec_x_y.size()/2,0);
  
  TCollection *coll = tnse2d.get_velocity_space().GetCollection();

  std::ostringstream velfilename;
  double ct = TDatabase::TimeDB->CURRENTTIME;
  //velfilename << "velocities." << ct << ".txt";
  velfilename << basename << ct << ".txt";
  std::ofstream velfile;
  velfile.open(velfilename.str().c_str());
  for(size_t i=0; i<cells.size(); ++i)
  {
    double temp[3];
    TBaseCell *c = cells.at(i);
    double x = vec_x_y.at(2*i);
    double y = vec_x_y.at(2*i+1);
    u1->FindGradientLocal(c, coll->GetIndex(c), x, y, temp);
    u1_at_xy.at(i) = temp[0];  
    u2->FindGradientLocal(c, coll->GetIndex(c), x, y, temp);
    u2_at_xy.at(i) = temp[0];
    /*Output::print("VELOCITIES: t ", std::scientific, ct,
      " x ", setw(8), vec_x_y.at(2*i), 
      " y ", setw(8), vec_x_y.at(2*i+1), 
      " u(x,y) ",setw(8), u_at_xy.at(i),
      " v(x,y) ",setw(8), v_at_xy.at(i));*/
    
    velfile << vec_x_y.at(2*i) << " "
            << vec_x_y.at(2*i+1) << " "
            << u1_at_xy.at(i) << " "
            << u2_at_xy.at(i) << endl;
    
  }
  velfile.close();
}

#endif
