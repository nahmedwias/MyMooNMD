#include <PrePost_Cylinder_Square.h>
#include <MooNMD_Io.h>
#include <FEDatabase3D.h>
#include <BoundFace.h>
#include <Database.h>

#include <FEFunction3D.h>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

void print_(const std::string name, std::vector<double> vec)
{
  int rank=0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if(rank == 0)
  {
    for(size_t i=0; i<vec.size(); i++)
      Output::print(name, ": ", i, "   ", std::scientific, vec.at(i));  
  }
  Output::print("------------------------\n");
}

void sort_front(size_t n_pres_nodes, int &count1, double* &press_cyl_x, double* &press_cyl_y,
                double* press_cyl_x_tmp, double* press_cyl_y_tmp)
{
  for(size_t i=0;i<n_pres_nodes;i++)
  {
    if (fabs(press_cyl_x_tmp[i]-0.45) < 1e-5)
    {
      press_cyl_x[count1] = press_cyl_x_tmp[i];
      press_cyl_y[count1] = press_cyl_y_tmp[i];
      // correct order
      for(int j=0;j<count1;j++)
      {
        if (press_cyl_y[count1] < press_cyl_y[j])
        {
          for(int k=count1-1;k>=j;k--)
            press_cyl_y[k+1] = press_cyl_y[k];
          press_cyl_y[j] = press_cyl_y_tmp[i];
          break;
        }
      }
      count1++;
    }
  }
}

void sort_back(size_t n_pres_nodes, int &count1, double* &press_cyl_x, double* &press_cyl_y,
               double *press_cyl_x_tmp, double *press_cyl_y_tmp)
{
  int count = count1;
  for(size_t i=0;i<n_pres_nodes;i++)
  {
    // back side
    if (fabs(press_cyl_x_tmp[i]-0.55) < 1e-5)
    {
      press_cyl_x[count1] = press_cyl_x_tmp[i];
      press_cyl_y[count1] = press_cyl_y_tmp[i];
      // correct order
      for(int j=count;j<count1;j++)
      {
        if (press_cyl_y[count1] > press_cyl_y[j])
        {
          for(int k=count1-1;k>=j;k--)
            press_cyl_y[k+1] = press_cyl_y[k];
          press_cyl_y[j] = press_cyl_y_tmp[i];
          break;
        }
      }
      count1++;
    }
  }
}

void sort_right(size_t n_pres_nodes, int &count1, double* &press_cyl_x, double* &press_cyl_y,
                double *press_cyl_x_tmp, double *press_cyl_y_tmp)
{
  int count = count1;
  for(size_t i=0;i<n_pres_nodes;i++)
  {
    if (fabs(press_cyl_y_tmp[i]-0.65) < 1e-5)
    {
      press_cyl_x[count1] = press_cyl_x_tmp[i];
      press_cyl_y[count1] = press_cyl_y_tmp[i];
      // correct order
      for(int j=count;j<count1;j++)
      {
        if (press_cyl_x[count1] > press_cyl_x[j])
        {
          for(int k=count1-1;k>=j;k--)
            press_cyl_x[k+1] = press_cyl_x[k];
          press_cyl_x[j] = press_cyl_x_tmp[i];
          break;
        }
      }
      count1++;
    }
  }
}

void sort_left(size_t n_pres_nodes, int &count1, double* &press_cyl_x, double* &press_cyl_y,
               double *press_cyl_x_tmp, double *press_cyl_y_tmp)
{
  int count = count1;
  for(size_t i=0;i<n_pres_nodes;i++)
  {
    if (fabs(press_cyl_y_tmp[i]-0.75) < 1e-5)
    {
      press_cyl_x[count1] = press_cyl_x_tmp[i];
      press_cyl_y[count1] = press_cyl_y_tmp[i];
      // correct order
      for(int j=count;j<count1;j++)
      {
        if (press_cyl_x[count1] < press_cyl_x[j])
        {
          for(int k=count1-1;k>=j;k--)
            press_cyl_x[k+1] = press_cyl_x[k];
          press_cyl_x[j] = press_cyl_x_tmp[i];
          break;
        }
      }
      count1++;
    }
  }
}

void sort_velo(size_t n_center_velo, double* center_velo, double* center_velo_tmp)
{
  int count=0;
  for(size_t i=0;i<n_center_velo;i++)
  {
    center_velo[count] = center_velo_tmp[i];
    // correct order
    for(int j=0;j<count;j++)
    {
      if (center_velo[count] < center_velo[j])
      {
        for(int k=count-1;k>=j;k--)
          center_velo[k+1] = center_velo[k];
        center_velo[j] = center_velo_tmp[i];
        break;
      }
    }
    count++;
  }
}
//===========================================================================
void velocit_at_cylinder(int &count, double* x )
{
  //   const int max_entries = 200;
  //   std::vector<double> center_velo_tmp(max_entries, -1);
  //   int found, found1;
  //   double xc;
}
//============================================================================

void Cylinder_Square::setParameters(ParameterDatabase& db_)
{
  double renr = db_["reynolds_number"];
  dimensionless_viscosity = 1./renr;
}

void Cylinder_Square::get_Drag_Lift(TFEFunction3D *u1, TFEFunction3D *u2,
                                    TFEFunction3D *u3, TFEFunction3D *p,
                                    TFEFunction3D *u1old, TFEFunction3D *u2old,
                                    double &cdrag, double &clift, const double dt)
{
  const double *pval = p->GetValues();
  double* valu1 = u1->GetValues();
  double* valu2 = u2->GetValues();
  double* valu3 = u3->GetValues();
  
  double* oldvalu1 = u1old->GetValues();
  double* oldvalu2 = u2old->GetValues();
  
  if ( dt < 1e-8)
  {
    OutPut("time step to small " << endl);
    exit(4711);
  }
  
  const TFESpace3D* vel_space = u1->GetFESpace3D();
  TCollection *coll = vel_space->GetCollection();
  size_t nCells=coll->GetN_Cells();
  TBaseCell* cell;
  TJoint* joint;
  TBoundFace* boundface;
  TBoundComp3D *BoundComp;
  
  int* UGlobalNumbers = vel_space->GetGlobalNumbers();
  int* UBeginIndex = vel_space->GetBeginIndex();
  std::vector<double> v(u1->GetLength(), 0.);
  
  for(size_t i=0; i<nCells; i++)
  {
    cell = coll->GetCell(i);
    // loop over all edges of cells
    for(int j=0; j<cell->GetN_Faces(); j++)
    {
      joint = cell->GetJoint(j);
      if ((joint->GetType() == BoundaryFace)) //(joint->GetType() == IsoBoundface)) // boundary edge 
      {
        boundface = (TBoundFace *)joint;  
        BoundComp = boundface->GetBoundComp();  // get boundary component
        int comp=BoundComp->GetID();              // boundary id 
        if ((comp>=4)&&(comp<=7)) 
        {
          FE3D FEEle = vel_space->GetFE3D(i,cell);   // finite element of cell
          TFE3D* eleCell =  TFEDatabase3D::GetFE3D(FEEle); 
          TFEDesc3D* FEDesc = eleCell->GetFEDesc3D();   // fe descriptor
          int N_DOF_Circ = FEDesc->GetN_JointDOF(); // # local dofs on joints
          int* DOF_Circ = FEDesc->GetJointDOF(j); // local dofs on joint j
          int* DOF = UGlobalNumbers + UBeginIndex[i]; // pointer to global dofs
          for(int k=0;k<N_DOF_Circ;k++)         // set fe on circle to 1 
            v[DOF[DOF_Circ[k]]] = 1.;
        }
      }
    }
  }
  cdrag = 0.; clift = 0.;
  FE3D LocalUsedElements[2];
  const TFESpace3D* pres_space = p->GetFESpace3D();;
  bool SecondDer[2] = { false, false};
  int N_Points;
  const double *weights, *xi, *eta, *zeta;
  std::vector<double>X(MaxN_QuadPoints_3D);
  std::vector<double>Y(MaxN_QuadPoints_3D);
  std::vector<double>Z(MaxN_QuadPoints_3D);
  std::vector<double>AbsDetjk(MaxN_QuadPoints_3D);
  
  
  BaseFunct3D* BaseFuncts = TFEDatabase3D::GetBaseFunct3D_IDFromFE3D();
  int* N_BaseFunct = TFEDatabase3D::GetN_BaseFunctFromFE3D();
  
  int* PGlobalNumbers = pres_space->GetGlobalNumbers();
  int* PBeginIndex = pres_space->GetBeginIndex();
  
  int tmp = MaxN_BaseFunctions3D;
  std::vector<double> FEFVal(tmp);
  std::vector<double> FEFVal1(tmp);
  std::vector<double> FEFVal2(tmp);
  std::vector<double> FEFVal3(tmp);
  std::vector<double> FEFVal4(tmp);
  std::vector<double> FEFVal5(tmp);
  std::vector<double> FEFVal6(tmp);
  
  double **OrigVals;
  
  double* aux = new double [MaxN_QuadPoints_3D*19];
  double *Derivatives[MaxN_QuadPoints_3D];
  
  int n_deriv=4;
  MultiIndex3D NeddedDerivatives[4]={D000, D100, D010, D001};
  for(int j=0;j<MaxN_QuadPoints_3D;j++)
    Derivatives[j] = aux + j*19;
  
  for(size_t i=0; i<nCells; i++)
  {
    cell = coll->GetCell(i);
#ifdef _MPI
    if(cell->IsHaloCell())
    {
      continue;
    }
#endif
    // ####################################################################
    // find local used elements on this cell
    // ####################################################################
    int N_LocalUsedElements = 2;
    LocalUsedElements[0] = vel_space->GetFE3D(i, cell);
    LocalUsedElements[1] = pres_space->GetFE3D(i, cell);
    
    TFEDatabase3D::GetOrig(N_LocalUsedElements, LocalUsedElements, coll, cell, SecondDer,
                           N_Points, xi, eta, zeta, weights, X.data(), Y.data(), Z.data(), AbsDetjk.data());
    FE3D CurrentElement = LocalUsedElements[1];
    BaseFunct3D BaseFunct = BaseFuncts[CurrentElement];
    int N_ = N_BaseFunct[CurrentElement];
    
    int* DOF = PGlobalNumbers + PBeginIndex[i];
    
    for(int j=0;j<N_;j++)
      FEFVal[j] = pval[DOF[j]];
    
    OrigVals = TFEDatabase3D::GetOrigElementValues(BaseFunct, D000);
    
    for(int j=0; j<N_Points; j++)
    {
      double *Orig = OrigVals[j];
      double value = 0.;
      for(int l=0; l<N_; l++)
        value += FEFVal[l] * Orig[l];
      
      Derivatives[j][0] = value;
    }
    // calculate all needed values of u1, u2, u3 
    CurrentElement = LocalUsedElements[0];
    BaseFunct = BaseFuncts[CurrentElement];
    N_ = N_BaseFunct[CurrentElement];
    DOF = UGlobalNumbers + UBeginIndex[i];
    for(int l=0;l<N_;l++)
    {
      FEFVal1[l] = valu1[DOF[l]];
      FEFVal2[l] = valu2[DOF[l]];
      FEFVal3[l] = valu3[DOF[l]];
      FEFVal4[l] = v[DOF[l]];
      FEFVal5[l] = oldvalu1[DOF[l]];
      FEFVal6[l] = oldvalu2[DOF[l]];
    }
    for(int k=0; k<n_deriv; k++)
    {
      OrigVals = TFEDatabase3D:: GetOrigElementValues(BaseFunct,NeddedDerivatives[k]);
      for(int j=0;j<N_Points;j++)
      {
        double *Orig = OrigVals[j];
        double value1 = 0.;
        double value2 = 0.;
        double value3 = 0.;
        double value4 = 0.;
        for(int l=0;l<N_;l++)
        {
          value1 += FEFVal1[l] * Orig[l];
          value2 += FEFVal2[l] * Orig[l];
          value3 += FEFVal3[l] * Orig[l];
          value4 += FEFVal4[l] * Orig[l];
        }
        Derivatives[j][k+1] = value1;
        Derivatives[j][k+5] = value2;
        Derivatives[j][k+9] = value3;
        Derivatives[j][k+13] = value4;
        if(k==0)
        {
          value1 = 0.;
          value2 = 0.;
          for(int l=0;l<N_;l++)
          {
            value1 += FEFVal5[l] * Orig[l];
            value2 += FEFVal6[l] * Orig[l];
          }// endfor l
          Derivatives[j][17] = value1;
          Derivatives[j][18] = value2;
        }//endif 
      }//endfor j
    }//endfor k
    // calculation
    for(int j=0;j<N_Points;j++)
    {
      double *Der = Derivatives[j];

      // nu * (u1_x*v_x+ u1_y*v_y + u1_z*v_z), v= (v,0,0)
      double value1  = (Der[1]-Der[17])*Der[13]/dt 
          + dimensionless_viscosity*(Der[2]*Der[14]+Der[3]*Der[15]+Der[4]*Der[16]);
      // (u1 * u1_x + u2* u1_y + u3* u1_z) * (1,0,0)
      value1 += (Der[1]*Der[2]+Der[5]*Der[3]+Der[9]*Der[4])*Der[13];
      // pressure times divergence of test function (1,0,0)
      value1 -= Der[0]*Der[14];

      double value2  = (Der[5]-Der[18])*Der[13]/dt  + 
          dimensionless_viscosity*(Der[6]*Der[14]+Der[7]*Der[15]+Der[8]*Der[16]);
      value2 += (Der[1]*Der[6]+Der[5]*Der[7]+Der[9]*Der[8])*Der[13];
      value2 -= Der[0]*Der[15];

      cdrag += AbsDetjk[j]*weights[j] * value1;
      clift += AbsDetjk[j]*weights[j] * value2;
    }
  }// for nCells
  
#ifdef _MPI
  //communicate the values of cd and cl and sum them up
  MPI_Comm comm = MPI_COMM_WORLD;

  double sendbuf[2] = {cdrag, clift};
  double recvbuf[2] = {0.0, 0.0};
  MPI_Allreduce(sendbuf, recvbuf, 2, MPI_DOUBLE, MPI_SUM, comm);
  cdrag = recvbuf[0];
  clift = recvbuf[1];
#endif
  
  cdrag *= -50;
  clift *= -50;
  
  delete u1;
  delete u2;
  delete u3;
  delete u1old;
  delete u2old;
}

double Cylinder_Square::get_p_diff(const std::array< double, int(3) >& point_A, 
                                   const std::array< double, int(3) >& point_B, const TFEFunction3D& p)
{
  std::vector<double> dP1(4,0.0);
  std::vector<double> dP2(4,0.0);

  double pdiff = 0;

#ifdef _MPI
  {//this whole scope is only of interest in mpi case
    MPI_Comm comm = MPI_COMM_WORLD;
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);

    //find out two processes which contain the points of interest
    int proc_A = p.GetFESpace3D()->GetCollection()->find_process_of_point(0.45, 0.2, 0.205);
    int proc_B = p.GetFESpace3D()->GetCollection()->find_process_of_point(0.55, 0.2, 0.205);

    //calculate value in point A
    if(my_rank==proc_A)
    {
      if(!p.FindGradient(point_A[0],point_A[1],point_A[2], dP1))
      {
        ErrThrow("Hey! I could not find point_A on the promised process!");
      }
    }

    //calculate value in point B
    if (my_rank==proc_B)
    {
      if(!p.FindGradient(point_B[0],point_B[1],point_B[2], dP2))
      {
        ErrThrow("Hey! I could not find point_B on the promised process!");
      }
    }

    // I'm using collective instead of point-to-point communication,
    // which is slower but right now much easier to implement...
    MPI_Bcast(&dP1.at(0), 4, MPI_DOUBLE, proc_A, comm);
    MPI_Bcast(&dP2.at(0), 4, MPI_DOUBLE, proc_B, comm);
    pdiff = dP1[0] - dP2[0];

  }
#else
  p.FindGradient(point_A[0],point_A[1],point_A[2], dP1);
  p.FindGradient(point_B[0],point_B[1],point_B[2], dP2);
  pdiff = dP1[0] - dP2[0];
#endif
  return pdiff;
}

void Cylinder_Square::compute_drag_lift_pdiff(TimeNavierStokes<3>& tnse3d)
{
#ifdef _MPI
  MPI_Comm comm = MPI_COMM_WORLD;
  int my_rank;
  MPI_Comm_rank(comm, &my_rank);
#else
  int my_rank = 0;
#endif
  //compute drag and lift and print them
  double drag, lift;

  const TFEVectFunct3D& u(tnse3d.get_velocity());
  TFEFunction3D& p(tnse3d.get_pressure()); //I want this const!!!
  
  const TFEVectFunct3D& uold (tnse3d.get_velocity_old());

  TFEFunction3D* u1 = u.GetComponent(0);
  TFEFunction3D* u2 = u.GetComponent(1);
  TFEFunction3D* u3 = u.GetComponent(2);
  
  TFEFunction3D* u1old = uold.GetComponent(0);
  TFEFunction3D* u2old = uold.GetComponent(1);

  const double dt = tnse3d.get_time_stepping_scheme().get_step_length();
  // compute drag lift
  get_Drag_Lift(u1, u2, u3, &p, 
                u1old, u2old, drag, lift, dt);


  std::array<double,3> point_A = {{0.45, 0.2, 0.205}};
  std::array<double,3> point_B = {{0.55, 0.2, 0.205}};

  double pdiff = get_p_diff(point_A, point_B, p);

  if(my_rank == 0)
  {
    double ct = tnse3d.get_time_stepping_scheme().current_time_;
    Output::print<1>(ct, " Drag: ", setprecision(16), drag, " ", lift, " ", pdiff);
  }
}


void Cylinder_Square::PrepareCenterlineVelocities(TCollection* coll)
{
  const int max_entries = 200;
  std::vector<double> center_velo_tmp(max_entries, -1);  
  size_t nCells = coll->GetN_Cells();
  TBaseCell *cell;
  int found, found1;
  double x[8], y[8], z[8];
  int count = 0;
  for(size_t i=0; i< nCells; i++)
  {
    cell = coll->GetCell(i);
    int n_vert = cell->GetN_Vertices();
    found1 = 0;
    double xc = 0.;
    //loop over vertices 
    for(int j=0; j<n_vert; j++)
    {
      cell->GetVertex(j)->GetCoords(x[j],y[j],z[j]);
      if (fabs(y[j]-0.7)<1e-5)
      {
        found = 0;
        found1 = 1;
        xc += x[j];
        for(int k=0;k<=count;k++)
        {
          if((fabs(center_velo_tmp[k]-x[j]) < 1e-5))
          {
            found++;
            break;
          }
        }
        // new entry
        if (!found)
        {
          center_velo_tmp[count] = x[j];
          count++;
          if (count >= max_entries)
            ErrThrow("PrepareCenterlineVelocities: max_entries too small ", max_entries);
        }
      }//endif 
    }//endfor n_vert
    // coordinate in the center
    if (found1)
    {
      found = 0;
      xc/=4;
      for(int k=0;k<=count;k++)
      {
        if ((fabs(center_velo_tmp[k]-xc)<1e-5))
        {
          found++;
          break;
        }
      }
      // new entry
      if (!found)
      {
        center_velo_tmp[count] = xc;
        count++;
        if (count >= max_entries)
          ErrThrow("PrepareCenterlineVelocities: max_entries too small ", max_entries);
      }
    }
  }

  // allocate array and initialize
  n_center_velo = count;
  center_velo.resize(3*count);
  sort_velo(n_center_velo, center_velo.data(), center_velo_tmp.data());
  // print_("a", center_velo);exit(0);
  
  counter_av = 0;
  // memory for friction velocities
  velo_friction.resize(4,0.);
  count_fric_vel = 0;
}

void Cylinder_Square::CenterlineVelocities(TimeNavierStokes<3>& tnse3d)
{
#ifdef _MPI
  int mpi_size, mpi_rank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &mpi_rank);
  MPI_Comm_size(comm, &mpi_size);
#endif


  TCollection * coll = tnse3d.get_velocity_space().GetCollection();
  size_t nCells = coll->GetN_Cells();
  TBaseCell* cell;
  double x[8], y[8], z[8];
  
  TFEFunction3D* u1 = tnse3d.get_velocity().GetComponent(0);
  TFEFunction3D* u2 = tnse3d.get_velocity().GetComponent(1);
  double val[5];
  
  std::vector<double> center_velo_x(2*n_center_velo,0.);
  double *center_velo_y =center_velo.data() + n_center_velo;
  std::vector<int> center_velo_num(n_center_velo, 0);
  
  for(size_t i=0; i<nCells; i++)
  {
    cell = coll->GetCell(i);
    int nv = cell->GetN_Vertices();
    int found=0;
    for(int j=0; j<nv; j++)
    {
      cell->GetVertex(j)->GetCoords(x[j], y[j], z[j]);
      if(fabs(y[j]-0.7) < 1e-5)
      {
        found = 1;
        for(size_t k=0; k<n_center_velo; k++)
        {
          if(fabs(center_velo[k] - x[j]) < 1e-5)
          {
            u1->FindGradientLocal(cell, i, x[j], y[j], z[j], val);
            center_velo_x[k] += val[0];
            u2->FindGradientLocal(cell,i,x[j],y[j],z[j],val+1);
            center_velo_y[k] += val[1];
            center_velo_num[k]++;
            // top and bottom twice, only half of number of mesh cells
            if ((fabs(z[j])<1e-5) || (fabs(z[j]- 0.4)<1e-5))
            {
              center_velo_x[k] += val[0];
              center_velo_y[k] += val[1];
              center_velo_num[k]++;
            }//endif
          }//endif
        }//endfor n_center_velo
      }//endif
    }//endfor nv
    // coordiante in the center
    if(found)
    {
      for(int j=0;j<nv;j++)
      {
        for(int j1=j+1;j1<nv;j1++)
        {
          if ((fabs(y[j]-0.7)<1e-5)&&(fabs(y[j1]-0.7)<1e-5) && (fabs(z[j]-z[j1])<1e-5))
          {
            double xc = (x[j]+x[j1])/2.0;
            for(size_t k=0;k<n_center_velo;k++)
            {
              if (fabs(center_velo[k]-xc)<1e-5)
              {
                u1->FindGradientLocal(cell,i,xc,y[j],z[j],val);
                center_velo_x[k] += val[0];
                u2->FindGradientLocal(cell,i,xc,y[j],z[j],val+1);
                center_velo_y[k] += val[1];
                center_velo_num[k]++;
                // top and bottom twice, only half of number of mesh cells
                if ((fabs(z[j])<1e-5) || (fabs(z[j]- 0.4)<1e-5))
                {
                  center_velo_x[k] += val[0];
                  center_velo_y[k] += val[1];
                  center_velo_num[k]++;
                }
              }
            }
          }
        }
      }
    }
  }
  double ct = tnse3d.get_time_stepping_scheme().current_time_;
#ifdef _SEQ
  for(size_t i=0;i<n_center_velo;i++)
  {
    center_velo_x[i]/=center_velo_num[i];
    center_velo_y[i]/=center_velo_num[i];
    center_velo[n_center_velo+i] = counter_av * center_velo[n_center_velo+i]/(counter_av+1)
                                       + center_velo_x[i]/(counter_av+1);
    center_velo[2*n_center_velo+i] = counter_av * center_velo[2*n_center_velo+i]/(counter_av+1)
                                       + center_velo_y[i]/(counter_av+1);
    Output::print(ct, " c_vel ", center_velo[i], " x ", center_velo_x[i], " ", center_velo[n_center_velo+i],
                  " y ", center_velo_y[i], " ", center_velo[2*n_center_velo+i]);
  }
#elif defined(_MPI)
  // assume here that each process holds a number of correct values of
  // center_velo_x, center_velo_y, and center_velo (integrated)
  // Those will be communicated to root and put out on root only

  int mroot = 0;
  std::vector<int> n_vals_from_proc{};
  std::vector<int> displacements{};
  if(mpi_rank == mroot)
  {
    n_vals_from_proc = std::vector<int>(mpi_size);
    displacements = std::vector<int>(mpi_size);
  }

  // Root gathers how many values to get from every process (per call)
  {
    int n_center_velo_int = n_center_velo; //implicit cast...
    int* sbuf = &n_center_velo_int;
    int scount = 1;
    int* rbuf = &n_vals_from_proc[0];
    int rcount = 1;
    MPI_Gather(sbuf, scount, MPI_INT, //send
               rbuf, rcount, MPI_INT, //receive
               mroot, comm);          //control
  }

  if(mpi_rank == mroot)
  {
    for(int p = 0 ; p < n_vals_from_proc.size() ; ++p)
    {
      if(p>0)
        displacements[p] = displacements[p-1] + n_vals_from_proc[p-1];
      Output::print("From processor ", p, " fetch ", n_vals_from_proc[p], " values.");
      Output::print("Displacements [", p, "] is ", displacements[p], ".");

    }
  }

  std::vector<double> center_velo_x_global{};
  std::vector<double> center_velo_y_global{};
  std::vector<double> center_velo_global{};
  int n_center_velo_global = 0;

  if(mpi_rank == mroot)
  {
    n_center_velo_global = std::accumulate(n_vals_from_proc.begin(),n_vals_from_proc.end(),0);
    center_velo_x_global = std::vector<double>(n_center_velo_global,0.0);
    center_velo_y_global = std::vector<double>(n_center_velo_global,0.0);
    center_velo_global   = std::vector<double>(3*n_center_velo_global,0.0);
  }

  // Root gathers the values which lie scattered on the other processes.
  {
    // gather center_velo_x
    double* sbuf = &center_velo_x[0];
    int scount = n_center_velo;
    double* rbuf = &center_velo_x_global[0];
    int* rcounts = &n_vals_from_proc[0];
    int* displs  = &displacements[0];
    MPI_Gatherv(sbuf, scount, MPI_DOUBLE,          //send
                rbuf, rcounts, displs, MPI_DOUBLE, //receive
                mroot, comm);                      //control


    // gather center_velo_y
    sbuf = &center_velo_y[0];
    rbuf = &center_velo_y_global[0];
    MPI_Gatherv(sbuf, scount, MPI_DOUBLE,          //send
                rbuf, rcounts, displs, MPI_DOUBLE, //receive
                mroot, comm);                      //control

    // gather center_velo (three times as long!)
    sbuf = &center_velo[0];
    rbuf = &center_velo_global[0];
    MPI_Gatherv(sbuf, scount, MPI_DOUBLE,          //send
                rbuf, rcounts, displs, MPI_DOUBLE, //receive
                mroot, comm);                      //control

    sbuf = &center_velo[n_center_velo];
    rbuf = &center_velo_global[n_center_velo_global];
    MPI_Gatherv(sbuf, scount, MPI_DOUBLE,          //send
                rbuf, rcounts, displs, MPI_DOUBLE, //receive
                mroot, comm);                      //control

    sbuf = &center_velo[2*n_center_velo];
    rbuf = &center_velo_global[2*n_center_velo_global];
    MPI_Gatherv(sbuf, scount, MPI_DOUBLE,          //send
                rbuf, rcounts, displs, MPI_DOUBLE, //receive
                mroot, comm);                      //control
  }

  // Root must somehow order (!) and output the values.
  if(mpi_rank == 0)
  {
    for(int cv =0; cv < n_center_velo_global ;++cv)
      //Output::print("Center velo global first line ", cv, " is ", center_velo_global[cv]);
      Output::print(ct, " c_vel ", center_velo[cv],
                    " x ", center_velo_x[cv], " ", center_velo[n_center_velo+cv],
                    " y ", center_velo_y[cv], " ", center_velo[2*n_center_velo+cv]);
  }
#endif
  counter_av++;
}

void Cylinder_Square::PrepareVelocityAtCylinder(TCollection* coll)
{
  const int max_entries = 200;
  std::vector<double> center_velo_tmp(max_entries, -1);  
  size_t nCells = coll->GetN_Cells();
  TBaseCell *cell;
  int found, found1;
  double x[8], y[8], z[8];
  int count = 0;
  for(size_t i=0; i< nCells; i++)
  {
    cell = coll->GetCell(i);
    int n_vert = cell->GetN_Vertices();
    found1 = 0;
    double yc = 0.;
    //loop over vertices 
    for(int j=0; j<n_vert; j++)
    {
      cell->GetVertex(j)->GetCoords(x[j],y[j],z[j]);
      if (fabs(x[j]-0.5)<1e-5)
      {
        found = 0;
        found1 = 1;
        yc += y[j];
        for(int k=0;k<=count;k++)
        {
          if((fabs(center_velo_tmp[k]-y[j]) < 1e-5))
          {
            found++;
            break;
          }
        }
        // new entry
        if (!found)
        {
          center_velo_tmp[count] = y[j];
          count++;
          if (count >= max_entries)
            ErrThrow("PrepareVelocityAtCylinder: max_entries too small ", max_entries);
        }
      }//endif 
    }//endfor n_vert
    // coordinate in the center
    if (found1)
    {
      found = 0;
      yc/=4;
      for(int k=0;k<=count;k++)
      {
        if ((fabs(center_velo_tmp[k]-yc)<1e-5))
        {
          found++;
          break;
        }
      }
      // new entry
      if (!found)
      {
        center_velo_tmp[count] = yc;
        count++;
        if (count >= max_entries)
          ErrThrow("PrepareCenterlineVelocities: max_entries too small ", max_entries);
      }
    }
  }

  // allocate array and initialize
  n_cyl_velo = count;
  cyl_velo.resize(2*count);
  sort_velo(n_cyl_velo, cyl_velo.data(), center_velo_tmp.data());
  
  counter_av_single = 0;
}

void Cylinder_Square::VelocityAtCylinder(TimeNavierStokes<3>& tnse3d)
{
  TCollection * coll = tnse3d.get_velocity_space().GetCollection();
  size_t nCells = coll->GetN_Cells();
  TBaseCell* cell;
  double x[8], y[8], z[8];
  
  TFEFunction3D* u1 = tnse3d.get_velocity().GetComponent(0);
  double val[4];
  
  std::vector<double> cyl_velo_x(n_cyl_velo,0.);
  std::vector<int> cyl_velo_num(n_cyl_velo, 0);
  
  for(size_t i=0; i<nCells; i++)
  {
    cell = coll->GetCell(i);
    int nv = cell->GetN_Vertices();
    int found=0;
    for(int j=0; j<nv; j++)
    {
      cell->GetVertex(j)->GetCoords(x[j], y[j], z[j]);
      if(fabs(x[j]-0.5) < 1e-5)
      {
        found = 1;
        for(size_t k=0; k<n_cyl_velo; k++)
        {
          if(fabs(cyl_velo[k] - y[j]) < 1e-5)
          {
            u1->FindGradientLocal(cell, i, x[j], y[j], z[j], val);
            cyl_velo_x[k] += val[0]; cyl_velo_num[k]++;
            // top and bottom twice, only half of number of mesh cells
            if ((fabs(z[j])<1e-5) || (fabs(z[j]- 0.4)<1e-5))
            {
              cyl_velo_x[k] += val[0]; cyl_velo_num[k]++;
            }//endif 
          }//endif
        }//endfor n_cyl_velo
      }//endif
    }//endfor nv
    // coordiante in the center
    if(found)
    {
      for(int j=0;j<nv;j++)
      {
        for(int j1=j+1;j1<nv;j1++)
        {
          if ((fabs(x[j]-0.5)<1e-5)&&(fabs(x[j1]-0.5)<1e-5) && (fabs(z[j]-z[j1])<1e-5))
          {
            double yc = (y[j]+y[j1])/2.0;
            for(size_t k=0;k<n_cyl_velo;k++)
            {
              if (fabs(cyl_velo[k]-yc)<1e-5)
              {
                u1->FindGradientLocal(cell,i,yc,y[j],z[j],val);
                cyl_velo_x[k] += val[0]; cyl_velo_num[k]++;
                // top and bottom twice, only half of number of mesh cells
                if ((fabs(z[j])<1e-5) || (fabs(z[j]- 0.4)<1e-5))
                {
                  cyl_velo_x[k] += val[0]; cyl_velo_num[k]++;
                }
              }
            }
          }
        }
      }
    }
  }
  double ct = tnse3d.get_time_stepping_scheme().current_time_;
  for(size_t i=0;i<n_cyl_velo;i++)
  {
    cyl_velo_x[i] /= cyl_velo_num[i];
    cyl_velo[n_cyl_velo+i] = counter_av_single * cyl_velo[n_cyl_velo+i]/(counter_av_single+1)
                                       + cyl_velo_x[i]/(counter_av_single+1);
    Output::print(ct, " cyl_v ", cyl_velo[i], " ", cyl_velo_x[i], " ", cyl_velo[n_cyl_velo+i]);
  }
  counter_av_single++;
}


void Cylinder_Square::PreparePressureAtCylinder(const TCollection* coll)
{
  const int max_entries = 257;
  int found;
  double press_cyl_x_tmp[max_entries], press_cyl_y_tmp[max_entries];
  TBaseCell *cell;
  TJoint *joint;
  TBoundFace *boundface;
  TBoundComp3D *BoundComp;
  TVertex *vertex;
  
  size_t N_Cells = coll->GetN_Cells();
  double sx, sy;
  int count = 0;
  for(size_t i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    int N_Faces=cell->GetN_Faces();
    for(int j=0;j<N_Faces;j++)              // loop over all edges of cell
    {
      joint=cell->GetJoint(j);
      if ((joint->GetType() == BoundaryFace))
      {
        boundface = (TBoundFace *)joint;  
        BoundComp = boundface->GetBoundComp();  // get boundary component
        int comp=BoundComp->GetID();              // boundary id 
        if ((comp>=4)&&(comp<=7)) 
        {
          int nv = cell->GetN_Vertices();
          double sx0=0.; double sy0 =0.; double sx1 =0.; double sy1 = 0.0;
          int nodesx = 0; int nodesy = 0;
          // find vertices on the boundary
          // compute barycenter of boundary faces
          for(int k=0;k<nv;k++)
          {
            double x, y, z;
            cell->GetVertex(k)->GetCoords(x,y,z);
            if ((fabs(x-0.45)<1e-5)||(fabs(x-0.55)<1e-5))
            {
              sx0 += x; sy0 += y; nodesx++;
            }
            if ((fabs(y-0.65)<1e-5)||(fabs(y-0.75)<1e-5))
            {
              sx1 += x; sy1 += y; nodesy++;
            }
          }
          if (nodesx==4)
          {
            sx = sx0/4; sy = sy0/4;
          }
          if (nodesy==4)
          {
            sx = sx1/4; sy = sy1/4;
          }
          found = 0;
          for(int k=0;k<=count;k++)
          {
            if ((fabs(press_cyl_x_tmp[k]-sx)<1e-5) && (fabs(press_cyl_y_tmp[k]-sy)<1e-5))
            {
              found++; break;
            }
          }
          // new entry
          if (!found)
          {
            press_cyl_x_tmp[count] = sx;
            press_cyl_y_tmp[count] = sy;
            count++;
            if (count >= max_entries)
              ErrThrow("PreparePressureAtCylinder: max_entries too small ", max_entries );
          }
        }
      }
    }
  }
  // allocate array and initialize
  n_pres_nodes = count;
  press_cyl.resize(3*count);// = new double[3*count];
  double* press_cyl_x = press_cyl.data();
  double* press_cyl_y = press_cyl.data() + count;
  memset(press_cyl_y+count, 0, count*sizeof(double));
  count = 0;
  // sort front 
  sort_front(n_pres_nodes, count, press_cyl_x, press_cyl_y, press_cyl_x_tmp, press_cyl_y_tmp);
  // sort left 
  sort_left(n_pres_nodes, count, press_cyl_x, press_cyl_y, press_cyl_x_tmp, press_cyl_y_tmp);
  // sort back 
  sort_back(n_pres_nodes, count, press_cyl_x, press_cyl_y, press_cyl_x_tmp, press_cyl_y_tmp);
  // sort right
  sort_right(n_pres_nodes, count, press_cyl_x, press_cyl_y, press_cyl_x_tmp, press_cyl_y_tmp);
  //print_("a", press_cyl);
  //exit(0);
  //print_("x", press_cyl);
  counter_av_pres = 0;
}

void Cylinder_Square::PressureAtCylinder(TimeNavierStokes<3>& tnse3d)
{
  TCollection * coll = tnse3d.get_velocity_space().GetCollection();
  size_t nCells = coll->GetN_Cells();
  TBaseCell* cell;
  TJoint* joint;
  TBoundFace* boundface;
  TBoundComp3D* boundcomp;
  double sx=0., sy=0., sz=0.;
  double* press_cyl_x = press_cyl.data();
  double* press_cyl_y = press_cyl.data() + n_pres_nodes;
  double* press_cyl_val = press_cyl_y + n_pres_nodes;
  
  std::vector<double> press_cyl_tmp(n_pres_nodes,0.);
  std::vector<int> count_pres_num(n_pres_nodes, 0);
  double val[4];
  for(size_t i=0; i<nCells; i++)
  {
    cell = coll->GetCell(i);
    int n_faces  = cell->GetN_Faces();
    for(int j=0; j<n_faces; j++)
    {
      joint = cell->GetJoint(j);
      if((joint->GetType() == BoundaryFace))
      {
        boundface = (TBoundFace *)joint;  
        boundcomp = boundface->GetBoundComp();  // get boundary component
        int comp=boundcomp->GetID();              // boundary id 
        if ((comp>=4)&&(comp<=7)) 
        {
          int nv = cell->GetN_Vertices();
          double sx0 = 0.; double sy0 =0.; double sz0=0.; 
          double sx1 = 0.; double sy1 =0.; double sz1=0.;
          int nodesx = 0; int nodesy = 0;
          // find vertices on the boundary, compute barycenter of boundary faces
          for(int k=0;k<nv;k++)
          {
            double x, y, z;
            cell->GetVertex(k)->GetCoords(x,y,z);
            if ((fabs(x-0.45)<1e-5)||(fabs(x-0.55)<1e-5))
            {
              sx0 += x; sy0 += y; sz0 += z; nodesx++;
            }
            if ((fabs(y-0.65)<1e-5)||(fabs(y-0.75)<1e-5))
            {
              sx1 += x; sy1 += y; sz1 += z; nodesy++;
            }
          }
          if (nodesx==4)
          {
            sx = sx0/4; sy = sy0/4; sz = sz0/4;
          }
          if (nodesy==4)
          {
            sx = sx1/4;sy = sy1/4; sz = sz1/4;
          }
          for(size_t k=0;k<=n_pres_nodes;k++)
          {
            if ((fabs(press_cyl_x[k]-sx)<1e-5) && (fabs(press_cyl_y[k]-sy)<1e-5))
            {
              tnse3d.get_pressure().FindGradientLocal(cell, i, sx, sy, sz, val);
              press_cyl_tmp[k] += val[0];
              count_pres_num[k]++;
              break;
            }
          }
        }//endif comp>=4
      }//endif joint type
    }
  }//endfor cells
  /*
   cout<<"no of " <<counter_av_pres<< "  " << count_pres_num.size()<<" " << n_pres_nodes<<endl;
  for(int i=0; i<n_pres_nodes;i++)
    cout<<"i " <<count_pres_num[i]<< " " << press_cyl_tmp[i] <<" " << press_cyl_val[i]<<endl;
   */
  for(size_t i=0; i<n_pres_nodes; i++)
    press_cyl_tmp[i] /= count_pres_num[i];
  double ct = tnse3d.get_time_stepping_scheme().current_time_;

  for(size_t i=0; i<n_pres_nodes; i++)
  {
    press_cyl_val[i] = counter_av_pres *(press_cyl_val[i]) /(counter_av_pres+1)
                               + press_cyl_tmp[i]/ (counter_av_pres+1);
    Output::print<1>(ct, " p_cyl ", press_cyl_tmp[i], " ", press_cyl_val[i]);
  }
  counter_av_pres++;
}


void Cylinder_Square::SetNoPenetrationValues(TimeNavierStokes<3>& tnse3d)
{
//   cout<<"setting No SetNoPenetrationValues"<<endl;
//   TFEFunction3D* u2 = tnse3d.get_velocity().GetComponent(1);
//   double* valu2 = u2->GetValues();
//   // diagonal matrices only A11, A22, A33
//   std::vector<std::shared_ptr<FEMatrix>> blocks 
//   = tnse3d.get_system_matrix_().get_blocks_uniquely();
// 
//   int *rowPtr = blocks.at(5)->GetRowPtr();
//   int *kCol = blocks.at(5)->GetKCol();
//   double *entries = blocks.at(5)->GetEntries();
//   int nrows = blocks.at(5)->GetN_Rows();
//   
//   const double dt = tnse3d.get_time_stepping_scheme().get_step_length();
//   TCollection * coll = tnse3d.get_velocity_space().GetCollection();
//   double hmin, hmax;
//   coll->GetHminHmax(&hmin, &hmax);
//   double comp = TDatabase::ParamDB->PENETRATION_CONSTANT * hmin * dt/8.;
// 
//   for(int i=0; i<nrows; i++)
//   {
//     for(int j=rowPtr[i]; j<rowPtr[i+1]; j++)
//     {
//       if(fabs(entries[j]) > comp)
//       {
//         valu2[kCol[j]] = 0.;
//       }
//     }
//   }
//   
//   TFEFunction3D* u3 = tnse3d.get_velocity().GetComponent(2);
//   double* valu3 = u3->GetValues();  
//   
//   rowPtr = blocks.at(10)->GetRowPtr();
//   kCol = blocks.at(10)->GetKCol();
//   entries = blocks.at(10)->GetEntries();
//   nrows = blocks.at(10)->GetN_Rows();
// 
//   for(int i=0; i<nrows; i++)
//   {
//     for(int j=rowPtr[i]; j<rowPtr[i+1]; j++)
//     {
//       if(fabs(entries[j]) > comp)
//       {
//         valu3[kCol[j]] = 0.;
//       }
//     }
//   }
}



void Cylinder_Square::ComputeFrictionVelocities(const TimeNavierStokes<3>& tnse3d)
{
  TCollection * coll = tnse3d.get_velocity_space().GetCollection();
  size_t nCells = coll->GetN_Cells();
  TBaseCell* cell;
  TJoint *joint;
  TBoundFace *boundface;
  TBoundComp3D *BoundComp;
  double x, y, z, values[8];
  double val[4], no[4];
  val[0]=val[1]=val[2]=val[3]=0.;
  no[0]=no[1]=no[2]=no[3]=0.;
  
  TFEFunction3D* u1 = tnse3d.get_velocity().GetComponent(0);
  TFEFunction3D* u2 = tnse3d.get_velocity().GetComponent(1);
  
  for(size_t i=0; i<nCells; i++)
  {
    cell=coll->GetCell(i);
    int n_faces = cell->GetN_Faces();
    for(int j=0; j<n_faces; j++)
    {
      joint = cell->GetJoint(j);
      if((joint->GetType() == BoundaryFace))
      { 
        boundface = (TBoundFace *)joint;
        BoundComp = boundface->GetBoundComp(); //get boundary component
        int comp=BoundComp->GetID(); // boundary id
        if((comp >=4) && (comp<=7))
        {
          int nv = cell->GetN_Vertices();
          //find vertex on the boundary
          for(int k=0; k<nv; k++)
          {
            cell->GetVertex(k)->GetCoords(x,y,z);
            if((fabs(x-0.45)<1e-5)||(fabs(x-0.55)<1e-5) 
                || (fabs(y-0.65)<1e-5)||(fabs(y-0.75)<1e-5))
            {
              switch(comp)
              {
                case 4:
                  u1->FindGradientLocal(cell, i, x, y, z, values);
                  val[0] -= values[1];
                  no[0] +=1;
                  break;
                case 5:
                  u2->FindGradientLocal(cell, i, x, y, z, values);
                  val[1] -= values[2];
                  no[1] +=1;
                  break;
                case 6:
                  u1->FindGradientLocal(cell, i, x, y, z, values);
                  val[2] -= values[1];
                  no[2] +=1;
                  break;
                case 7:
                  u2->FindGradientLocal(cell, i, x, y, z, values);
                  val[3] -= values[2];
                  no[3] +=1;
                  break;
              }
            }//
          }
        }//endif comp >=4
      }//endif boudaryface
    }
  }
  
  if(count_fric_vel == 0)
  {
    for(size_t i=0; i<velo_friction.size(); i++)
      velo_friction[0] = val[i]/no[i];
  }
  else
  {
    for(size_t i=0; i<velo_friction.size(); i++)
      velo_friction[i] = count_fric_vel *(velo_friction[i]) /(count_fric_vel+1)
      + (val[i]/no[i])/ (count_fric_vel+1);
  }
  
  double ct = tnse3d.get_time_stepping_scheme().current_time_;
  Output::print<1>(ct, " fric velo ", velo_friction[0], " ", velo_friction[1], 
                   " ", velo_friction[2], " ", velo_friction[3]);
  
  Output::print<1>(ct, " y+ ", sqrt(fabs(velo_friction[0])*TDatabase::ParamDB->RE_NR), " ", 
                   sqrt(fabs(velo_friction[1])*TDatabase::ParamDB->RE_NR), " ",
                   sqrt(fabs(velo_friction[2])*TDatabase::ParamDB->RE_NR), " ",
                   sqrt(fabs(velo_friction[3])*TDatabase::ParamDB->RE_NR));

  count_fric_vel += 1;
}
