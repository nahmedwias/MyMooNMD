#include <ChannelFlowRoutines.h>
#include <BaseCell.h>
#include <Collection.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <Enumerations.h>
#include <PeriodicJoint.h>
#include <FEDatabase3D.h>
#include <algorithm>    // std::stable_sort
#include <functional>
#include <numeric>

#ifdef _MPI
#include <ParFECommunicator3D.h>
#endif

void print(const std::string name, std::vector<double> vec)
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

void ChannelTau180::computeAverageVelocity(
        std::array<std::vector<double>, 9> velo, const TimeNavierStokes<3> &tnse3d)
{
#ifdef _MPI
  ErrThrow("MPI communications is not supported for this",
           "function src/Geometry/ChannelFlowRoutines.C");
#endif
  // finite element space
  std::shared_ptr<const TFESpace3D> space = tnse3d.get_velocity_space();
  size_t ndofs = space->GetN_DegreesOfFreedom();
  /// resize the array of vectors
  std::fill(velo.begin(), velo.end(), std::vector<double>(ndofs));  
  // collection
  const TCollection* coll = space->GetCollection();  
  // number of mesh cells
  size_t nCells = coll->GetN_Cells();
  // finite element function from the class itself
  const TFEVectFunct3D& u = tnse3d.get_velocity();
  TFEFunction3D* u1 =u.GetComponent(0);
  TFEFunction3D* u2 =u.GetComponent(1);
  TFEFunction3D* u3 =u.GetComponent(2);
  // values which are the return from the function "FindGradientLocal"
  double v1[4], v2[4], v3[4];
  int *dofs;
  std::vector<double> areas(ndofs);
  int *N_BasFunct=TFEDatabase3D::GetN_BaseFunctFromFE3D();
  for(size_t i=0; i<nCells; ++i)
  {
    TBaseCell *cell = space->GetCollection()->GetCell(i);
    double area = cell->GetMeasure();
    // mapping of the local to global dofs
    dofs =  space->GetGlobalNumbers() + space->GetBeginIndex()[i];
    // computation of barycenter
    double sx=0, sy=0, sz=0;
    int ldof;
    nBasisFunction=N_BasFunct[space->GetFE3D(i,cell)];
    for(size_t j=0; j<nBasisFunction; ++j)
    {
      ldof = dofs[j]; // local dofs
      sx += xDofs[ldof];
      sy += yDofs[ldof];
      sz += zDofs[ldof];
    }
    sx /=nBasisFunction;  sy /=nBasisFunction; sz /=nBasisFunction;
    // compute the derivatives
    double x, y, z;
    for(size_t j=0; j<nBasisFunction; ++j)
    {
      ldof = dofs[j];
      // coordinates
      x=xDofs[ldof]; y=yDofs[ldof]; z=zDofs[ldof];
      // save areas for computing average 
      areas.at(ldof) += area;
      if(reynolds_number==180)
      {
        if ((sx >0)&& (x<-6))
          x = -x;
        if((sy >0)&& (y<-2))
          y = -y;
      }
      //cout<<"bary: " << x << "  " << y << "  "<< z <<endl;
      // compute the local gradients
      u1->FindGradientLocal(cell, i, x, y, z, v1);
      u2->FindGradientLocal(cell, i, x, y, z, v2);
      u3->FindGradientLocal(cell, i, x, y, z, v3);
      // cout << "xyz: " << v1[1]<<"  " << v1[2]<<"  " << v1[0]<<endl;
      // fill the vectors "velo"
      velo[0].at(ldof) += v1[1]*area; /*u1x*/ velo[1].at(ldof) += v1[2]*area; // u1y
      velo[2].at(ldof) += v1[3]*area; // u1z
      velo[3].at(ldof) += v2[1]*area; /*u2x*/ velo[4].at(ldof) += v2[2]*area; // u2y
      velo[5].at(ldof) += v2[3]*area; // u2z
      velo[6].at(ldof) += v3[1]*area; /*u3x*/ velo[7].at(ldof) += v3[2]*area; // u3y
      velo[8].at(ldof) += v3[3]*area; // u3z
    }
  }
  // divide the vector for each component by the area
  for(size_t i=0; i<velo.size(); i++)
  {
    std::transform(velo[i].begin(), velo[i].begin()+ndofs, areas.begin(), 
                   velo[i].begin(), std::divides<double>());
  }
}

void ChannelTau180::setParameters(ParameterDatabase &db)
{
  // example Channel Flow with Reynolds Number 180 and 395
  if(db["example"].is(8))
  {
    TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY=180;
    TDatabase::ParamDB->INTERNAL_MEAN_COMPUTATION=0;
    //
    if (TDatabase::ParamDB->RE_NR==180)
    {
      TDatabase::ParamDB->INTERNAL_BULK_MEAN = 15.6803;
      TDatabase::ParamDB->INTERNAL_BULK_SIMULATION = 15.6803;
    }
    if (TDatabase::ParamDB->RE_NR==395)
    {
      TDatabase::ParamDB->INTERNAL_BULK_MEAN = 17.5452;
      TDatabase::ParamDB->INTERNAL_BULK_SIMULATION = 17.5452;
    }
    if (TDatabase::ParamDB->RE_NR==590)
    {
      TDatabase::ParamDB->INTERNAL_BULK_MEAN = 18.6544005283276;
      TDatabase::ParamDB->INTERNAL_BULK_SIMULATION = 18.6544005283276;
    }
    if (TDatabase::ParamDB->CELL_MEASURE==0)
    {
      TDatabase::ParamDB->CELL_MEASURE = 2;
      OutPut("CELL_MEASURE changed to " <<
      TDatabase::ParamDB->CELL_MEASURE << endl);
    }
    //parameter = 2 for using Gauss3
    //
    TDatabase::ParamDB->INTERNAL_QUAD_RULE=2;

    Output::print("INTERNAL_QUAD_RULE: ",
                  TDatabase::ParamDB->INTERNAL_QUAD_RULE);
  }  
  double renr = db["reynolds_number"];
  reynolds_number = renr;
}

void ChannelTau180::setZCoordinates(TCollection* Coll, int level)
{
  int i, j, k, nCells, nLayers, nVetrices, layer_int, grid_type;
  TBaseCell *cell;
  TVertex *vertex;
  double x,y,z, layer,gamma=TDatabase::ParamDB->CHANNEL_GRID_STRETCH;

  nCells=Coll->GetN_Cells();
  // number of layers on initial grid * 2^level
  nLayers=TDatabase::ParamDB->N_CELL_LAYERS*(int)(pow(2.0,level));
  //TODO Fix me N_CELL_LAYERS
  // nLayers = n_layers*(int)(pow(2.0, level));
  grid_type=TDatabase::ParamDB->GRID_TYPE;
  // cout<<grid_type<<" is grid type "<< endl;

  // first loop over cells
  for(i=0;i<nCells;i++)
  {
    cell=Coll->GetCell(i);
    nVetrices=cell->GetN_Vertices();
    for (j=0;j<nVetrices;j++)
    {
      // read coordinates
      vertex=cell->GetVertex(j);
      vertex->GetCoords(x, y, z);
      // check if z coordinate fits to the prescribed distribution
      // compute inverse mapping
      double t1, t2;
      t1 = (atanh(tanh(gamma)*(z-1))/gamma+1)*nLayers/2.0; // on grid 0
      t2 = acos(1-z)*nLayers/Pi; // on grid 1
      layer = (grid_type==0) ? t1 : t2;

      layer_int=(int)(layer+1e-7);
      // not the correct position
      double temp=layer-layer_int;
      if(fabs(temp) > 1e-5)
      {
        //OutPut(fabs(layer-layer_int) << " ");
        // find closest odd integer
        double a=fabs(layer_int/2.0 - (int)(layer_int/2.0+1e-7));
        k = (a > 0.1) ? layer_int : layer_int+1;
        // cout<<" k " << k << " ";
        // compute new z coordinate
        t1=1 + tanh(gamma*(2.0*k/nLayers -1))/tanh(gamma);
        t2=1-cos(Pi*k/nLayers);
        z = (grid_type==0) ? t1 : t2;
        // set new z coordinate
        vertex->SetCoords(x,y,z);
      }
    }
  }
}
void ChannelTau180::checkZCoordinates(TCollection* Coll, int level)
{
  int i, j, k, nCells, nLayers, nVertices, found, grid_type;
  TBaseCell *cell;
  TVertex *vertex;
  double x,y,z, *zcoor, gamma=TDatabase::ParamDB->CHANNEL_GRID_STRETCH;

  nCells=Coll->GetN_Cells();
  grid_type=TDatabase::ParamDB->GRID_TYPE;
  // number of layers on initial grid * 2^level
  nLayers=TDatabase::ParamDB->N_CELL_LAYERS*(int)(pow(2.0,level));
  Output::print("number of layers on grid level: ", level, "  " , nLayers, "  " ,TDatabase::ParamDB->N_CELL_LAYERS );

  zcoor=new double[nLayers+1];
  for (i=0;i<=nLayers;i++)
  {
    double t1, t2;
    t1 = 1 + tanh(gamma*(2.0*i/nLayers -1))/tanh(gamma);
    t2 = 1-cos(Pi*i/nLayers);
    zcoor[i] = (grid_type==0) ? t1 : t2;
  }
  // first loop over cells
  for(i=0;i<nCells;i++)
  {
    cell=Coll->GetCell(i);
        nVertices=cell->GetN_Vertices();
    for (j=0;j<nVertices;j++)
    {
      // read coordinates
      vertex=cell->GetVertex(j);
      vertex->GetCoords(x, y, z);
      // check if z coordinate fits to the prescribed 
      // distribution compute inverse mapping
      found=0;
      for (k=0;k<=nLayers;k++)
      {
        if (fabs(zcoor[k] - z)<1e-6)
        {
          found++;
          break;
        }
      }
      if (found)
        continue;
      OutPut("coordinate " << z << " not in list !!!");
      exit(4711);
    }
  }
}

void ChannelTau180::setRefineDesc(TCollection* coll)
{
  int nCells=coll->GetN_Cells();
  // reset clipboards to -1
  TBaseCell *cell, *neigh;
  int nJoints;
  for(int i=0; i<nCells; ++i)
  {
    cell= coll->GetCell(i);
    nJoints=cell->GetN_Joints();
    for(int j=0;j<nJoints;++j)
    {
      neigh=cell->GetJoint(j)->GetNeighbour(cell);
      if(neigh) neigh->SetClipBoard(-1);
    }
    cell->SetClipBoard(-1);
  }
  // set the clip board for each cell collection
  for(int i=0; i<nCells; ++i)
    coll->GetCell(i)->SetClipBoard(i);
  
  // set the refinement descriptor
  for(int i=0; i<nCells; ++i)
  {
    cell= coll->GetCell(i);
    nJoints=cell->GetN_Joints();
    for(int j=0;j<nJoints;++j)
    {
      neigh=cell->GetJoint(j)->GetNeighbour(cell);
      if( neigh && neigh->GetClipBoard()==-1)
      {
        (nJoints==6) ? 
            neigh->SetRefDesc(TDatabase::RefDescDB[Hexahedron]) :
            neigh->SetRefDesc(TDatabase::RefDescDB[Tetrahedron]);  
      }
    }
  }
  Output::print("Refinement Descriptor is set for Channel flow ");
}

void ChannelTau180::setPeriodicFaceJoints(TCollection* Coll)
{
  int i, j, N_Cells, N_Faces, l1, jj, ii, N_Faces1;
  double Y[MaxN_QuadPoints_3D];
  double Z[MaxN_QuadPoints_3D];
  double X1[MaxN_QuadPoints_3D], Y1[MaxN_QuadPoints_3D];
  double Z1[MaxN_QuadPoints_3D];
  std::vector<double> X(MaxN_QuadPoints_3D);
  double x,y,z,x1,y1,z1;
  TBaseCell *cell, *cell1;
  TJoint *joint, *joint1;
  const int *TmpFV, *TmpLen;
  int MaxLen, MaxLen1, x_vert_m1, x_vert_1, y_vert_0, y_vert_2;
  int found, x1_vert_m1, x1_vert_1, y1_vert_0, y1_vert_2;
  double y_bottom, y_top, x_bottom, x_top;

  if (reynolds_number==180)
  {
    y_bottom=-2*Pi/3;
    y_top=2*Pi/3;
    x_bottom=-2*Pi;
    x_top=2*Pi;
  }
  else
  {
    y_bottom=-Pi/2;
    y_top=Pi/2;
    x_bottom=-Pi;
    x_top=Pi;
  }

  N_Cells=Coll->GetN_Cells();

  Output::print("SetPeriodicJoints ");
  // first loop over cells
  for(i=0;i<N_Cells;i++)
  {
    cell=Coll->GetCell(i);
    N_Faces=cell->GetN_Faces();
    // check if face on boundary
    for (j=0;j<N_Faces;j++)
    {
      joint=cell->GetJoint(j);
      //not on boundary
      if ((!(joint->GetType() == BoundaryFace)) &&
        (!(joint->GetType() == IsoBoundFace)))
        continue;
      // find vertices on the face
      cell->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen);
      // Output::print("vertices on face ", TmpLen[j]);
      y_vert_0= y_vert_2=x_vert_m1=x_vert_1=0;
      // compute coordinates of vertices
      for (l1=0;l1<TmpLen[j];l1++)
      {
        cell->GetVertex(TmpFV[j*MaxLen+l1])->GetCoords(X[l1], Y[l1], Z[l1]);
        // check if vertex on y=y_bottom
        if (fabs(Y[l1]-y_bottom)<1e-5)
          y_vert_0++;
        // check if vertex on y=y_top
        if (fabs(Y[l1]-y_top)<1e-5)
          y_vert_2++;
        // check if vertex on x=x_bottom
        if (fabs(X[l1]-x_bottom)<1e-5)
          x_vert_m1++;
        // check if vertex on x=x_top
        if (fabs(X[l1]-x_top)<1e-5)
          x_vert_1++;
      }
      found=0;
      if (y_vert_0==TmpLen[j])
        found++;
      if (y_vert_2==TmpLen[j])
        found++;
      if (x_vert_m1==TmpLen[j])
        found++;
      if (x_vert_1==TmpLen[j])
        found++;
      // face not at the interesting boundaries
      if (!found)
        continue;
      // compute barycentre
      x=y=z=0;
      for (l1=0;l1<TmpLen[j];l1++)
      {
        x+=X[l1];
        y+=Y[l1];
        z+=Z[l1];
      }
      x /= TmpLen[j];
      y /= TmpLen[j];
      z /= TmpLen[j];
      //Output::print("bary ", x, " ", y, " ", z);
     /*
      for (int k=0;k<TmpLen[j];k++)
      {
        Output::print("face ", X[k], " ", Y[k], " ", Z[k]);
      }
      */
      // inner loop over the cells
      for(ii=i+1;ii<N_Cells;ii++)
      {
        cell1=Coll->GetCell(ii);
        N_Faces1=cell1->GetN_Faces();
        // check if face on boundary
        for (jj=0;jj<N_Faces1;jj++)
        {
          joint1=cell1->GetJoint(jj);
          //not on boundary
          if (!(joint1->GetType() == BoundaryFace))
            continue;
          // find vertices on the face
          cell1->GetShapeDesc()->GetFaceVertex(TmpFV, TmpLen, MaxLen1);

          x1=y1=z1=0;
          y1_vert_0= y1_vert_2=x1_vert_m1=x1_vert_1=0;
          // compute coordinates of vertices
          // count number of vertices which fulfill the criteria
          for (l1=0;l1<TmpLen[jj];l1++)
          {
            cell1->GetVertex(TmpFV[jj*MaxLen1+l1])->
              GetCoords(X1[l1], Y1[l1], Z1[l1]);
            // check if vertex on y=y_bottom
            if (fabs(Y1[l1]-y_bottom)<1e-5)
              y1_vert_0++;
            // check if vertex on y=y_top
            if (fabs(Y1[l1]-y_top)<1e-5)
              y1_vert_2++;
            // check if vertex on x=x_bottom
            if (fabs(X1[l1]-x_bottom)<1e-5)
              x1_vert_m1++;
            // check if vertex on x=x_top
            if (fabs(X1[l1]-x_top)<1e-5)
              x1_vert_1++;
            x1+=X1[l1];
            y1+=Y1[l1];
            z1+=Z1[l1];
          }

          // bary center
          x1 /= TmpLen[jj];
          y1 /= TmpLen[jj];
          z1 /= TmpLen[jj];
          found=0;
          // Output::print("baryopp ", x1, " ", y1, " ", z1);
          // all vertices of original face are on y_bottom
          // and all vertices of second face are on y_top
          if ((y_vert_0==TmpLen[j])&&(y1_vert_2==TmpLen[jj]))
          {
            // check if the x,z - coordinates of the barycenters are the same
            if((fabs(x-x1)>1e-5)||(fabs(z-z1)>1e-5))
              continue;
            // the face match
            found++;
          }
          // do the same vice versa
          if ((y_vert_2==TmpLen[j])&&(y1_vert_0==TmpLen[jj]))
          {
            if  ((fabs(x-x1)>1e-5)||(fabs(z-z1)>1e-5))
              continue;
            found++;
          }
          // all vertices of original face are on x_bottom
          // and all vertices of second face are on x_top
          if ((x_vert_m1==TmpLen[j])&&(x1_vert_1==TmpLen[jj]))
          {
            // check if the y,z - coordinates of the barycenters are the same
            if  ((fabs(z-z1)>1e-5)||(fabs(y-y1)>1e-5))
              continue;
            found++;
          }
          // do the same vice versa
          if ((x_vert_1==TmpLen[j])&&(x1_vert_m1==TmpLen[jj]))
          {
            if  ((fabs(z-z1)>1e-5)||(fabs(y-y1)>1e-5))
              continue;
            found++;
          }
          if (!found)
            continue;
          /*for (int k=0;k<TmpLen[jj];k++)
          {
            Output::print("opp ", X1[k], " ", Y1[k], " ", Z1[k]);
          }*/
          // delete old joints
          delete joint;
          delete joint1;
          // make new joint
          joint=new TPeriodicJoint(cell, cell1);
          // set joint
          cell->SetJoint(j,joint);
          cell1->SetJoint(jj,joint);
          // find opposite vertex to local vertex zero of face
          if (((y_vert_0==TmpLen[j])&&(y1_vert_2==TmpLen[jj]))
            || ((y_vert_2==TmpLen[j])&&(y1_vert_0==TmpLen[jj])))
          {
            found=-1;
            for (l1=0;l1<TmpLen[jj];l1++)
            {
              if ((fabs(X[0]-X1[l1])<1e-5)&& (fabs(Z[0]-Z1[l1])<1e-5))
              {
                found=l1;
                break;
              }
            }
          }
          if (((x_vert_m1==TmpLen[j])&&(x1_vert_1==TmpLen[jj]))
            || ((x_vert_1==TmpLen[j])&&(x1_vert_m1==TmpLen[jj])))
          {
            found=-1;
            for (l1=0;l1<TmpLen[jj];l1++)
            {
              if ((fabs(Z[0]-Z1[l1])<1e-5)&& (fabs(Y[0]-Y1[l1])<1e-5))
              {
                found=l1;
                break;
              }
            }
          }
          //OutPut("opposite to zero vertex " << found << endl);
          joint->SetMapType(found);
        }
      }
    }
  }
}

void ChannelTau180::GetCoordinatesOfDof(const TimeNavierStokes<3>& tnse3d)
{
  int rank = 0;
#ifdef _MPI
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif
  std::shared_ptr<const TFESpace3D> feSpace=tnse3d.get_velocity_space();
  /// memory for the dofs of x,y and z coordinates
  size_t nDofs=feSpace->GetN_DegreesOfFreedom();
  xDofs.resize(nDofs); yDofs.resize(nDofs); zDofs.resize(nDofs);
  
  for(size_t i=0; i<nDofs; ++i)
    feSpace->GetDOFPosition(i, xDofs.at(i), yDofs.at(i), zDofs.at(i));
  /// tried the above options (I && II) but nothing works
  /// restricted to the other version "localDofs()"
  size_t nCells=feSpace->GetCollection()->GetN_Cells();

  // initialize the layer with zero
  zLayers.resize(1000); nZLayers=0;
  std::fill(zLayers.begin(),zLayers.end(),-4711.);
  
    // loop over all cells - some checks on correct elements
  for(size_t i=0; i<nCells; ++i)
  {
    TBaseCell *cell=feSpace->GetCollection()->GetCell(i);
    FE3D cE=feSpace->GetFE3D(i, cell);

    if((cE != C_Q2_3D_H_A && cE != C_Q2_3D_H_M) && (cE != C_Q1_3D_H_A && cE != C_Q1_3D_H_M))
    {
      ErrThrow("coordinates of dofs are not tested for the"
               "fini element ", cE, " yet");
    }
   }// endfor i<N_Cells
  
  // loop over ndofs
  for(size_t i=0; i<nDofs; i++)
  {
#ifdef _MPI
    const int* masters = tnse3d.get_velocity_space()->get_communicator().GetMaster();
    if(masters[i] != rank)
      continue;
#endif
    if(    fabs(xDofs.at(i)) < 1e-6 
        && fabs(yDofs.at(i)) < 1e-6 )
    {
       zLayers.at(nZLayers) = zDofs.at(i);
      ++nZLayers;
    }
  }
#ifdef _MPI
  {
    int nZLayers_glob = 0;
    int* recievebuf = new int[size]; 
    // gather the nzlayers to root
    MPI_Gather(&nZLayers, 1, MPI_INT, // send
               recievebuf, 1, MPI_INT, // recieve
               0, MPI_COMM_WORLD); // control
    
    // count the global number of zlayers on root
    int* displs = new int[size];
    if(rank==0)
    {
      for(int i =0; i< size;i++)
      {
        displs[i] = nZLayers_glob;
        nZLayers_glob += recievebuf[i];
      }
    }
    
    double* sendzlayer = &zLayers.at(0);
    double* recievezlayer = new double[nZLayers_glob];
    // gather the zlayers from all processes to root
    MPI_Gatherv(sendzlayer, nZLayers, MPI_DOUBLE, // send 
                recievezlayer, recievebuf, displs, // recieve
                MPI_DOUBLE, 0, MPI_COMM_WORLD);  // control
    
    if(rank==0)
    {
      std::sort(recievezlayer, recievezlayer+nZLayers_glob);
      zLayers.resize(nZLayers_glob);
      std::copy(recievezlayer,recievezlayer + nZLayers_glob, zLayers.begin());
      /*for(int i =0; i< nZLayers_glob;i++)
      {
        Output::print(zLayers[i]);
      }*/
    }
    
    // give the information back to all processes
    MPI_Bcast(&nZLayers_glob, 1, MPI_INT, 0, MPI_COMM_WORLD);
    nZLayers = nZLayers_glob;
    
    if(rank != 0)
      zLayers.resize(nZLayers);
    
    MPI_Bcast(&zLayers.at(0) , nZLayers, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /*for(int i =0; i< nZLayers;i++)
    {
      Output::print(zLayers.at(i));
    }*/
    delete [] recievebuf; delete [] recievezlayer; delete [] displs;
  }
#else
// sort the zlayers
std::sort(zLayers.begin(), zLayers.begin()+nZLayers);
#endif

  if(rank == 0)
  {
    for(size_t i=0;i<nZLayers;i++)
      Output::print(i, " ", zLayers[i]);
  }
}

void ChannelTau180::set_up_memory()
{
  std::vector<double> temp(nZLayers);
  for(int i=0; i<3; i++)
    MeanVelocity.push_back(temp);
  
  for(int i=0; i<6; i++)
    ReynoldsStress.push_back(temp);
  
  DerivmeanVelo.resize(nZLayers);
}


void ChannelTau180::computeMeanVelocity(const TimeNavierStokes<3>& tnse3d)
{
  const TFEVectFunct3D& U = tnse3d.get_velocity();
  size_t nuDofs = U.GetComponent(0)->GetLength();
  // component of solution vector
  const std::vector<double> u1(U.GetComponent(0)->GetValues(), 
                           U.GetComponent(0)->GetValues()+nuDofs);
  const std::vector<double> u2(U.GetComponent(1)->GetValues(), 
                           U.GetComponent(1)->GetValues()+nuDofs);
  const std::vector<double> u3(U.GetComponent(2)->GetValues(), 
                           U.GetComponent(2)->GetValues()+nuDofs);
  std::deque<std::vector<double>> u;
  u.push_back(u1);  u.push_back(u2);  u.push_back(u3);

  // spatial mean average=>compute the summation of all values per layer
  // and then the average 
  std::deque<std::vector<double>> spatialMean;
  std::shared_ptr<const TFESpace3D> space = tnse3d.get_velocity_space();
  for(auto & it : u)
  {
    spatialMean.push_back(compute_sum_of_velocity(it, space));
  }
  // communicate the meanvelocity and sum up
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for(size_t dim = 0; dim < spatialMean.size();++dim)
  {
    int n_elems = nZLayers;
    double* sbuf = &spatialMean.at(dim).at(0);// send buffer
    double* rbuf = new double[n_elems]; // recieve buffer
    
    MPI_Allreduce(sbuf, rbuf, n_elems,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for(int i =0; i < n_elems; ++i)
    {
      spatialMean.at(dim).at(i) = rbuf[i];
      // if(rank == 0)
      //  Output::print("dim: ", dim, " i: ", i, " ", rbuf[i]);
    }    
    delete[] rbuf;
  }
  //MPI_Finalize();
#endif

  // compute the summation of layers dof for averaging
  std::vector<int> sum_layer_dofs(nZLayers);
  count_dofs_per_layer(sum_layer_dofs, space);  
  // communicate the summation of layers dofs 
#ifdef _MPI
    int n_elems = nZLayers;
    int* sbuf = &sum_layer_dofs.at(0);
    int* rbuf = new int[n_elems];
    
    MPI_Allreduce(sbuf, rbuf, n_elems,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    for(int i =0; i < n_elems; ++i)
    {
      sum_layer_dofs.at(i) = rbuf[i];
      // Output::print("Summation i: ", i, " ", sum_layer_dofs[i]);
    }    
    delete[] rbuf;
    // MPI_Finalize();
#endif

  // compute average
  for(auto &v : spatialMean)
    std::transform(v.begin(), v.begin()+nZLayers, sum_layer_dofs.begin(), 
                   v.begin(), std::divides<double>());

  double T=TDatabase::TimeDB->CURRENTTIME;
  double t0=TDatabase::TimeDB->T0;
  if(T >= t0)
  {
    int i=0;
    for(auto &it : MeanVelocity)
    {
      temporalMean(spatialMean.at(i),it);
      i++;
    }
  }
  else
  {
    temporalMean(spatialMean.at(0), MeanVelocity.front());
  }  
  /// compute the eddy viscosity model
  // std::array<std::vector<double>, 6> eddyviscosity;
  // eddy_viscosity(eddyviscosity, tnse3d);
  // compute the ReynoldsStress
  // first compute ui*uj, i=1,2,3 and the spatial mean velocities
  std::deque<std::vector<double>>uiuj;
  std::vector<double> temp(nZLayers, 0);
  for(size_t i=0; i<6; i++)
    uiuj.push_back(temp);
  for(size_t i=0; i<nuDofs; i++)
  {
#ifdef _MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int* masters = tnse3d.get_velocity_space()->get_communicator().GetMaster();
    if(masters[i] != rank)
      continue;
#endif 
    for(size_t j=0; j<nZLayers; j++)
    {
      if(fabs(zDofs[i]-zLayers[j]) < 1e-6)
      {
        uiuj[0][j]=uiuj[0][j]+u1.at(i)*u1.at(i);
        uiuj[1][j]=uiuj[1][j]+u2.at(i)*u2.at(i);
        uiuj[2][j]=uiuj[2][j]+u3.at(i)*u3.at(i);
        uiuj[3][j]=uiuj[3][j]+u1.at(i)*u2.at(i);
        uiuj[4][j]=uiuj[4][j]+u1.at(i)*u3.at(i);
        uiuj[5][j]=uiuj[5][j]+u2.at(i)*u3.at(i);
        break;
      }
    }
  }
#ifdef _MPI  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  for(size_t dim = 0; dim < 6;++dim)
  {
    int n_elems = nZLayers;
    double *sbuf = &uiuj.at(dim).at(0); // send buffer
    double* rbuf = new double[n_elems]; // recieve buffer
    
    MPI_Allreduce(sbuf, rbuf, n_elems,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    for(int i =0; i < n_elems; ++i)
    {
      uiuj.at(dim).at(i) = rbuf[i];
      // if(rank == 0)
      //  Output::print("dim: ", dim, " i: ", i, " ", rbuf[i]);
    }    
    delete[] rbuf;      
  }
#endif
  //compute average
  for(auto &it : uiuj)
    std::transform(it.begin(), it.begin()+nZLayers, sum_layer_dofs.begin(), 
                       it.begin(), std::divides<double>());
  if(T >= t0)
  {
    int i=0;
    for(auto &it : ReynoldsStress) // R_ij: in volker's paper
    {
      temporalMean(uiuj.at(i),it);
      i++;
    }
  }
  /// compute root mean square (formula 11 Volker's paper)
  std::deque<std::vector<double>>rmsInstnsities;
  if(T>=t0)
    rmsInstnsities=getrms(ReynoldsStress, MeanVelocity);
  // print and save the data
  print_quantity_of_interest(MeanVelocity, rmsInstnsities, ReynoldsStress);
}

void ChannelTau180::temporalMean(std::vector< double > spatial_mean, 
                                std::vector< double >& temporal_mean)
{
  double step_length=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double currentTime = TDatabase::TimeDB->CURRENTTIME;
  double tstart = TDatabase::TimeDB->T0;
  
  double factor;
  if(currentTime >= tstart)
    factor= step_length/(currentTime - tstart + step_length);
  else
    factor = step_length/(currentTime+step_length);
  for(size_t i=0; i<spatial_mean.size(); i++)
  {
    temporal_mean.at(i) = temporal_mean.at(i)
              + factor*(spatial_mean.at(i) - temporal_mean.at(i));
  }
}

void ChannelTau180::count_dofs_per_layer(std::vector<int> &summ, std::shared_ptr<const TFESpace3D> fesp)
{
  for(size_t i=0; i<zDofs.size(); ++i)
  {
#ifdef _MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(fesp->get_communicator().GetMaster()[i] != rank)
      continue;
#endif 
    for(size_t j=0; j<nZLayers; ++j)
    {
      if(fabs(zDofs[i]-zLayers[j]) < 1e-6)
      {
        summ[j]++;
        break;
      }
    }
  }
}

std::vector< double >
  ChannelTau180::compute_sum_of_velocity(std::vector< double > in, std::shared_ptr<const TFESpace3D> fesp)
{
  // summation in the computation of mean velocities  
  std::vector<double>temp(nZLayers);  
  for(size_t i=0; i<in.size(); ++i)
  {
#ifdef _MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(fesp->get_communicator().GetMaster()[i] != rank)
      continue;
#endif 
    for(size_t j=0; j<temp.size(); ++j)
    {
      if(fabs(zDofs[i]-zLayers[j]) < 1e-6)
      {
        temp.at(j) = temp.at(j) + in.at(i);
        break;
      }
    }
  }
  return temp;
}

void ChannelTau180::eddy_viscosity(
  std::array< std::vector< double>, int(6)> eddy,
  const TimeNavierStokes<3>& tnse3d)
{
  std::shared_ptr<const TFESpace3D> space = tnse3d.get_velocity_space(); 
  size_t nuDofs = space->GetN_DegreesOfFreedom();
      /// fill with zeros
  std::fill(eddy.begin(), eddy.end(), std::vector<double>(nuDofs));
  if(!TDatabase::ParamDB->CHANNEL_STATISTICS2_WITH_MODEL)
    return;
  else
  {
    ErrThrow("ChannelFlowRoutines::ChannelTau180:: ",
             "Eddy viscosity model is not implemented yet");
  }
  /// else compute the contribution from eddy viscosity model
  std::array<std::vector<double>, 9> gradu;
  double area;
  size_t nCells=space->GetCollection()->GetN_Cells();
  size_t nCellsPerLayer;
  double hxhy;
  ///TODO: DISCTYPE is no longer used in the Database 
  /**switch(TDatabase::ParamDB->DISCTYPE)
  {
    case GALERKIN:
      break;
    case SMAGORINSKY:
      ChannelTau180::computeAverageVelocity(gradu, tnse3d);
      // Reynolds number only: 395 or 180 are used      
      area = (reynolds_number == 395) ? 4.*Pi*Pi : 32.*Pi*Pi/3;
      nCellsPerLayer = 2*nCells/(nZLayers-1);
      hxhy = area/(2*nCellsPerLayer);
      Output::print("nCells per layer:  ", nCellsPerLayer, " hxhy: ", hxhy);
      break;
    default:
      ErrThrow("mean velocity for DISCTYPE ", 
               TDatabase::ParamDB->DISCTYPE," is not implemented yet");
  }*/
}

double ChannelTau180::getFrictionVelocity(std::vector< double > vec, 
                                          std::vector<double> &meanDeriv)
{
  double temp;
  /// computing the friction velocity
  /// first compute the derivative of mean velocity
  for(size_t i=0; i<nZLayers; ++i)
  {
    if(i==0)
    {
      meanDeriv.push_back(vec.at(i+1)/zLayers.at(i+1));
    }
    else if(i==nZLayers-1)
    {
      double t = (vec.at(i-1) - vec.at(i))/
                 (zLayers.at(i-1) - zLayers.at(i));
      meanDeriv.push_back(t);
    }
    else
    {
      double D = (zLayers.at(i)-zLayers.at(i-1)) * (zLayers.at(i+1)-zLayers.at(i-1))*
                 (zLayers.at(i+1)-zLayers.at(i));

      double b = ( (vec.at(i)-vec.at(i-1))* (zLayers[i+1]*zLayers[i+1]
                    -zLayers[i-1]*zLayers[i-1])
               - ( vec.at(i+1)-vec.at(i-1))*(zLayers[i]*zLayers[i]
                   -zLayers[i-1]*zLayers[i-1]) )/D;

        double a = ( (zLayers[i]-zLayers[i-1])*(vec.at(i+1)-vec.at(i-1))
                  - (zLayers[i+1]-zLayers[i-1])*(vec.at(i)-vec.at(i-1)) )/D;
      meanDeriv.push_back(2*a*zLayers[i]+b);
    }
  }
  temporalMean(meanDeriv, DerivmeanVelo);
  temp = 1./reynolds_number*(DerivmeanVelo.at(0)-
                            DerivmeanVelo.at(nZLayers-1))*0.5;  
  return temp;
}

std::deque<std::vector< double>> 
ChannelTau180::getrms(std::deque<std::vector<double>> reynoldStress, 
                            std::deque<std::vector<double>> meanvelo)
{
  std::deque<std::vector<double>> temp(3);
  for(size_t i=0; i<nZLayers; ++i)
  { 
    double val =  2./3.*(reynoldStress.at(0)[i] - pow(meanvelo.at(0)[i],2));
    val -= 1./3.*(reynoldStress.at(1)[i] - pow(meanvelo.at(1)[i],2));
    val -= 1./3.*(reynoldStress.at(2)[i] - pow(meanvelo.at(2)[i],2));
    temp.at(0).push_back(val);
    /// 
    val = 2./3.*(reynoldStress.at(1)[i] - pow(meanvelo.at(1)[i],2));
    val -=  1./3.*(reynoldStress.at(0)[i] - pow(meanvelo.at(0)[i],2));
    val -= 1./3.*(reynoldStress.at(2)[i] - pow(meanvelo.at(2)[i],2));
    temp.at(1).push_back(val);
    /// 
    val = 2./3.*(reynoldStress.at(2)[i] - pow(meanvelo.at(2)[i],2));
    val -=  1./3.*(reynoldStress.at(0)[i] - pow(meanvelo.at(0)[i],2));
    val -= 1./3.*(reynoldStress.at(1)[i] - pow(meanvelo.at(1)[i],2));
    temp.at(2).push_back(val);
  }
  return temp;
}

void ChannelTau180::print_quantity_of_interest(std::deque< std::vector< double > > meanvelocity, 
                                 std::deque< std::vector< double > > rms,
                                 std::deque< std::vector< double > > R)
{
  int rank=0;
#ifdef _MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  double t = TDatabase::TimeDB->CURRENTTIME; 
  double bulkVelocity=0.;
  for(size_t i=0; i<nZLayers-1; ++i)
    bulkVelocity += 0.5*(meanvelocity.at(0)[i+1] + meanvelocity.at(0)[i])*(zLayers[i+1]-zLayers[i]);
  bulkVelocity /=2.;
  if(rank==0)
    Output::print("bulk velocity : t ", t, " ", setw(8), bulkVelocity, " tau ", 
                (TDatabase::ParamDB->INTERNAL_BULK_MEAN-bulkVelocity)+1);
  // set the bulk velocity
  TDatabase::ParamDB->INTERNAL_BULK_SIMULATION = bulkVelocity;

  if(!(t>= TDatabase::TimeDB->T0))
  {
    return;
  }
  // frictionvelocity and derivative of mean velocity
  std::vector<double> dmu;  
  double u_tau = getFrictionVelocity(meanvelocity.front(), dmu);
  if(rank ==0)
  {
    for(size_t i=0; i<nZLayers; ++i)
    {
      Output::print("t ", std::scientific, t, " ", setw(8), zLayers.at(i), " ",setw(8), 
                    reynolds_number*(1-fabs(1-zLayers[i])), " mu ", setw(8), meanvelocity.at(0)[i],
                    " dmu ", setw(8), dmu[i], 
                    " mv ", setw(8), meanvelocity.at(2)[i], " mw ", setw(8), meanvelocity.at(1)[i]);
    }
  }
  //
  double rmsu, rmsv, rmsw;
  double R12, R13, R23;
  double R12_abs=0., R13_abs=0., R23_abs=0.;
  if(rank ==0)
  {
    for(size_t i=0; i<nZLayers; ++i)
    {
      // formula 11 Volker's paper
      rmsu = sqrt(fabs(rms.at(0)[i]));
      rmsv = sqrt(fabs(rms.at(1)[i]));
      rmsw = sqrt(fabs(rms.at(2)[i]));
      // formula 10 Volker's paper
      R12 = R.at(3)[i] - meanvelocity.at(0)[i] * meanvelocity.at(1)[i]; 
      R12_abs += fabs(R12);
      R13 = R.at(4)[i] - meanvelocity.at(0)[i] * meanvelocity.at(2)[i]; 
      R13_abs += fabs(R13);
      R23 = R.at(5)[i] - meanvelocity.at(1)[i] * meanvelocity.at(2)[i]; 
      R23_abs += fabs(R23);
      
      Output::print("t ", std::scientific, t, " ", setw(8), zLayers.at(i), " ", setw(8),
                    reynolds_number*(1-fabs(1-zLayers[i])), " rms_u ", setw(8), rmsu, 
                    " rms_v ", setw(8), rmsw, " rms_w ", setw(8), rmsv,
                    " R_uv ", setw(8), R13, " R_uw ", setw(8), R12, 
                    " R_vw ", setw(8), R23);
    }
    Output::print("u_tau ", u_tau, " zero statistics ", setw(8), R12_abs/nZLayers, 
                  " ", setw(8), R13_abs/nZLayers, " ", setw(8), R23_abs/nZLayers);
  }
}

