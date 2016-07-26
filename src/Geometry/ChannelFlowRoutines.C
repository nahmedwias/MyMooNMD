#include <ChannelFlowRoutines.h>
#include <BaseCell.h>
#include <Collection.h>
#include <TNSE3D_Routines.h>
#include <Database.h>
#include <MooNMD_Io.h>
#include <Enumerations.h>
#include <PeriodicJoint.h>
#include <FEDatabase3D.h>
#include <algorithm>    // std::stable_sort
#include <functional>
#include <numeric>

void print(const std::string name, std::vector<double> vec)
{
  for(size_t i=0; i<vec.size(); i++)
    Output::print(name, ": ", i, "   ", vec.at(i));  
  Output::print("------------------------\n");
}

void localDofs(int i, double* x, double* y, double* z,
               double &locx, double &locy, double &locz)
{
  switch(i)
  {
    case 0:
      locx=x[0]; locy=y[0]; locz=z[0];
      break;
    case 1:
      locx=0.5*(x[0]+x[1]); locy =0.5*(y[0]+y[1]); locz=0.5*(z[0]+z[1]);
      break;
    case 2:
      locx=x[1]; locy=y[1]; locz=z[1];
      break;
    case 3:
      locx=0.5*(x[0]+x[3]); locy=0.5*(y[0]+y[3]); locz=0.5*(z[0]+z[3]);
      break;
    case 4:
      locx=0.25*(x[0]+x[1]+x[2]+x[3]); locy=0.25*(y[0]+y[1]+y[2]+y[3]); locz=0.25*(z[0]+z[1]+z[2]+z[3]);
      break;
    case 5:
      locx=0.5*(x[1]+x[2]); locy=0.5*(y[1]+y[2]); locz=0.5*(z[1]+z[2]);
      break;
    case 6:
      locx=x[3]; locy=y[3]; locz=z[3];
      break;
    case 7:
      locx=0.5*(x[2]+x[3]); locy=0.5*(y[2]+y[3]); locz=0.5*(z[2]+z[3]);
      break;
    case 8:
      locx=x[2]; locy=y[2]; locz=z[2];
      break;
    case 9:
      locx=0.5*(x[0]+x[4]); locy=0.5*(y[0]+y[4]); locz=0.5*(z[0]+z[4]);
      break;
    case 10:
      locx=0.25*(x[0]+x[1]+x[4]+x[5]); locy=0.25*(y[0]+y[1]+y[4]+y[5]); locz=0.25*(z[0]+z[1]+z[4]+z[5]);
      break;
    case 11:
      locx=0.5*(x[1]+x[5]); locy=0.5*(y[1]+y[5]); locz=0.5*(z[1]+z[5]);
      break;
    case 12:
      locx=0.25*(x[0]+x[3]+x[4]+x[7]); locy=0.25*(y[0]+y[3]+y[4]+y[7]); locz=0.25*(z[0]+z[3]+z[4]+z[7]);
      break;
    case 13:
      locx=0.125*(x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]);
      locy=0.125*(y[0]+y[1]+y[2]+y[3]+y[4]+y[5]+y[6]+y[7]);
      locz=0.125*(z[0]+z[1]+z[2]+z[3]+z[4]+z[5]+z[6]+z[7]);
      break;
    case 14:
      locx=0.25*(x[1]+x[2]+x[5]+x[6]); locy=0.25*(y[1]+y[2]+y[5]+y[6]); locz=0.25*(z[1]+z[2]+z[5]+z[6]);
      break;
    case 15:
      locx=0.5*(x[3]+x[7]); locy=0.5*(y[3]+y[7]); locz=0.5*(z[3]+z[7]);
      break;
    case 16:
      locx=0.25*(x[2]+x[3]+x[6]+x[7]); locy=0.25*(y[2]+y[3]+y[6]+y[7]); locz=0.25*(z[2]+z[3]+z[6]+z[7]);
      break;
    case 17:
      locx=0.5*(x[2]+x[6]); locy=0.5*(y[2]+y[6]); locz=0.5*(z[2]+z[6]);
      break;
    case 18:
      locx=x[4]; locy=y[4]; locz=z[4];
      break;
    case 19:
      locx=0.5*(x[4]+x[5]); locy=0.5*(y[4]+y[5]); locz=0.5*(z[4]+z[5]);
      break;
    case 20:
      locx=x[5]; locy=y[5]; locz=z[5];
      break;
    case 21:
      locx=0.5*(x[4]+x[7]); locy=0.5*(y[4]+y[7]); locz=0.5*(z[4]+z[7]);
      break;
    case 22:
      locx=0.25*(x[4]+x[5]+x[6]+x[7]); locy=0.25*(y[4]+y[5]+y[6]+y[7]); locz=0.25*(z[4]+z[5]+z[6]+z[7]);
      break;
    case 23:
      locx=0.5*(x[5]+x[6]); locy=0.5*(y[5]+y[6]); locz=0.5*(z[5]+z[6]);
      break;
    case 24:
      locx=x[7]; locy=y[7]; locz=z[7];
      break;
    case 25:
      locx=0.5*(x[6]+x[7]); locy=0.5*(y[6]+y[7]); locz=0.5*(z[6]+z[7]);
      break;
    case 26:
      locx=x[6]; locy=y[6]; locz=z[6];
      break;
    default:
      ErrThrow("C_Q2_3D_H_A has 26 degrees of freedoms only" );
      break;
  }
}


void ChannelFlowRoutines::computeAverageVelocity(
        std::array<std::vector<double>, 9> velo,
        int ndofs, Time_NSE3D &tnse3d)
{
  /// resize the array of vectors
  std::fill(velo.begin(), velo.end(), std::vector<double>(ndofs));
  // finite element space
  const TFESpace3D& space = tnse3d.get_velocity_space();
  // collection
  TCollection* coll = space.GetCollection();  
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
  double reynoldsNumber = TDatabase::ParamDB->RE_NR;
  std::vector<double> areas(ndofs);
  for(size_t i=0; i<nCells; ++i)
  {
    TBaseCell *cell = space.GetCollection()->GetCell(i);
    double area = cell->GetMeasure();
    // mapping of the local to global dofs
    dofs =  space.GetGlobalNumbers() + space.GetBeginIndex()[i];
    // computation of barycenter
    double sx=0, sy=0, sz=0;
    int ldof;
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
      if(reynoldsNumber==180)
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

  for(size_t i=0; i<velo.size(); i++)
  {
    std::transform(velo[i].begin(), velo[i].begin()+ndofs, areas.begin(), 
                   velo[i].begin(), std::divides<double>());
  }
}

void ChannelFlowRoutines::setParameters(ParameterDatabase &db)
{
  // example Channel Flow with Reynolds Number 180 and 395
  if(db["example"].is(7))
  {
    TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY=180;
    TDatabase::ParamDB->INTERNAL_MEAN_COMPUTATION=0;
    TDatabase::ParamDB->INTERNAL_QUAD_RULE=2;

    Output::print("INTERNAL_QUAD_RULE: ",
                  TDatabase::ParamDB->INTERNAL_QUAD_RULE);
  }  
}

void ChannelFlowRoutines::setZCoordinates(TCollection* Coll, int level)
{
  int i, j, k, nCells, nLayers, nVetrices, layer_int, grid_type;
  TBaseCell *cell;
  TVertex *vertex;
  double x,y,z, layer,gamma=TDatabase::ParamDB->CHANNEL_GRID_STRETCH;

  nCells=Coll->GetN_Cells();
  // number of layers on initial grid * 2^level
  nLayers=TDatabase::ParamDB->N_CELL_LAYERS*(int)(pow(2.0,level));
  grid_type=TDatabase::ParamDB->GRID_TYPE;
  // cout<<grid_type<<" is grid type "<< endl;

  Output::print("number of layers on grid level: ", level, "  " , nLayers );

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
void ChannelFlowRoutines::checkZCoordinates(TCollection* Coll, int level)
{
  int i, j, k, nCells, nLayers, nVertices, found, grid_type;
  TBaseCell *cell;
  TVertex *vertex;
  double x,y,z, *zcoor, gamma=TDatabase::ParamDB->CHANNEL_GRID_STRETCH;

  nCells=Coll->GetN_Cells();
  grid_type=TDatabase::ParamDB->GRID_TYPE;
  // number of layers on initial grid * 2^level
  nLayers=TDatabase::ParamDB->N_CELL_LAYERS*(int)(pow(2.0,level));
  Output::print("number of layers on grid level: ", level, "  " , nLayers );
  
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

void ChannelFlowRoutines::setRefineDesc(TCollection* coll)
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

void ChannelFlowRoutines::setPeriodicFaceJoints(TCollection* Coll)
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
  double RE=TDatabase::ParamDB->RE_NR;
  double y_bottom, y_top, x_bottom, x_top;

  if (RE==180)
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

void ChannelFlowRoutines::GetCoordinatesOfDof(const Time_NSE3D& tnse3d)
{
  const double reynoldsNumber=TDatabase::ParamDB->RE_NR;
  double per_x, per_y;
  if(reynoldsNumber==180)
  {
    per_x=2.*Pi; per_y=2.*Pi/3.;
  }
  else
  {
    per_x=Pi; per_y=Pi/2.;
  }
  size_t Max=1000;  
  const TFESpace3D& feSpace=tnse3d.get_velocity_space();
  /// memory for the dofs of x,y and z coordinates
  size_t nDofs=tnse3d.get_size();
  xDofs.resize(nDofs); yDofs.resize(nDofs); zDofs.resize(nDofs);
  ///TODO: try to fix the function "localDofs();" with the
  /// already implemented function which returns dofs of position
  /// OutCome: some values have different -ve signs for the
  /// periodic joints: may be or may not be fixed????
  /// option: I
  // // feSpace.GetDOFPosition(xDofs.data(), yDofs.data(),zDofs.data());
  // // cout<<xDofs[0]<< "  " << yDofs[0] << "  " << zDofs[0] << endl;
  /// option: II
  /*for(int i=0; i<tnse3d.get_velocity_space().GetN_DegreesOfFreedom();++i)
  {
    feSpace.GetDOFPosition(i, localX, localY, localZ);
    xDofs.at(i)=localX; yDofs.at(i)=localY; zDofs.at(i)=localZ;
  }*/
  /// tried the above options (I && II) but nothing works
  /// restricted to the other version "localDofs()"
  size_t nCells=feSpace.GetCollection()->GetN_Cells();
  int *N_BaseFunct=TFEDatabase3D::GetN_BaseFunctFromFE3D();
  int globalDof;
  double x[8], y[8], z[8], localX, localY, localZ;

  // initialize the layer with zero
  zLayers.resize(Max); nZLayers=0;
  std::fill(zLayers.begin(),zLayers.end(),-4711.);
  // loop over all cells
  for(size_t i=0; i<nCells; ++i)
  {
    TBaseCell *cell=feSpace.GetCollection()->GetCell(i);
    FE3D cE=feSpace.GetFE3D(i, cell);
    nBasisFunction=N_BaseFunct[cE];

    int nVertices=cell->GetN_Vertices();
    for(int j=0; j<nVertices; ++j)
    {
      TVertex* vertex=cell->GetVertex(j);
      vertex->GetCoords(x[j],y[j],z[j]);
    }
    if(cE != C_Q2_3D_H_A && cE != C_Q2_3D_H_M)
    {
      ErrThrow("coordinates of dofs are not implemented for the"
               "fini element ", cE, " yet");
    }
    /// mapping from local to global degrees of freedom
    int *dof=feSpace.GetGlobalNumbers()+feSpace.GetBeginIndex()[i];
    for(size_t j=0; j<nBasisFunction;++j)
    {
      /// compute the local dofs
      localDofs(j, x, y, z, localX, localY, localZ);
      globalDof=dof[j];
      /// peridoic b.c. are set onlocy from one side
      if((fabs(localX-per_x)<1e-6) ||(fabs(localY-per_y)<1e-6))
          continue;
      // Output::print(i ," loc " ,j ," glob " ,globalDof ," x " 
      // ,localX ," y " ,localY ," z " ,localZ );
      // save the (xyz)Dofs for later use
      xDofs[globalDof]=localX;
      yDofs[globalDof]=localY;
      zDofs[globalDof]=localZ;

      for(size_t k=0;k<Max; ++k)
      {
        if(fabs(zLayers[k]-localZ) < 1e-6)
           break;
        if(zLayers[k] == -4711.)
        {
          zLayers[k]=localZ;
          nZLayers++;
          break;
        }
      }
    }// endfor j<nBasisFunction
  }// endfor i<N_Cells
  // sort the z layer
  std::sort(zLayers.begin(), zLayers.begin()+nZLayers);

  for(size_t i=0;i<nZLayers;i++)
    cout<<i << " " << zLayers[i] << endl;

}

void ChannelFlowRoutines::computeMeanVelocity(Time_NSE3D& tnse3d)
{
  const TFEVectFunct3D& u = tnse3d.get_velocity();
  TFEFunction3D* u1 =u.GetComponent(0);
  TFEFunction3D* u2 =u.GetComponent(1);
  TFEFunction3D* u3 =u.GetComponent(2);
  // number of velocity dofs for each component
  size_t nuDofs = u1->GetLength();
  // velocity vectors
  std::vector<double> velocity_u1(u1->GetValues(),
                                  u1->GetValues()+nuDofs);
  std::vector<double> velocity_u2(u2->GetValues(),
                                    u2->GetValues()+nuDofs);
  std::vector<double> velocity_u3(u3->GetValues(),
                                    u3->GetValues()+nuDofs);

  // compute the summation of all values per layer
  std::list<std::vector<double>> velocity;
  velocity.push_back(velocity_u1);
  velocity.push_back(velocity_u2);
  velocity.push_back(velocity_u3);
  
  /// vector for spatial mean velocity
  std::list<std::vector<double>> spatialmeanAverage;  
  // loop over all velocity components
  for(auto & v : velocity)
  {
    spatialmeanAverage.push_back(spatialMean(v));
  }
  //FIXME update the summation to compute the average
  summation(velocity_u1.size());
  //TODO put in extra function
  for(auto &v : spatialmeanAverage)
  {
    std::transform(v.begin(), v.begin()+nZLayers, sum_layer_dofs.begin(), 
                   v.begin(), std::divides<double>());
  }
  // compute the mean velocity. average in time
  std::deque<std::vector<double>> meanVelocityTimeAverage;
  double T=TDatabase::TimeDB->CURRENTTIME;
  double t0=TDatabase::TimeDB->T0;
  if(T >= t0)
  {
    for(auto v : spatialmeanAverage)
    {
      meanVelocityTimeAverage.push_back(meanVelocity(v));
    }
  }
  else
  {
    //TODO check needed for all components or only for one
    meanVelocityTimeAverage.push_back(meanVelocity(spatialmeanAverage.front()));
  }
  /// compute the eddy viscosity model
  std::array<std::vector<double>, 6> eddyviscosity;
  eddy_viscosity(eddyviscosity, tnse3d);
  /// now computing the mean velocity of the product
  /// u_i*u_j, i,j=1,2,3
  /// first compute the product and then the mean velocities
  std::list<std::vector<double>> ui_uj;
  ui_uj.push_back(product(velocity_u1,velocity_u1)); // u1*u1
  ui_uj.push_back(product(velocity_u2,velocity_u2)); // u2*u2
  ui_uj.push_back(product(velocity_u3,velocity_u3)); // u3*u3
  ui_uj.push_back(product(velocity_u1,velocity_u2)); // u1*u2
  ui_uj.push_back(product(velocity_u1,velocity_u3)); // u1*u3
  ui_uj.push_back(product(velocity_u2,velocity_u3)); // u2*u3
  
  std::list<std::vector<double>> uiujmean;
  for(auto &v : ui_uj)
  {
    uiujmean.push_back(spatialMean(v));
  }  
  // averaging
  for(auto &v : uiujmean)
  {
    std::transform(v.begin(), v.begin()+nZLayers, sum_layer_dofs.begin(), 
                   v.begin(), std::divides<double>());
  }
  // Reynolds stress tensor
  std::deque<std::vector<double>> ReynoldsStress; // R_ij: in volker's paper
  if(T >= t0)
  {
    for(auto v : uiujmean)
    {
      ReynoldsStress.push_back(meanVelocity(v));
    }
  }
  /// compute root mean square (formula 11 Volker's paper)
  std::deque<std::vector<double>> rmsInstnsities;
  cout<<"here it is : " << endl;
  cout<<TDatabase::TimeDB->CURRENTTIME<<endl;
  if(T >= t0)
    rmsInstnsities=getrms(ReynoldsStress,meanVelocityTimeAverage);
  /// save the output into the files
  saveData(meanVelocityTimeAverage, rmsInstnsities, ReynoldsStress);  
}

void ChannelFlowRoutines::summation(size_t length)
{
  sum_layer_dofs.resize(nZLayers);  
  for(size_t i=0; i<length; ++i)
  {
    for(size_t j=0; j<nZLayers; ++j)
    {
      if(fabs(zDofs[i]-zLayers[j]) < 1e-6)
      {
        sum_layer_dofs.at(j)++;
      }
    }
  }
}

std::vector< double >
  ChannelFlowRoutines::spatialMean(std::vector< double > in)
{
  // summation in the computation of mean velocities  
  std::vector<double>temp(nZLayers);
  // cout<<temp.size()<<"  " << in.size() << endl;
  for(size_t i=0; i<in.size(); ++i)
  {
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

std::vector< double > ChannelFlowRoutines::meanVelocity(std::vector<double> vecin)
{
  double step_length=TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double currentTime = TDatabase::TimeDB->CURRENTTIME;
  double tstart = TDatabase::TimeDB->T0;
  
  std::vector<double> temp(nZLayers);
  
  double factor;
  if(currentTime >= tstart)
    factor= step_length/(currentTime - tstart + step_length);
  else
    factor = step_length/(currentTime+step_length);
  
  for(size_t i=0; i<vecin.size(); ++i)
    temp.at(i) = temp.at(i) + factor*(vecin.at(i) - temp.at(i));
  
  return temp;
}

std::vector< double > ChannelFlowRoutines::product(std::vector< double > vec1, 
                                                   std::vector< double > vec2)
{
  std::vector<double> res(vec1.size());
  std::transform(vec1.begin(), vec1.end(), vec2.begin(), res.begin(),
                 std::multiplies<double>());
  return res;
}


void ChannelFlowRoutines::eddy_viscosity(
  std::array< std::vector< double>, int(6)> eddy,
  Time_NSE3D& tnse3d)
{
  size_t nuDofs = tnse3d.get_size();
      /// fill with zeros
  std::fill(eddy.begin(), eddy.end(), std::vector<double>(nuDofs));
  if(!TDatabase::ParamDB->CHANNEL_STATISTICS2_WITH_MODEL)
    return;
  /// else compute the contribution from eddy viscosity model
  std::array<std::vector<double>, 9> gradu;
  double area;
  double re_nr = TDatabase::ParamDB->RE_NR;
  const TFESpace3D& space = tnse3d.get_velocity_space();
  size_t nCells=space.GetCollection()->GetN_Cells();
  size_t nCellsPerLayer;
  double hxhy;
  /// 
  switch(TDatabase::ParamDB->DISCTYPE)
  {
    case GALERKIN:
      break;
    case SMAGORINSKY:
      ChannelFlowRoutines::computeAverageVelocity(gradu, nuDofs, tnse3d);
      // Reynolds number only: 395 or 180 are used
      area = (re_nr == 395) ? 4.*Pi*Pi : 32.*Pi*Pi/3;
      nCellsPerLayer = 2*nCells/(nZLayers-1);
      hxhy = area/(2*nCellsPerLayer);
      Output::print("nCells per layer:  ", nCellsPerLayer, " hxhy: ", hxhy);
      break;
    default:
      ErrThrow("mean velocity for DISCTYPE ", 
               TDatabase::ParamDB->DISCTYPE," is not implemented yet");
  }
  const TFEVectFunct3D& uvec = tnse3d.get_velocity();
  TFEFunction3D* u1 =uvec.GetComponent(0);
  TFEFunction3D* u2 =uvec.GetComponent(1);
  TFEFunction3D* u3 =uvec.GetComponent(2);
  
  std::vector<double> velocity_u1(u1->GetValues(),
                                  u1->GetValues()+nuDofs);
  std::vector<double> velocity_u2(u2->GetValues(),
                                    u2->GetValues()+nuDofs);
  std::vector<double> velocity_u3(u3->GetValues(),
                                    u3->GetValues()+nuDofs);
/*
  double u[3];
  double delta;
  for(size_t i=0; i<nuDofs; ++i)
  {
    for(size_t j=0; j<nZLayers; ++j)
    {
      if(fabs(zDofs[i]-zLayers[j]) < 1e-6)
      {
        // filter width
        switch(TDatabase::ParamDB->DISCTYPE)
        {
          case SMAGORINSKY:
            // only Q2 elements are tested
            if(TDatabase::ParamDB->VELOCITY_SPACE != 2 || 
               TDatabase::ParamDB->VELOCITY_SPACE != 12)
            {
              Output::print("change velocity space to 2 or 12");
              ErrThrow("space ", TDatabase::ParamDB->VELOCITY_SPACE, " is not tested yet");
            }
            /// compute delta (given by the smallest edge), 
            /// take average values between the layers
            if(j==0)
              delta = 2.*zLayers.at(1);
            else
            {
              if(j==nZLayers-1)
                delta = 4-2*zLayers.at(nZLayers-2);
              else
                delta = zLayers.at(j+1) - zLayers.at(j-1);
            }
            break;
          default:
            ErrThrow("DISCTYPE ", TDatabase::ParamDB->DISCTYPE, " is not tested yet");
        }
        u[0] = velocity_u1.at(i);
        u[1] = velocity_u2.at(i);
        u[2] = velocity_u3.at(i);
        //TODO fill the vectors with right entries
        ErrThrow("not completed yet: ");
      }
    }
  }
*/
}

double ChannelFlowRoutines::getFrictionVelocity(std::vector< double > vec)
{
  double temp;
  /// computing the friction velocity
  /// first compute the derivative of mean velocity
  //TODO not clear completely
  std::vector<double> meanDeriv;
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
  // print("mean", meanDeriv);
  std::vector<double> meanVeloDerivTime;
  meanVeloDerivTime.resize(nZLayers);
  //
  meanVeloDerivTime=meanVelocity(meanDeriv);
  double re_nr = TDatabase::ParamDB->RE_NR;
  temp = 1/re_nr*(meanVeloDerivTime.at(0)-meanVeloDerivTime.at(nZLayers-1))*0.5;  
  return temp;
}

std::deque<std::vector< double>> ChannelFlowRoutines::getrms(std::deque<std::vector<double>> reynoldStress, 
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

void ChannelFlowRoutines::saveData(std::deque< std::vector< double > > m, 
                                 std::deque< std::vector< double > > rms,
                                 std::deque< std::vector< double > > R)
{
  double t = TDatabase::TimeDB->CURRENTTIME; 
  double bulkVelocity=0.;
  for(size_t i=0; i<nZLayers; ++i)
  {
    bulkVelocity += 0.5*(m.at(0)[i+1] + m.at(0)[i])*(zLayers[i+1]-zLayers[i]);
  }
  bulkVelocity /=2.;
  Output::print("bulk velocity : t ", t, " ", setw(8), bulkVelocity, " tau ", 
                (TDatabase::ParamDB->INTERNAL_BULK_MEAN-bulkVelocity)+1);
  // set the bulk velocity
  TDatabase::ParamDB->INTERNAL_BULK_SIMULATION = bulkVelocity;

  if(!(t>= TDatabase::TimeDB->T0))
  {
    return;
  }
  // frictionvelocity
  double u_tau = getFrictionVelocity(m.front());

  double re=TDatabase::ParamDB->RE_NR;
  for(size_t i=0; i<nZLayers; ++i)
  {
    Output::print("t ", t, " ", setw(8), zLayers.at(i), " ",setw(8), 
                         re*(1-fabs(1-zLayers[i])), " mu ", setw(8), m.at(0)[i], 
                       " mv ", setw(8), m.at(2)[i], " mw ", setw(8), m.at(1)[i]);
  }
  //
  double rmsu, rmsv, rmsw;
  double R12, R13, R23;
  double R12_abs=0., R13_abs=0., R23_abs=0.;
  for(size_t i=0; i<nZLayers; ++i)
  {
    // formula 11 Volker's paper
    rmsu = sqrt(fabs(rms.at(0)[i]))/u_tau;
    rmsv = sqrt(fabs(rms.at(1)[i]))/u_tau;
    rmsw = sqrt(fabs(rms.at(2)[i]))/u_tau;
    // formula 10 Volker's paper
    R12 = R.at(3)[i] - m.at(0)[i] * m.at(1)[i]; R12 /= (u_tau*u_tau);
    R12_abs += fabs(R12);
    R13 = R.at(4)[i] - m.at(0)[i] * m.at(2)[i]; R13 /= (u_tau*u_tau);
    R13_abs += fabs(R13);
    R23 = R.at(5)[i] - m.at(1)[i] * m.at(2)[i]; R23 /= (u_tau*u_tau);
    R23_abs += fabs(R23);

    Output::print("t ", t, " ", setw(8), zLayers.at(i), " ", setw(8),
                  re*(1-fabs(1-zLayers[i])), " rms_u* ", setw(8), rmsu, 
                  " rms_v* ", setw(8), rmsw, " rms_w* ", setw(8), rmsv,
                  " R_uv ", setw(8), R13, " R_uw ", setw(8), R12, 
                  " R_vw ", setw(8), R23);
  }
  Output::print("u_tau ", u_tau, " zero statistics ", setw(8), R12_abs/nZLayers, 
                " ", setw(8), R13_abs/nZLayers, " ", setw(8), R23_abs/nZLayers);
  // Output::print("u_tau: ", );
}

