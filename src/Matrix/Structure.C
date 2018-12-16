#include <Structure.h>
#include <MooNMD_Io.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <FEDatabase3D.h>
#include <FESpace1D.h>
#include <FESpace2D.h>
#include <FESpace3D.h>

#include <cstring>
#include <vector>
#include <stdlib.h>
#include <algorithm>


/* generate the matrix structure, both space are 1D */
TStructure::TStructure( const TFESpace1D * space )
: TStructure()
{
  int i,j,k,l,n,N_, n1,n2, m;
  int *GlobalNumbers;
  int *BeginIndex;
  int *Numbers, *nieb_Numbers;
  int N_Hanging, Offset;
  int N_Inner, NE, nieb_i, nieb_e, nieb_n1;
  int *KColAux;
  int index, oldindex, EdgeIntegrals=space->IsDGSpace();

  TCollection *Coll;

  TBaseCell *cell, *neigh;
  FE1D CurrentElement, CurrentNeighbour;

  // all dof are treated as unknowns !!!
  // no boundary description is used so far!!!
  N_Inner=space->GetN_Inner();
  ActiveBound = N_Inner;
  nRows = N_Inner;
  nColumns = nRows;
  N_Hanging = 0;
  nHangingEntries=0;
  // AuxPtr[i] will contain an upper bound for the number of ...
  // matrix entries in row i
  l=nRows+1;
  std::vector<int> AuxPtr(l, 0);

  GlobalNumbers=space->GetGlobalNumbers();
  BeginIndex=space->GetBeginIndex();

  // loop over all elements
  N_= space->GetN_Cells();
  Coll = space->GetCollection();
  
  // associate each cell with her number in the collection
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
  }

  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    CurrentElement = space->GetFE1D(i, cell);
    n1 = TFEDatabase2D::GetFE1D(CurrentElement)->GetN_DOF();
    for(j=0;j<n1;j++)
    {
      k=Numbers[j];
      AuxPtr[k]+=n1;

     //DG additional entries necessary if integrals over the edges(points in 1D)
     // appear in the discretization
     if(EdgeIntegrals)
     {
      l=cell->GetN_Edges();                   // # edges
      for(m=0;m<l;m++)                        // for all edges
      {
       neigh = (cell->GetJoint(m))->GetNeighbour(cell);

       if(neigh)
       {
        n = neigh->GetClipBoard();
        CurrentNeighbour = space->GetFE1D(n, neigh);
        n2 = TFEDatabase2D::GetFE1D(CurrentNeighbour)->GetN_DOF();
        AuxPtr[k]+=n2;
       }//if(neigh)
      }
     }//if (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)

    }                                             // for(j=0;j<n1;
  }                                               //  for(i=0;i<N

  nEntries = 0;

//   for(i=0;i<=N_;i++)
//     cout << i << "   " << AuxPtr[i] << endl;
//   cout << endl;
// exit(0);
  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  N_=nRows;
  l=AuxPtr[0];
  AuxPtr[0]=0;
  //   cout << l << "  " << AuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=AuxPtr[i+1];
    AuxPtr[i+1]=AuxPtr[i]+l;
    l=k;
    //     cout << AuxPtr[i+1] << endl;
  }

  /*
  cout << endl;
  cout << "AuxPtr array" << endl;
  for(i=0;i<=N_;i++)
    cout << i << "   " << AuxPtr[i] << endl;
  cout << endl;
  */
// exit(0);

  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[nRows];                               // upper bound for number of matrix entries
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  N_= space->GetN_Cells();
  Coll = space->GetCollection();
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    CurrentElement = space->GetFE1D(i, cell);
    n1 = TFEDatabase2D::GetFE1D(CurrentElement)->GetN_DOF();

    for(j=0;j<n1;j++)                             // test space
    {
      n=Numbers[j];
      for(k=0;k<n1;k++)                           // ansatz space
      {
        m=Numbers[k];
        // this  node is a real node (inner or Neumann)
        index=AuxPtr[n];
        l=KColAux[index];
        // check whether this column is already in this row
        while(l!=-1 && l!=m)
        {
          index++;
          l=KColAux[index];
        }
        if(l==-1)
        {
          // this is a new column for this row
          KColAux[index]=m;
          nEntries++;
        }
      }                                           // endfor k
      
     // DG part
     if(EdgeIntegrals)
     {
      NE=cell->GetN_Edges();                   // # edges
      for(nieb_e=0;nieb_e<NE;nieb_e++)                        // for all edges
      {
       neigh = (cell->GetJoint(nieb_e))->GetNeighbour(cell);

       if(neigh)
       {
        nieb_i = neigh->GetClipBoard();
        CurrentNeighbour = space->GetFE1D(nieb_i, neigh);
        nieb_Numbers=GlobalNumbers+BeginIndex[nieb_i];
        nieb_n1 = TFEDatabase2D::GetFE1D(CurrentNeighbour)->GetN_DOF();
 
        for(k=0;k<nieb_n1;k++)                           // ansatz space
         {
          m=nieb_Numbers[k];
          // this  node is a real node (inner or Neumann)
          index=AuxPtr[n];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
           index++;
           l=KColAux[index];
          }
         if(l==-1)
          {
           // this is a new column for this row
           KColAux[index]=m;
           nEntries++;
          }   
         } // for(k=0;k<nieb_n1;k++)  
       }  //if(neigh)
       
      } //  for(nieb_e=0;nieb_e<NE;nieb_e++)
     }//if (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)

      
    }                                             // endfor j
  }                                               // for(i=0;

  
  /*
  // check
  cout << endl;
  cout << "check" << endl;
  N_=ActiveBound;
  N_=nRows;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=AuxPtr[i+1];
    for(j=AuxPtr[i];j<index;j++) //  && KColAux[j]!=-1
  cout << setw(4) << KColAux[j];
  cout << endl;
  }
  exit(0);
  */

  //   cout << "Number of matrix entries: ";
  //   cout << nEntries << endl;
  //   cout << endl;

  // compress KColAux array to KCol by deleting all -1's
  // build the RowPtr array
  N_=ActiveBound;
  columns.resize(nEntries);
  rows=AuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=AuxPtr[i+1];
    for(j=AuxPtr[i];j<m && KColAux[j]!=-1;j++)
    {
      columns[index]=KColAux[j];
      index++;
    }                                             // endfor j
    rows[i]=oldindex;
    //     cout << setw(4) << i << " rows[i]: " << rows[i] << endl;
  }                                               // endfor i

  // cout << "index: " << index << endl;
  // cout << "rows[N_]: " << rows[N_] << endl;
  Offset=index-rows[N_];
  for(i=0,j=ActiveBound;i<=N_Hanging;i++,j++)
  {
    // cout << "HangingRows[i]: " << HangingRows[i] << endl;
    rows[j]+=Offset;
    // cout << setw(4) << j << " rows[j]: " << rows[j] << endl;
  }

  /*
    // print out the whole matrix structure
    cout << endl;
    N_=nRows;
    for(i=0;i<N_;i++)
    {
      cout << rows[i] << "---" << rows[i+1]-1 << endl;
      cout << "Rows: " << setw(4) << i << ": ";
      end=rows[i+1];
      for(j=rows[i];j<end;j++)
        cout << setw(4) << columns[j];
  cout << endl;
  }
//   */

  // free KColAux
  delete [] KColAux;
#ifndef _MPI
  this->info();
#endif
  this->Sort();
}


/* generate the matrix structure, both space are 2D */
TStructure::TStructure( const TFESpace2D* space )
: TStructure()
{
  int i,j,k,l,n,N_, n1,n2,m,p,q,r;
  int *GlobalNumbers;
  int *BeginIndex;
  int *Numbers, *NumbersNeighbour;

  int N_Dirichlet;
  int N_Inner;
  int N_Hanging;
  int *BoundNodeBounds, N_NonDiri, N_BoundaryNodeTypes;
  int HangingBound;

  int *KColAux, *HangingKColAux;
  int index, oldindex;

  int Offset, *DOF, end;
  THangingNode **HangingNodes;
  THangingNode *hn;

  TCollection *Coll;

  TBaseCell *cell, *neigh;
  FE2D CurrentElement, CurrentNeighbour;

  ActiveBound = space->GetN_ActiveDegrees();
  HangingBound = space->GetHangingBound();

  N_BoundaryNodeTypes = space->GetN_DiffBoundaryNodeTypes();
  BoundNodeBounds = space->GetN_BoundaryNodes();
  N_NonDiri = 0;
  for(i=0;i<N_BoundaryNodeTypes;i++)
    N_NonDiri += BoundNodeBounds[i];

  N_Dirichlet = space->GetN_Dirichlet();
  N_Inner = space->GetN_Inner();
  N_Hanging = space->GetN_Hanging();

  nRows = ActiveBound+N_Hanging+N_Dirichlet;
  nColumns = nRows;
  // assembles matrices without shorter rows for non-active dof
   if (TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE)
  {
    N_NonDiri += N_Dirichlet;
    N_Dirichlet = 0;
    ActiveBound = N_Inner + N_NonDiri;
  }
  ColOrder = 0;
  // AuxPtr[i] will contain an upper bound for the number of
  // matrix entries in row i
  l=nRows+1;
  std::vector<int> AuxPtr(l, 0);

  l=N_Hanging+1;
  std::vector<int> HangingAuxPtr(l, 0);

  GlobalNumbers = space->GetGlobalNumbers();

  BeginIndex = space->GetBeginIndex();

  Offset=ActiveBound;

  // loop over all elements
  N_ = space->GetN_Cells();
  Coll = space->GetCollection();

  // associate each cell with her number in the collection
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);                        // set clipboard to number of the cell in collection
  }

  // loop over the mesh cells
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    CurrentElement = space->GetFE2D(i, cell);
    n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

    // loop over the local degrees of freedom
    for(j=0;j<n1;j++)
    {
      k=Numbers[j];
      if(k<ActiveBound)
      {
        AuxPtr[k]+=n1;

        // additional entries necessary if integrals over the edges
        // appear in the discretization
        if (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)
        {
          l=cell->GetN_Edges();                   // # edges
          for(m=0;m<l;m++)                        // for all edges
          {
                                                  // neighbour cell
            neigh = cell->GetJoint(m)->GetNeighbour(cell);
            if(neigh)
            {
              n = neigh->GetClipBoard();
              CurrentNeighbour = space->GetFE2D(n, neigh);
              n2 = TFEDatabase2D::GetFE2D(CurrentNeighbour)->GetN_DOF();
              AuxPtr[k]+=n2;
            }                                     //endif
          }                                       //endfor m
        }
      }
      else
      {
                                                  // not adjusted for edge stabilization because hanging nodes are no more used
        if(k<HangingBound) HangingAuxPtr[k-Offset]+=n1;
      }                                           // endif
    }                                             // endfor j
  }                                               // endfor i

  // add rows for Dirichlet nodes in  space
  nEntries=N_Dirichlet;
  Offset=ActiveBound+N_Hanging;
  for(i=0,j=Offset;i<N_Dirichlet;i++,j++)
  {
    AuxPtr[j]=1;
  }

  // add couplings for hanging nodes of  space
  HangingNodes = space->GetHangingNodes();
  Offset=ActiveBound;
  for(i=0,j=Offset;i<N_Hanging;i++,j++)
  {
    // cout << "i: " << i << " j: " << j << endl;
    hn=HangingNodes[i];
    // cout << hn;
    // there is the additional entry in diagonal
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes() + 1;
    AuxPtr[j]=n;
    nEntries+=n;
    // cout << "AuxPtr[j]: " << AuxPtr[j] << endl;
  }

  // additional space for storing the columns caused by the
  // hanging nodes of  space => some new columns
  Offset=ActiveBound;
  HangingNodes = space->GetHangingNodes();
  // cout << "N_Hanging: " << N_Hanging << endl;
  for(i=0;i<N_Hanging;i++)
  {
    hn=HangingNodes[i];
    m=HangingAuxPtr[i];
    // cout << "HangingAuxPtr[i]: " << m << endl;
    DOF=hn->GetDOF();
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
    for(j=0;j<n;j++)
    {
      k=DOF[j];
      if(k<ActiveBound)
        AuxPtr[k] += m;
    }
  }

  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  N_=nRows;
  l=AuxPtr[0];
  AuxPtr[0]=0;

  for(i=0;i<N_;i++)
  {
    k=AuxPtr[i+1];
    AuxPtr[i+1]=AuxPtr[i]+l;
    l=k;
  }

  /*
  cout << endl;
  cout << "AuxPtr array" << endl;
  for(i=0;i<=N_;i++)
    cout << i << "   " << AuxPtr[i] << endl;
  cout << endl;
  */

  // cout << "Upper bound: " << AuxPtr[nRows] << endl;

  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[nRows];                               // upper bound for number of matrix entries
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  // sum up the array HangingAuxPtr, HangingAuxPtr[i] will now contain
  // the index for Hanging KColAux array, the column numbers for row i
  // are in the intervall // [ HangingAuxPtr[i], HangingAuxPtr[i+1] )
  N_=N_Hanging;
  l=HangingAuxPtr[0];
  HangingAuxPtr[0]=0;
  // cout << HangingAuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=HangingAuxPtr[i+1];
    HangingAuxPtr[i+1]=HangingAuxPtr[i]+l;
    l=k;
    // cout << HangingAuxPtr[i+1] << endl;
  }

  /*
  cout << endl;
  cout << "HangingAuxPtr array" << endl;
  for(i=0;i<=N_;i++)
    cout << i << "   " << HangingAuxPtr[i] << endl;
  cout << endl;
  */

  // cout << "Upper bound for hanging nodes: ";
  // cout << HangingAuxPtr[N_Hanging] << endl;

  // get memory for HangingKColAux array, initialize it with -1
  l=HangingAuxPtr[N_Hanging];                     //upper bound for number of matrix entries
  // cout << "l= " << l << endl;
  HangingKColAux=new int[l];
  memset(HangingKColAux, -1, sizeof(int)*l);
  nHangingEntries=0;

  N_ = space->GetN_Cells();
  Coll = space->GetCollection();
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers=GlobalNumbers+BeginIndex[i];

    CurrentElement = space->GetFE2D(i, cell);
    n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

    for(j=0;j<n1;j++)
    {
      for(k=0;k<n1;k++)
      {
        m=Numbers[k];
        n=Numbers[j];
        if(n<ActiveBound)
        {
          // this  node is a real node (inner or Neumann)
          index=AuxPtr[n];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
            index++; l=KColAux[index];
          }
          if(l==-1)
          {
            // this is a new column for this row
            KColAux[index]=m;
            nEntries++;
          }
        }
        else
        {
          if(n<HangingBound)
          {
            // this node is a hanging node in  space
            index=HangingAuxPtr[n-ActiveBound];
            l=HangingKColAux[index];
            // check whether this column is already in this row
            while(l!=-1 && l!=m)
            {
              index++; l=HangingKColAux[index];
            }
            if(l==-1)
            {
              // this is a new column for this row
              HangingKColAux[index]=m;
              nHangingEntries++;
            }
          }                                       // endif
        }                                         // endif
      }                                           // endfor k

      // additional entries necessary if integrals over the edges
      // appear in the discretization
      if (TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS)
      {
        r=cell->GetN_Edges();
        for(p=0;p<r;p++)                          // for all edges
        {
                                                  // neighbour cell
          neigh=cell->GetJoint(p)->GetNeighbour(cell);
          if(neigh)
          {
            q = neigh->GetClipBoard();
            NumbersNeighbour=GlobalNumbers+BeginIndex[q];

            CurrentNeighbour = space->GetFE2D(q, neigh);
            n2 = TFEDatabase2D::GetFE2D(CurrentNeighbour)->GetN_DOF();

            for(k=0;k<n2;k++)
            {
              m=NumbersNeighbour[k];
              n=Numbers[j];
              if(n<ActiveBound)
              {
                // this  node is a real node (inner or Neumann)
                index=AuxPtr[n];
                l=KColAux[index];
                // check whether this column is already in this row
                while(l!=-1 && l!=m)
                {
                  index++; l=KColAux[index];
                }
                if(l==-1)
                {
                  // this is a new column for this row
                  KColAux[index]=m;
                  nEntries++;
                }
              }
              else
              {
                if(n<HangingBound)
                {
                  // this node is a hanging node in  space
                  index=HangingAuxPtr[n-ActiveBound];
                  l=HangingKColAux[index];
                  // check whether this column is already in this row
                  while(l!=-1 && l!=m)
                  {
                    index++; l=HangingKColAux[index];
                  }
                  if(l==-1)
                  {
                    // this is a new column for this row
                    HangingKColAux[index]=m;
                    nHangingEntries++;
                  }
                }                                 // endif
              }                                   // endif
            }                                     // endfor k
          }                                       // endif
        }                                         // endfor p
      }
    }                                             // endfor j
  }                                               // endfor i

  // for(i=0;i<AuxPtr[nRows]_;i++)
  // {
  // cout << "KColAux: " << KColAux[i] << endl;
  // }

  /*
  // check hanging node data
  cout << endl;
  cout << "check hanging node data" << endl;
  N_=N_Hanging;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=HangingAuxPtr[i+1];
    for(j=HangingAuxPtr[i];j<index && HangingKColAux[j]!=-1;j++)
      cout << setw(4) << HangingKColAux[j];
  cout << endl;
  }
  */

  /*
  cout << "Number of matrix entries (hanging nodes): ";
  cout << nHangingEntries << endl;
  cout << endl;
  */

  // compress HangingKColAux array to hangingColums by deleting all -1's
  // build the HangingRows array
  N_=N_Hanging;
  hangingColums.resize(nHangingEntries, 0);
  HangingRows=HangingAuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=HangingAuxPtr[i+1];
    for(j=HangingAuxPtr[i];j<m && HangingKColAux[j]!=-1;j++)
    {
      hangingColums[index]=HangingKColAux[j];
      index++;
    }                                             // endfor j
    HangingRows[i]=oldindex;
    // cout << HangingRows[i] << endl;
  }                                               // endfor i
  HangingRows[N_]=index;

  // free HangingKColAux
  delete [] HangingKColAux;

  // add the additional columns from hanging nodes to other nodes
  Offset=ActiveBound;
  HangingNodes = space->GetHangingNodes();
  for(i=0;i<N_Hanging;i++)
  {
    // cout << "hanging node: " << i << endl;
    hn=HangingNodes[i];
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
    DOF=hn->GetDOF();

    for(j=0;j<n;j++)                              // loop over all nodes in coupling
    {
      k=DOF[j];
      // cout << "k= " << k << endl;
      if(k<ActiveBound)
      {
        // node is either inner or Neumann node

        end=HangingAuxPtr[i+1];
        for(oldindex=HangingAuxPtr[i];oldindex<end;oldindex++)
        {
          m=hangingColums[oldindex];
          // cout << "m= " << m << endl;
          index=AuxPtr[k];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
            index++; l=KColAux[index];
          }
          if(l==-1)
          {
            // this is a new column for this row
            KColAux[index]=m;
            nEntries++;
          }
        }                                         // endfor
      }                                           // endif
    }                                             // endfor j
  }                                               // endfor i

  /*
  // check
  cout << endl;
  cout << "check" << endl;
  N_=ActiveBound;
  N_=nRows;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=AuxPtr[i+1];
    for(j=AuxPtr[i];j<index && KColAux[j]!=-1;j++)
  cout << setw(4) << KColAux[j];
  cout << endl;
  }
  */

  /*
  cout << "Number of matrix entries: ";
  cout << nEntries << endl;
  cout << endl;
  */

  // compress KColAux array to KCol by deleting all -1's
  // build the rows array
  N_=ActiveBound;
  columns.resize(nEntries);
  rows=AuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=AuxPtr[i+1];
    for(j=AuxPtr[i];j<m && KColAux[j]!=-1;j++)
    {
      columns[index]=KColAux[j];
      index++;
    }                                             // endfor j
    rows[i]=oldindex;
    // cout << setw(4) << i << " rows[i]: " << rows[i] << endl;
  }                                               // endfor i

  // cout << "index: " << index << endl;
  // cout << "rows[N_]: " << rows[N_] << endl;
  Offset=index-rows[N_];
  for(i=0,j=ActiveBound;i<=N_Hanging;i++,j++)
  {
    // cout << "HangingRows[i]: " << HangingRows[i] << endl;
    rows[j]+=Offset;
    // cout << setw(4) << j << " rows[j]: " << rows[j] << endl;
  }

  j=ActiveBound+N_Hanging;
  Offset=rows[j];
  for(i=0;i<N_Dirichlet;i++,j++)
  {
    rows[j+1]=rows[j]+1;
    // cout << setw(4) << j+1 << " rows[j+1]: " << rows[j+1] << endl;
  }

  // add information for hanging and Dirichlet nodes into matrix
  HangingNodes = space->GetHangingNodes();
  Offset=ActiveBound;
  m=ActiveBound;
  for(i=0,j=Offset;i<N_Hanging;i++,j++)
  {
    // cout << "i: " << i << " j: " << j << endl;
    hn=HangingNodes[i];
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
    DOF=hn->GetDOF();
    index=AuxPtr[j];
    columns[index]=m;
    index++;
    m++;
    for(k=0;k<n;k++)
    {
      // cout << "index: " << index << " DOF[k]" << DOF[k] << endl;
      columns[index]=DOF[k];
      index++;
    }
  }

  // add Dirichlet rows
  j=HangingBound;
  index=rows[ActiveBound+N_Hanging];
  for(i=0;i<N_Dirichlet;i++)
  {
    // cout << "index: " << index << endl;
    columns[index]=j;
    j++;
    index++;
  }

  /*    // print out the whole matrix structure
      cout << endl;
      N_=nRows;
      for(i=0;i<N_;i++)
      {
        cout << rows[i] << "---" << rows[i+1]-1 << endl;
        cout << "Rows: " << setw(4) << i << ": ";
        end=rows[i+1];
        for(j=rows[i];j<end;j++)
          cout << setw(4) << columns[j];
        cout << endl;
  }
  */

  // free KColAux
  delete [] KColAux;

#ifndef _MPI
  this->info();
#endif
  this->Sort();
}

#ifdef __3D__
/* generate the matrix structure, both spaces are 3D */
TStructure::TStructure( const TFESpace3D *space )
{
  int i,j,k,l,n,N_, n1, m; 
  int *Numbers;

  int N_Dirichlet;
  int N_Inner;
  int N_Hanging;
  int *BoundNodeBounds, N_NonDiri, N_BoundaryNodeTypes;
  int HangingBound;

  int *KColAux, *HangingKColAux;
  int index, oldindex;

  int Offset, *DOF, end;
  THangingNode **HangingNodes;
  THangingNode *hn;

  TCollection *Coll;

  TBaseCell *cell;
  FE3D CurrentElement;

  ActiveBound = space->GetN_ActiveDegrees();
  HangingBound = space->GetHangingBound();

  N_BoundaryNodeTypes = space->GetN_DiffBoundaryNodeTypes();
  BoundNodeBounds = space->GetN_BoundaryNodes();
  N_NonDiri = 0;
  for(i=0;i<N_BoundaryNodeTypes;i++)
    N_NonDiri += BoundNodeBounds[i];

  N_Dirichlet = space->GetN_Dirichlet();
  N_Inner = space->GetN_Inner();
  N_Hanging = space->GetN_Hanging();

  nRows = ActiveBound+N_Hanging+N_Dirichlet;
  nColumns = nRows;
  // assembles matrices without shorter rows for non-active dof
  if (TDatabase::ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE)
  {
    N_NonDiri += N_Dirichlet;
    N_Dirichlet = 0;
    ActiveBound = N_Inner + N_NonDiri;
  }
  ColOrder = 0;
  // AuxPtr[i] will contain an upper bound for the number of 
  // matrix entries in row i
  l=nRows+1;
  std::vector<int> AuxPtr(l, 0);

  l=N_Hanging+1;
  std::vector<int> HangingAuxPtr(l, 0);

  Offset=ActiveBound;

  // loop over all elements 
  N_ = space->GetN_Cells();
  Coll = space->GetCollection();
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers = space->GetGlobalDOF(i);

    CurrentElement = space->GetFE3D(i, cell);
    n1 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

    for(j=0;j<n1;j++)
    {
      k=Numbers[j];
      if(k<ActiveBound) 
      {
        AuxPtr[k]+=n1;
      }
      else
      {
        if(k<HangingBound) HangingAuxPtr[k-Offset]+=n1;
      } // endif
    } // endfor j
  } // endfor i

  // add rows for Dirichlet nodes in  space
  nEntries=N_Dirichlet;
  Offset=ActiveBound+N_Hanging;
  for(i=0,j=Offset;i<N_Dirichlet;i++,j++)
  {
    AuxPtr[j]=1;
  }

// #ifdef __3D__
  // add couplings for hanging nodes of  space
  HangingNodes = space->GetHangingNodes();
  Offset=ActiveBound;
  for(i=0,j=Offset;i<N_Hanging;i++,j++)
  {
    // cout << "i: " << i << " j: " << j << endl;
    hn=HangingNodes[i];
    // there is the additional entry in diagonal
    n=TFEDatabase3D::GetHNDesc3D(hn->GetType())->GetN_Nodes() + 1;
    AuxPtr[j]=n;
    nEntries+=n;
    // cout << "AuxPtr[j]: " << AuxPtr[j] << endl;
  }

  // additional space for storing the columns caused by the 
  // hanging nodes of  space => some new columns
  Offset=ActiveBound;
  HangingNodes = space->GetHangingNodes();
  // cout << "N_Hanging: " << N_Hanging << endl;
  for(i=0;i<N_Hanging;i++)
  {
    hn=HangingNodes[i];
    m=HangingAuxPtr[i];
    // cout << "HangingAuxPtr[i]: " << m << endl;
    DOF=hn->GetDOF();
    n=TFEDatabase3D::GetHNDesc3D(hn->GetType())->GetN_Nodes();
    for(j=0;j<n;j++)
    {
      k=DOF[j];
      if(k<ActiveBound)
        AuxPtr[k] += m;
    }
  }
// #endif

  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  N_=nRows;
  l=AuxPtr[0];
  AuxPtr[0]=0;
  // cout << AuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=AuxPtr[i+1];
    AuxPtr[i+1]=AuxPtr[i]+l;
    l=k;
    // cout << AuxPtr[i+1] << endl;
  }

  /*
  cout << endl;
  cout << "AuxPtr array" << endl;
  for(i=0;i<=N_;i++)
    cout << i << "   " << AuxPtr[i] << endl;
  cout << endl;
  */

  // cout << "Upper bound: " << AuxPtr[nRows] << endl;

  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[nRows]; // upper bound for number of matrix entries 
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

// #ifdef __3D__
  // sum up the array HangingAuxPtr, HangingAuxPtr[i] will now contain 
  // the index for Hanging KColAux array, the column numbers for row i 
  // are in the intervall // [ HangingAuxPtr[i], HangingAuxPtr[i+1] )
  N_=N_Hanging;
  l=HangingAuxPtr[0];
  HangingAuxPtr[0]=0;
  // cout << HangingAuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=HangingAuxPtr[i+1];
    HangingAuxPtr[i+1]=HangingAuxPtr[i]+l;
    l=k;
    // cout << HangingAuxPtr[i+1] << endl;
  }

  /*
  cout << endl;
  cout << "HangingAuxPtr array" << endl;
  for(i=0;i<=N_;i++)
    cout << i << "   " << HangingAuxPtr[i] << endl;
  cout << endl;
  */

  // cout << "Upper bound for hanging nodes: ";
  // cout << HangingAuxPtr[N_Hanging] << endl;

  // get memory for HangingKColAux array, initialize it with -1
  l=HangingAuxPtr[N_Hanging]; //upper bound for number of matrix entries 
  // cout << "l= " << l << endl;
  HangingKColAux=new int[l];
  memset(HangingKColAux, -1, sizeof(int)*l);
// #endif
  nHangingEntries=0;

  N_ = space->GetN_Cells();
  Coll = space->GetCollection();
  for(i=0;i<N_;i++)
  {
    cell = Coll->GetCell(i);
    Numbers = space->GetGlobalDOF(i);

    CurrentElement = space->GetFE3D(i, cell);
    n1 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

    for(j=0;j<n1;j++)
    {
      for(k=0;k<n1;k++)
      {
        m=Numbers[k];
        n=Numbers[j];
        if(n<ActiveBound)
        {
          // this  node is a real node (inner or Neumann)
          index=AuxPtr[n];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
            index++; l=KColAux[index];
          }
          if(l==-1)
          {
            // this is a new column for this row
            KColAux[index]=m;
            nEntries++;
          }
        }
        else
        {
// #ifdef __3D__
          if(n<HangingBound)
          {
            // this node is a hanging node in  space
            index=HangingAuxPtr[n-ActiveBound];
            l=HangingKColAux[index];
            // check whether this column is already in this row
            while(l!=-1 && l!=m)
            {
              index++; l=HangingKColAux[index];
            }
            if(l==-1)
            {
              // this is a new column for this row
              HangingKColAux[index]=m;
              nHangingEntries++;
            }
          } // endif
// #endif
        } // endif
      } // endfor k
    } // endfor j
  } // endfor i

// #ifdef __3D__
  // check hanging node data
/*
  cout << endl;
  cout << "check hanging node data" << endl;
  N_=N_Hanging;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=HangingAuxPtr[i+1];
    for(j=HangingAuxPtr[i];j<index && HangingKColAux[j]!=-1;j++) 
      cout << setw(4) << HangingKColAux[j];
    cout << endl;
  }
*/
  
  // cout << "Number of matrix entries (hanging nodes): ";
  // cout << nHangingEntries << endl;
  // cout << endl;

  // compress HangingKColAux array to hangingColums by deleting all -1's
  // build the HangingRows array
  N_=N_Hanging;
  hangingColums.resize(nHangingEntries, 0);
  HangingRows=HangingAuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=HangingAuxPtr[i+1];
    for(j=HangingAuxPtr[i];j<m && HangingKColAux[j]!=-1;j++)
    {
      hangingColums[index]=HangingKColAux[j];
      index++;
    } // endfor j
    HangingRows[i]=oldindex;
    // cout << HangingRows[i] << endl;
  } // endfor i
  HangingRows[N_]=index;

  // free HangingKColAux
  delete[] HangingKColAux;

  // add the additional columns from hanging nodes to other nodes
  Offset=ActiveBound;
  HangingNodes = space->GetHangingNodes();
  for(i=0;i<N_Hanging;i++)
  {
    // cout << "hanging node: " << i << endl;
    hn=HangingNodes[i];
    n=TFEDatabase3D::GetHNDesc3D(hn->GetType())->GetN_Nodes();
    DOF=hn->GetDOF();

    for(j=0;j<n;j++) // loop over all nodes in coupling
    {
      k=DOF[j];
      // cout << "k= " << k << endl;
      if(k<ActiveBound)
      {
        // node is either inner or Neumann node

        end=HangingAuxPtr[i+1];
        for(oldindex=HangingAuxPtr[i];oldindex<end;oldindex++)
        {
          m=hangingColums[oldindex];
          // cout << "m= " << m << endl;
          index=AuxPtr[k];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
            index++; l=KColAux[index];
          }
          if(l==-1)
          {
            // this is a new column for this row
            KColAux[index]=m;
            nEntries++;
          }
        } // endfor
      } // endif
    } // endfor j
  } // endfor i
// #endif

  /*
  // check
  cout << endl;
  cout << "check" << endl;
  N_=ActiveBound;
  N_=nRows;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=AuxPtr[i+1];
    for(j=AuxPtr[i];j<index && KColAux[j]!=-1;j++) 
      cout << setw(4) << KColAux[j];
    cout << endl;
  }
  */
  
  // cout << "Number of matrix entries: ";
  // cout << nEntries << endl;
  // cout << endl;

  // compress KColAux array to KCol by deleting all -1's
  // build the rows array
  N_=ActiveBound;
  columns.resize(nEntries);
  rows=AuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=AuxPtr[i+1];
    for(j=AuxPtr[i];j<m && KColAux[j]!=-1;j++)
    {
      columns[index]=KColAux[j];
      index++;
    } // endfor j
    rows[i]=oldindex;
    // cout << setw(4) << i << " rows[i]: " << rows[i] << endl;
  } // endfor i

  // cout << "index: " << index << endl;
  // cout << "rows[N_]: " << rows[N_] << endl;
  Offset=index-rows[N_];
  for(i=0,j=ActiveBound;i<=N_Hanging;i++,j++)
  {
    // cout << "HangingRows[i]: " << HangingRows[i] << endl;
    rows[j]+=Offset;
    // cout << setw(4) << j << " rows[j]: " << rows[j] << endl;
  }

  j=ActiveBound+N_Hanging;
  Offset=rows[j];
  for(i=0;i<N_Dirichlet;i++,j++)
  {
    rows[j+1]=rows[j]+1;
    // cout << setw(4) << j+1 << " rows[j+1]: " << rows[j+1] << endl;
  }

// #ifdef __3D__
  // add information for hanging and Dirichlet nodes into matrix
  HangingNodes = space->GetHangingNodes();
  Offset=ActiveBound;
  m=ActiveBound;
  for(i=0,j=Offset;i<N_Hanging;i++,j++)
  {
    // cout << "i: " << i << " j: " << j << endl;
    hn=HangingNodes[i];
    n=TFEDatabase3D::GetHNDesc3D(hn->GetType())->GetN_Nodes();
    DOF=hn->GetDOF();
    index=AuxPtr[j];
    columns[index]=m;
    index++;
    m++;
    for(k=0;k<n;k++)
    {
      // cout << "index: " << index << " DOF[k]" << DOF[k] << endl;
      columns[index]=DOF[k];
      index++;
    }
  }
// #endif

  // add Dirichlet rows
  j=HangingBound;
  index=rows[ActiveBound+N_Hanging];
  for(i=0;i<N_Dirichlet;i++)
  {
    // cout << "index: " << index << endl;
    columns[index]=j;
    j++;
    index++;
  }

/*
  // print out the whole matrix structure
  cout << endl;
  N_=nRows;
  for(i=0;i<N_;i++)
  {
    cout << rows[i] << "---" << rows[i+1]-1 << endl;
    cout << "Rows: " << setw(4) << i << ": ";
    end=rows[i+1];
    for(j=rows[i];j<end;j++)
      cout << setw(4) << columns[j];
    cout << endl;
  }
*/

  // free KColAux
  delete[] KColAux;


#ifndef _MPI
  this->info();
#endif
  this->Sort();
} 

#endif // 3D



/* generate the matrix structure, both spaces are 2D */
TStructure::TStructure(const TFESpace2D* testspace,
                       const TFESpace2D* ansatzspace, bool is_empty)
 : TStructure()
{
  nRows = testspace->GetN_DegreesOfFreedom();
  nColumns = ansatzspace->GetN_DegreesOfFreedom();
  ActiveBound = testspace->GetN_ActiveDegrees();
  rows = std::vector<int>(nRows+1, 0);

  if (is_empty)
  {//no need to create anything else...
    return;
  }

  TCollection *coll;
  TBaseCell *cell;
  int i,j,k,l,m,n, N_, n1, n2, index, oldindex;
  int TestN_BoundNodeTypes, AnsatzN_BoundNodeTypes;
  int *TestN_BoundNodes, *AnsatzN_BoundNodes;
  int TestSumBoundNodes, AnsatzSumBoundNodes;
  int *TestGlobalNumbers, *AnsatzGlobalNumbers;
  int *TestNumbers, *AnsatzNumbers;
  int *TestBeginIndex, *AnsatzBeginIndex;
  int *KColAux, end;
  FE2D CurrentElement; 

  int TestN_Hanging, AnsatzN_Hanging;
  int TestActiveBound, TestHangingBound;
  int AnsatzActiveBound, AnsatzHangingBound;
  int Offset, *DOF;
  THangingNode **TestHangingNodes;
  THangingNode **AnsatzHangingNodes;
  THangingNode *hn;
  int *HangingKColAux;

  int N, N1, N2;

  // test if both spaces are defined on the same triangulation
  if(testspace->GetCollection() != ansatzspace->GetCollection())
  {
    ErrThrow("Structure2D.C : grid for test and ansatz space is not the same!");
  }

  // get collection of mesh cells
  coll = testspace->GetCollection();
  
  // test space and ansatz space differ
  // get information from the spaces
  AnsatzN_BoundNodeTypes = ansatzspace->GetN_DiffBoundaryNodeTypes();
  TestN_BoundNodeTypes = testspace->GetN_DiffBoundaryNodeTypes();

  //int AnsatzN_Dirichlet = ansatzspace->GetN_Dirichlet();
  //int TestN_Dirichlet = testspace->GetN_Dirichlet();

  //int AnsatzN_Inner = ansatzspace->GetN_Inner();
  //int TestN_Inner = testspace->GetN_Inner();

  AnsatzN_BoundNodes = ansatzspace->GetN_BoundaryNodes();
  TestN_BoundNodes = testspace->GetN_BoundaryNodes();

  AnsatzSumBoundNodes = 0;
  for(i=0;i<AnsatzN_BoundNodeTypes;i++)
    AnsatzSumBoundNodes += AnsatzN_BoundNodes[i];

  TestSumBoundNodes = 0;
  for(i=0;i<TestN_BoundNodeTypes;i++)
    TestSumBoundNodes += TestN_BoundNodes[i];

  AnsatzGlobalNumbers = ansatzspace->GetGlobalNumbers();
  TestGlobalNumbers = testspace->GetGlobalNumbers();

  AnsatzBeginIndex = ansatzspace->GetBeginIndex();
  TestBeginIndex = testspace->GetBeginIndex();

  TestN_Hanging = testspace->GetN_Hanging();
  TestActiveBound = testspace->GetN_ActiveDegrees();
  TestHangingBound= testspace->GetHangingBound();
  TestHangingNodes = testspace->GetHangingNodes();

  AnsatzN_Hanging = ansatzspace->GetN_Hanging();
  AnsatzActiveBound = ansatzspace->GetN_ActiveDegrees();
  AnsatzHangingBound= ansatzspace->GetHangingBound();
  AnsatzHangingNodes = ansatzspace->GetHangingNodes();

  // AuxPtr[i] will contain an upper bound for the number of 
  // matrix entries in row i
  l = nRows + 1;
  std::vector<int> AuxPtr(l);

  l=TestN_Hanging+1;
  std::vector<int> HangingAuxPtr(l, 0);

  Offset = TestActiveBound;

  // number of mesh cells
  N_ = coll->GetN_Cells();
  // loop over all mesh cells
  for(i=0;i<N_;i++)
  {
    // get cell i 
    cell = coll->GetCell(i);
    // get global number of d.o.f. of test and ansatz function 
    // which live on this mesh cell
    TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
    AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[i];
    
    // get fe spaces on the mesh cell
    CurrentElement = testspace->GetFE2D(i, cell);
    n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

    CurrentElement = ansatzspace->GetFE2D(i, cell);
    n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
    // update vector which stores information on the length of the rows
    for(j=0;j<n1;j++)
    {
      k = TestNumbers[j];
      // real dof or Dirichlet node
      if(k<TestActiveBound || k>=TestHangingBound)
      {
        AuxPtr[k] += n2;
      }
      else
      {
        HangingAuxPtr[k-Offset]+=n2;
      } // endif

      if(AnsatzN_Hanging)
      {
        for(l=0;l<n2;l++)
        {
          m = AnsatzNumbers[l];
          if(m>=AnsatzActiveBound && m<AnsatzHangingBound)
          {
            // hanging node in ansatz space
            hn = AnsatzHangingNodes[m-AnsatzActiveBound];
            n = TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
            AuxPtr[k] += n;
          }
        }
      } // endif AnsatzN_Hanging
    } // endfor j
  } // endfor i

  // additional space for storing the columns caused by the 
  // hanging nodes of  space => some new columns
  Offset=TestActiveBound;
  // cout << "TestN_Hanging: " << TestN_Hanging << endl;
  for(i=0;i<TestN_Hanging;i++)
  {
    hn=TestHangingNodes[i];
    m=HangingAuxPtr[i];
    // cout << "HangingAuxPtr[i]: " << m << endl;
    DOF=hn->GetDOF();
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
    for(j=0;j<n;j++)
    {
      k=DOF[j];
      if(k<TestActiveBound)
        AuxPtr[k] += m;
    }
  }

  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  N_=nRows;
  l=AuxPtr[0];
  AuxPtr[0]=0;
  // cout << AuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=AuxPtr[i+1];
    // cout << "i: " << i << "  l: " << l << endl;
    AuxPtr[i+1]=AuxPtr[i]+l;
    l=k;
    // cout << AuxPtr[i+1] << endl;
  }

  // cout << "Upper bound: " << AuxPtr[nRows] << endl;
  
  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[nRows]; // upper bound for number of matrix entries 
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  // sum up the array HangingAuxPtr, HangingAuxPtr[i] will now contain 
  // the index for Hanging KColAux array, the column numbers for row i 
  // are in the intervall // [ HangingAuxPtr[i], HangingAuxPtr[i+1] )
  N_=TestN_Hanging;
  l=HangingAuxPtr[0];
  HangingAuxPtr[0]=0;
//   cout << HangingAuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=HangingAuxPtr[i+1];
    HangingAuxPtr[i+1]=HangingAuxPtr[i]+l;
    l=k;
    // cout << HangingAuxPtr[i+1] << endl;
  }

  // get memory for HangingKColAux array, initialize it with -1
  l=HangingAuxPtr[TestN_Hanging]; //upper bound for number of matrix entries 
  // cout << "l= " << l << endl;
  HangingKColAux=new int[l];
  memset(HangingKColAux, -1, sizeof(int)*l);
  nHangingEntries=0;

  nEntries=0;
  N_ = coll->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
    AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[i];

    cell = coll->GetCell(i);

    CurrentElement = testspace->GetFE2D(i, cell);
    n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

    CurrentElement = ansatzspace->GetFE2D(i, cell);
    n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

    // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;
    
    for(j=0;j<n1;j++)
    {
      for(k=0;k<n2;k++)
      {
        m=AnsatzNumbers[k];
        n=TestNumbers[j];

        if(n<TestActiveBound || n>=TestHangingBound)
        {
          index=AuxPtr[n];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m && index<AuxPtr[n+1])
          {
            index++; l=KColAux[index];
          }
          if(l==-1)
          {
            // this is a new column for this row
            KColAux[index]=m;
            nEntries++;
          }
        }
        else
        {
          // this node is a hanging node in  space
          index=HangingAuxPtr[n-TestActiveBound];
          l=HangingKColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
            index++; l=HangingKColAux[index];
          }
          if(l==-1)
          {
            // this is a new column for this row
            HangingKColAux[index]=m;
            nHangingEntries++;
          }
        } // endif

        // hanging nodes in ansatz space
        if(m>=AnsatzActiveBound && m<AnsatzHangingBound)
        {
          hn = AnsatzHangingNodes[m-AnsatzActiveBound];
          N2 = TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
          DOF = hn->GetDOF();
          for(N1=0;N1<N2;N1++)
          {
            index=AuxPtr[n];
            l = KColAux[index];
            N = DOF[N1];
            while(l!=-1 && l!=N && index<AuxPtr[n+1])
            {
              index++; l=KColAux[index];
            }
            if(l==-1)
            {
              // this is a new column for this row
              KColAux[index]=N;
              nEntries++;
            }
          } // endfor N1
        } // endif
      } // endfor k
    } // endfor j
  } // endfor i
  
/*
  // check
  cout << endl;
  cout << "check" << endl;
  N_=nRows;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=AuxPtr[i+1];
    for(j=AuxPtr[i];j<index && KColAux[j]!=-1;j++) 
      cout << setw(4) << KColAux[j];
    cout << endl;
  }
*/
  
  /*
  cout << "Number of matrix entries: ";
  cout << nEntries << endl;
  cout << endl;
  */

  // compress HangingKColAux array to hangingColums by deleting all -1's
  // build the HangingRows array
  N_=TestN_Hanging;
  hangingColums.resize(nHangingEntries);
  HangingRows=HangingAuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=HangingAuxPtr[i+1];
    for(j=HangingAuxPtr[i];j<m && HangingKColAux[j]!=-1;j++)
    {
      hangingColums[index]=HangingKColAux[j];
      index++;
    } // endfor j
    HangingRows[i]=oldindex;
    // cout << HangingRows[i] << endl;
  } // endfor i
  HangingRows[N_]=index;

  // free HangingKColAux
  delete [] HangingKColAux;

  // add the additional columns from hanging nodes to other nodes
  Offset = TestActiveBound;
  for(i=0;i<TestN_Hanging;i++)
  {
    // cout << "hanging node: " << i << endl;
    hn=TestHangingNodes[i];
    n=TFEDatabase2D::GetHNDesc2D(hn->GetType())->GetN_Nodes();
    DOF=hn->GetDOF();

    for(j=0;j<n;j++) // loop over all nodes in coupling
    {
      k=DOF[j];
      // cout << "k= " << k << endl;
      if(k<TestActiveBound)
      {
        // node is either inner or Neumann node

        end=HangingAuxPtr[i+1];
        for(oldindex=HangingAuxPtr[i];oldindex<end;oldindex++)
        {
          m=hangingColums[oldindex];
          // cout << "m= " << m << endl;
          index=AuxPtr[k];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m)
          {
            index++; l=KColAux[index];
          }
          if(l==-1)
          {
            // this is a new column for this row
            KColAux[index]=m;
            nEntries++;
          }
        } // endfor
      } // endif
    } // endfor j
  } // endfor i

  // compress KColAux array to KCol by deleting all -1's
  // build the rows array
  N_=nRows;
  columns.resize(nEntries);
  rows=AuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=AuxPtr[i+1];
    for(j=AuxPtr[i];j<m && KColAux[j]!=-1;j++)
    {
      columns[index]=KColAux[j];
      index++;
    } // endfor j
    rows[i]=oldindex;
    // cout << setw(4) << i << " rows[i]: " << rows[i] << endl;
  } // endfor i
  rows[nRows]=nEntries;

  Sort();
  
/*
  // print out the whole matrix structure
  cout << endl;
  N_=nRows;
  for(i=0;i<N_;i++)
  {
    cout << rows[i] << "---" << rows[i+1]-1 << endl;
    cout << "Rows: " << setw(4) << i << ": ";
    end=rows[i+1];
    for(j=rows[i];j<end;j++)
      cout << setw(4) << columns[j];
    cout << endl;
  }
*/
  
  // free KColAux
  delete [] KColAux;

#ifndef _MPI
  this->info();
#endif
  this->Sort();
}

#ifdef __3D__
/* generate the matrix structure, both spaces are 3D */
TStructure::TStructure(const TFESpace3D *testspace, const TFESpace3D *ansatzspace, bool is_empty)
 : TStructure()
{
  nRows = testspace->GetN_DegreesOfFreedom();
  nColumns = ansatzspace->GetN_DegreesOfFreedom();
  ActiveBound = testspace->GetN_ActiveDegrees();
  rows = std::vector<int>(nRows+1, 0);

  if (is_empty)
  {//no need to create anything else...
    return;
  }

  TCollection *coll;
  TBaseCell *cell;
  int i,j,k,l,m,n, N_, n1, n2, index, oldindex;
  int TestN_BoundNodeTypes, AnsatzN_BoundNodeTypes;
  int *TestN_BoundNodes, *AnsatzN_BoundNodes;
  int TestSumBoundNodes, AnsatzSumBoundNodes;
  int *TestGlobalNumbers, *AnsatzGlobalNumbers;
  int *TestNumbers, *AnsatzNumbers;
  int *TestBeginIndex, *AnsatzBeginIndex;
  int *KColAux;
  FE3D CurrentElement; 

  if(testspace->GetCollection() != ansatzspace->GetCollection())
  {
    return;
  }

  coll = testspace->GetCollection();
  
  ActiveBound = testspace->GetN_ActiveDegrees();

  AnsatzN_BoundNodeTypes = ansatzspace->GetN_DiffBoundaryNodeTypes();
  TestN_BoundNodeTypes = testspace->GetN_DiffBoundaryNodeTypes();

  //int AnsatzN_Dirichlet = ansatzspace->GetN_Dirichlet();
  //int TestN_Dirichlet = testspace->GetN_Dirichlet();

  //int AnsatzN_Inner = ansatzspace->GetN_Inner();
  //int TestN_Inner = testspace->GetN_Inner();

  AnsatzN_BoundNodes = ansatzspace->GetN_BoundaryNodes();
  TestN_BoundNodes = testspace->GetN_BoundaryNodes();

  AnsatzSumBoundNodes = 0;
  for(i=0;i<AnsatzN_BoundNodeTypes;i++)
    AnsatzSumBoundNodes += AnsatzN_BoundNodes[i];

  TestSumBoundNodes = 0;
  for(i=0;i<TestN_BoundNodeTypes;i++)
    TestSumBoundNodes += TestN_BoundNodes[i];

  AnsatzGlobalNumbers = ansatzspace->GetGlobalNumbers();
  TestGlobalNumbers = testspace->GetGlobalNumbers();

  AnsatzBeginIndex = ansatzspace->GetBeginIndex();
  TestBeginIndex = testspace->GetBeginIndex();

  nRows = testspace->GetN_DegreesOfFreedom();
  nColumns = ansatzspace->GetN_DegreesOfFreedom();

  l = nRows + 1;
  std::vector<int> AuxPtr(l, 0);

  N_ = coll->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    cell = coll->GetCell(i);
    TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
    AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[i];

    CurrentElement = testspace->GetFE3D(i, cell);
    n1 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

    CurrentElement = ansatzspace->GetFE3D(i, cell);
    n2 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

    for(j=0;j<n1;j++)
    {
      k = TestNumbers[j];
      AuxPtr[k] += n2;
    } // endfor j
  } // endfor i

  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  N_=nRows;
  l=AuxPtr[0];
  AuxPtr[0]=0;
  // cout << AuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=AuxPtr[i+1];
    // cout << "i: " << i << "  l: " << l << endl;
    AuxPtr[i+1]=AuxPtr[i]+l;
    l=k;
    // cout << AuxPtr[i+1] << endl;
  }

  // cout << "Upper bound: " << AuxPtr[nRows] << endl;
  
  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[nRows]; // upper bound for number of matrix entries 
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  nEntries=0;
  N_ = coll->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
    AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[i];

    cell = coll->GetCell(i);

    CurrentElement = testspace->GetFE3D(i, cell);
    n1 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

    CurrentElement = ansatzspace->GetFE3D(i, cell);
    n2 = TFEDatabase3D::GetFE3D(CurrentElement)->GetN_DOF();

    // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;
    
    for(j=0;j<n1;j++)
    {
      for(k=0;k<n2;k++)
      {
        m=AnsatzNumbers[k];
        n=TestNumbers[j];

        index=AuxPtr[n];
        l=KColAux[index];
        // check whether this column is already in this row
        while(l!=-1 && l!=m && index<AuxPtr[n+1])
        {
          index++; l=KColAux[index];
        }
        if(l==-1)
        {
          // this is a new column for this row
          KColAux[index]=m;
          nEntries++;
        }
      } // endfor k
    } // endfor j
  } // endfor i
  
/*
  // check
  cout << endl;
  cout << "check" << endl;
  N_=nRows;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=AuxPtr[i+1];
    for(j=AuxPtr[i];j<index && KColAux[j]!=-1;j++) 
      cout << setw(4) << KColAux[j];
    cout << endl;
  }
*/
  
  // cout << "Number of matrix entries: ";
  // cout << nEntries << endl;
  // cout << endl;

  // compress KColAux array to KCol by deleting all -1's
  // build the rows array
  N_=nRows;
  columns.resize(nEntries);
  rows=AuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=AuxPtr[i+1];
    for(j=AuxPtr[i];j<m && KColAux[j]!=-1;j++)
    {
      columns[index]=KColAux[j];
      index++;
    } // endfor j
    rows[i]=oldindex;
    // cout << setw(4) << i << " rows[i]: " << rows[i] << endl;
  } // endfor i
  rows[nRows]=nEntries;
  
  Sort();
  
/*
  // print out the whole matrix structure
  cout << endl;
  N_=nRows;
  for(i=0;i<N_;i++)
  {
    cout << rows[i] << "---" << rows[i+1]-1 << endl;
    cout << "Rows: " << setw(4) << i << ": ";
    end=rows[i+1];
    for(j=RowPtr[i];j<end;j++)
      cout << setw(4) << columns[j];
    cout << endl;
  }
*/
  
  // free KColAux
  delete[] KColAux;

#ifndef _MPI
  this->info();
#endif
  this->Sort();
}
#endif



/* generate the matrix structure, both spaces are 2D */
/* both spaces are defined on different grids */
TStructure::TStructure(const TFESpace2D * testspace, int test_level,
                       const TFESpace2D * ansatzspace, int ansatz_level)
{
  TCollection *coll_coarse, *coll_fine;
  const TBaseCell *cell_coarse, *cell_fine;
  const TBaseCell *cell_child_1, *cell_child_2, *cell_child_3, *cell_child_4, *cell_child_5;
  const TBaseCell *cell_parent, *cell_tmp, *cell_child_6;
  int i,j,k,l,m,n, N_, n1, n2, index, oldindex;
  int TestN_BoundNodeTypes, AnsatzN_BoundNodeTypes;
  int *TestN_BoundNodes, *AnsatzN_BoundNodes;
  int TestSumBoundNodes, AnsatzSumBoundNodes;
  int *TestGlobalNumbers, *AnsatzGlobalNumbers;
  int *TestNumbers, *AnsatzNumbers;
  int *TestBeginIndex, *AnsatzBeginIndex;
  int *KColAux;
  int  N_child,N_child_1, N_child_2, N_child_3, N_child_4, N_child_5;
  int N_child_6;
  int ii, j1, j2, j3, j4, j5, j6, level_diff;
  FE2D CurrentElement; 

  level_diff = abs(test_level-ansatz_level);
  if (level_diff>6)
  {
     ErrThrow("level difference greater than 5 not implemented !!!");
  }

  if (test_level < ansatz_level)
  {
    // get collection of mesh cells
    coll_coarse = testspace->GetCollection();
    coll_fine = ansatzspace->GetCollection();
  }
  else
  {
    // get collection of mesh cells
    coll_fine = testspace->GetCollection();
    coll_coarse = ansatzspace->GetCollection();
  }

 // get information from the spaces
  AnsatzN_BoundNodeTypes = ansatzspace->GetN_DiffBoundaryNodeTypes();
  TestN_BoundNodeTypes = testspace->GetN_DiffBoundaryNodeTypes();

  //int AnsatzN_Dirichlet = ansatzspace->GetN_Dirichlet();
  //int TestN_Dirichlet = testspace->GetN_Dirichlet();

  //int AnsatzN_Inner = ansatzspace->GetN_Inner();
  //int TestN_Inner = testspace->GetN_Inner();

  AnsatzN_BoundNodes = ansatzspace->GetN_BoundaryNodes();
  TestN_BoundNodes = testspace->GetN_BoundaryNodes();

  AnsatzSumBoundNodes = 0;
  for(i=0;i<AnsatzN_BoundNodeTypes;i++)
    AnsatzSumBoundNodes += AnsatzN_BoundNodes[i];

  TestSumBoundNodes = 0;
  for(i=0;i<TestN_BoundNodeTypes;i++)
    TestSumBoundNodes += TestN_BoundNodes[i];

  AnsatzGlobalNumbers = ansatzspace->GetGlobalNumbers();
  TestGlobalNumbers = testspace->GetGlobalNumbers();

  AnsatzBeginIndex = ansatzspace->GetBeginIndex();
  TestBeginIndex = testspace->GetBeginIndex();

  nRows = testspace->GetN_DegreesOfFreedom();
  nColumns = ansatzspace->GetN_DegreesOfFreedom();

  l = nRows + 1;
  std::vector<int> AuxPtr(l, 0);

  // set numeration of the cells
  N_ = coll_coarse->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    cell_coarse = coll_coarse->GetCell(i);
    cell_coarse->SetClipBoard(i);
  }
  N_ = coll_fine->GetN_Cells();
  for(i=0;i<N_;i++)
  {
    cell_fine = coll_fine->GetCell(i);
    cell_fine->SetClipBoard(i);
  }

  if (test_level < ansatz_level)
  {
    // number of mesh cells
    N_ = coll_coarse->GetN_Cells();
    // loop over all mesh cells of the test space
    for(i=0;i<N_;i++)
    {
      // get cell i 
       cell_coarse = coll_coarse->GetCell(i);
       // get global number of d.o.f. of test and ansatz function 
       // which live on this mesh cell
       TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
       // get fe spaces on the mesh cell
       CurrentElement = testspace->GetFE2D(i, cell_coarse);
       n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
       
       // find all children of the coarse mesh cell 
       N_child = cell_coarse->GetN_Children();
       for (j1=0;j1<N_child;j1++)
       {
         // first refinement level 
         cell_child_1 =  cell_coarse->GetChild(j1);
         N_child_1 = cell_child_1->GetN_Children();
         // finest level reached 
         if (N_child_1==0)
         {
           ii = cell_child_1->GetClipBoard();
           AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

           CurrentElement = ansatzspace->GetFE2D(ii, cell_child_1);
           n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
           // update vector which stores information on the length of the rows
           for(j=0;j<n1;j++)
           {
              k = TestNumbers[j];
              AuxPtr[k] += n2;
           } // endfor j
         }
         else
         {
            // second refinement level 
           for (j2=0;j2<N_child_1;j2++)
           {
             cell_child_2 =  cell_child_1->GetChild(j2);
             N_child_2 = cell_child_2->GetN_Children();
             // finest level reached
             if (N_child_2==0)
             {
                ii = cell_child_2->GetClipBoard();
                AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

                CurrentElement = ansatzspace->GetFE2D(ii, cell_child_2);
                n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                // update vector which stores information on the length of the rows
                for(j=0;j<n1;j++)
                {
                   k = TestNumbers[j];
                   AuxPtr[k] += n2;
                } // endfor j
             }
             else
             {
               // third refinement level 
               for (j3=0;j3<N_child_2;j3++)
               {
                 cell_child_3 =  cell_child_2->GetChild(j3);
                 N_child_3 = cell_child_3->GetN_Children();
                 // finest level reached
                 if (N_child_3==0)
                 {
                   ii = cell_child_3->GetClipBoard();
                   AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

                   CurrentElement = ansatzspace->GetFE2D(ii, cell_child_3);
                   n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                   // update vector which stores information on the length of the rows
                   for(j=0;j<n1;j++)
                   {
                      k = TestNumbers[j];
                      AuxPtr[k] += n2;
                   } // endfor j
                 }
                 else
                 {                
                   // fourth refinement level 
                   for (j4=0;j4<N_child_3;j4++)
                   {
                     cell_child_4 =  cell_child_3->GetChild(j4);
                     N_child_4 = cell_child_4->GetN_Children();
                     // finest level reached
                     if (N_child_4==0)
                     {
                       ii = cell_child_4->GetClipBoard();
                       AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

                       CurrentElement = ansatzspace->GetFE2D(ii, cell_child_4);
                       n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                       // update vector which stores information on the length of the rows
                       for(j=0;j<n1;j++)
                       {
                         k = TestNumbers[j];
                         AuxPtr[k] += n2;
                       } // endfor j
                     }
                     else
                     {
                       // fifth refinement level 
                       for (j5=0;j5<N_child_4;j5++)
                       {
                          cell_child_5 =  cell_child_4->GetChild(j5);
                          N_child_5 = cell_child_5->GetN_Children();
                          // finest level reached
                          if (N_child_5==0)
                          {
                             ii = cell_child_5->GetClipBoard();
                             AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

                             CurrentElement = ansatzspace->GetFE2D(ii, cell_child_5);
                             n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                             // update vector which stores information on the length of the rows
                             for(j=0;j<n1;j++)
                             {
                               k = TestNumbers[j];
                               AuxPtr[k] += n2;
                             } // endfor j
                          }
        else 
        {
           // sixth refinement level 
           for (j6=0;j6<N_child_5;j6++)
           {
        cell_child_6 =  cell_child_5->GetChild(j6);
        N_child_6 = cell_child_6->GetN_Children();
        // finest level reached
        if (N_child_6==0)
        {
            ii = cell_child_6->GetClipBoard();
            AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
            
            CurrentElement = ansatzspace->GetFE2D(ii, cell_child_6);
            n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
            // update vector which stores information on the length of the rows
            for(j=0;j<n1;j++)
            {
          k = TestNumbers[j];
          AuxPtr[k] += n2;
            } // endfor j
        }
        else
        {
            ErrThrow("some error");
        }
           } // endfor j6
        }  // endfor else (j5) 
                       } // endfor j5
                     } // endfor else (j4)
                   } // endfor j4
                 } // endfor else (j3)
               } // endfor j3
             } // endfor else (j2)
           } // endfor j2
         } //endfor else (j1)
       } // endfor j1   
    } // endfor i
  }
  else
  {
    // test_level >= ansatz_level
    // number of mesh cells
    N_ = coll_fine->GetN_Cells();
    // loop over all mesh cells of the test space
    for(i=0;i<N_;i++)
    {
      // get cell i 
      cell_fine = coll_fine->GetCell(i);
      // get global number of d.o.f. of test and ansatz function 
      // which live on this mesh cell
      TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
      // get fe spaces on the mesh cell
      CurrentElement = testspace->GetFE2D(i, cell_fine);
      n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
      // find parent cell on coarse grid
      cell_parent = cell_fine;
      for (j1=0;j1<level_diff;j1++)
      {
        cell_tmp = cell_parent->GetParent();
        cell_parent = cell_tmp;
      }
      ii = cell_parent->GetClipBoard();
      AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

      CurrentElement = ansatzspace->GetFE2D(ii, cell_parent);
      n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
      // update vector which stores information on the length of the rows
      for(j=0;j<n1;j++)
      {
        k = TestNumbers[j];
        AuxPtr[k] += n2;
      } // endfor j
    }
  }

  // sum up the array AuxPtr, AuxPtr[i] will now contain the index for
  // KColAux array, the column numbers for row i are in the intervall
  // [ AuxPtr[i],AuxPtr[i+1] )
  N_=nRows;
  l=AuxPtr[0];
  AuxPtr[0]=0;
  // cout << AuxPtr[0] << endl;
  for(i=0;i<N_;i++)
  {
    k=AuxPtr[i+1];
    // cout << "i: " << i << "  l: " << l << endl;
    AuxPtr[i+1]=AuxPtr[i]+l;
    l=k;
    // cout << AuxPtr[i+1] << endl;
  }

  // cout << "Upper bound: " << AuxPtr[nRows] << endl;

  // get memory for KColAux array, initialize it with -1
  l=AuxPtr[nRows]; // upper bound for number of matrix entries 
  KColAux=new int[l];
  memset(KColAux, -1, sizeof(int)*l);

  // compute column indices
  nEntries=0;
  if (test_level < ansatz_level)
  {
    // number of mesh cells
    N_ = coll_coarse->GetN_Cells();
    for(i=0;i<N_;i++)
    {
      cell_coarse = coll_coarse->GetCell(i);
      TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
      CurrentElement = testspace->GetFE2D(i, cell_coarse);
      n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();

      // find all children of the coarse mesh cell 
      N_child = cell_coarse->GetN_Children();
      // first level of refinement
      for (j1=0;j1<N_child;j1++)
      {
        cell_child_1 =  cell_coarse->GetChild(j1);
        N_child_1 = cell_child_1->GetN_Children();
        // finest level reached
        if (N_child_1==0)
        {
          ii = cell_child_1->GetClipBoard();
          AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
          CurrentElement = ansatzspace->GetFE2D(ii, cell_child_1);
          n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
          // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
          for(j=0;j<n1;j++)
          {
            for(k=0;k<n2;k++)
            {
              m=AnsatzNumbers[k];
              n=TestNumbers[j];
              index=AuxPtr[n];
              l=KColAux[index];
              // check whether this column is already in this row
              while(l!=-1 && l!=m && index<AuxPtr[n+1])
              {
                index++; 
                l=KColAux[index];
              }
              if(l==-1)
              {
                 // this is a new column for this row
                KColAux[index]=m;
                nEntries++;
              }
            } // endfor k
          } // endfor j
        }
        else
        {
          // second level of refinement
          for (j2=0;j2<N_child_1;j2++)
          {
            cell_child_2 =  cell_child_1->GetChild(j2);
            N_child_2 = cell_child_2->GetN_Children();
            // finest level reached
            if (N_child_2==0)
            {
              ii = cell_child_2->GetClipBoard();
              AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
              CurrentElement = ansatzspace->GetFE2D(ii, cell_child_2);
              n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
              // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
              for(j=0;j<n1;j++)
              {
                for(k=0;k<n2;k++)
                {
                  m=AnsatzNumbers[k];
                  n=TestNumbers[j];
                  index=AuxPtr[n];
                  l=KColAux[index];
                  // check whether this column is already in this row
                  while(l!=-1 && l!=m && index<AuxPtr[n+1])
                  {
                    index++; 
                    l=KColAux[index];
                  }
                  if(l==-1)
                  {
                    // this is a new column for this row
                    KColAux[index]=m;
                    nEntries++;
                  }
                } // endfor k
              } // endfor j
            }
            else
            {
              // third level of refinement
              for (j3=0;j3<N_child_2;j3++)
              {
                cell_child_3 =  cell_child_2->GetChild(j3);
                N_child_3 = cell_child_3->GetN_Children();
                // finest level reached
                if (N_child_3==0)
                {
                  ii = cell_child_3->GetClipBoard();
                  AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
                  CurrentElement = ansatzspace->GetFE2D(ii, cell_child_3);
                  n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                  // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
                  for(j=0;j<n1;j++)
                  {
                    for(k=0;k<n2;k++)
                    {
                      m=AnsatzNumbers[k];
                      n=TestNumbers[j];
                      index=AuxPtr[n];
                      l=KColAux[index];
                      // check whether this column is already in this row
                      while(l!=-1 && l!=m && index<AuxPtr[n+1])
                      {
                        index++; 
                        l=KColAux[index];
                      }
                      if(l==-1)
                      {
                        // this is a new column for this row
                        KColAux[index]=m;
                        nEntries++;
                      }
                    } // endfor k
                  } // endfor j
                }
                else
                {
                  // fourth level of refinement
                  for (j4=0;j4<N_child_3;j4++)
                  {
                    cell_child_4 =  cell_child_3->GetChild(j4);
                    N_child_4 = cell_child_4->GetN_Children();
                    // finest level reached
                    if (N_child_4==0)
                    {
                      ii = cell_child_4->GetClipBoard();
                      AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
                      CurrentElement = ansatzspace->GetFE2D(ii, cell_child_4);
                      n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                      // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
                      for(j=0;j<n1;j++)
                      {
                        for(k=0;k<n2;k++)
                        {
                          m=AnsatzNumbers[k];
                          n=TestNumbers[j];
                          index=AuxPtr[n];
                          l=KColAux[index];
                          // check whether this column is already in this row
                          while(l!=-1 && l!=m && index<AuxPtr[n+1])
                          {
                            index++; 
                            l=KColAux[index];
                          }
                          if(l==-1)
                          {
                            // this is a new column for this row
                            KColAux[index]=m;
                            nEntries++;
                          }
                        } // endfor k
                      } // endfor j
                    }
                    else
                    {
                      // fifth level of refinement
                      for (j5=0;j5<N_child_4;j5++)
                      {
                        cell_child_5 =  cell_child_4->GetChild(j5);
                        N_child_5 = cell_child_5->GetN_Children();
                        // finest level reached
                        if (N_child_5==0)
                        {
                          ii = cell_child_5->GetClipBoard();
                          AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
                          CurrentElement = ansatzspace->GetFE2D(ii, cell_child_5);
                          n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
                          // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
                          for(j=0;j<n1;j++)
                          {
                            for(k=0;k<n2;k++)
                            {
                              m=AnsatzNumbers[k];
                              n=TestNumbers[j];
                              index=AuxPtr[n];
                              l=KColAux[index];
                              // check whether this column is already in this row
                              while(l!=-1 && l!=m && index<AuxPtr[n+1])
                              {
                                index++; 
                                l=KColAux[index];
                              }
                              if(l==-1)
                              {
                                // this is a new column for this row
                                KColAux[index]=m;
                                nEntries++;
                              }
                            } // endfor k
                          } // endfor j
                        }
                        else
                        {
                          ErrThrow("4711");
                        } // end else (j5)
                      } // endfor j5                                        
                    } // end else (j4)
                  } // endfor j4                  
                } // end else (j3)
              } // endfor j3
            } // end else (j2)
          } // endfor j2
        } // end else (j1)
      } // endfor j1       
    } // endfor i
  } // endfor test_level < ansatz_level
  else
  {
    // test_level >= ansatz_level
    // number of mesh cells
    N_ = coll_fine->GetN_Cells();
    // loop over all mesh cells of the test space
    for(i=0;i<N_;i++)
    {
      // get cell i 
      cell_fine = coll_fine->GetCell(i);
      // get global number of d.o.f. of test and ansatz function 
      // which live on this mesh cell
      TestNumbers = TestGlobalNumbers + TestBeginIndex[i];
      // get fe spaces on the mesh cell
      CurrentElement = testspace->GetFE2D(i, cell_fine);
      n1 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
      // find parent cell on coarse grid
      cell_parent = cell_fine;
      for (j1=0;j1<level_diff;j1++)
      {
        cell_tmp = cell_parent->GetParent();
        cell_parent = cell_tmp;
      }
      ii = cell_parent->GetClipBoard();
      AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];

      ii = cell_parent->GetClipBoard();
      AnsatzNumbers = AnsatzGlobalNumbers + AnsatzBeginIndex[ii];
      CurrentElement = ansatzspace->GetFE2D(ii, cell_parent);
      n2 = TFEDatabase2D::GetFE2D(CurrentElement)->GetN_DOF();
      // cout << "i: " << i << " n1: " << n1 << " n2: " << n2 << endl;    
      for(j=0;j<n1;j++)
      {
        for(k=0;k<n2;k++)
        {
          m=AnsatzNumbers[k];
          n=TestNumbers[j];
          index=AuxPtr[n];
          l=KColAux[index];
          // check whether this column is already in this row
          while(l!=-1 && l!=m && index<AuxPtr[n+1])
          {
            index++; 
            l=KColAux[index];
          }
          if(l==-1)
          {
                                // this is a new column for this row
            KColAux[index]=m;
            nEntries++;
          }
        } // endfor k
      } // endfor j
    } // endfor i
  } // end else

/*
  // check
  cout << endl;
  cout << "check" << endl;
  N_=nRows;
  for(i=0;i<N_;i++)
  {
    cout << "Row: " << setw(4) << i << ": ";
    index=AuxPtr[i+1];
    for(j=AuxPtr[i];j<index && KColAux[j]!=-1;j++) 
      cout << setw(4) << KColAux[j];
    cout << endl;
  }
*/
#ifdef _MPI
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0)
#endif 
  { 
  cout << "Number of matrix entries: ";
  cout << nEntries << endl;
  cout << endl;
  }
  // compress KColAux array to KCol by deleting all -1's
  // build the rows array
  N_=nRows;
  columns.resize(nEntries);
  rows=AuxPtr;

  index=0;
  for(i=0;i<N_;i++)
  {
    oldindex=index;
    m=AuxPtr[i+1];
    for(j=AuxPtr[i];j<m && KColAux[j]!=-1;j++)
    {
      columns[index]=KColAux[j];
      index++;
    } // endfor j
    rows[i]=oldindex;
    // cout << setw(4) << irowsrows[i]: " << rows[i] << endl;
  } // endfor i
  rows[nRows]=nEntries;
  
/*
  // print out the whole matrix structure
  OutPut(endl);
  N_=nRows;
  for(i=0;i<N_;i++)
  {
     OutPut( rows[i] << "---" << rows[i+1]-1 << endl);
     OutPut("Rows: " << setw(4) << i << ": ");
     end=rows[i+1];
     for(j=rows[i];j<end;j++)
        OutPut( setw(4) << KCol[j]);
     OutPut(endl);
  }
*/
  
  // free KColAux
  delete [] KColAux;
  
#ifndef _MPI
  this->info();
#endif
  this->Sort();
}


TStructure::TStructure() 
 : nRows(0), nColumns(0), nEntries(0), columns(), rows(),
   ActiveBound(0), ColOrder(0), nHangingEntries(0), hangingColums(),
   HangingRows()
{
}

TStructure::TStructure(int n, int nEntries, int *col_ptr, int *row_ptr)
 : TStructure(n, n, nEntries, col_ptr, row_ptr)
{
}

/* generate the matrix structure, all arrays are already defined */
TStructure::TStructure(int nRows, int nCols, int nEntries, int *col_ptr, 
                       int *row_ptr)
 : TStructure(nRows, nCols, nRows, nEntries, col_ptr, row_ptr)
{
  
}

TStructure::TStructure(int nRows, int nCols, int nActive, int nEntries,
                       int *col_ptr, int *row_ptr)
 : nRows(nRows), nColumns(nCols), nEntries(nEntries), columns(nEntries, 0), 
   rows(nRows+1), ActiveBound(nActive), ColOrder(0), nHangingEntries(0),
   hangingColums(), HangingRows()
{
  std::copy(col_ptr, col_ptr + this->GetN_Entries(), this->columns.begin());
  std::copy(row_ptr, row_ptr + this->GetN_Rows() + 1, this->rows.begin());
  this->Sort();
}

TStructure::TStructure(int nRows, int nCols)
 : nRows(nRows), nColumns(nCols), nEntries(0), columns(), 
   rows(nRows+1, 0), ActiveBound(nRows), ColOrder(0), nHangingEntries(0),
   hangingColums(), HangingRows()
{
}


TStructure::TStructure(int n) : TStructure(n, n)
{
}

/* sort column numbers: diag is first element, other numbers are increasing */
void TStructure::SortDiagFirst()
{
  int i,j,k;
  int end, begin;

  end = 0;
  for(i=0;i<ActiveBound;i++)
  {
    begin = end;
    end = rows[i+1];
    k = columns[begin];
    for(j=rows[i];j<end;j++)
    {
      if(columns[j] == i)
      {
        // diag entry
        columns[begin] = i;
        columns[j] = k;
        break;
      } // endif
    } // endfor j
    SortRow(&columns[0]+begin+1, &columns[0]+end);
  } // endfor i
  
  this->ColOrder = 2;
}

/* sort one row [BeginPtr, AfterEndPtr) */
void TStructure::SortRow(int *BeginPtr, int *AfterEndPtr)
{
  int *IPtr, *JPtr, T;

  for(IPtr=BeginPtr;IPtr<AfterEndPtr;IPtr++)
  {
    for(JPtr=IPtr+1;JPtr<AfterEndPtr;JPtr++)
    {
      if( *IPtr > *JPtr )
      {
        T = *IPtr;
        *IPtr = *JPtr;
        *JPtr = T;
      }
    } // endfor JPtr
  } // endfor IPtr
  
  this->ColOrder = 1;
}

/* sort numbers within each row */
void TStructure::Sort()
{
  int end, begin;

  end = 0;
  for(int i=0; i<nRows; i++)
  {
    begin = end;
    end = rows[i+1];
    SortRow(&columns[0]+begin, &columns[0]+end);
  } // endfor i
}

size_t TStructure::get_n_entries_in_row(size_t row_index) const
{
  if((int) row_index > nRows)
  {
    ErrThrow("unable to return the number of entries in row ", row_index,
             ". There are only ", this->nRows, " rows.");
  }
  auto n = this->rows[row_index+1]-this->rows[row_index];
  if(n < 0)
    ErrThrow("it seems this structure is in an invalid state");
  return (size_t)n;
}


void TStructure::reset_n_entries()
{
  //throw if number of rows changed
  if ((int)this->rows.size() - 1 != this->nRows)
    ErrThrow("TStructure: nRows != rows.size() - 1 ",nRows, " != ", rows.size() + 1);

  //throw if there are hanging entries
  if(this->nHangingEntries != 0)
    ErrThrow("TStructure has hanging node entries. "
            "TStructure::reset_n_entries will not yet reset their "
            "number and arrays correctly.");

  //reset number of entries to last entry of row ptr
  this->nEntries = this->rows.back();

  //columns array must be resized  to nEntries - everything behind that is erased
  this->columns.resize(this->nEntries);
}

int TStructure::index_of_entry(const int i, const int j) const
{
  if(i < 0 || i >= this->GetN_Rows())
  {
    ErrThrow("row index is too large or too small");
  }
  if(j < 0 || j > this->GetN_Columns())
  {

//      cout << "Hier" << endl;
//      cout << j << endl;
//      cout << this->GetN_Columns() << endl ;
//      cout << this->GetN_Rows() << endl ;
      
    ErrThrow("column index is too large or too small");
  }
  
  for (int m=rows[i];m < rows[i+1]; m++) 
  {
    if (columns[m]== j) 
    {
      // index found in sparsity pattern
      return m;
    }
  }
  // index not in the sparsity pattern
  return -1;
}

unsigned int TStructure::getNActiveEntries() const
{
  return this->rows[this->ActiveBound];
}

/* return a new structure for a transposed matrix */
std::shared_ptr<TStructure> TStructure::GetTransposed() const
{
  if(nHangingEntries!=0)
  {
    ErrThrow("TStructure::GetTransposed(): Hanging entries not supported!");
  }
  // variables for transposed structure:
  int nRowsT = nColumns;
  int nColsT = nRows;
  // number of entries does not change
  int *rowsT = new int[nRowsT+1];  memset(rowsT, 0, (nRowsT+1)*SizeOfInt);
  int *colsT = new int[nEntries]; memset(colsT, 0, nEntries *SizeOfInt);
  
  int *ColB_count = new int[nColumns]; 
  memset(ColB_count, 0, nColumns*SizeOfInt);
  
  // count number of entries per column (will be number of entries in each row)
  for(int i=0;i<rows[nRows];i++)
    rowsT[columns[i]]++;
  // change to increasing numbering as in rows
  for(int i=0,k=0;i<=nRowsT;i++)
  {
    int j = rowsT[i];
    rowsT[i] = k;
    k += j;
  }
  
  // fill 'colsT'
  // loop over (non-transposed) rows
  for(int i=0; i<nRows; i++)
  {
    // loop over all entries in this row
    for(int j=rows[i]; j<rows[i+1]; j++)
    {
      int col = columns[j]; // (non-transposed) column = transposed row
      int l  = rowsT[col];
      int offset = ColB_count[col];
      colsT[l+offset] = i;
      ColB_count[col]++;
    }
  }
  delete []ColB_count;
  
  std::shared_ptr<TStructure> structureT(new TStructure(nRowsT, nColsT, nEntries, colsT, rowsT)); 

  //the new TStructure makes its own copy of these
  delete[] colsT;
  delete[] rowsT;

  return structureT;
}

std::shared_ptr<TStructure> get_product_structure(
    TStructure const & strucA, TStructure const & strucB)
{
  const int n_A_rows = strucA.GetN_Rows();   // = n_C_rows
  const int n_A_cols = strucA.GetN_Columns();
  const int n_B_rows = strucB.GetN_Rows();
  const int n_B_cols = strucB.GetN_Columns();   // = n_C_cols
  
  if(n_A_cols != n_B_rows)
  {
    ErrThrow("dimension mismatch during matrix-matrix multiplication");
  }
  const int * const a_rows = strucA.GetRowPtr();
  const int * const a_cols = strucA.GetKCol();
  
  // everything needed to call the constructor of TStructure later on:
  int n_c_entries = 0; // number of entries in product structure C
  int * c_rows = new int[n_A_rows + 1]; // row pointer
  memset(c_rows,0.0, (n_A_rows + 1)*SizeOfInt);
  int * c_cols; // columns pointer
  std::vector<std::vector<int> > dofs(n_A_rows);
  
  // loop over all rows of C
  for(int row = 0; row < n_A_rows; row++)
  {
    // loop over all columns of C
    for(int col = 0; col < n_B_cols; col++)
    {
      // check whether 'this row of A' x 'this column of B' would give an entry
      int n_a_entries_in_row = a_rows[row+1] - a_rows[row];
      // loop over all entries in this row in A
      for(int i = 0; i < n_a_entries_in_row; i++)
      {
        if(strucB.index_of_entry(a_cols[i+a_rows[row]], col) != -1)
        {
          dofs[row].push_back(col);
          break;
        }
      }
    }
    n_c_entries += dofs[row].size();
    c_rows[row+1] = n_c_entries; // c_rows[0] has been set to 0 already
  }
  
  // now fill the array c_cols
  c_cols = new int[n_c_entries];
  // loop over all rows of C
  for(int row = 0; row < n_A_rows; row++)
  {
    // loop over all columns of C
    int nEntries_in_this_row = c_rows[row+1] - c_rows[row];
    for(int col = 0; col < nEntries_in_this_row; col++)
    {
      c_cols[c_rows[row] + col] = dofs[row].at(col);
    }
  }  
  return std::make_shared<TStructure>(n_A_rows, n_B_cols, n_c_entries, c_cols,
                                      c_rows);
}


TStructure* TStructure::get_structure_of_product_with_transpose_from_right() 
const
{
  if(this->ColOrder != 1)
    ErrThrow("TStructure::get_structure_of_product_with_transpose_from_right ",
             "only works for structures which are ordered. ", this->ColOrder);
  
  int nProductEntries = 0; // number of entries in product structure
  int nProductRows = nRows;
  int nProductColumns = nRows;

  std::vector<int> productRowPtr(nProductRows + 1, 0); // row pointer

  // this is a temporary storing structure "gridPlaces"
  std::vector<std::vector<int> > gridPlaces(nRows); 
  // gridPlaces[row] stores, in which columns in this row there are non-zero 
  // entries

  // loop over all rows of the product
  for(int row = 0; row < nProductRows; row++)
  {
    // loop over all columns of the product
    for(int col = 0; col < nProductColumns; col++)
    {
      // check whether 'this row of A' x 'this column of A^T' would give an 
      // entry
      // this boils down to checking 'row1 of A' x 'row2 of A' with row1 = 
      // row, row2 = col
      int row1 = row;
      int row2 = col; //two definitions to fix ideas

      // work on two segments of column array
      const int* row1ColBegin = &columns[rows[row1]];
      const int row1ColSize = rows[row1+1] - rows[row1];
      const int* row2ColBegin = &columns[rows[row2]];
      const int row2ColSize = rows[row2+1] - rows[row2];

      // find out if the two arryas contain a common value
      // exploit the fact that both are sorted
      size_t index1 = 0;
      size_t index2 = 0;
      while ((int)index1 < row1ColSize && (int)index2 < row2ColSize)
      {
        if (row1ColBegin[index1] > row2ColBegin[index2])
        {
          index2++;
        }
        else if (row2ColBegin[index2] > row1ColBegin[index1])
        {
          index1++;
        }
        else
        {
          //we found a pair of indices with equal entry in KCol
          gridPlaces[row].push_back(col);
          break;
        }
      }

    } //end for columns of the product

    // update the total number of entries in the product
    nProductEntries += gridPlaces[row].size();
    // productRowPtr[0] has been set to 0 already
    productRowPtr[row+1] = nProductEntries;
  }//end for rows of the product

  // FROM HERE IT'S ONLY TRANSFORMING THE TEMPORARY DATA STRUCTURE TO THE 
  // REQUIRED ONE
  // now fill the array productColumnsPtr
   //initialise the product's columns pointer
  std::vector<int> productColumnPtr(nProductEntries, 0);

  // loop over all rows of C
  for(int row = 0; row < nProductRows; row++)
  {
    // loop over all columns of C
    int n_entries_in_this_row = productRowPtr[row+1] - productRowPtr[row];
    for(int col = 0; col < n_entries_in_this_row; col++)
    {
      productColumnPtr[productRowPtr[row] + col] = gridPlaces[row].at(col);
    }
  }
  // hand over a pointer to the product's structure.
  return new TStructure(nProductRows, nProductColumns, nProductEntries, 
                        &productColumnPtr[0], &productRowPtr[0]);
}


TStructure* TStructure::get_structure_of_product_with_transpose_from_right(
  const TStructure& B) const
{
  if(this->ColOrder != 1 || B.GetColOrder() != 1)
    ErrThrow("TStructure::get_structure_of_product_with_transpose_from_right ",
             "only works for structures which are ordered. ", this->ColOrder,
             "  ", B.GetColOrder());
  if(this->nColumns != B.GetN_Rows() || B.GetN_Columns() != this->nColumns)
    ErrThrow("dimension mismatch, inner matrix has wrong dimensions, ",
             B.GetN_Rows(), "  ", B.GetN_Columns());
  
  int nProductEntries = 0; // number of entries in product structure
  int nProductRows = this->nRows;
  int nProductColumns = this->nRows;
  
  // lambda funcion to check if the product of 'B*A^T' has an entry (i,j)
  auto check_entry = [this, B](int i, int j) -> bool
    {
      // work on two segments of column array
      const int* row1ColBegin = &B.GetKCol()[B.GetRowPtr()[i]];
      const int row1ColSize = B.GetRowPtr()[i+1] - B.GetRowPtr()[i];
      const int* row2ColBegin = &this->GetKCol()[this->GetRowPtr()[j]];
      const int row2ColSize = this->GetRowPtr()[j+1] 
                              - this->GetRowPtr()[j];

      // find out if the two arryas contain a common value
      // exploit the fact that both are sorted
      size_t index1 = 0;
      size_t index2 = 0;
      while ((int)index1 < row1ColSize && (int)index2 < row2ColSize)
      {
        if (row1ColBegin[index1] > row2ColBegin[index2])
        {
          index2++;
        }
        else if (row2ColBegin[index2] > row1ColBegin[index1])
        {
          index1++;
        }
        else
        {
          return true;
        }
      }
      return false;
   };
  
  std::vector<int> productRowPtr(nProductRows + 1, 0); // row pointer
  // this is a temporary storing structure "gridPlaces"
  std::vector<std::vector<int>> gridPlaces(this->nRows);
  // gridPlaces[row] stores, in which columns in this row there are non-zero 
  // entries
  // loop over all rows of the product
  for(int row = 0; row < nProductRows; row++)
  {
    // loop over all columns of the product
    for(int col = 0; col < nProductColumns; col++)
    {
      //Output::print("row ", row, "  column ", col);
      // check whether 'this row of A' x 'this column of B*A^T' would give an 
      // entry
      // this boils down to checking 'row1 of A' x 'row2 of B*A' with row1 = 
      // row, row2 = col
      int row1 = row;
      //int row2 = col; //two definitions to fix ideas

      // work on two segments of column array
      const int* row1ColBegin = &columns[rows[row1]];
      const int row1ColSize = rows[row1+1] - rows[row1];
      //const int* row2ColBegin = &columns[rows[row2]];
      //const int row2ColSize = rows[row2+1] - rows[row2];

      // find out if the two arryas contain a common value
      // exploit the fact that both are sorted
      size_t index1 = 0;
      size_t index2 = 0;
      while ((int)index1 < row1ColSize && (int)index2 < this->nColumns)
      {
        //Output::print("  index1 ", index1, "  index2 ", index2);
        if (row1ColBegin[index1] > (int)index2)
        {
          index2++;
        }
        else if ((int)index2 > row1ColBegin[index1])
        {
          index1++;
          // I think we should never be here (Ulrich)
        }
        else
        {
          // now we need to check if 'B*A^T' has an entry at (index2,col)
          if(check_entry(index2, col))
          {
            //we found a pair of indices with equal entry in KCol
            gridPlaces[row].push_back(col);
            break;
          }
          ++index1;
          ++index2;
        }
      }

    } //end for columns of the product

    // update the total number of entries in the product
    nProductEntries += gridPlaces[row].size();
    // productRowPtr[0] has been set to 0 already
    productRowPtr[row+1] = nProductEntries;
  }//end for rows of the product

  // FROM HERE IT'S ONLY TRANSFORMING THE TEMPORARY DATA STRUCTURE TO THE 
  // REQUIRED ONE
  // now fill the array productColumnsPtr
   //initialise the product's columns pointer
  std::vector<int> productColumnPtr(nProductEntries, 0);

  // loop over all rows of C
  for(int row = 0; row < nProductRows; row++)
  {
    // loop over all columns of C
    int n_entries_in_this_row = productRowPtr[row+1] - productRowPtr[row];
    for(int col = 0; col < n_entries_in_this_row; col++)
    {
      productColumnPtr[productRowPtr[row] + col] = gridPlaces[row].at(col);
    }
  }
  // hand over a pointer to the product's structure.
  return new TStructure(nProductRows, nProductColumns, nProductEntries, 
                        &productColumnPtr[0], &productRowPtr[0]);
}


void TStructure::fortran_shift()
{
  // the first entry in the this->rows is always zero (when indices start with 
  // zero), so it is used as an indicator if the entries in this->rows and 
  // this->columns are in fortran or c++ style.
  int s = 1; // forward shift
  if(this->is_fortran_shifted())
  {
    // backward shift (indices start with 1)
    s = -1;
  }
  // else // forward shift (indices start with 0)
  std::for_each(this->rows.begin(), this->rows.end(),       [s](int& i){i+=s;});
  std::for_each(this->columns.begin(), this->columns.end(), [s](int& i){i+=s;});
}

bool TStructure::is_fortran_shifted() const
{
  if(this->rows[0] == 1)
  {
    return true;
  }
  else if(this->rows[0] == 0)
  {
    return false;
  }
  else
    ErrThrow("broken matrix structure, first index in row vector is neither ",
             "0 nor 1, but ", this->rows[0]);
}

void TStructure::info() const
{
  Output::info<3>("TStructure","Information on the stored matrix structure");
  Output::dash<3>("Number of rows: ", nRows);
  Output::dash<3>("Number of columns: ", nColumns);
  Output::dash<3>("Number of matrix entries: ", nEntries);
}

void TStructure::draw(std::string filename) const
{
  std::ofstream out_stream(filename);
  if(!out_stream)
  {
    ErrMsg("cannot open '" << filename << "' for output");
    return;
  }
  
  double scale = 3;
  int BX = (int) (scale * this->nColumns); // width of picture
  int BY = (int) (scale * this->nRows);    // height of picture
  int offset = 2 * scale; // picture is a bit away from picture boundary
  
  out_stream << "%!PS-Adobe-3.0\n";
  out_stream << "%%Creator: ParMooN (Ulrich Wilbrandt)\n";
  out_stream << "%%DocumentFonts: Helvetica\n";
  out_stream << "%%BoundingBox: 0 0 " << 2*offset + BX << " " << 2*offset + BY;
  out_stream << endl;
  out_stream << "%%Pages: 1\n";
  out_stream << "%%EndComments\n";
  out_stream << "%%EndProlog\n";
  out_stream << "%%Page: 1 1\n";
  out_stream << "%% n_rows " << this->nRows << ",  n_columns " 
             << this->nColumns << ",  n_entries " << this->nEntries << endl;
  out_stream << "/M {" << scale-1 << " " << scale-1 << " rectfill} def\n";
  for(int row = 0; row < this->nRows; ++row)
  {
    for(int col = this->rows[row]; col < this->rows[row+1]; ++col)
    {
      out_stream << (offset + scale * this->columns[col]) << " " 
                 << (offset + scale * (this->nRows - row)) << " M\n";
    }
  }
  out_stream << "stroke" << endl;
  out_stream << "showpage" << endl;
  out_stream << "%%Trailer" << endl;
  out_stream << "%%Pages: 1" << endl;
  out_stream.close();
  Output::print<2>("postscript picture of matrix structure drawn in file ",
                   filename);
}


bool operator==(const TStructure &lhs, const TStructure &rhs)
{
  return lhs.nRows == rhs.nRows 
      && lhs.nColumns == rhs.nColumns
      && lhs.nEntries == rhs.nEntries
      && lhs.nHangingEntries == rhs.nHangingEntries;
}

bool operator!=(const TStructure &lhs, const TStructure &rhs)
{
  return !(rhs == lhs);
}

