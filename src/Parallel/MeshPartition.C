// =======================================================================
// @(#)MeshPartition.C
// 
// Purpose:     partition the domain into "npart" parts for parallel computing
// 
// Author:      Sashikumaar Ganesan
// History:      start of implementation  07/09/09 (Sashikumaar Ganesan)
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <Database.h>
#include <Domain.h>
#include <Output2D.h>
#include <string.h>
#include <sstream>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <malloc.h>
#include <Joint.h>
#include <SubDomainJoint.h>
#include <SubDomainHaloJoint.h>
#include <Vertex.h>
#include <BaseCell.h>
#include <Quadrangle.h>
#include <MacroCell.h>
#include <Edge.h>

#ifdef __2D__
  #include <FEDatabase2D.h>
#endif

#ifdef _MPI
extern "C"
{
  #include <metis.h>
  #include <parmetis.h>
}
#endif


static void Sort(int *Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, k, *rr, len, s;
  int Mid, Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr= new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;
    
  } while (i<r);

  delete [] rr;

}
 
static void Sort(TVertex **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, k, *rr, len, s;
  TVertex *Mid, *Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr= new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete [] rr;

}



static void Sort(TEdge **Array, int length)
{
  int n=0, l=0, r=length-1, m, max=0;
  int i, j, k, *rr, len, s;
  TEdge *Mid, *Temp;
  double lend = length;

  len=(int)(4*log(lend)/log(2.0) +2);
  

  
  rr= new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
	if(max<n) max=n;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

    
//   printf("%d :Time taken for Domain Decomposition is %d\n", max, len);
  
  delete [] rr;

}



static void Sort(TBaseCell **Array, int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, k, *rr, len, s;
  TBaseCell *Mid, *Temp;
  double lend = length;

  len=(int)(2*log(lend)/log((double) 2.0)+2);
  rr= new int[len];

  do
  {
    do
    {
      i=l;
      j=r;

      m=(l+r)/2;
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) i++;

        while(Array[j] < Mid) j--;

        if (i<=j)
        {
          Temp=Array[i];
          Array[i]=Array[j];
          Array[j]=Temp;
          i++; j--;
        }
      } while (i<=j);

      if (l<j)
      {
        rr[++n]=r;
        r=j;
      }
    } while (l<j);

    if (n>0) r=rr[n--];

    if (i<r) l=i;

  } while (i<r);

  delete [] rr;

}




static int GetIndex(TVertex **Array, int Length, TVertex *Element)
{
  int l=0, r=Length, m=(r+l)/2;
  TVertex *Mid;

  Mid=Array[m];
  while(Mid!=Element)
  {
    if(Mid>Element)
    {
      l=m;
    }
    else
    {
      r=m;
    }
    m=(r+l)/2;
    Mid=Array[m];
  }

  return m;
}

#ifdef _MPI

#ifdef __2D__


#else // 3D
void Partition_Mesh3D_old(MPI_Comm comm, TDomain *Domain, int &MaxRankPerV)
{
 idxtype *Cell_Rank, *Vert_Rank; 
 int i, j, k, m, M, n, rank, size, N_Cells, N, ID, Neib_ID, out_rank= TDatabase::ParamDB->Par_P0;
 int *GlobalCellIndex, N_RootVertices, N_VertInCell, N_JointsInCell, N_AllLocVert;
 int *VertexNumbers, *PointNeighb, maptype, MaxCpV;
 int MaxLen, N_Edges, N_FaceEdges, ii, jj, GblCellNr;
 int N_EdgeDel=0, N_VertexDel=0, N_CellDel=0;
 const int *TmpFV, *TmpLen, *TmpVE, *EdgeVertex, *NeibEdgeVertex;
 
 int VertexCellNo, VertexcellID, N_CellsInThisVert, *VertNeibRanks;
 int N_LocalCells, N_OwnCells, N_OwnIncidentCells, N_NeibIncidentCells;
 int a, b, N_SubDomInThisVert, *HaloCellIndex;
 int M1, M2, N1, N2, N_CellsIn_a, N_CellsIn_b;
 int N_CrossNeibs, *CrossNeibsRank, *HaloCellGlobalNo, *HaloCellLocVertNo;
 int N_SubDomIn_a, N_SubDomIn_b, *Temp, kk, GlobCellNo, test_b, EdgeCellID;

 bool UPDATE, UPDATE_1;

 TVertex *Vert_a, **VertexDel, *Last;
 TCollection *coll;
 TBaseCell *cell, *neib_cell, *Vertexcell, **SubDomainCells, *cell_a, *cell_b, **IncCells; 
 TBaseCell **CellDel, *LastCell;
 TJoint *Joint, *NewJoint;
 TEdge *edge, **EdgeDel, *LastEdge;
 TShapeDesc *ShapeDesc, *NeibShapeDesc;
 TVertex *CurrVert, *NeibVert_a; 
 
 MPI_Status status, status1, status2, status3;
 MPI_Request request, request1;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // root will not take part in computations
  MaxCpV = 0;
  MaxRankPerV = -1;

//   if(rank==0)
//    MaxRankPerV = 1;

  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  cell = coll->GetCell(0);
  N_VertInCell = cell->GetN_Vertices();
  N_JointsInCell = cell->GetN_Joints();
  N_AllLocVert = N_VertInCell*N_Cells;
  VertexNumbers= new int[N_AllLocVert];
  Vert_Rank = new idxtype[N_AllLocVert];  
  Cell_Rank = new idxtype[N_Cells];  

  VertexDel = new TVertex*[N_AllLocVert];
  EdgeDel = new TEdge*[N_Cells*12]; // max 12 in Hexa
  CellDel = new TBaseCell*[N_Cells];

  if(rank==0)
   printf("Number of  ranks: %d\n",  size);

  if(N_Cells<size)
   {
     printf("Number of cells less than number of processors !!!, %d\n",  N_Cells);
     MPI_Finalize();
     exit(0);
   }

 if(size==1)
   {
    cout <<  "Total number of process should be grater than 1 (or 2 if root is not involved in computation), but you have " << size<<endl;
     MPI_Finalize();
     exit(0);
   }

 //check, all cell have same number of vertices  
  for(i=1;i<N_Cells;i++)
   if(((coll->GetCell(i))->GetN_Vertices() != N_VertInCell) || ( (coll->GetCell(i))->GetN_Joints() != N_JointsInCell) )
    {
     cout << N_JointsInCell << "Mesh partition for heterogeneous cells are not yet implemented " << N_VertInCell<<endl;
     MPI_Finalize();
     exit(0);
    }


// partition the mesh in the root and send info to all processors
 if(rank==0)
  {
   cout <<  "Total Cells " << N_Cells<<endl;

   int  m1, edgecut=0, *NumberVertex;
   int etype, numflag=0; // c-style numbering (array starting with 0)
   int type=TDatabase::ParamDB->Par_P2;

   idxtype *MetisVertexNumbers; 

   double t1, t2;

   TVertex **Vertices;    

   
    MetisVertexNumbers= new idxtype[N_AllLocVert];
    NumberVertex=new int[N_AllLocVert];
    Vertices=new TVertex*[N_AllLocVert];
    
    N = 0;    
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
     
      for(j=0;j<N_VertInCell;j++)
       {
        Vertices[N]=cell->GetVertex(j);
        N++;
       }
     } // for i=0

    // sort the Vertices array
    Sort(Vertices, N);

    Last=NULL;
    N_RootVertices=-1;
    for(i=0;i<N_AllLocVert;i++)
     {
      if((CurrVert=Vertices[i])!=Last)
       {
        N_RootVertices++;
        Last=CurrVert;
       }
      NumberVertex[i]=N_RootVertices;
     }
    N_RootVertices++;


    m=0;
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);

      for(j=0;j<N_VertInCell;j++)
       {
        CurrVert=cell->GetVertex(j);
        N=GetIndex(Vertices, N_AllLocVert, CurrVert);
        VertexNumbers[m]=NumberVertex[N];
        m++;
       } // endfor j
     } //endfor i

   memcpy(MetisVertexNumbers, VertexNumbers, N_AllLocVert*SizeOfInt);

   cout << "Total Vertices " <<N_RootVertices<<endl;
   //cout << "Vertices N " <<N_AllLocVert<<endl; 

  for(i=1;i<size;i++)
   MPI_Send(&N_RootVertices, 1, MPI_INT, i, 75, comm);
  
  //   for (i=0; i<N_LocVertices; i++) 
 //    cout<< i <<"  " << "elements " << MetisVertexNumbers[i] <<endl;

    if(N_VertInCell==4)
     { etype = 2; } // tetrahedral
    else if(N_VertInCell==8)
     { etype = 3; } // Hexahedral
    else
     {
      cout<<" Error only  Tetra or Hexa mesh can be partitioned !!" <<endl;
      MPI_Finalize();
      exit(0);
     }
 
  for(i=1;i<size;i++)
    MPI_Send(VertexNumbers, N_AllLocVert, MPI_INT, i, 80, comm);  
  
  t1 = MPI_Wtime();
   if(type == 0)
    METIS_PartMeshNodal(&N_Cells, &N_RootVertices, MetisVertexNumbers, &etype, &numflag,
                       &size, &edgecut, Cell_Rank, Vert_Rank);
   else if(type == 1)
    METIS_PartMeshDual(&N_Cells, &N_RootVertices, MetisVertexNumbers, &etype, &numflag,
                       &size, &edgecut, Cell_Rank, Vert_Rank);
   else
    {
     cout<<" Error METIS_PartMesh implemented for Par_P2 = 0 or 1 !!" <<endl;
     MPI_Abort(comm, 0);
    }

  t2 = MPI_Wtime();
  OutPut( "Time taken for METIS mesh partinioning "<< t2-t1<< " sec"<<endl); 


  // collect the info for subdomain mesh manupulation
   PointNeighb = new int[N_RootVertices];
   memset(PointNeighb, 0, N_RootVertices*SizeOfInt);

   // find cells per vertex 
   for (i=0;i<N_AllLocVert;i++)
    PointNeighb[VertexNumbers[i]]++;
 
   // find maximum cells per vertex ( maxCpV )
   for (i=0;i<N_RootVertices;i++)
    if (PointNeighb[i] > MaxCpV) MaxCpV = PointNeighb[i];

//       printf("Max number of cells per vertex %d \n", MaxCpV);

  /** PointNeighb's first column contains number of neib cells associated with each vertex*/
  /** further columns contain the cell numbers associated with this vertex*/
   MaxCpV++;
  for(i=1;i<size;i++)
   MPI_Send(&MaxCpV, 1, MPI_INT, i, 85, comm);

   delete [] PointNeighb;
   PointNeighb = new int[N_RootVertices*MaxCpV];
   memset(PointNeighb, 0, (N_RootVertices*MaxCpV)*SizeOfInt);

   for(i=0;i<N_Cells;i++)
    {
     cell = coll->GetCell(i);

     for(j=0;j<N_VertInCell;j++)
      {
       M = VertexNumbers[i*N_VertInCell + j] *MaxCpV ;
       PointNeighb[M]++;
       PointNeighb[M + PointNeighb[M]  ] = i;
      } // for(j=0;j<k;j++) 
    } //  for(i=0;i<N_Cells;i++) 

   for(i=1;i<size;i++)
    MPI_Send(PointNeighb, N_RootVertices*MaxCpV, MPI_INT, i, 90, comm); 

   delete [] MetisVertexNumbers;
   delete [] NumberVertex;
   delete [] Vertices;
  }
 else
  {
    MPI_Recv(&N_RootVertices, 1, MPI_INT, 0, 75, comm, &status);
    // printf("%d SubDomain_N_Cells in rank test 1, %d  \n", rank, N_RootVertices );

    MPI_Recv(VertexNumbers, N_AllLocVert, MPI_INT, 0, 80, comm, &status);

    MPI_Recv(&MaxCpV, 1, MPI_INT, 0, 85, comm, &status);

    PointNeighb = new int[N_RootVertices*MaxCpV];

    MPI_Recv(PointNeighb, N_RootVertices*MaxCpV, MPI_INT, 0, 90, comm, &status);     
  } //else if(rank==0)


    /**Metis partition done!! code for all processors */  
    MPI_Bcast(Cell_Rank, N_Cells, MPI_INT, 0, comm);
    MPI_Bcast(Vert_Rank, N_RootVertices, MPI_INT, 0, comm);

//     if(rank==1)
//     for(i=0;i<N_RootVertices;i++)
//      printf("%d   mesh partition %d Vert_Rank %d \n", rank, i, Vert_Rank[i]);

    /** put the rank ID and global cell No in each cell */
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      cell->SetSubDomainNo(Cell_Rank[i]);
      cell->SetGlobalCellNo(i);
//       if(i == 111615)
// 	 printf("rank %d Mesh Partition N_CellsIn_b %d\n", rank, Cell_Rank[i]);
//       
     }
//      MPI_Finalize();
//      exit(0); 
//   if(rank!=0) // rank will not take part in computations
   {
 /** set all subdomain vertices */  
     // set the clipboard
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);

      if(cell->GetSubDomainNo()==rank) 
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

    /** fill the VertNeibRanks info */
    /** first column contains how many ranks contain this vertex */
    /** further colums contain the rank ID of the subdomains */
    VertNeibRanks = new int[N_RootVertices*MaxCpV];
    memset(VertNeibRanks, 0, N_RootVertices*MaxCpV*SizeOfInt);

    HaloCellIndex = new int[size];
    for(i=0;i<size;i++)
     HaloCellIndex[i] = -1;

    HaloCellGlobalNo = new int[MaxCpV];
    HaloCellLocVertNo = new int[MaxCpV];

    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      ID = cell->GetSubDomainNo();

        /**run only through own cells */
      if(ID==rank) 
       {
        cell->SetAsOwnCell();

        //set SubDomainVert if any vert in this cell is so
        //needed for moving meshes and setting cross vertex
        for(j=0;j<N_VertInCell;j++)
         {
          CurrVert = cell->GetVertex(j);
	  
//           if(  i == 99172 && j==3  ) //111615
// 	  {
//  	    N = VertexNumbers[i*N_VertInCell + j]*MaxCpV ;
// 	    N_CellsInThisVert =  PointNeighb[N];
// 	       printf("Raaaaaaaank %d vert %d N_CellsInThisVert  %d\n", rank, j, N_CellsInThisVert);
// 	       
//          for(ii=1;ii<=N_CellsInThisVert;ii++)
//            {
//             VertexCellNo = PointNeighb[N + ii];
//             Vertexcell = coll->GetCell(VertexCellNo);
//             VertexcellID = Vertexcell->GetSubDomainNo();
//             
// 	    if(VertexcellID==2)
//              printf("Raaank %d Neib No %d Neibcell GlobalNo %d NeibID  %d\n", rank, ii, Vertexcell->GetGlobalCellNo(), VertexcellID); 
// 	    
// 	    
// 	   }
// 	       
// 	  }
// 	  
//           if(  i == 44171 && j==1  )
// 	  {
//  	    N = VertexNumbers[i*N_VertInCell + j]*MaxCpV ;
// 	    N_CellsInThisVert =  PointNeighb[N];
// 	       printf("Rbbbbbbbbbbbbk %d vert %d N_CellsInThisVert  %d\n", rank, j, N_CellsInThisVert);
// 	       
//          for(ii=1;ii<=N_CellsInThisVert;ii++)
//            {
//             VertexCellNo = PointNeighb[N + ii];
//             Vertexcell = coll->GetCell(VertexCellNo);
//             VertexcellID = Vertexcell->GetSubDomainNo();
//             
// 	    if(VertexcellID==0)
//              printf("Rbbbnk %d Neib No %d Neibcell GlobalNo %d NeibID  %d\n", rank, ii, Vertexcell->GetGlobalCellNo(), VertexcellID); 
// 	    
// 	    
// 	   }
// 	       
// 	  }
  

          if(CurrVert->GetClipBoard() != -1)
           continue;

          CurrVert->SetClipBoard(5);
          N = VertexNumbers[i*N_VertInCell + j]*MaxCpV ;

          N_CellsInThisVert =  PointNeighb[N];
          N_OwnIncidentCells = 0;
          N_NeibIncidentCells = 0;

          // check! any subdomain cell containg this vert has to be added
          // Setting all halo cells
          UPDATE = TRUE;
          for(ii=1;ii<=N_CellsInThisVert;ii++)
           {
            VertexCellNo = PointNeighb[N + ii];
            Vertexcell = coll->GetCell(VertexCellNo);
            VertexcellID = Vertexcell->GetSubDomainNo();

            if(VertexcellID != rank)
             {
              //vertex belong to diff subdomain 
              if(UPDATE)
               {
                cell->SetAsDependentCell();
                CurrVert->SetAsSubDomainVert();
                HaloCellIndex[VertexcellID]= VertexCellNo;
                UPDATE = FALSE;
               }
              Vertexcell->SetAsHaloCell();

//               if(HaloCellIndex[VertexcellID]< VertexCellNo)
//                 HaloCellIndex[VertexcellID]= VertexCellNo;

               /** changed on 4 Feb 2012 - by Sashi */
               if(HaloCellIndex[VertexcellID] == -1)
                HaloCellIndex[VertexcellID]= VertexCellNo;

             }
            } //  for(ii=1;ii< 

//           if(rank==2 && i==44171)            
//             for(ii=0;ii<size;ii++)
//                printf(" HaloCellIndex rank %d HaloCellIndex %d \n", rank, HaloCellIndex[ii] );
     


           if(CurrVert->IsSubDomainVert())
             for(ii=1;ii<=N_CellsInThisVert;ii++)
              {
               VertexCellNo = PointNeighb[N + ii];
               Vertexcell = coll->GetCell(VertexCellNo);
               VertexcellID = Vertexcell->GetSubDomainNo();

               // all cells associated with this vert are DepCells
               Vertexcell->SetAsDependentCell();

               if(VertexcellID==rank)
                continue;

               N_CellsIn_b = VertNeibRanks[N];
               UPDATE = TRUE;

               for(jj=1;jj<=N_CellsIn_b;jj++)
                if( VertNeibRanks[N + jj]==VertexcellID)
                 {
                  UPDATE = FALSE;
                  break;
                 }

               if(UPDATE)
                {
                 N2 = HaloCellIndex[VertexcellID];
                 HaloCellGlobalNo[VertNeibRanks[N]] = N2;

                 //find the local index of this vertex in the neib cell
                 jj=0;
                 while(CurrVert != (coll->GetCell(N2))->GetVertex(jj)) jj++;
                 HaloCellLocVertNo[VertNeibRanks[N]] = jj;

                 VertNeibRanks[N]++;
                 VertNeibRanks[N + VertNeibRanks[N]] = VertexcellID;

                 if(MaxRankPerV<VertNeibRanks[N])
                   MaxRankPerV=VertNeibRanks[N];
                } //   if(UPDATE)
              } //  for(ii=1;ii<=N_CellsInT

           N_SubDomIn_a = VertNeibRanks[N];
           Temp = VertNeibRanks+(N+1);
   
 
//         if(rank==2 && i==44171 && j==1)
//          {
//  
// 	    
//           for(jj=0;jj<N_SubDomIn_a;jj++)       
//             printf("SetSubDomainInfo rank %d vert %d Neib %d GCell  %d neibvert %d\n", rank, j, Temp[jj], HaloCellGlobalNo[jj], HaloCellLocVertNo[jj]);
// 	 }
	 
/*        if(rank==0 && i==111615 && j==1)
         {
 
	    
          for(jj=0;jj<N_SubDomIn_a;jj++)       
            printf("SetSubDomainInfo rank %d vert %d Neib %d GCell  %d neibvert %d\n", rank, j, Temp[jj], HaloCellGlobalNo[jj], HaloCellLocVertNo[jj]);
	 }*/
    
           /** set only one cell (lowest index cell) from each subdomain as the Halo cell for this vertex */
           CurrVert->SetSubDomainInfo(N_SubDomIn_a,  Temp, HaloCellGlobalNo, HaloCellLocVertNo);

            //reset
            for(jj=0;jj<size;jj++)
             HaloCellIndex[jj] = -1;

           } // for(j=0;j<N_V
       } // if(ID==ra
     } //  for(i=0;i<N_Cel 

    // including own rank
    MaxRankPerV++;
    
//     MPI_Finalize();
//     exit(0);  
 
     
    delete [] HaloCellIndex;
    delete [] HaloCellGlobalNo;
    delete [] HaloCellLocVertNo;

     // get own SubDomain collection of cells 
     N_LocalCells = 0;
     N_OwnCells = 0;

     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
       ID = cell->GetSubDomainNo();

       if(ID==rank || cell->IsHaloCell() )
        N_LocalCells++;

       if(ID==rank)
        N_OwnCells++;
      } //  for (i=0;i<N_Cells;i++)

      
     printf("Rank: %d N_own Cells %d\n ", rank,  N_OwnCells);   
      
      
     if(N_LocalCells)
      {
       // collect the own collection
       SubDomainCells = new TBaseCell*[N_LocalCells];
       GlobalCellIndex = new int[N_LocalCells];
      }
     else
      {
       SubDomainCells = NULL;
       GlobalCellIndex = NULL;
      }


     N=0;
     M=0;
     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
       ID = cell->GetSubDomainNo();

       // inthe collection, first own cells then Halo cells
       if(ID==rank)
        {
         SubDomainCells[N] =  cell;
         GlobalCellIndex[N] =  i;
         N++;
        }
      else if(cell->IsHaloCell() )
        {
         SubDomainCells[N_OwnCells + M] =  cell;
         GlobalCellIndex[N_OwnCells + M] =  i;
         M++;
        }
      else //cell no longer needed for this rank, delete!!!!!
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
          EdgeDel[N_EdgeDel++] = cell->GetEdge(j);

        for(j=0;j<N_VertInCell;j++)
         VertexDel[N_VertexDel++] = cell->GetVertex(j);

         CellDel[N_CellDel++] = cell;
       }
      } //  for (i=0;i<N_Cells;i++)

/*      if(rank ==0)
        {
	  CellDel = new TBaseCell*[50];
	  
         cell =  SubDomainCells[20554];
	 j = 4;
         edge = cell->GetEdge(j);
         
//           edge->GetNeibSubDomainRanks(N, CrossNeibsRank);
          
          edge->GetNeibs(N, CellDel);
	  
	              printf("j %d N  %d\n ", j, N); 
	  
          for(j=0;j<N;j++)
             printf("j %d Neibcell %d\n ", j, CellDel[j]->GetSubDomainNo());
        }        
     
     MPI_Finalize();
     exit(0);     */   

//   if(rank==out_rank)
//      printf("Mesh Partition 1 rank %d \n" , rank);     

    // delete other processors cells from own collection, excluding halo cells
    if(N_VertexDel)
     Sort(VertexDel, N_VertexDel); 
    if(N_EdgeDel)
     Sort(EdgeDel, N_EdgeDel);
    if(N_CellDel)
     Sort(CellDel, N_CellDel);

    Last=NULL;
    for(i=0;i<N_VertexDel;i++)
     {
      if( VertexDel[i] != Last)
       {
        delete VertexDel[i];
        Last = VertexDel[i];
       }
     }
   delete [] VertexDel;


    LastEdge=NULL;
    for(i=0;i<N_EdgeDel;i++)
     {
      if( EdgeDel[i] != LastEdge)
       {
        delete EdgeDel[i];
        LastEdge = EdgeDel[i];
       }
     }
   delete [] EdgeDel;

    LastCell=NULL;
    for(i=0;i<N_CellDel;i++)
     {
      if( CellDel[i] != LastCell)
       {
        delete CellDel[i];
        LastCell = CellDel[i];
       }
     }
   delete [] CellDel;

    if((N+M)!=N_LocalCells)
     {
      printf("Error in mesh partition %d \n", rank);
      MPI_Abort(comm, 0);
     } 
         
//      if(rank==out_rank)
//      printf("Mesh Partition 2 rank %d \n" , rank);           
         
         
//     MPI_Finalize();
//     exit(0);  


 /** edges will be flaged as SubDomainEdges and subdomain infos are filled */
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         {
          edge = cell->GetEdge(j);
          if(edge->GetClipBoard()==-1)
           {
            edge->SetClipBoard(5);
            edge->InitSubDomainInfo(rank);
           } // if(edge->GetCli
         } // for(j=0;j<N_Edges;
       } // if(cell->IsDependentCell
     } // for(i=0;i<N_Cells;i++

 /** change all subdomain joint from JointEqN into TSubDomainJoint */
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      /**run only through Dependent cells */
      if(cell->IsDependentCell()) 
       {
        GblCellNr = GlobalCellIndex[i];
        ID = cell->GetSubDomainNo();

	
	//if( GblCellNr == 111615)
          // printf("rank %d local cell no %d\n", rank, i);
	
        /** change the faces, also edges in 3D */
        for(j=0;j<N_JointsInCell;j++)
         {
          Joint = cell->GetJoint(j);

          if(Joint->GetType() == JointEqN)
           {
            neib_cell = Joint->GetNeighbour(cell);
            Neib_ID = neib_cell->GetSubDomainNo();


            if(ID!=Neib_ID)
             {
              GlobCellNo =  neib_cell->GetGlobalCellNo();

              // this joint belongs to two SubDomains
              cell->SetAsSubDomainInterfaceCell();

              // find neib cell local face number
              m=0;
              while(Joint!=neib_cell->GetJoint(m)) m++;

              delete Joint;
              NewJoint = new TSubDomainJoint(cell, neib_cell, Neib_ID, GlobCellNo, m);
              cell->SetJoint(j, NewJoint);
              neib_cell->SetJoint(m, NewJoint);
              NewJoint->SetMapType();

              // set all edges in this face as SubDomainEdges
              ShapeDesc = cell->GetShapeDesc();
              ShapeDesc->GetFaceEdge(TmpFV, TmpLen, MaxLen);
              N_FaceEdges = TmpLen[j];
              for(n=0;n<N_FaceEdges;n++)
               {
                edge = cell->GetEdge(TmpFV[j*MaxLen+n]);
                edge->SetAsNotCrossEdgeFor(rank, Neib_ID);
               } // for(n=0;n<N_Edge

              // remove the Neib_ID from VertNeibRanks
              // these face vertices cannot be cross vertices for processors ID and Neib_ID
              ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
              for(ii=0;ii<TmpLen[j];ii++)
              {
               jj = TmpFV[j*MaxLen + ii];
               M = VertexNumbers[GblCellNr*N_VertInCell + jj] *MaxCpV;
               N = VertNeibRanks[M];

               for(jj=1;jj<=N;jj++)
                if(VertNeibRanks[M + jj] == Neib_ID )
                 {
                  VertNeibRanks[M + jj] = VertNeibRanks[M + N];
                  VertNeibRanks[M]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)
              } // for(ii=0;ii<TmpL
             }//if(ID!=Neib_ID
           } // if(Joint->GetT
         }// for (j=0;j< ;
       } // if(cell->IsDependentCell() 
     }// for(i=0;i<N_Cel
    
     
    /** find cross edges (if any) */
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        M = GlobalCellIndex[i];
        N_Edges=cell->GetN_Edges();
        ShapeDesc= cell->GetShapeDesc();
        ShapeDesc->GetEdgeVertex(EdgeVertex);

        for(j=0;j<N_Edges;j++)
         {
          edge = cell->GetEdge(j);

          if(edge->GetClipBoard()==-1 && edge->IsSubDomainEdge())
           {
            edge->SetClipBoard(5);
            N_CrossNeibs = edge->GetN_CrossNeibs();

            if(N_CrossNeibs==0)
             continue;

            cell->SetAsCrossEdgeCell();

            a = EdgeVertex[2*j];
            b = EdgeVertex[2*j+ 1];

            Vert_a = cell->GetVertex(a);
            edge->SetCrossNeibInfo(Vert_a);
            edge->GetCrossEdgeNeibs(N_CrossNeibs, CrossNeibsRank);

            M1 = VertexNumbers[M*N_VertInCell + a]*MaxCpV;
            N1 = VertNeibRanks[M1];

            M2 = VertexNumbers[M*N_VertInCell + b]*MaxCpV;
            N2 = VertNeibRanks[M2];

            for(jj=0;jj<N_CrossNeibs;jj++)
             {
              Neib_ID = CrossNeibsRank[jj];

               for(kk=1;kk<=N1;kk++)
                if(VertNeibRanks[M1 + kk] == Neib_ID )
                 {
                  VertNeibRanks[M1 + kk] = VertNeibRanks[M1 + N1];
                  VertNeibRanks[M1]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)

               for(kk=1;kk<=N2;kk++)
                if(VertNeibRanks[M2 + kk] == Neib_ID )
                 {
                  VertNeibRanks[M2 + kk] = VertNeibRanks[M2 + N2];
                  VertNeibRanks[M2]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)
             } // for(jj=0;jj<N_CrossNei
           } //   if(edge->GetClipBoard()==-1 && edge->IsSubDomainEdge())
          } // for(j=0;j<N_Edge
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel
         
    /** find cross vertex (if any) */
    // set the clipboard
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

    /** set incident cell list for all vertices */
    IncCells = new TBaseCell*[MaxCpV];
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      M = GlobalCellIndex[i];    
     
      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        {
         // set the vertexcells for all vertices in dep. cells
         CurrVert = cell->GetVertex(j);
         if(CurrVert->GetClipBoard() != -1)
          continue;

         CurrVert->SetClipBoard(5);

         N = VertexNumbers[M*N_VertInCell + j] *MaxCpV;
         N_CellsInThisVert =  PointNeighb[N];

         for(ii=1;ii<=N_CellsInThisVert;ii++)
          IncCells[ii-1] = coll->GetCell(PointNeighb[N + ii]);

         CurrVert->SetVertexCells(N_CellsInThisVert, IncCells);
        } //  for (j=0;j<N_VertInCell;
     } //  for(i=0;i<N_OwnCell
    delete [] IncCells;


     // test 
//        i = 111615; j=1;
//         if(rank==0 )
//          {
//           N = VertexNumbers[i*N_VertInCell + j]*MaxCpV ;
//           N_CellsIn_b = VertNeibRanks[N];
//           printf("test rank %d Mesh Partition N_CellsIn_b %d\n", rank, N_CellsIn_b);
// 	    
//           for(jj=1;jj<=N_CellsIn_b;jj++)       
//             printf("test rank %d vert %d Neib %d VertNeibRanks %d\n", rank, j, jj, VertNeibRanks[N + jj]);
//          }
//          
//          i = 44171;        j=1;
//         if(rank==2 )
//          {
//           N = VertexNumbers[i*N_VertInCell + j]*MaxCpV ;
//           N_CellsIn_b = VertNeibRanks[N];
//           printf("rank %d Mesh Partition N_CellsIn_b %d\n", rank, N_CellsIn_b);
// 	    
//           for(jj=1;jj<=N_CellsIn_b;jj++)       
//             printf("rank %d vert %d Neib %d VertNeibRanks %d\n", rank, j, jj, VertNeibRanks[N + jj]);
//          }       
//      MPI_Finalize();
//      exit(0); 

     // set the clipboard
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        M = GlobalCellIndex[i];

        for(j=0;j<N_VertInCell;j++)
         {
          CurrVert = cell->GetVertex(j);

          // continue if the vert is already handled or not a subdomain vert
          if( (CurrVert->GetClipBoard() != -1) || !(CurrVert->IsSubDomainVert()))
           continue;

//         if(rank==0 && M==99172 && j==3)
//          {
//           N = VertexNumbers[M*N_VertInCell + j]*MaxCpV ;
//           N_CellsIn_b = VertNeibRanks[N];
//           printf("AddCrossNeib rank %d local cell no %d N_CellsIn_b %d\n", rank, i, N_CellsIn_b);
//   
//           for(jj=1;jj<=N_CellsIn_b;jj++)       
//             printf("AddCrossNeib rank %d vert %d Neib %d VertNeibRanks %d\n", rank, j, jj, VertNeibRanks[N + jj]);
// 
//           TDatabase::ParamDB->Par_P5 =1;
// 	  if( cell->IsCrossVertexCell())
// 	    printf("rank %dvert %d Is SubDomainVert\n", rank,  j);
// // 	    
//          }

          CurrVert->SetClipBoard(5);

          M1 = VertexNumbers[M*N_VertInCell + j] *MaxCpV ;
          N1 = VertNeibRanks[M1];

          if(N1!=0)
           cell->SetAsCrossVertexCell();

          for(ii=1;ii<=N1;ii++)
            CurrVert->AddCrossNeib(VertNeibRanks[M1+ii]);

//            TDatabase::ParamDB->Par_P5 =0;
          } // for(j=0;j<N_VertInCell;j
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel


    if(N_LocalCells)
     { Domain->ReplaceTreeInfo(N_LocalCells, SubDomainCells, GlobalCellIndex, N_OwnCells); }
    else
     { Domain->SetN_OwnCells(N_OwnCells); }

        // initialize iterators
    TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
    TDatabase::IteratorDB[It_LE]->SetParam(Domain);
    TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
    TDatabase::IteratorDB[It_Between]->SetParam(Domain);
    TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

    delete [] VertNeibRanks;
   } // 

//         M=18482;
//         j=3;
//         if(rank==0)	  
//          {
//           cell = SubDomainCells[M];
// 	  CurrVert = cell->GetVertex(j);
//           CurrVert->GetCrossNeibsInfo(N,  Temp, HaloCellGlobalNo, HaloCellLocVertNo);
//  
// //                printf("AAAA rank %d  N %d\n", rank, N);  
//           for(jj=0;jj<N;jj++)       
//             printf("rank %d Finalvert %d cross Neib rank %d GCellNo %d LocvertNo %d\n", rank, j, Temp[jj],  HaloCellGlobalNo[jj], HaloCellLocVertNo[jj] );  
// 
//             if(cell->IsCrossVertexCell())
//               printf("rank %d Finallocal Cell No %d Global cell No %d\n", rank, M, cell->GetGlobalCellNo());  
//          }
// 
//         M=9068;
//         j=1;
//         if(rank==2)	  
//          {
//           cell = SubDomainCells[M];
// 	  CurrVert = cell->GetVertex(j);
//           CurrVert->GetCrossNeibsInfo(N,  Temp, HaloCellGlobalNo, HaloCellLocVertNo);
//  
//  
//           for(jj=0;jj<N;jj++)       
//             printf("rank %d Finalvert %d cross Neib rank %d GCellNo %d LocvertNo %d\n", rank, j, Temp[jj],  HaloCellGlobalNo[jj], HaloCellLocVertNo[jj] );  
// 
//           if(cell->IsCrossVertexCell())
//            printf("rank %d Finallocal Cell No %d Global cell No %d\n", rank, M, cell->GetGlobalCellNo());  
//          }          
//    

 
//    delete [] Vert_Rank; // for more no. processor showing double free or corruption 
   delete [] VertexNumbers;
   delete [] Cell_Rank;
   delete [] PointNeighb;

   
   //Barrier is needed, before calling FECommunicator, since neib process must have all info
   MPI_Barrier(comm);
   
  if(rank==out_rank)
     printf("Mesh Partition end rank %d \n" , rank);

//      MPI_Finalize();
//     exit(0);  



} // Partition_Mesh3D()


void Partition_Mesh3D(MPI_Comm comm, TDomain *Domain, int &MaxRankPerV)
{
 idxtype *Cell_Rank, *Vert_Rank; 
 int i, j, k, m, M, n, rank, size, N_Cells, N, ID, Neib_ID, out_rank= TDatabase::ParamDB->Par_P0;
 int *GlobalCellIndex, N_RootVertices, N_VertInCell, N_JointsInCell, N_AllLocVert;
 int *VertexNumbers, *PointNeighb, maptype, MaxCpV;
 int MaxLen, N_Edges, N_FaceEdges, ii, jj, GblCellNr;
 int N_EdgeDel=0, N_VertexDel=0, N_CellDel=0;
 const int *TmpFV, *TmpLen, *TmpVE, *EdgeVertex, *NeibEdgeVertex;
 
 int VertexCellNo, VertexcellID, N_CellsInThisVert, *VertNeibRanks;
 int N_LocalCells, N_OwnCells, N_OwnIncidentCells, N_NeibIncidentCells;
 int a, b, N_SubDomInThisVert, *HaloCellIndex;
 int M1, M2, N1, N2, N_CellsIn_a, N_CellsIn_b;
 int N_CrossNeibs, *CrossNeibsRank, *HaloCellGlobalNo, *HaloCellLocVertNo;
 int N_SubDomIn_a, N_SubDomIn_b, *Temp, kk, GlobCellNo, test_b, EdgeCellID;

 bool UPDATE, UPDATE_1;

 TVertex *Vert_a, **VertexDel, *Last;
 TCollection *coll;
 TBaseCell *cell, *neib_cell, *Vertexcell, **SubDomainCells, *cell_a, *cell_b, **IncCells; 
 TBaseCell **CellDel, *LastCell;
 TJoint *Joint, *NewJoint;
 TEdge *edge, **EdgeDel, *LastEdge;
 TShapeDesc *ShapeDesc, *NeibShapeDesc;
 TVertex *CurrVert, *NeibVert_a; 
 
 MPI_Status status, status1, status2, status3;
 MPI_Request request, request1;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  // root will not take part in computations
  MaxCpV = 0;
  MaxRankPerV = -1;

//   if(rank==0)
//    MaxRankPerV = 1;
  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  cell = coll->GetCell(0);
  N_VertInCell = cell->GetN_Vertices();
  N_JointsInCell = cell->GetN_Joints();
  N_AllLocVert = N_VertInCell*N_Cells;
  VertexNumbers= new int[N_AllLocVert];
  Vert_Rank = new idxtype[N_AllLocVert];  
  Cell_Rank = new idxtype[N_Cells];  

 // if(rank==0)
 //  printf("Number of  ranks: %d\n",  size);

  if(N_Cells<size)
   {
     printf("Number of cells less than number of processors !!!, %d\n",  N_Cells);
     MPI_Finalize();
     exit(0);
   }

 if(size==1)
   {
    cout <<  "Total number of process should be grater than 1 (or 2 if root is not involved in computation), but you have " << size<<endl;
     MPI_Finalize();
     exit(0);
   }

 //check, all cell have same number of vertices  
  for(i=1;i<N_Cells;i++)
   if(((coll->GetCell(i))->GetN_Vertices() != N_VertInCell) || ( (coll->GetCell(i))->GetN_Joints() != N_JointsInCell) )
    {
     cout << N_JointsInCell << "Mesh partition for heterogeneous cells are not yet implemented " << N_VertInCell<<endl;
     MPI_Finalize();
     exit(0);
    }


// partition the mesh in the root and send info to all processors
 if(rank==0)
  {
  // cout <<  "Total Cells " << N_Cells<<endl;

   int  m1, edgecut=0, *NumberVertex;
   int etype, numflag=0; // c-style numbering (array starting with 0)
   int type=TDatabase::ParamDB->Par_P2;

   idxtype *MetisVertexNumbers; 

   double t1, t2;

   TVertex **Vertices;    

   
    MetisVertexNumbers= new idxtype[N_AllLocVert];
    NumberVertex=new int[N_AllLocVert];
    Vertices=new TVertex*[N_AllLocVert];
   
  /** *********************************************/
  /** STEP 1 : STORE VERTEX POINTERS CELL-WISE */
  /** *********************************************/
    N = 0;    
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      for(j=0;j<N_VertInCell;j++)
       {
        Vertices[N]=cell->GetVertex(j);
        N++;
       }
     } // for i=0

  /** *********************************************/
  /** STEP 2 : SORT THE VERTICES ARRAY */ 
  /** *********************************************/
    Sort(Vertices, N);
    
  /** ***************************************************/
  /**STEP 3: STORE THE SORTED POINTER ARRAY AS INDICES */
  /** ***************************************************/
    Last=NULL;
    N_RootVertices=-1;
    for(i=0;i<N_AllLocVert;i++)
     {
      if((CurrVert=Vertices[i])!=Last)
       {
        N_RootVertices++;
        Last=CurrVert;
       }
      NumberVertex[i]=N_RootVertices;
     }
    N_RootVertices++;


  /** *********************************************/   
  /** STEP 4 : STORE THE INDICES CELL-WISE */
  /** *********************************************/ 
    m=0;
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);

      for(j=0;j<N_VertInCell;j++)
       {
        CurrVert=cell->GetVertex(j);
        N=GetIndex(Vertices, N_AllLocVert, CurrVert);
        VertexNumbers[m]=NumberVertex[N];
        m++;
       } // endfor j
     } //endfor i

   memcpy(MetisVertexNumbers, VertexNumbers, N_AllLocVert*SizeOfInt);

 //  cout << "Total Vertices " <<N_RootVertices<<endl;
   //cout << "Vertices N " <<N_AllLocVert<<endl;
   
  /** 1. SEND INFO TO ALL PROCESSORS */
  for(i=1;i<size;i++)
   MPI_Send(&N_RootVertices, 1, MPI_INT, i, 75, comm);
  
  //   for (i=0; i<N_LocVertices; i++) 
 //    cout<< i <<"  " << "elements " << MetisVertexNumbers[i] <<endl;

    if(N_VertInCell==4)
     { etype = 2; } // tetrahedral
    else if(N_VertInCell==8)
     { etype = 3; } // Hexahedral
    else
     {
      cout<<" Error only  Tetra or Hexa mesh can be partitioned !!" <<endl;
      MPI_Finalize();
      exit(0);
     }
     
  /** 2. SEND INFO TO ALL PROCESSORS */    
  for(i=1;i<size;i++)
    MPI_Send(VertexNumbers, N_AllLocVert, MPI_INT, i, 80, comm);  
  
  /** *********************************************/
  /** STEP 5 : MESH PARTITION */
  /** *********************************************/
  t1 = MPI_Wtime();
   if(type == 0)
    METIS_PartMeshNodal(&N_Cells, &N_RootVertices, MetisVertexNumbers, &etype, &numflag,
                       &size, &edgecut, Cell_Rank, Vert_Rank);
   else if(type == 1)
    METIS_PartMeshDual(&N_Cells, &N_RootVertices, MetisVertexNumbers, &etype, &numflag,
                       &size, &edgecut, Cell_Rank, Vert_Rank);
   else
    {
     cout<<" Error METIS_PartMesh implemented for Par_P2 = 0 or 1 !!" <<endl;
     MPI_Abort(comm, 0);
    }

  t2 = MPI_Wtime();
  OutPut( "Time taken for METIS mesh partinioning "<< t2-t1<< " sec"<<endl); 


  /** *********************************************/
  /** STEP 6 : TO FIND MAXIMUM CELLS PER VERTEX */
  /** *********************************************/ 
   PointNeighb = new int[N_RootVertices];
   memset(PointNeighb, 0, N_RootVertices*SizeOfInt);

  
  for (i=0;i<N_AllLocVert;i++)
    PointNeighb[VertexNumbers[i]]++;
 
   // find maximum cells per vertex ( maxCpV )
   for (i=0;i<N_RootVertices;i++)
    if (PointNeighb[i] > MaxCpV) MaxCpV = PointNeighb[i];

//       printf("Max number of cells per vertex %d \n", MaxCpV);
   MaxCpV++;
  for(i=1;i<size;i++)
   MPI_Send(&MaxCpV, 1, MPI_INT, i, 85, comm);

   delete [] PointNeighb;
  
  /** ********************************************************************************************/
  /** STEP 7 : CREATE AN ARRAY CONTAINING INFO REGARDING THE NEIGHBOURS OF VERTICES (CELL INDEX) */
  /** PointNeighb's first column contains number of neib cells associated with each vertex*/
  /** further columns contain the cell numbers associated with this vertex*/
  /** ********************************************************************************************/
   PointNeighb = new int[N_RootVertices*MaxCpV];
   memset(PointNeighb, 0, (N_RootVertices*MaxCpV)*SizeOfInt);

   for(i=0;i<N_Cells;i++)
    {
     cell = coll->GetCell(i);

     for(j=0;j<N_VertInCell;j++)
      {
       M = VertexNumbers[i*N_VertInCell + j] *MaxCpV ;
       PointNeighb[M]++;
       PointNeighb[M + PointNeighb[M]  ] = i;
      } // for(j=0;j<k;j++) 
    } //  for(i=0;i<N_Cells;i++) 

   for(i=1;i<size;i++)
    MPI_Send(PointNeighb, N_RootVertices*MaxCpV, MPI_INT, i, 90, comm); 

   delete [] MetisVertexNumbers;
   delete [] NumberVertex;
   delete [] Vertices;
  }
 else
  {
    MPI_Recv(&N_RootVertices, 1, MPI_INT, 0, 75, comm, &status);
    // printf("%d SubDomain_N_Cells in rank test 1, %d  \n", rank, N_RootVertices );

    MPI_Recv(VertexNumbers, N_AllLocVert, MPI_INT, 0, 80, comm, &status);

    MPI_Recv(&MaxCpV, 1, MPI_INT, 0, 85, comm, &status);

    PointNeighb = new int[N_RootVertices*MaxCpV];

    MPI_Recv(PointNeighb, N_RootVertices*MaxCpV, MPI_INT, 0, 90, comm, &status);     
  } //else if(rank==0)


    /**Metis partition done!! code for all processors */  
    MPI_Bcast(Cell_Rank, N_Cells, MPI_INT, 0, comm);
    MPI_Bcast(Vert_Rank, N_RootVertices, MPI_INT, 0, comm);
  
  /** ********************************************************/
  /** STEP 8.1 : SETTING SUBDOMAIN NUMBER AND GLOBAL CELL NO */
  /** ********************************************************/
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      cell->SetSubDomainNo(Cell_Rank[i]);
      cell->SetGlobalCellNo(i); 
     }

  /** *********************************************/
  /** STEP 8.2 : SET ALL VERTICES TO -1 */
  /** *********************************************/ 
    
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);

      if(cell->GetSubDomainNo()==rank) 
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

  /** *********************************************/
  /** STEP 9 : Fill the VertNeibRanks info        */
  /** first column contains how many ranks contain this vertex further columns contain the rank ID of the subdomains */
  /** *********************************************/ 
    VertNeibRanks = new int[N_RootVertices*MaxCpV];
    memset(VertNeibRanks, 0, N_RootVertices*MaxCpV*SizeOfInt);

    HaloCellIndex = new int[size];
    for(i=0;i<size;i++)
      HaloCellIndex[i] = -1;

    HaloCellGlobalNo = new int[MaxCpV];
    HaloCellLocVertNo = new int[MaxCpV];

    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      ID = cell->GetSubDomainNo();

  /**run only through own cells */
      if(ID==rank) 
       {
        cell->SetAsOwnCell();

        //set SubDomainVert if any vert in this cell is so needed for moving meshes and setting cross vertex
        for(j=0;j<N_VertInCell;j++)
         {
          CurrVert = cell->GetVertex(j);
	  
	  if(CurrVert->GetClipBoard() != -1)
           continue;

          CurrVert->SetClipBoard(5);
          N = VertexNumbers[i*N_VertInCell + j]*MaxCpV ;

          N_CellsInThisVert =  PointNeighb[N];
          N_OwnIncidentCells = 0;
          N_NeibIncidentCells = 0;

          // check! any subdomain cell containg this vert has to be added
          // Setting all halo cells
          UPDATE = TRUE;
          for(ii=1;ii<=N_CellsInThisVert;ii++)
           {
            VertexCellNo = PointNeighb[N + ii];
            Vertexcell = coll->GetCell(VertexCellNo);
            VertexcellID = Vertexcell->GetSubDomainNo();

            if(VertexcellID != rank)
             {
              //vertex belong to diff subdomain 
              if(UPDATE)
               {
                cell->SetAsDependentCell();
                CurrVert->SetAsSubDomainVert();
                HaloCellIndex[VertexcellID]= VertexCellNo;
                UPDATE = FALSE;
               }
              Vertexcell->SetAsHaloCell();

               /** changed on 4 Feb 2012 - by Sashi */
               if(HaloCellIndex[VertexcellID] == -1)
                HaloCellIndex[VertexcellID]= VertexCellNo;

             }
            } //  for(ii=1;ii< 

           if(CurrVert->IsSubDomainVert())
             for(ii=1;ii<=N_CellsInThisVert;ii++)
              {
               VertexCellNo = PointNeighb[N + ii];
               Vertexcell = coll->GetCell(VertexCellNo);
               VertexcellID = Vertexcell->GetSubDomainNo();

               // all cells associated with this vert are DepCells
               Vertexcell->SetAsDependentCell();

               if(VertexcellID==rank)
                continue;

               N_CellsIn_b = VertNeibRanks[N];
               UPDATE = TRUE;

               for(jj=1;jj<=N_CellsIn_b;jj++)
                if( VertNeibRanks[N + jj]==VertexcellID)
                 {
                  UPDATE = FALSE;
                  break;
                 }

               if(UPDATE)
                {
                 N2 = HaloCellIndex[VertexcellID];
                 HaloCellGlobalNo[VertNeibRanks[N]] = N2;

                 //find the local index of this vertex in the neib cell
                 jj=0;
                 while(CurrVert != (coll->GetCell(N2))->GetVertex(jj)) jj++;
                 HaloCellLocVertNo[VertNeibRanks[N]] = jj;

                 VertNeibRanks[N]++;
                 VertNeibRanks[N + VertNeibRanks[N]] = VertexcellID;

                 if(MaxRankPerV<VertNeibRanks[N])
                   MaxRankPerV=VertNeibRanks[N];
                } //   if(UPDATE)
              } //  for(ii=1;ii<=N_CellsInT

           N_SubDomIn_a = VertNeibRanks[N];
           Temp = VertNeibRanks+(N+1);
 
    
           /** set only one cell (lowest index cell) from each subdomain as the Halo cell for this vertex */
           CurrVert->SetSubDomainInfo(N_SubDomIn_a,  Temp, HaloCellGlobalNo, HaloCellLocVertNo);

            //reset
            for(jj=0;jj<size;jj++)
             HaloCellIndex[jj] = -1;

           } // for(j=0;j<N_V
       } // if(ID==ra
     } //  for(i=0;i<N_Cel 
    
    MaxRankPerV++;    // including own rank
     
    delete [] HaloCellIndex;
    delete [] HaloCellGlobalNo;
    delete [] HaloCellLocVertNo;

  /** *************************************************/
  /** STEP 10 : Get own SubDomain collection of cells */
  /** *************************************************/ 
     N_LocalCells = 0;
     N_OwnCells = 0;

     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
       ID = cell->GetSubDomainNo();

       if(ID==rank || cell->IsHaloCell() )
        N_LocalCells++;

       if(ID==rank)
        N_OwnCells++;
      } //  for (i=0;i<N_Cells;i++)

      
   //  printf("Rank: %d N_own Cells %d\n ", rank,  N_OwnCells);   
      
      
     if(N_LocalCells)
      {
       // collect the own collection
       SubDomainCells = new TBaseCell*[N_LocalCells];
       GlobalCellIndex = new int[N_LocalCells];
      }
     else
      {
       SubDomainCells = NULL;
       GlobalCellIndex = NULL;
      }


      
  /** *************************************************************/
  /** STEP 11 : Fill subdomain own cells and their global numbers */
  /** *************************************************************/ 
  VertexDel = new TVertex*[N_AllLocVert];
  EdgeDel = new TEdge*[N_Cells*12]; // max 12 in Hexa
  CellDel = new TBaseCell*[N_Cells];
     N=0;
     M=0;
     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
       ID = cell->GetSubDomainNo();

       // in the collection, first own cells then Halo cells
       if(ID==rank)
        {
         SubDomainCells[N] =  cell;
         GlobalCellIndex[N] =  i;
         N++;
        }
      else if(cell->IsHaloCell() )
        {
         SubDomainCells[N_OwnCells + M] =  cell;
         GlobalCellIndex[N_OwnCells + M] =  i;
         M++;
        }
      else //cell no longer needed for this rank, delete!!!!!
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
          EdgeDel[N_EdgeDel++] = cell->GetEdge(j);

        for(j=0;j<N_VertInCell;j++)
         VertexDel[N_VertexDel++] = cell->GetVertex(j);

         CellDel[N_CellDel++] = cell;
       }
      } //  for (i=0;i<N_Cells;i++)

  /** ***********************************************************************************/
  /** STEP 12 : Delete other processors cells from own collection, excluding halo cells */
  /** ***********************************************************************************/ 
    // 
/*    if(N_VertexDel)
     Sort(VertexDel, N_VertexDel); 
    if(N_EdgeDel)
     Sort(EdgeDel, N_EdgeDel);*/
//     if(N_CellDel)
//      Sort(CellDel, N_CellDel);

  /*  Last=NULL;
    for(i=0;i<N_VertexDel;i++)
     {
      if( VertexDel[i] != Last)
       {
        delete VertexDel[i];
        Last = VertexDel[i];
       }
     }
   delete [] VertexDel;


    LastEdge=NULL;
    for(i=0;i<N_EdgeDel;i++)
     {
      if( EdgeDel[i] != LastEdge)
       {
        delete EdgeDel[i];
        LastEdge = EdgeDel[i];
       }
     }
   delete [] EdgeDel;

    LastCell=NULL;
    for(i=0;i<N_CellDel;i++)
     {
      if( CellDel[i] != LastCell)
       {
        delete CellDel[i];
        LastCell = CellDel[i];
       }
     }
   delete [] CellDel;

    if((N+M)!=N_LocalCells)
     {
      printf("Error in mesh partition %d \n", rank);
      MPI_Abort(comm, 0);
     } 
     
   /** ***********************************************************************************/
  /** STEP 13 : Edges will be flaged as SubDomainEdges and subdomain infos are filled    */
  /** ************************************************************************************/ 
  
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         {
          edge = cell->GetEdge(j);
          if(edge->GetClipBoard()==-1)
           {
            edge->SetClipBoard(5);
            edge->InitSubDomainInfo(rank);
           } // if(edge->GetCli
         } // for(j=0;j<N_Edges;
       } // if(cell->IsDependentCell
     } // for(i=0;i<N_Cells;i++

              
  
 /** change all subdomain joint from JointEqN into TSubDomainJoint */
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      /**run only through Dependent cells */
      if(cell->IsDependentCell()) 
       {
        GblCellNr = GlobalCellIndex[i];
        ID = cell->GetSubDomainNo();

	//if( GblCellNr == 111615)
          // printf("rank %d local cell no %d\n", rank, i);
	
        /** change the faces, also edges in 3D */
        for(j=0;j<N_JointsInCell;j++)
         {
          Joint = cell->GetJoint(j);

          if(Joint->GetType() == JointEqN)
           {
            neib_cell = Joint->GetNeighbour(cell);
            Neib_ID = neib_cell->GetSubDomainNo();


            if(ID!=Neib_ID)
             {
              GlobCellNo =  neib_cell->GetGlobalCellNo();

              // this joint belongs to two SubDomains
              cell->SetAsSubDomainInterfaceCell();

              // find neib cell local face number
              m=0;
              while(Joint!=neib_cell->GetJoint(m)) m++;

              delete Joint;
              NewJoint = new TSubDomainJoint(cell, neib_cell, Neib_ID, GlobCellNo, m);
              cell->SetJoint(j, NewJoint);
              neib_cell->SetJoint(m, NewJoint);
              NewJoint->SetMapType();

              // set all edges in this face as SubDomainEdges
              ShapeDesc = cell->GetShapeDesc();
              ShapeDesc->GetFaceEdge(TmpFV, TmpLen, MaxLen);
              N_FaceEdges = TmpLen[j];
              for(n=0;n<N_FaceEdges;n++)
               {
                edge = cell->GetEdge(TmpFV[j*MaxLen+n]);
                edge->SetAsNotCrossEdgeFor(rank, Neib_ID);
               } // for(n=0;n<N_Edge

              // remove the Neib_ID from VertNeibRanks
              // these face vertices cannot be cross vertices for processors ID and Neib_ID
              ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
              for(ii=0;ii<TmpLen[j];ii++)
              {
               jj = TmpFV[j*MaxLen + ii];
               M = VertexNumbers[GblCellNr*N_VertInCell + jj] *MaxCpV;
               N = VertNeibRanks[M];

               for(jj=1;jj<=N;jj++)
                if(VertNeibRanks[M + jj] == Neib_ID )
                 {
                  VertNeibRanks[M + jj] = VertNeibRanks[M + N];
                  VertNeibRanks[M]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)
              } // for(ii=0;ii<TmpL
             }//if(ID!=Neib_ID
           } // if(Joint->GetT
         }// for (j=0;j< ;
       } // if(cell->IsDependentCell() 
     }// for(i=0;i<N_Cel

     
    /** find cross edges (if any) */
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        M = GlobalCellIndex[i];
        N_Edges=cell->GetN_Edges();
        ShapeDesc= cell->GetShapeDesc();
        ShapeDesc->GetEdgeVertex(EdgeVertex);

        for(j=0;j<N_Edges;j++)
         {
          edge = cell->GetEdge(j);

          if(edge->GetClipBoard()==-1 && edge->IsSubDomainEdge())
           {
            edge->SetClipBoard(5);
            N_CrossNeibs = edge->GetN_CrossNeibs();

            if(N_CrossNeibs==0)
             continue;

            cell->SetAsCrossEdgeCell();

            a = EdgeVertex[2*j];
            b = EdgeVertex[2*j+ 1];

            Vert_a = cell->GetVertex(a);
            edge->SetCrossNeibInfo(Vert_a);
            edge->GetCrossEdgeNeibs(N_CrossNeibs, CrossNeibsRank);

            M1 = VertexNumbers[M*N_VertInCell + a]*MaxCpV;
            N1 = VertNeibRanks[M1];

            M2 = VertexNumbers[M*N_VertInCell + b]*MaxCpV;
            N2 = VertNeibRanks[M2];

            for(jj=0;jj<N_CrossNeibs;jj++)
             {
              Neib_ID = CrossNeibsRank[jj];

               for(kk=1;kk<=N1;kk++)
                if(VertNeibRanks[M1 + kk] == Neib_ID )
                 {
                  VertNeibRanks[M1 + kk] = VertNeibRanks[M1 + N1];
                  VertNeibRanks[M1]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)

               for(kk=1;kk<=N2;kk++)
                if(VertNeibRanks[M2 + kk] == Neib_ID )
                 {
                  VertNeibRanks[M2 + kk] = VertNeibRanks[M2 + N2];
                  VertNeibRanks[M2]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)
             } // for(jj=0;jj<N_CrossNei
           } //   if(edge->GetClipBoard()==-1 && edge->IsSubDomainEdge())
          } // for(j=0;j<N_Edge
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel
         
    /** find cross vertex (if any) */
    // set the clipboard
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

    /** set incident cell list for all vertices */
    IncCells = new TBaseCell*[MaxCpV];
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      M = GlobalCellIndex[i];    
     
      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        {
         // set the vertexcells for all vertices in dep. cells
         CurrVert = cell->GetVertex(j);
         if(CurrVert->GetClipBoard() != -1)
          continue;

         CurrVert->SetClipBoard(5);

         N = VertexNumbers[M*N_VertInCell + j] *MaxCpV;
         N_CellsInThisVert =  PointNeighb[N];

         for(ii=1;ii<=N_CellsInThisVert;ii++)
          IncCells[ii-1] = coll->GetCell(PointNeighb[N + ii]);

         CurrVert->SetVertexCells(N_CellsInThisVert, IncCells);
        } //  for (j=0;j<N_VertInCell;
     } //  for(i=0;i<N_OwnCell
    delete [] IncCells;

     // set the clipboard
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        M = GlobalCellIndex[i];

        for(j=0;j<N_VertInCell;j++)
         {
          CurrVert = cell->GetVertex(j);

          // continue if the vert is already handled or not a subdomain vert
          if( (CurrVert->GetClipBoard() != -1) || !(CurrVert->IsSubDomainVert()))
           continue;


          CurrVert->SetClipBoard(5);

          M1 = VertexNumbers[M*N_VertInCell + j] *MaxCpV ;
          N1 = VertNeibRanks[M1];

          if(N1!=0)
           cell->SetAsCrossVertexCell();

          for(ii=1;ii<=N1;ii++)
            CurrVert->AddCrossNeib(VertNeibRanks[M1+ii]);

//            TDatabase::ParamDB->Par_P5 =0;
          } // for(j=0;j<N_VertInCell;j
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel


    if(N_LocalCells)
     { Domain->ReplaceTreeInfo(N_LocalCells, SubDomainCells, GlobalCellIndex, N_OwnCells); }
    else
     { Domain->SetN_OwnCells(N_OwnCells); }

        // initialize iterators
    TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
    TDatabase::IteratorDB[It_LE]->SetParam(Domain);
    TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
    TDatabase::IteratorDB[It_Between]->SetParam(Domain);
    TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   delete [] VertNeibRanks;

//    delete [] Vert_Rank; // for more no. processor showing double free or corruption 
   delete [] VertexNumbers;
   delete [] Cell_Rank;
   delete [] PointNeighb;

   
   //Barrier is needed, before calling FECommunicator, since neib process must have all info
   MPI_Barrier(comm);
   
  if(rank==out_rank)
     printf("Mesh Partition end rank %d \n" , rank);

//      MPI_Finalize();
//     exit(0);  
}
#endif // 3D

void Domain_Crop(MPI_Comm comm, TDomain *Domain)
{
  
#if 1
  //Variable list
  
 idxtype *Cell_Rank, *Vert_Rank;
 int i, j, k, m, M, n, rank, size, N_Cells, N, ID, Neib_ID, out_rank= TDatabase::ParamDB->Par_P0;
 int *GlobalCellIndex, N_RootVertices, N_VertInCell, N_JointsInCell, N_AllLocVert;
 int *VertexNumbers, *PointNeighb, maptype, MaxCpV;
 int MaxLen, N_Edges, N_FaceEdges, ii, jj, GblCellNr;
 int N_EdgeDel=0, N_VertexDel=0, N_CellDel=0;
 const int *TmpFV, *TmpLen, *TmpVE, *EdgeVertex, *NeibEdgeVertex;

 int MaxRankPerV,VertexCellNo, VertexcellID, N_CellsInThisVert, *VertNeibRanks, *NumberVertex;
 int N_LocalCells, N_OwnCells, N_OwnIncidentCells, N_NeibIncidentCells;
 int a, b, N_SubDomInThisVert, *HaloCellIndex;
 int M1, M2, N1, N2, N_CellsIn_a, N_CellsIn_b;
 int N_CrossNeibs, *CrossNeibsRank, *HaloCellGlobalNo, *HaloCellLocVertNo;
 int N_SubDomIn_a, N_SubDomIn_b, *Temp, kk, GlobCellNo, test_b, EdgeCellID;

 bool UPDATE, UPDATE_1;

 TVertex *Vert_a, **VertexDel, *Last;
 TCollection *coll;
 TBaseCell *cell, *neib_cell, *Vertexcell, **SubDomainCells,**OwnCells, *cell_a, *cell_b, **IncCells;
 TBaseCell **CellDel, *LastCell;
 TJoint *Joint, *NewJoint;
 TEdge *edge, **EdgeDel, *LastEdge;
 TShapeDesc *ShapeDesc, *NeibShapeDesc;
 TVertex *CurrVert, *NeibVert_a;
  
 MPI_Status status, status1, status2, status3;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
 
 // root will not take part in computations
  MaxCpV = 0;
  MaxRankPerV = -1;

  coll = Domain->GetCollection(It_Finest, 0);
  N_Cells = coll->GetN_Cells();
  cell = coll->GetCell(0);
  N_VertInCell = cell->GetN_Vertices();
  N_JointsInCell = cell->GetN_Joints();
  N_AllLocVert = N_VertInCell*N_Cells;
  VertexNumbers= new int[N_AllLocVert];
  Vert_Rank = new idxtype[N_AllLocVert];
  Cell_Rank = new idxtype[N_Cells];

  VertexDel = new TVertex*[N_AllLocVert];
  EdgeDel = new TEdge*[N_Cells*12]; // max 12 in Hexa
  CellDel = new TBaseCell*[N_Cells];

   TVertex **Vertices;
   
#endif
   
  
   
   int  m1, edgecut=0;
   int etype, numflag=0; // c-style numbering (array starting with 0)
   int type=TDatabase::ParamDB->Par_P2;

   double t1, t2;

    NumberVertex=new int[N_AllLocVert];
    Vertices=new TVertex*[N_AllLocVert];

    N = 0;
   /** *********************************************/
    /** STEP 1 : STORE VERTEX POINTERS CELL-WISE */
    /** *********************************************/
    
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      for(j=0;j<N_VertInCell;j++)
       {
        Vertices[N]=cell->GetVertex(j);
        N++;
       }
     } 
  /** *********************************************/
   /** STEP 2 : SORT THE VERTICES ARRAY */ 
  /** *********************************************/
    Sort(Vertices, N);

  /** ***************************************************/
    /**STEP 3: STORE THE SORTED POINTER ARRAY AS INDICES */
  /** ***************************************************/
    Last=NULL;
    N_RootVertices=-1;
    for(i=0;i<N_AllLocVert;i++)
     {
      if((CurrVert=Vertices[i])!=Last)
       {
        N_RootVertices++;
        Last=CurrVert;
       }
      NumberVertex[i]=N_RootVertices;
     }
    N_RootVertices++;

 /** *********************************************/   
  /** STEP 4 : STORE THE INDICES CELL-WISE */
 /** *********************************************/ 
    m=0;
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);

      for(j=0;j<N_VertInCell;j++)
       {
        CurrVert=cell->GetVertex(j);
        N=GetIndex(Vertices, N_AllLocVert, CurrVert);
        VertexNumbers[m]=NumberVertex[N];
        m++;
       } 
     } 


  /** *********************************************/
   /** STEP 6 : TO FIND MAXIMUM CELLS PER VERTEX */
  /** *********************************************/ 
  
  PointNeighb = new int[N_RootVertices];
  memset(PointNeighb, 0, N_RootVertices*SizeOfInt);
  
  for (i=0;i<N_AllLocVert;i++)
    PointNeighb[VertexNumbers[i]]++;
  
  for (i=0;i<N_RootVertices;i++)
    if (PointNeighb[i] > MaxCpV) MaxCpV = PointNeighb[i];
  
   MaxRankPerV = MaxCpV;
   MaxCpV++; // ACCOUNTING FOR AN EXTRA COLUMN GIVING INFO OF NUMBER OF CELLS SHARING THE VERTEX
   
   delete [] PointNeighb;
  
  /** ********************************************************************************************/
  /** STEP 7 : CREATE AN ARRAY CONTAINING INFO REGARDING THE NEIGHBOURS OF VERTICES (CELL INDEX) */
  /** ********************************************************************************************/
   
   PointNeighb = new int[N_RootVertices*MaxCpV] ;  
   memset(PointNeighb, 0, (N_RootVertices*MaxCpV)*SizeOfInt);

   for(i=0;i<N_Cells;i++)
    {
     for(j=0;j<N_VertInCell;j++)
      {
       M = VertexNumbers[i*N_VertInCell + j] *MaxCpV ;
       PointNeighb[M]++;
       PointNeighb[M + PointNeighb[M]] = i;
      }  
    }
    
   delete [] NumberVertex;
   delete [] Vertices;
   


    /**Metis partition done!! code for all processors */
    //MPI_Bcast(Cell_Rank, N_Cells, MPI_INT, 0, comm);
    

/********** VARIABLES FOR MULTIGRID ****************/	

  int Nchildren,childn,parentglobalno;

 
  Nchildren = coll->GetCell(0)->GetParent()->GetN_Children();

  for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);	
      cell->SetSubDomainNo(cell->GetParent()->GetSubDomainNo());
      childn = cell->GetParent()->GetChildNumber(cell);	
      parentglobalno = cell->GetParent()->GetGlobalCellNo();	
      cell->SetGlobalCellNo(parentglobalno*Nchildren+childn);
      cell->SetLocalCellNo(i);
      Cell_Rank[i] = cell->GetSubDomainNo();
     }
  
  /** *********************************************/
  /** STEP 8.2 : SET ALL VERTICES TO -1 */
  /** *********************************************/ 
    
    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);

      if(cell->GetSubDomainNo()==rank) 
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

  /** *********************************************/
  /** STEP 9 : Fill the VertNeibRanks info        */
  /** first column contains how many ranks contain this vertex further columns contain the rank ID of the subdomains */
  /** *********************************************/ 
    VertNeibRanks = new int[N_RootVertices*MaxCpV];
    memset(VertNeibRanks, 0, N_RootVertices*MaxCpV*SizeOfInt);

    HaloCellIndex = new int[size];
    for(i=0;i<size;i++)
      HaloCellIndex[i] = -1;

    HaloCellGlobalNo = new int[MaxCpV];
    HaloCellLocVertNo = new int[MaxCpV];

    for(i=0;i<N_Cells;i++)
     {
      cell = coll->GetCell(i);
      ID = cell->GetSubDomainNo();

  /**run only through own cells */
      if(ID==rank) 
       {
        cell->SetAsOwnCell();
      //set SubDomainVert if any vert in this cell is so needed for moving meshes and setting cross vertex
        for(j=0;j<N_VertInCell;j++)
         {
          CurrVert = cell->GetVertex(j);
	  
	  if(CurrVert->GetClipBoard() != -1)
           continue;

          CurrVert->SetClipBoard(5);
          N = VertexNumbers[i*N_VertInCell + j]*MaxCpV ;

          N_CellsInThisVert =  PointNeighb[N];
          N_OwnIncidentCells = 0;
          N_NeibIncidentCells = 0;

          // check! any subdomain cell containg this vert has to be added
          // Setting all halo cells
          UPDATE = TRUE;
          for(ii=1;ii<=N_CellsInThisVert;ii++)
           {
            VertexCellNo = PointNeighb[N + ii];
            Vertexcell = coll->GetCell(VertexCellNo);
            VertexcellID = Vertexcell->GetSubDomainNo();

            if(VertexcellID != rank)
             {
              //vertex belong to diff subdomain 
              if(UPDATE)
               {
                cell->SetAsDependentCell();
		
                CurrVert->SetAsSubDomainVert();
                HaloCellIndex[VertexcellID]= VertexCellNo;
                UPDATE = FALSE;
               }
              Vertexcell->SetAsHaloCell();

               /** changed on 4 Feb 2012 - by Sashi */
               if(HaloCellIndex[VertexcellID] == -1)
                HaloCellIndex[VertexcellID]= VertexCellNo;

             }
            } //  for(ii=1;ii< 

           if(CurrVert->IsSubDomainVert())
             for(ii=1;ii<=N_CellsInThisVert;ii++)
              {
               VertexCellNo = PointNeighb[N + ii];
               Vertexcell = coll->GetCell(VertexCellNo);
               VertexcellID = Vertexcell->GetSubDomainNo();

               // all cells associated with this vert are DepCells
               Vertexcell->SetAsDependentCell();

               if(VertexcellID==rank)
                continue;

               N_CellsIn_b = VertNeibRanks[N];
               UPDATE = TRUE;

               for(jj=1;jj<=N_CellsIn_b;jj++)
                if( VertNeibRanks[N + jj]==VertexcellID)
                 {
                  UPDATE = FALSE;
                  break;
                 }

               if(UPDATE)
                {
                 N2 = HaloCellIndex[VertexcellID];
                 HaloCellGlobalNo[VertNeibRanks[N]] = N2;

                 //find the local index of this vertex in the neib cell
                 jj=0;
                 while(CurrVert != (coll->GetCell(N2))->GetVertex(jj)) jj++;
                 HaloCellLocVertNo[VertNeibRanks[N]] = jj;

                 VertNeibRanks[N]++;
                 VertNeibRanks[N + VertNeibRanks[N]] = VertexcellID;

                 if(MaxRankPerV<VertNeibRanks[N])
                   MaxRankPerV=VertNeibRanks[N];
                } //   if(UPDATE)
              } //  for(ii=1;ii<=N_CellsInT

           N_SubDomIn_a = VertNeibRanks[N];
           Temp = VertNeibRanks+(N+1);
 
    
           /** set only one cell (lowest index cell) from each subdomain as the Halo cell for this vertex */
           CurrVert->SetSubDomainInfo(N_SubDomIn_a,  Temp, HaloCellGlobalNo, HaloCellLocVertNo);

            //reset
            for(jj=0;jj<size;jj++)
             HaloCellIndex[jj] = -1;

           } // for(j=0;j<N_V
       } // if(ID==ra
     } //  for(i=0;i<N_Cel 
    
    MaxRankPerV++;    // including own rank

    delete [] HaloCellIndex;
    delete [] HaloCellGlobalNo;
    delete [] HaloCellLocVertNo;

  /** *************************************************/
  /** STEP 10 : Get own SubDomain collection of cells */
  /** *************************************************/ 
     N_LocalCells = 0;
     N_OwnCells = 0;

     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
       ID = cell->GetSubDomainNo();

       if(ID==rank || cell->IsHaloCell() )
        N_LocalCells++;

       if(ID==rank)
        N_OwnCells++;
      } //  for (i=0;i<N_Cells;i++)

      
    // printf("Rank: %d N_own Cells %d\n ", rank,  N_OwnCells);   
      

     if(N_LocalCells)
      {
       // collect the own collection
       SubDomainCells = new TBaseCell*[N_LocalCells];
       GlobalCellIndex = new int[N_LocalCells];
      }
     else
      {
       SubDomainCells = NULL;
       GlobalCellIndex = NULL;
      }


      
  /** *************************************************************/
  /** STEP 11 : Fill subdomain own cells and their global numbers */
  /** *************************************************************/ 
  VertexDel = new TVertex*[N_AllLocVert];
  EdgeDel = new TEdge*[N_Cells*12]; // max 12 in Hexa
  CellDel = new TBaseCell*[N_Cells];
     N=0;
     M=0;
     for(i=0;i<N_Cells;i++)
      {
       cell = coll->GetCell(i);
       ID = cell->GetSubDomainNo();

       // in the collection, first own cells then Halo cells
       if(ID==rank)
        {
         SubDomainCells[N] =  cell;
         GlobalCellIndex[N] =  i;
         N++;
        }
      else if(cell->IsHaloCell() )
        {
         SubDomainCells[N_OwnCells + M] =  cell;
         GlobalCellIndex[N_OwnCells + M] =  i;
         M++;
        }
      else //cell no longer needed for this rank, delete!!!!!
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
          EdgeDel[N_EdgeDel++] = cell->GetEdge(j);

        for(j=0;j<N_VertInCell;j++)
         VertexDel[N_VertexDel++] = cell->GetVertex(j);

         CellDel[N_CellDel++] = cell;
       }
      } //  for (i=0;i<N_Cells;i++)

      /** ***********************************************************************************/
  /** STEP 12 : Delete other processors cells from own collection, excluding halo cells */
  /** ***********************************************************************************/ 
    // 
/*    if(N_VertexDel)
     Sort(VertexDel, N_VertexDel); 
    if(N_EdgeDel)
     Sort(EdgeDel, N_EdgeDel);*/
//     if(N_CellDel)
//      Sort(CellDel, N_CellDel);

  /*  Last=NULL;
    for(i=0;i<N_VertexDel;i++)
     {
      if( VertexDel[i] != Last)
       {
        delete VertexDel[i];
        Last = VertexDel[i];
       }
     }
   delete [] VertexDel;


    LastEdge=NULL;
    for(i=0;i<N_EdgeDel;i++)
     {
      if( EdgeDel[i] != LastEdge)
       {
        delete EdgeDel[i];
        LastEdge = EdgeDel[i];
       }
     }
   delete [] EdgeDel;

    LastCell=NULL;
    for(i=0;i<N_CellDel;i++)
     {
      if( CellDel[i] != LastCell)
       {
        delete CellDel[i];
        LastCell = CellDel[i];
       }
     }
   delete [] CellDel;

    if((N+M)!=N_LocalCells)
     {
      printf("Error in mesh partition %d \n", rank);
      MPI_Abort(comm, 0);
     } 
     
   /** ***********************************************************************************/
  /** STEP 13 : Edges will be flaged as SubDomainEdges and subdomain infos are filled    */
  /** ************************************************************************************/ 
  
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         {
          edge = cell->GetEdge(j);
          if(edge->GetClipBoard()==-1)
           {
            edge->SetClipBoard(5);
            edge->InitSubDomainInfo(rank);
           } // if(edge->GetCli
         } // for(j=0;j<N_Edges;
       } // if(cell->IsDependentCell
     } // for(i=0;i<N_Cells;i++

              
  
 /** change all subdomain joint from JointEqN into TSubDomainJoint */
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      /**run only through Dependent cells */
      if(cell->IsDependentCell()) 
       {
        GblCellNr = GlobalCellIndex[i];
        ID = cell->GetSubDomainNo();

	//if( GblCellNr == 111615)
          // printf("rank %d local cell no %d\n", rank, i);
	
        /** change the faces, also edges in 3D */
        for(j=0;j<N_JointsInCell;j++)
         {
          Joint = cell->GetJoint(j);

          if(Joint->GetType() == JointEqN)
           {
            neib_cell = Joint->GetNeighbour(cell);
            Neib_ID = neib_cell->GetSubDomainNo();


            if(ID!=Neib_ID)
             {
              GlobCellNo =  neib_cell->GetGlobalCellNo();

              // this joint belongs to two SubDomains
              cell->SetAsSubDomainInterfaceCell();

              // find neib cell local face number
              m=0;
              while(Joint!=neib_cell->GetJoint(m)) m++;

              delete Joint;
              NewJoint = new TSubDomainJoint(cell, neib_cell, Neib_ID, GlobCellNo, m);
              cell->SetJoint(j, NewJoint);
              neib_cell->SetJoint(m, NewJoint);
              NewJoint->SetMapType();

              // set all edges in this face as SubDomainEdges
              ShapeDesc = cell->GetShapeDesc();
              ShapeDesc->GetFaceEdge(TmpFV, TmpLen, MaxLen);
              N_FaceEdges = TmpLen[j];
              for(n=0;n<N_FaceEdges;n++)
               {
                edge = cell->GetEdge(TmpFV[j*MaxLen+n]);
                edge->SetAsNotCrossEdgeFor(rank, Neib_ID);
               } // for(n=0;n<N_Edge

              // remove the Neib_ID from VertNeibRanks
              // these face vertices cannot be cross vertices for processors ID and Neib_ID
              ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
              for(ii=0;ii<TmpLen[j];ii++)
              {
               jj = TmpFV[j*MaxLen + ii];
               M = VertexNumbers[GblCellNr*N_VertInCell + jj] *MaxCpV;
               N = VertNeibRanks[M];

               for(jj=1;jj<=N;jj++)
                if(VertNeibRanks[M + jj] == Neib_ID )
                 {
                  VertNeibRanks[M + jj] = VertNeibRanks[M + N];
                  VertNeibRanks[M]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)
              } // for(ii=0;ii<TmpL
             }//if(ID!=Neib_ID
           } // if(Joint->GetT
         }// for (j=0;j< ;
       } // if(cell->IsDependentCell() 
     }// for(i=0;i<N_Cel

     
    /** find cross edges (if any) */
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        M = GlobalCellIndex[i];
        N_Edges=cell->GetN_Edges();
        ShapeDesc= cell->GetShapeDesc();
        ShapeDesc->GetEdgeVertex(EdgeVertex);

        for(j=0;j<N_Edges;j++)
         {
          edge = cell->GetEdge(j);

          if(edge->GetClipBoard()==-1 && edge->IsSubDomainEdge())
           {
            edge->SetClipBoard(5);
            N_CrossNeibs = edge->GetN_CrossNeibs();

            if(N_CrossNeibs==0)
             continue;

            cell->SetAsCrossEdgeCell();

            a = EdgeVertex[2*j];
            b = EdgeVertex[2*j+ 1];

            Vert_a = cell->GetVertex(a);
            edge->SetCrossNeibInfo(Vert_a);
            edge->GetCrossEdgeNeibs(N_CrossNeibs, CrossNeibsRank);

            M1 = VertexNumbers[M*N_VertInCell + a]*MaxCpV;
            N1 = VertNeibRanks[M1];

            M2 = VertexNumbers[M*N_VertInCell + b]*MaxCpV;
            N2 = VertNeibRanks[M2];

            for(jj=0;jj<N_CrossNeibs;jj++)
             {
              Neib_ID = CrossNeibsRank[jj];

               for(kk=1;kk<=N1;kk++)
                if(VertNeibRanks[M1 + kk] == Neib_ID )
                 {
                  VertNeibRanks[M1 + kk] = VertNeibRanks[M1 + N1];
                  VertNeibRanks[M1]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)

               for(kk=1;kk<=N2;kk++)
                if(VertNeibRanks[M2 + kk] == Neib_ID )
                 {
                  VertNeibRanks[M2 + kk] = VertNeibRanks[M2 + N2];
                  VertNeibRanks[M2]--;
                  break;
                 } // for(jj=1;jj<=N;jj++)
             } // for(jj=0;jj<N_CrossNei
           } //   if(edge->GetClipBoard()==-1 && edge->IsSubDomainEdge())
          } // for(j=0;j<N_Edge
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel
         
    /** find cross vertex (if any) */
    // set the clipboard
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

    /** set incident cell list for all vertices */
    IncCells = new TBaseCell*[MaxCpV];
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      M = GlobalCellIndex[i];    
     
      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        {
         // set the vertexcells for all vertices in dep. cells
         CurrVert = cell->GetVertex(j);
         if(CurrVert->GetClipBoard() != -1)
          continue;

         CurrVert->SetClipBoard(5);

         N = VertexNumbers[M*N_VertInCell + j] *MaxCpV;
         N_CellsInThisVert =  PointNeighb[N];

         for(ii=1;ii<=N_CellsInThisVert;ii++)
          IncCells[ii-1] = coll->GetCell(PointNeighb[N + ii]);

         CurrVert->SetVertexCells(N_CellsInThisVert, IncCells);
        } //  for (j=0;j<N_VertInCell;
     } //  for(i=0;i<N_OwnCell
    delete [] IncCells;

     // set the clipboard
    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];

      if(cell->IsDependentCell())
       for (j=0;j<N_VertInCell;j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
     }

    for(i=0;i<N_OwnCells;i++)
     {
      cell = SubDomainCells[i];
      if(cell->IsDependentCell())
       {
        M = GlobalCellIndex[i];

        for(j=0;j<N_VertInCell;j++)
         {
          CurrVert = cell->GetVertex(j);

          // continue if the vert is already handled or not a subdomain vert
          if( (CurrVert->GetClipBoard() != -1) || !(CurrVert->IsSubDomainVert()))
           continue;


          CurrVert->SetClipBoard(5);

          M1 = VertexNumbers[M*N_VertInCell + j] *MaxCpV ;
          N1 = VertNeibRanks[M1];

          if(N1!=0)
           cell->SetAsCrossVertexCell();

          for(ii=1;ii<=N1;ii++)
            CurrVert->AddCrossNeib(VertNeibRanks[M1+ii]);

//            TDatabase::ParamDB->Par_P5 =0;
          } // for(j=0;j<N_VertInCell;j
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

    if(N_LocalCells)
     { Domain->ReplaceTreeInfo(N_LocalCells, SubDomainCells, GlobalCellIndex, N_OwnCells); }
    else
     { Domain->SetN_OwnCells(N_OwnCells); }

        // initialize iterators
    TDatabase::IteratorDB[It_EQ]->SetParam(Domain);
    TDatabase::IteratorDB[It_LE]->SetParam(Domain);
    TDatabase::IteratorDB[It_Finest]->SetParam(Domain);
    TDatabase::IteratorDB[It_Between]->SetParam(Domain);
    TDatabase::IteratorDB[It_OCAF]->SetParam(Domain);

   delete [] VertNeibRanks;

coll = Domain->GetCollection(It_Finest, 0);
 for(i=0;i<N_LocalCells;i++)
   {  
     cell = coll->GetCell(i);
     cell->SetCellIndex(i);
   }
   
 
  delete [] VertexNumbers;
  //delete [] Vert_Rank;
  delete [] Cell_Rank;
  delete [] PointNeighb;
 
   //Barrier is needed, before calling FECommunicator, since neib process must have all info
   MPI_Barrier(comm);
} // Partition_Mesh2D()

#endif //mpi
