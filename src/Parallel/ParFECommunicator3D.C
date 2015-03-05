// =======================================================================
// @(#)ParFECommunicator3D.C
//
// Class:      TParFECommunicator3D
// Purpose:    Class containing all info needed for communication between subdomains
//
// Author:     Sashikumaar Ganesan (01.10.09)
//
// History:    Start of implementation 01.10.09 (Sashikumaar Ganesan)
//
// ======================================================================= 
#ifdef _MPI

#include "mpi.h"
#include <ParFECommunicator3D.h>
#include <FEDatabase3D.h>
#include <Database.h>
#include <SubDomainJoint.h>
#include <Edge.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>


#define GLOBAL_NO 0
#define DOF_NO 1

#define HALOCELL 0
#define NONHALO  1

TParFECommunicator3D::TParFECommunicator3D(MPI_Comm comm, TFESpace3D *fespace, TSquareStructure3D* Sqstruct)
{
 int N_U;

 Comm    = comm;
 FESpace = fespace;
 sqstruct = Sqstruct;
 MaxSubDomainPerDof = fespace->GetMaxSubDomainPerDof();

 if(MaxSubDomainPerDof<0)
  {
   printf("Error: SetMaxSubDomainPerDof in FeSpace before calling ParFECommunicator2D \n");
   MPI_Finalize();
   exit(0);
  }

  N_U = FESpace->GetN_DegreesOfFreedom();

  MPI_Allreduce(&N_U, &MaxN_LocalDofAllRank, 1, MPI_INT, MPI_MAX, Comm);

 // ConstructDofRankIndex();
  
 
   if(TDatabase::ParamDB->SC_SMOOTHER_SCALAR==5)
   {
    ConstructDofMapRe();
#ifdef _HYBRID
    /**    */
    ColorAndReorder(0,N_Master,N_CMaster,ptrCMaster);
    ColorAndReorder(N_Master,N_Master+N_Int,N_CInt,ptrCInt);
    ColorAndReorder(N_Master+N_Int,N_Master+N_Int+N_Dept,N_CDept,ptrCDept);
#endif
    /**    */
    SetNewGN();
   }

   else
      ConstructDofMap();
     
  SetFENeibCommunicationSteps();
 

 //ScheduleParFEComm3D();

 // MapDofFromNeib3D();

 // ConstructGlobalDofFromNeib3D();
}

   int partition( int a[],int b[], int l, int r);

static void quickSort( int a[], int b[],int l, int r)
{
   int j;

   if( l < r ) 
   {
       // divide and conquer
       j = partition( a,b, l, r);
       quickSort( a, b,l, j-1);
       quickSort( a,b, j+1, r);
   }
}



int partition( int a[], int b[], int l, int r) {
   int pivot, i, j, t,t2;
   pivot = a[l];
   i = l; j = r+1;
		
   while( 1)
   {
   	do ++i; while( a[i] <= pivot && i <= r );
   	do --j; while( a[j] > pivot );
   	if( i >= j ) break;
   	t = a[i]; a[i] = a[j]; a[j] = t;
	t2=b[i]; b[i]=b[j]; b[j]=t2;
   }
   t = a[l]; a[l] = a[j]; a[j] = t;
   t2=b[l]; b[l]=b[j]; b[j]=t2;
   return j;
}

static void SortPos(int *Array, int *Array2,int length)
{
  int n=0, l=0, r=length-1, m;
  int i, j, k, *rr, len, s;
  int Mid, Temp, Temp2;
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
      if(m>length || m<0) printf("accessing out of bounds\n");
      Mid=Array[m];

      do
      {
        while(Array[i] > Mid) {i++;
               if(i>length || i<0) printf("accessing out of bounds\n");
	}
        while(Array[j] < Mid) {j--;
	   if(j>length || j<0) printf("accessing out of bounds\n");
	}

        if (i<=j)
        {
	   if(i>length || i<0) printf("accessing out of bounds\n");
	     if(j>length || j<0) printf("accessing out of bounds\n");
          Temp=Array[i];  Temp2=Array2[i];
          Array[i]=Array[j]; Array2[i]=Array2[j];
          Array[j]=Temp; Array2[j]=Temp2;
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

static int GetLocalCellIndex(int N_DependentCells, int *DeptCellGlobalNo, int *DependentCellIndex, int GlobalCellNo)
{
  
  int l=0, r=N_DependentCells, m=(r+l)/2;
  int Mid;

  Mid= DeptCellGlobalNo[m];

  while(Mid!=GlobalCellNo)
  {
    if(Mid>GlobalCellNo)
    {  r=m; }
    else
    { l=m; }

    m=(r+l)/2;
    Mid=DeptCellGlobalNo[m];
  }

    return DependentCellIndex[m];
}

static int GetDeptIndex(int N, int *array, int val)
{

  int l=0, r=N, m=(r+l)/2;
  int Mid;

  Mid= array[m];

  while(Mid != val)
  {
    if(Mid>val)
    {  r=m; }
    else
    { l=m; }

    m=(r+l)/2;
    Mid=array[m];
  }

  return m;
}


/** Construct DofRankIndex based on used FESpace */
void TParFECommunicator3D::ConstructDofRankIndex_old()
{
 int rank, size, i, j, k, l, m, m1, P, N_Cells, N_U, N_LocDof, ID;
 int M, N, N_Vert, N_Joints, N_JointDOF, Neib_ID;
 int *DOF, *GlobalNumbers, *BeginIndex, *JointDof, Disp, N_EdgeDOF, *EdgeDof;
 int N_CrossEdgeNeibs, *CrossEdgeNeibsRank, N_Edges, N_VertInCell;
 int N_VertCrossNeibs, *VertCrossNeibs, VertDof;
 int test, N_Dof;

 double x,y,z;

 bool UPDATE;

 TCollection *Coll;
 TBaseCell *cell, *Neib_cell;
 TVertex *Vert;
 TJoint *Joint;
 FE3D FeId;
 TFEDesc3D *FeDesc;
 TEdge *edge;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_U = FESpace->GetN_DegreesOfFreedom();

  BeginIndex = FESpace->GetBeginIndex();
  GlobalNumbers = FESpace->GetGlobalNumbers();
  N_DofRankIndex = new int[N_U];
  DofRankIndex = new int[MaxSubDomainPerDof*N_U];
  memset(N_DofRankIndex, 0, N_U*SizeOfInt);

  Disp = N_U*MaxSubDomainPerDof;
  for(i=0; i<Disp; i++)
   DofRankIndex[i] = -1;       
  
  N_DependentCells=0;
  /** own dofs */
  // find how many SubDomain contain each dof and their corresponding IDs
  for(i=0; i<N_Cells; i++)
   {
    cell = Coll->GetCell(i);
    ID = cell->GetSubDomainNo();

//     if(rank!=ID)
//      continue;

    DOF = GlobalNumbers + BeginIndex[i];
    N_LocDof = BeginIndex[i+1] - BeginIndex[i];

    for(j=0; j<N_LocDof; j++)
     {
      N = DOF[j];
      M = N_DofRankIndex[N];

      /** adding own process ID */
      UPDATE = TRUE;
      for(k=0; k<M; k++)
      if(ID==DofRankIndex[N*MaxSubDomainPerDof + k ])
       {
        UPDATE = FALSE;
        break;
       } // for(k=0; k<M; k++)

      if(UPDATE)
       {
        DofRankIndex[N*MaxSubDomainPerDof + M ] = ID;
        N_DofRankIndex[N]++;
	
//    if(rank==0 && N==72)
//      printf("%d index %d own N_neibs  %d \n", rank, N,  N_DofRankIndex[N]);
// 
//     if(rank==3  && N==5566)
//      printf("%d index %d own N_neibs %d \n", rank, 5566,  N_DofRankIndex[5566]);  	
// 	
       }
      }  // for(j=0; j<N_DOF; j++)
      
 
    /** adding the neib process ID */
    if(cell->IsDependentCell())
     {
      N_DependentCells++;
      N_Vert= cell->GetN_Vertices();

      for(j=0; j<N_Vert; j++)
       {
        Vert = cell->GetVertex(j);
        Vert->SetClipBoard(-1);
       } // for(j=0; j<N_Vert
      } // if(cell->IsDependen
   } //  for(i=0; i<N_Cells

   
 
   
   
/*    if(rank==2 )
     printf("%d index %d  N_neibs after  %d \n", rank, 5981,  N_DofRankIndex[5981]);

    if(rank==3 )
     printf("%d index %d N_neibs after %d \n", rank, 5566,  N_DofRankIndex[5566]);  	
	*/   
   
  if(N_DependentCells)
   {
    DependentCellIndex = new int[N_DependentCells];
    DeptCellGlobalNo = new int[N_DependentCells];
   }

    /** for cross edge cells set clioboard (if any) */
    N_DependentCells=0;
    for(i=0;i<N_Cells;i++)
     {
      cell = Coll->GetCell(i);
      ID = cell->GetSubDomainNo();

//       if(rank!=ID)
//        continue;

      if(cell->IsDependentCell())
       {
        DependentCellIndex[N_DependentCells] = i;
        DeptCellGlobalNo[N_DependentCells] = cell->GetGlobalCellNo();
        N_DependentCells++;
       }

      if(cell->IsCrossEdgeCell())
       {
        N_Edges=cell->GetN_Edges();

        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);

       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

  /** dofs on faces/joints */
  test=0;
  for(i=0; i<N_Cells; i++)
   {
    cell = Coll->GetCell(i);
    ID = cell->GetSubDomainNo();

//    if(rank==1 &&  cell->GetGlobalCellNo()==39 )
//     for(j=0; j<4; j++)
//      {
//       (cell->GetVertex(j))->GetCoords(x,y,z);
//        printf("ConstructDofRankIndex vert %d j %d x %f y %f z %f\n", i, j, x,y,z);
//      }
    
    
//     if(rank!=ID)
//      continue;

// //     if(!(  i==18|| i==17|| i==16|| i==13|| i==12|| i==11|| i==18))
// //      continue;

//     if(!( i==0 || i==14 || i==8 || i==16 || i==10|| i==18|| i==15|| i==7|| i==18|| i==5))
//      continue;

//  if(!(  i==0 || i==9 || i==10  || i==4 || i==2 || i==1 || i==14  || i==3 || i==11 ))
//      continue;

    if(cell->IsDependentCell())
     {
      DOF = GlobalNumbers + BeginIndex[i];
      N_Joints = cell->GetN_Joints();
      FeId = FESpace->GetFE3D(i, cell);
      FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
      N_JointDOF = FeDesc->GetN_JointDOF();
      N_Dof = FeDesc->GetN_DOF();

//       if(rank==8 && i==5)
//       for(j=0; j<N_Dof; j++)
//        printf("ConstructDofRankIndex  j %d dof %d \n", j, DOF[j]);

//    if(rank==0 &&  cell->GetGlobalCellNo()==35 )
//     for(j=0; j<4; j++)
//      {
//       (cell->GetVertex(j))->GetCoords(x,y,z);
//        printf("ConstructDofRankIndex vert %d j %d x %f y %f z %f\n", i, j, x,y,z);
//      }



     // Discspace (assumed all cells haveing same FEID)
     if(N_JointDOF<=0)
      break;

      for(j=0; j<N_Joints; j++)
       {
        Joint = cell->GetJoint(j);

//  if((i==5 && j==5))
//         if( (i==5 && j==3))
//  if( (i==5 && j==2))
//  if((i==10 && j==4) || (i==4 && j==4))
//         continue;

         // put neib joint subdomain ID in all DOF in the joint
        if(Joint->GetType() == SubDomainJoint)
         {
          Neib_ID = ((TSubDomainJoint *)Joint)->GetNeibRank();
          JointDof = FeDesc->GetJointDOF(j);

            Neib_cell = Joint->GetNeighbour(cell);
//             Neib_ID = Neib_cell->GetSubDomainNo();


          for(k=0; k<N_JointDOF; k++)
           {
            N = DOF[ JointDof[k] ];
            M = N_DofRankIndex[N];

//            if(rank==0 && cell->GetGlobalCellNo()==35 && j== 0)
//            printf("Rank %d DOF %d \n", rank,  N );

// 
//           if(rank==0 && i==169 && j==2)
//            printf("Rank %d DOF %d \n", rank,  N,  Neib_ID );



            UPDATE = TRUE;
            for(l=0; l<M; l++)
            if(Neib_ID==DofRankIndex[N*MaxSubDomainPerDof + l ])
             {

//           if(rank==TDatabase::ParamDB->Par_P5 && Neib_ID==TDatabase::ParamDB->Par_P6 )
//              printf("DofRankIndex cell %d joint %d JointDof  %d Dof %d M %d !!!\n",
//                       i, j, JointDof[k], N, M);

//           if(rank==TDatabase::ParamDB->Par_P5 && Neib_ID==TDatabase::ParamDB->Par_P6 )
//            test++;

              UPDATE = FALSE;
              break;
             } // for(l=0; l<M; l++)

            if(UPDATE)
             {
              DofRankIndex[N*MaxSubDomainPerDof + M ] = Neib_ID;
              N_DofRankIndex[N]++;

//                if(rank==0 && N==66)
//                 printf("%d Glob_CellNo %d face N_DOFNeibs  %d localjointno  %d Neib_ID  %d Neib GCellNo %d\n", rank, cell->GetGlobalCellNo(),  N_DofRankIndex[N], j, Neib_ID, Neib_cell->GetGlobalCellNo());
// 	       
//                if(rank==1 && Neib_cell->GetGlobalCellNo()==39)
//                 printf("%d Glob_CellNo %d face N_DOFNeibs  %d localjointno  %d Neib_ID  %d Neib GCellNo %d\n", rank, i,  N_DofRankIndex[N], j, Neib_ID, Neib_cell->GetGlobalCellNo());
//                if(rank==2  && N==3745)
//                 printf("%d Glob_CellNo %d face N_DOFNeibs  %d localjointno  %d  Neib_ID  %d, Neib GCellNo %d\n", rank, i,  N_DofRankIndex[3745], j, Neib_ID, Neib_cell->GetGlobalCellNo());  	
//      
//               if(rank==TDatabase::ParamDB->Par_P5 && N==19)
//                 printf("Rank %d i %d   j %d  Neib_ID %d \n", rank, i,  j, DofRankIndex[N*MaxSubDomainPerDof + M ]);
             }
           } //for(k=0; k<N_JointDOF; 
//            printf("%d SubDomainJoint Neib_ID %d \n", rank, Neib_ID);
         }// if(Joint->GetT  
       } // for(j=0; j<N_Joints; j+
     }//if(cell->IsDep
   } //  for(i=0; i<N_Cells


//    N=72;
//    if(rank==0)
//       printf("%d index %d N_neibs  %d \n", rank, N,  N_DofRankIndex[N]);
   
/*  MPI_Finalize();
 exit(0);   */   
   
//   M=0;
//   if(rank==TDatabase::ParamDB->Par_P5)
//     for(i=0; i<N_U; i++)
//      for(j=0; j<N_DofRankIndex[i]; j++)
//       if(DofRankIndex[i*MaxSubDomainPerDof + j ]==TDatabase::ParamDB->Par_P6)
//         M++;
// 
//   if(rank==TDatabase::ParamDB->Par_P5)
//    printf("Face  Test %d M %d!!!\n", test, M);
// 
//   M=0;
//   if(rank==TDatabase::ParamDB->Par_P6)
//     for(i=0; i<N_U; i++)
//      for(j=0; j<N_DofRankIndex[i]; j++)
//       if(DofRankIndex[i*MaxSubDomainPerDof + j ]==TDatabase::ParamDB->Par_P5)
//         M++;
// 

//     if(rank==1 )
//      printf("%d index %d  N_neibs face  %d \n", rank, 4038,  N_DofRankIndex[4038]);
// 
//     if(rank==2 )
//      printf("%d index %d N_neibs face   %d \n", rank, 3745,  N_DofRankIndex[3745]);  	
// 	   


  /** dofs on edges */
   // put neib cross edge subdomain ID in all DOF in the edge 
   for(i=0; i<N_Cells; i++)
    {
     cell = Coll->GetCell(i);
     ID = cell->GetSubDomainNo();

//      if(rank!=ID)
//       continue;

     if(cell->IsCrossEdgeCell())
       { 
        DOF = GlobalNumbers + BeginIndex[i];
        N_Edges=cell->GetN_Edges();
        FeId = FESpace->GetFE3D(i, cell);
        FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
        
        if(!(FeDesc->IsEdgeVertData_Filled()) )
        {
           printf("Rank %d, Error! Edge and vertes data are not set in FEdesc3D for this FE : %d\n", rank, FeId);
           MPI_Abort(MPI_COMM_WORLD,  0);  
        }
        
        N_EdgeDOF = FeDesc->GetN_EdgeDOF();

       // Discspace or non-conforming space
       if(N_EdgeDOF<=0)
         break;

        for(j=0;j<N_Edges;j++)
         {
          edge = cell->GetEdge(j);

          if( edge->IsSubDomainCrossEdge() &&  (edge->GetClipBoard()==-1) )
           {
            edge->SetClipBoard(5);
            edge->GetCrossEdgeNeibs(N_CrossEdgeNeibs, CrossEdgeNeibsRank);
            EdgeDof = FeDesc->GetEdgeDOF(j);

            for(k=0; k<N_CrossEdgeNeibs; k++)
             {
              Neib_ID = CrossEdgeNeibsRank[k];

              for(l=0; l<N_EdgeDOF; l++) 
               {
                N = DOF[ EdgeDof[l] ];
                M = N_DofRankIndex[N];

               UPDATE = TRUE;
               for(m=0; m<M; m++)
                if(Neib_ID== DofRankIndex[N*MaxSubDomainPerDof + m])
                 {

//          if(rank==TDatabase::ParamDB->Par_P6 && Neib_ID==TDatabase::ParamDB->Par_P5 )
//             printf("Cross edge Update FALSE rank %d edge %d M %d dofa % dofb %d NeibID %d\n", rank, j , M, N, DOF[ EdgeDof[1] ], Neib_ID);
// 

                  UPDATE = FALSE;
                  break;
                 }

               if(UPDATE)
                {
                 DofRankIndex[N*MaxSubDomainPerDof + M] = Neib_ID;
                 N_DofRankIndex[N]++;
                }
               } // for(l=0; l<N_EdgeDOF; l++)

//                M = N_DofRankIndex[N];
// 
//           if(rank==TDatabase::ParamDB->Par_P5 && Neib_ID==TDatabase::ParamDB->Par_P6 )
//             printf("Cross edge !!!\n");
// 
//                UPDATE = TRUE;
//                for(l=0; l<M; l++)
//                 if(Neib_ID== DofRankIndex[N*MaxSubDomainPerDof + l ])
//                  {
// 
//          if(rank==TDatabase::ParamDB->Par_P6 && Neib_ID==TDatabase::ParamDB->Par_P5 )
//             printf("Cross edge Update FALSE rank %d edge %d M %d dofa % dofb %d NeibID %d\n", rank, j , M, N, DOF[ EdgeDof[1] ], Neib_ID);
// 
// 
//                   UPDATE = FALSE;
//                   break;
//                  }
// 
// 
// 
//                if(UPDATE)
//                 {
//                  for(l=0; l<N_EdgeDOF; l++) 
//                   {
//                    P = DOF[ EdgeDof[l] ];
//                    DofRankIndex[P*MaxSubDomainPerDof + N_DofRankIndex[P] ] = Neib_ID;
//                    N_DofRankIndex[P]++;
// 
//           if(rank==TDatabase::ParamDB->Par_P5 && Neib_ID==TDatabase::ParamDB->Par_P6 )
//              printf("DofRankIndex Rank %d cell %d joint %d Edgeindex  %d Dof %d M %d !!!\n",
//                       rank, i, j,  EdgeDof[l], P, N_DofRankIndex[P]);
//                   }//for(l=0; l<N_Ed
//                  } // if(UPDATE)
              } // for(k=0; k<N_CrossEd
           } // if( edge->IsSubDomainCrossEdge() &&  (edge
         } // for(j=0;j<N_Edges
       } // if(cell->IsCrossEdgeCe
   } //  for(i=0; i<N_Cells

//   M=0;
//   if(rank==TDatabase::ParamDB->Par_P5)
//     for(i=0; i<N_U; i++)
//      for(j=0; j<N_DofRankIndex[i]; j++)
//       if(DofRankIndex[i*MaxSubDomainPerDof + j ]==TDatabase::ParamDB->Par_P6)
//         M++;
// 
//   if(rank==TDatabase::ParamDB->Par_P5)
//    printf("Edge  Test rank %d Neib %d  M %d!!!\n", rank,TDatabase::ParamDB->Par_P6, M);
// 
//   M=0;
//   if(rank==TDatabase::ParamDB->Par_P6)
//     for(i=0; i<N_U; i++)
//      for(j=0; j<N_DofRankIndex[i]; j++)
//       if(DofRankIndex[i*MaxSubDomainPerDof + j ]==TDatabase::ParamDB->Par_P5)
//         M++;
// 
// 
//   if(rank==TDatabase::ParamDB->Par_P6)
//    printf("Edge  Test rank %d Neib %d M %d!!!\n", rank,TDatabase::ParamDB->Par_P5, M);

//  MPI_Finalize();
//  exit(0);


  /** dofs on vertices */
 // put neib cross vertex subdomain ID in all DOF in the vertex 
   for(i=0; i<N_Cells; i++)
    {
     cell = Coll->GetCell(i);
     ID = cell->GetSubDomainNo();

//      if(rank!=ID)
//       continue;

     if(cell->IsCrossVertexCell())
       { 
        DOF = GlobalNumbers + BeginIndex[i];
        N_VertInCell = cell->GetN_Vertices();
        FeId = FESpace->GetFE3D(i, cell);
        FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);

       // Discspace or non-conforming space
       if(FeDesc->GetN_VertDOF()<=0)
        break;

        for(j=0;j<N_VertInCell;j++)
         {
          Vert=cell->GetVertex(j);

          if( Vert->IsCrossVert() &&  (Vert->GetClipBoard()==-1) )
           {
            Vert->SetClipBoard(5);
            Vert->GetCrossNeibs(N_VertCrossNeibs, VertCrossNeibs);
            VertDof = FeDesc->GetVertDOF(j);
            N = DOF[VertDof];

             for(k=0; k<N_VertCrossNeibs; k++)
              {
               Neib_ID = VertCrossNeibs[k];
               M = N_DofRankIndex[N];

//           if(rank==TDatabase::ParamDB->Par_P5 && Neib_ID==TDatabase::ParamDB->Par_P6 )
//              printf("DofRankIndex cell %d joint %d VertDof  %d Dof %d M %d !!!\n",
//                       i, j, VertDof, N, M);


               UPDATE = TRUE;
               for(l=0; l<M; l++)
                if(Neib_ID== DofRankIndex[N*MaxSubDomainPerDof + l ])
                 {
                  UPDATE = FALSE;
                  break;
                 }

               if(UPDATE)
                {
                 DofRankIndex[N*MaxSubDomainPerDof + M ] = Neib_ID;
                 N_DofRankIndex[N]++;

                } // if(UPDATE)
              } // for(k=0; k<N_VertCrossNeibs
           } //  if( Vert->IsCrossVert() &&  (V
         } //  for(j=0;j<N_VertIn
       } // if(cell->IsCrossVertexCell())
   } //  for(i=0; i<N_Cells

//   M=0;
//   if(rank==TDatabase::ParamDB->Par_P5)
//     for(i=0; i<N_U; i++)
//      for(j=0; j<N_DofRankIndex[i]; j++)
//       if(DofRankIndex[i*MaxSubDomainPerDof + j ]==TDatabase::ParamDB->Par_P6)
//         M++;
// 
// 
//   if(rank==TDatabase::ParamDB->Par_P5)
//    printf("Rank %d Vert  Test %d M %d!!!\n", rank, test, M);
// 
//   M=0;
//   if(rank==TDatabase::ParamDB->Par_P6)
//     for(i=0; i<N_U; i++)
//      for(j=0; j<N_DofRankIndex[i]; j++)
//       if(DofRankIndex[i*MaxSubDomainPerDof + j ]==TDatabase::ParamDB->Par_P5)
//         M++;
// 
// 
//   if(rank==TDatabase::ParamDB->Par_P6)
//    printf("Rank %d Vert  Test %d M %d!!!\n", rank, test, M);
// 

//  MPI_Finalize();
//  exit(0);


  N_Neibs = 0;
  for(i=0;i<N_U;i++)
    if(N_Neibs <N_DofRankIndex[i])
     {
      N_Neibs=N_DofRankIndex[i]; 
      break;
     }

  N = 0;
  if(N_Neibs)
   {
    NeibsRank =  new int[size];

    for(i=0;i<N_U;i++)
     {
      M = N_DofRankIndex[i];

      for(j=0; j<M; j++)
       {
        ID = DofRankIndex[i*MaxSubDomainPerDof + j];

        UPDATE = TRUE;
        for(k=0; k<N; k++)
         if(ID== NeibsRank[k])  
	  {
           UPDATE = FALSE;
	   break;
	  }

	if(UPDATE)
	 {
	  NeibsRank[N] = ID;
	  N++;
	 } // if(UPDATE)
       } //  for(j=0; j<M; 
     } // or(i=0;i<N_U
   } // if(N_Neibs)

  N_Neibs = N;

  N_DeptDofs = 0;
  for(i=0;i<N_U;i++)
    if(N_DofRankIndex[i]>1)
     N_DeptDofs++;

  MPI_Allreduce(&N_DeptDofs, &MaxN_DeptDofs_All, 1, MPI_INT, MPI_MAX, Comm); 
//   if(rank==0)
//    MaxN_DeptDofs_All = 0;

  N=0;
  if(N_DeptDofs)
   {
    DeptDofs = new int [N_DeptDofs];
    N_DeptDofNeibs = new int [N_DeptDofs];
    DeptDofMaster = new int[N_DeptDofs];
    DeptDofNeibRanks = new int [N_DeptDofs*MaxSubDomainPerDof];
    DeptDofNeibRanksLocalDOF = new int [N_DeptDofs*MaxSubDomainPerDof];
    DeptDofLocalIndexOfNeibRanks = new int [N_DeptDofs*MaxSubDomainPerDof];
    DepDofIndexOfLocDof = new int [N_U];
 
    memset(DepDofIndexOfLocDof, -1, N_U*SizeOfInt);    
//     for(i=0;i<N_DeptDofs;i++)
//      DeptDofMaster[i] = -1;

    for(i=0;i<N_DeptDofs*MaxSubDomainPerDof;i++)
     DeptDofNeibRanksLocalDOF[i] = -1;

    N=0;
    for(i=0;i<N_U;i++)
     {
     if(N_DofRankIndex[i]>1)
      {
       DepDofIndexOfLocDof[i] = N;
       DeptDofs[N++] = i;
      }
     }      
   }//if(N_DeptDofs)

    //sort DeptCellGlobalNo in assending order, useful in finding local cell no
    for(i=0;i<N_DependentCells-1;i++)
     for(j=i+1;j<N_DependentCells;j++)
      if(DeptCellGlobalNo[i] > DeptCellGlobalNo[j])
       {
        M =  DeptCellGlobalNo[i];
        DeptCellGlobalNo[i] = DeptCellGlobalNo[j];
        DeptCellGlobalNo[j] = M;

         M = DependentCellIndex[i];
         DependentCellIndex[i] = DependentCellIndex[j];
         DependentCellIndex[j] = M;
       }

//    N = 0;
//    if(rank==TDatabase::ParamDB->Par_P5)
//     for(k=0; k<N_U; k++)
//      if(N_DofRankIndex[k]>1)
//      printf("rank %d sl %d index %d SubDomainJoint Neib_ID %d \n", rank, N++, k,  DeptDofs[N]);

// if(rank==TDatabase::ParamDB->Par_P5)
//   printf("Rank %d ConstructDofRankIndex   N_DependentCells info %d \n", rank, N_DependentCells);
//   printf("Rank %d ConstructDofRankIndex %d \n", rank, N_Neibs-1);

//   printf("rank %d Final N_Neibs  %d   \n", rank,   N_Neibs);

//    if(rank==TDatabase::ParamDB->Par_P5)
//     for(k=0; k<N_Neibs; k++)
//      printf("Rank %d ConstructDofRankIndex %d \n", rank, NeibsRank[k]);

//   if(rank==TDatabase::ParamDB->Par_P5)
//    for(i=0; i<N_DependentCells; i++)
//      printf("%d ConstructDofRankIndex done %d\n", i, DependentCellIndex[i]);


// //                if(rank==0 )
//                 printf("Final %d  face N_DOFNeibs  %d   \n", rank,   N_DofRankIndex[4038]);
// // 
// //                if(rank==1 )
// //                 printf("Final %d  face N_DOFNeibs  %d    \n", rank,   N_DofRankIndex[3745]);  	
// //      

//   MPI_Finalize();
//   exit(0);
  
  if(TDatabase::ParamDB->SC_VERBOSE>4)
  if(rank==TDatabase::ParamDB->Par_P0)
   printf("ConstructDofRankIndex done !!!\n");

//   MPI_Finalize();
//   exit(0);
}


static int GetLocalIndex(int N, int *array, int val)
{
 int m=0;
 while(array[m] != val)
    m++;
 return m;
}


void TParFECommunicator3D::ConstructDofRankIndex()
{
 
 int rank, size, i,jj,ii, j, k, l, m, m1, P, N_Cells, N_U, N_LocDof, ID;
 int M, N, N_Vert, N_Joints, N_JointDOF, Neib_ID;
 int *DOF, *GlobalNumbers, *BeginIndex, *JointDof, Disp, N_EdgeDOF, *EdgeDof,*LocalIndex;
 int N_CrossEdgeNeibs, *CrossEdgeNeibsRank, N_Edges, N_VertInCell;
 int N_VertCrossNeibs, *VertCrossNeibs, VertDof;
 int test, N_Dof;

 double x,y,z;

 bool UPDATE;

 TCollection *Coll;
 TBaseCell *cell, *Neib_cell;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_U = FESpace->GetN_DegreesOfFreedom();
  
  printf("N_U is %d----Rank %d\n",N_U,rank);

  /** *************************************************************************************/
  /** ARRAY CONTAINING GLOBAL DOF numbers (GLOBAL means LOCALLY GLOBAL for the subdomain) */ 
  /** BeginIndex: START ID IN GlobalNumbers ARRAY Cellwise*********************************/
  BeginIndex = FESpace->GetBeginIndex();		
  GlobalNumbers = FESpace->GetGlobalNumbers();          
   
  /** Array containing number of ranks(process ID) associated (surrounding) with each DOF (including halo DOF) */ 
  N_DofRankIndex = new int[N_U];
  /** Array containing the rank ID surrounding each DOF (includeing HALO) */
  DofRankIndex = new int[MaxSubDomainPerDof*N_U];
  memset(N_DofRankIndex, 0, N_U*SizeOfInt);
  /** Array containing global number of all Local cells */
  LocalIndex = new int[N_Cells]; 

 
  Disp = N_U*MaxSubDomainPerDof;
  
  for(i=0; i<Disp; i++)       /** Initializing DofRankIndex Array */                   
   DofRankIndex[i] = -1;       
  
  N_DependentCells=0;
    
  
  /** Array contains a flag denoting which type of cells (Halo and Other) have been visited to determine its Master */ 
  int **celltypevisited;	
  int **MappingData;
  
  celltypevisited = new int*[2];
  for(i=0;i<2;i++)
  {
    celltypevisited[i] = new int[N_U];
    memset(celltypevisited[i], 0, N_U*SizeOfInt);
  }
  
  MappingData = new int*[2];
  for(i=0;i<2;i++)
  {
    MappingData[i] = new int[N_U];
    memset(MappingData[i], 0, N_U*SizeOfInt);
  }
  
 /** START ---> [ DofRankIndex ------- N_DofRankIndex -------- N_DependentCells ] **/ 
 
  for(i=0; i<N_Cells; i++)
   {
      cell = Coll->GetCell(i);
      ID = cell->GetSubDomainNo();
      LocalIndex[i] = cell->GetGlobalCellNo();
      DOF = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];

      for(j=0; j<N_LocDof; j++)
      {
	N = DOF[j];
	M = N_DofRankIndex[N];
	/** adding own process ID */
	
	UPDATE = TRUE;
        /**======================================================================
	  * FOLLOWING 'for' LOOP CHECKS WHETHER THE CURRENT RANK ID (i.e) 
          * VARIABLE ID HAS ALREADY BEEN INCLUDED IN THE LIST OF THE PRESENT DOF
          * =====================================================================*/
	for(k=0; k<M; k++)
	  if(ID==DofRankIndex[N*MaxSubDomainPerDof + k])
          {
           UPDATE = FALSE;
	   break;
	  } // for(k=0; k<M; k++)

        if(UPDATE)
	{
	  DofRankIndex[N*MaxSubDomainPerDof + M ] = ID;
	  N_DofRankIndex[N]++;
		  
	  /** Keeping a note of which type of Cells contain this dof */
	 if(cell->IsHaloCell())
	 {
	   MappingData[GLOBAL_NO][N]    = cell->GetGlobalCellNo();
	   MappingData[DOF_NO][N]       = j; 
	   celltypevisited[HALOCELL][N] = 1;
	 }
	 else
	   celltypevisited[NONHALO][N]  = 1;
	}
       }  // for(j=0; j<N_DOF; j++)
    } //  for(i=0; i<N_Cells
   
  /** END ---> [ DofRankIndex ------- N_DofRankIndex -------- N_DependentCells ] **/ 
    
  /** START ---> [ TemporaryMaster ------------ Verify] **/ 
  
    int *Verify;
    int temp,temp_globalno,temp_dofno;
   
    
    Master = new int[N_U];
    Verify = new int[N_U];
	    
    memset(Verify, 0, N_U*SizeOfInt);
	    
    for(i=0;i<N_U;i++)
    {
        temp = size ;
	for(j=0;j<N_DofRankIndex[i];j++)
	{
	  if(temp > DofRankIndex[i*MaxSubDomainPerDof + j])
	  temp = DofRankIndex[i*MaxSubDomainPerDof + j];
	}
	Master[i] = temp;
	if( celltypevisited[0][i] == 1 && celltypevisited[1][i] == 0)
	  Verify[i] = 1;			      
   }
    
  
   /** END ---> [ Master ------------ Verify] **/ 
  
    int *N_DOFtobeverified;	  /** Array containing how many DOF's need to be verified from other processors(ranks) */
    int *N_DOFverifiedbythisrank; /** Array containing how many DOF's need to be verified by this rank for other processors(ranks) */
    int VerifiedDOF = 0;	  /** Total no of DOF's that need to be verified from other ranks */
    int Verifiedbythisrank =0;    /** Total no of DOF verified by this rank */
    int **sendbuf,**recvbuf,*verrecvbuf,*versendbuf;
     
    /** PROTECTED VARIABLES  * int *sdispl,*rdispl;*/
    int *temp_arr;
	
    sendbuf 	= new int*[2];
    recvbuf 	= new int*[2];
    sdispl  	= new int[size];
    rdispl  	= new int[size];
    temp_arr	= new int[size];
	
    N_DOFtobeverified = new int[size];
    memset(N_DOFtobeverified,0,size*SizeOfInt);
    N_DOFverifiedbythisrank = new int[size];
    memset(N_DOFtobeverified,0,size*SizeOfInt);
	
    for(i=0;i<N_U;i++)
    {
      if(Verify[i] == 1)
      {
	N_DOFtobeverified[Master[i]]++;
	VerifiedDOF++;
      }
    }
		
    MPI_Alltoall(N_DOFtobeverified,1, MPI_INT,N_DOFverifiedbythisrank ,1, MPI_INT, Comm);
					
    for(i=0;i<size;i++)
      Verifiedbythisrank +=  N_DOFverifiedbythisrank[i];
		
    for(i=0;i<2;i++)
    {
      sendbuf[i] = new int[VerifiedDOF];
      recvbuf[i] = new int[Verifiedbythisrank];
    }
		
    //printf("Rank F %d----------%d\n",rank,Verifiedbythisrank);
    sdispl[0] = 0;
    for(i=1;i<size;i++)
    {
      sdispl[i] = N_DOFtobeverified[i-1] + sdispl[i-1];
      if(rank == 0);
      //printf(" %d \n",sdispl[i]);
    }
					
    rdispl[0] = 0;
    for(i=1;i<size;i++)
      rdispl[i] = N_DOFverifiedbythisrank[i-1] + rdispl[i-1];
					
					
    memcpy (temp_arr,sdispl, size*SizeOfInt );
					
    for(i=0;i<N_U;i++)
     {
      if(Verify[i] == 1)
      {
	//printf("%d ",temp_arr[Master[i]]);
	sendbuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
	sendbuf[DOF_NO][temp_arr[Master[i]]] = MappingData[DOF_NO][i];
	temp_arr[Master[i]]++;
      }
    }

    MPI_Alltoallv(sendbuf[GLOBAL_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[GLOBAL_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
    MPI_Alltoallv(sendbuf[DOF_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[DOF_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
	      
    verrecvbuf  	= new int[Verifiedbythisrank];
    versendbuf  	= new int[VerifiedDOF];

    for(i=0;i<Verifiedbythisrank;i++)
    {
      temp = GetLocalIndex(N_Cells,LocalIndex,recvbuf[GLOBAL_NO][i]);
      verrecvbuf[i] = Master[GlobalNumbers[temp * N_LocDof + recvbuf[DOF_NO][i]]];
    }
	
    MPI_Alltoallv(verrecvbuf, N_DOFverifiedbythisrank, rdispl, MPI_INT,versendbuf,N_DOFtobeverified, sdispl, MPI_INT, Comm);
					
    for(i=0;i<VerifiedDOF;i++)
    {
      temp_globalno = sendbuf[GLOBAL_NO][i];
      temp_dofno    = sendbuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
      if(Verify[temp] != 1)
	printf("Error : This degree of Freedom (%d) didn't require verification\n",temp);
      Master[temp] = versendbuf[i];
    }
					
/**===================================== END : MASTER FOR ALL DEGREES OF FREEDOM HAS BEEN DECIDED ========================*/

    delete [] verrecvbuf;
    delete [] versendbuf;
    delete [] N_DOFtobeverified;
    delete [] N_DOFverifiedbythisrank;
    delete [] temp_arr;

		
    for(i=0;i<2;i++)
    {
      delete [] sendbuf[i];
      delete [] recvbuf[i];
    }

    delete [] sendbuf;
    delete [] recvbuf;
		
					
/**==================================START : NOW FIRST GATHERING INFORMATION REGARDING WHICH 
 *                                    DOF's have to be received from OTHER PROCESSORS==================*/
					
    int **SlaveBuf,**MasterBuf,SizeDofSend = 0,*Test_Recv,*Test_Send;
    /** PROTECTED VARIABLES * int *N_DofSend,*N_DofRecv,*DofSend,*DofRecv;*/
    N_Slave = 0;
    N_OwnDof = 0;
					
    N_DofSend = new int[size];
    N_DofRecv = new int[size];
    temp_arr  = new int[size];
    SlaveBuf  = new int*[2];
    MasterBuf = new int*[2];
	
		
    memset(N_DofRecv,0,size*SizeOfInt);
			
    for(i=0;i<N_U;i++)
    {
     if(Master[i] != rank)
     {  
      N_Slave++;
      N_DofRecv[Master[i]]++;
     }
     else
       N_OwnDof++;
    }
    printf("Rank %d      OwnDofs is %d---------\n",rank,N_OwnDof);
    OwnDofs = new int[N_OwnDof];
    rdispl[0] = 0; 

    MPI_Alltoall(N_DofRecv,1, MPI_INT,N_DofSend ,1, MPI_INT, Comm);

    for(i=1;i<size;i++)
      rdispl[i] = rdispl[i-1] + N_DofRecv[i-1];
	
    sdispl[0] = 0;
    for(i=0;i<size;i++)
    {
      SizeDofSend += N_DofSend[i];
      sdispl[i] = sdispl[i-1] + N_DofSend[i-1];
    }


    DofSend  = new int[SizeDofSend];
    Test_Send = new int[SizeDofSend];
    DofRecv  = new int[N_Slave];
    Test_Recv = new int[N_Slave];
    
    N_SendDof = SizeDofSend;
	
    for(i=0;i<2;i++)
    {
      SlaveBuf[i]  = new int[N_Slave];
      MasterBuf[i] = new int[SizeDofSend];
    }
	
    memcpy (temp_arr,rdispl, size*SizeOfInt );
	
    m=0;
    for(i=0;i<N_U;i++)
    {
      if(Master[i] != rank)
      {
	SlaveBuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
	SlaveBuf[DOF_NO]   [temp_arr[Master[i]]] = MappingData[DOF_NO]   [i];
	DofRecv[temp_arr[Master[i]]] = i;
	temp_arr[Master[i]]++;
      }
      else
      {
	OwnDofs[m]=i;
	m++;
      }
    }
			
    MPI_Alltoallv(SlaveBuf[GLOBAL_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[GLOBAL_NO], N_DofSend, sdispl, MPI_INT, Comm);
    MPI_Alltoallv(SlaveBuf[DOF_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[DOF_NO], N_DofSend, sdispl, MPI_INT, Comm);
		
    for(i=0;i<SizeDofSend;i++)
    {
      temp_globalno = MasterBuf[GLOBAL_NO][i];
      temp_dofno    = MasterBuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
      DofSend[i] = temp;  
    }
    if(TDatabase::ParamDB->SC_VERBOSE>1)   
	printf(" Rank %d ------ NUMBER OF DOF's to be sent = %d -------- NUMBER OF DOF's to be recv = %d\n",rank,N_SendDof,N_Slave); 

    if(TDatabase::ParamDB->SC_VERBOSE>2)
      if(rank==TDatabase::ParamDB->Par_P0)
	printf("ConstructDofRankIndex done !!!\n");
}

void TParFECommunicator3D::ConstructDofMapRe()
{
   
 int rank, size, i,jj,ii, j, k, l, m, m1, P, N_Cells, N_U, N_Active, N_LocDof, ID;
 int M, N, N_Vert, N_Joints, N_JointDOF, Neib_ID;
 int *DOF, *GlobalNumbers, *BeginIndex, *JointDof, Disp, N_EdgeDOF, *EdgeDof,*LocalIndex,*Verify;
 int N_CrossEdgeNeibs, *CrossEdgeNeibsRank, N_Edges, N_VertInCell,N_OwnCells;
 int N_VertCrossNeibs, *VertCrossNeibs, VertDof;
 int test, N_Dof;

 double x,y,z;

 bool UPDATE;

 TCollection *Coll;
 TBaseCell *cell, *Neib_cell;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_Active = FESpace->GetN_ActiveDegrees();
  N_U = FESpace->GetN_DegreesOfFreedom();
  N_OwnCells = Coll->GetN_OwnCells();
  
 // printf("N_U is %d----Rank %d\n",N_U,rank);
  
  for(i=0;i<N_OwnCells;i++)
  {
    cell = Coll->GetCell(i);
    if(cell->IsHaloCell()) printf("This shudnt happen--------\n");
  }
  for(i=N_OwnCells;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(!cell->IsHaloCell()) printf("This shudnt happen--------\n");
  }

  /** *************************************************************************************/
  /** ARRAY CONTAINING GLOBAL DOF numbers (GLOBAL means LOCALLY GLOBAL for the subdomain) */ 
  /** BeginIndex: START ID IN GlobalNumbers ARRAY Cellwise*********************************/
  BeginIndex = FESpace->GetBeginIndex();		
  GlobalNumbers = FESpace->GetGlobalNumbers();          
   
  /** Array containing number of ranks(process ID) associated (surrounding) with each DOF (including halo DOF) */ 
  N_DofRankIndex = new int[N_U];
 
  memset(N_DofRankIndex, 1, N_U*SizeOfInt);
  /** Array containing global number of all Local cells */
  LocalIndex = new int[N_Cells]; 
  
  int **MappingData;
  MappingData = new int*[2];
  for(i=0;i<2;i++)
  {
    MappingData[i] = new int[N_U];
    memset(MappingData[i], 0, N_U*SizeOfInt);
  }
  
  Master = new int[N_U];
  for(i=0;i<N_U;i++) Master[i]=rank;
    
  Verify = new int[N_U];
  memset(Verify, 0, N_U*SizeOfInt);
  
 /** START ---> [ DofRankIndex ------- N_DofRankIndex -------- N_DependentCells ] **/ 
 for(i=0; i<N_Cells; i++)
  {
     cell = Coll->GetCell(i);
     LocalIndex[i] = cell->GetGlobalCellNo();
     
     if(cell->IsDependentCell() && !cell->IsHaloCell())
     {
      DOF = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];
      for(j=0; j<N_LocDof; j++)
      {
	N = DOF[j];
	Verify[N] = -1;
      }  // for(j=0; j<N_DOF; j++)
      
     }
  } //  for(i=0; i<N_Cells
    
 for(i=N_OwnCells;i<N_Cells;i++)
 {
      cell = Coll->GetCell(i);
      ID = cell->GetSubDomainNo();
      DOF = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];

      for(j=0; j<N_LocDof; j++)
      {
	N = DOF[j];
	//if(N>=N_Active) continue;
	if(Verify[N] == 0) 
	{
	  Verify[N] = 1;
	  Master[N] = ID;
	}
	if(ID < Master[N]) Master[N] = ID;
	MappingData[GLOBAL_NO][N]    = cell->GetGlobalCellNo();
	MappingData[DOF_NO][N]       = j; 
      }
      
  }
 
   
    int temp,temp_globalno,temp_dofno;
  
  
   /** END ---> [ Master ------------ Verify] **/ 
  
    int *N_DOFtobeverified;	  /** Array containing how many DOF's need to be verified from other processors(ranks) */
    int *N_DOFverifiedbythisrank; /** Array containing how many DOF's need to be verified by this rank for other processors(ranks) */
    int VerifiedDOF = 0;	  /** Total no of DOF's that need to be verified from other ranks */
    int Verifiedbythisrank =0;    /** Total no of DOF verified by this rank */
    int **sendbuf,**recvbuf,*verrecvbuf,*versendbuf;
     
    /** PROTECTED VARIABLES  * int *sdispl,*rdispl;*/
    int *temp_arr;
	
    sendbuf 	= new int*[2];
    recvbuf 	= new int*[2];
    sdispl  	= new int[size];
    rdispl  	= new int[size];
    temp_arr	= new int[size];
	
    N_DOFtobeverified = new int[size];
    memset(N_DOFtobeverified,0,size*SizeOfInt);
    N_DOFverifiedbythisrank = new int[size];
    memset(N_DOFtobeverified,0,size*SizeOfInt);
	
    for(i=0;i<N_U;i++)
    {
      if(Verify[i] == 1)
      {
	N_DOFtobeverified[Master[i]]++;
	VerifiedDOF++;
      }
    }
		
    MPI_Alltoall(N_DOFtobeverified,1, MPI_INT,N_DOFverifiedbythisrank ,1, MPI_INT, Comm);
					
    for(i=0;i<size;i++)
      Verifiedbythisrank +=  N_DOFverifiedbythisrank[i];
		
    for(i=0;i<2;i++)
    {
      sendbuf[i] = new int[VerifiedDOF];
      recvbuf[i] = new int[Verifiedbythisrank];
    }
		
    //printf("Rank F %d----------%d\n",rank,Verifiedbythisrank);
    sdispl[0] = 0;
    for(i=1;i<size;i++)
    {
      sdispl[i] = N_DOFtobeverified[i-1] + sdispl[i-1];
      if(rank == 0);
      //printf(" %d \n",sdispl[i]);
    }
					
    rdispl[0] = 0;
    for(i=1;i<size;i++)
      rdispl[i] = N_DOFverifiedbythisrank[i-1] + rdispl[i-1];
					
					
    memcpy (temp_arr,sdispl, size*SizeOfInt );
					
    for(i=0;i<N_U;i++)
     {
      if(Verify[i] == 1)
      {
	//printf("%d ",temp_arr[Master[i]]);
	sendbuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
	sendbuf[DOF_NO][temp_arr[Master[i]]] = MappingData[DOF_NO][i];
	temp_arr[Master[i]]++;
      }
    }

    MPI_Alltoallv(sendbuf[GLOBAL_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[GLOBAL_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
    MPI_Alltoallv(sendbuf[DOF_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[DOF_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
	      
    verrecvbuf  	= new int[Verifiedbythisrank];
    versendbuf  	= new int[VerifiedDOF];

    for(i=0;i<Verifiedbythisrank;i++)
    {
      temp = GetLocalIndex(N_Cells,LocalIndex,recvbuf[GLOBAL_NO][i]);
      verrecvbuf[i] = Master[GlobalNumbers[temp * N_LocDof + recvbuf[DOF_NO][i]]];
    }
	
    MPI_Alltoallv(verrecvbuf, N_DOFverifiedbythisrank, rdispl, MPI_INT,versendbuf,N_DOFtobeverified, sdispl, MPI_INT, Comm);
					
    for(i=0;i<VerifiedDOF;i++)
    {
      temp_globalno = sendbuf[GLOBAL_NO][i];
      temp_dofno    = sendbuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
      if(Verify[temp] != 1)
	printf("Error : This degree of Freedom (%d) didn't require verification\n",temp);
      Master[temp] = versendbuf[i];
    }
					
/**===================================== END : MASTER FOR ALL DEGREES OF FREEDOM HAS BEEN DECIDED ========================*/

    delete [] verrecvbuf;
    delete [] versendbuf;
    delete [] N_DOFtobeverified;
    delete [] N_DOFverifiedbythisrank;
    delete [] temp_arr;

		
    for(i=0;i<2;i++)
    {
      delete [] sendbuf[i];
      delete [] recvbuf[i];
    }

    delete [] sendbuf;
    delete [] recvbuf;
		
		
/**==================================START : NOW FIRST GATHERING INFORMATION REGARDING WHICH 
 *                                    DOF's have to be received from OTHER PROCESSORS==================*/
  int *Verify2;
  bool dept=FALSE;
   Verify2 = new int[N_U];
  memset(Verify2, 0, N_U*SizeOfInt);
  for(i=N_OwnCells;i<N_Cells;i++)
  {
   cell = Coll->GetCell(i);
   DOF = GlobalNumbers + BeginIndex[i];
   N_LocDof = BeginIndex[i+1] - BeginIndex[i];

    for(j=0; j<N_LocDof; j++)
     {
      N = DOF[j];
      Verify2[N]=-1;
     }  // for(j=0; j<N_DOF; j++)
  }
  
  for(i=0; i<N_Cells; i++)
  {
     cell = Coll->GetCell(i);

     if(cell->IsDependentCell() && !cell->IsHaloCell())
     {
      DOF = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];

      for(j=0; j<N_LocDof; j++)
      {
	N = DOF[j];
	if(Master[N] !=rank)
	    dept = TRUE;
      }  // for(j=0; j<N_DOF; j++) 
      
      if(dept)
      {
       for(j=0; j<N_LocDof; j++)
       {
	 N = DOF[j];
	 if(Verify2[N]!=-1 && Master[N]==rank)
	    Verify[N]=-2;
       }  // for(j=0; j<N_DOF; j++) 
      }
	
     }
  } //  for(i=0; i<N_Cells
   delete [] Verify2;

    int **SlaveBuf,**MasterBuf,**SlaveBufH,**MasterBufH,SizeDofSend = 0;
    int *temp_arr2,*temp_arr3,*sdisp, *rdisp ,*N_DofSendTemp, *N_DofRecvTemp;
    // PROTECTED VARIABLES * int *N_DofSend,*N_DofRecv,*DofSend,*DofRecv;
   
    N_OwnDof = 0;
    N_Slave = 0;  N_Halo = 0; N_Master = 0; N_Int=0; N_Dept = 0;
    N_SendDof = 0; N_SendDofH = 0; N_ActiveSlave=0; N_ActiveHalo=0;
    
    N_DofSend = new int[size]; N_DofSendH = new int[size];
    N_DofRecv = new int[size]; N_DofRecvH = new int[size];
    temp_arr  = new int[size];
    SlaveBuf  = new int*[2];
    MasterBuf = new int*[2];

    memset(N_DofRecv,0,size*SizeOfInt);
    memset(N_DofRecvH,0,size*SizeOfInt);
   
    for(i=0;i<N_U;i++)
    {
     if(Master[i] != rank)
     {
       if(Verify[i] == -1)
       {
	 N_Slave++; 
	 N_DofRecv[Master[i]]++;
	 if(i<N_Active) N_ActiveSlave++;
       }
       else
       {
        N_Halo++;
        N_DofRecvH[Master[i]]++;
	if(i<N_Active) N_ActiveHalo++;
       }
     }
     else
     {
       N_OwnDof++;
       if(i<N_Active)
       {
        if(Verify[i] == -1)
         N_Master++; 
        else if(Verify[i] == -2)
	 N_Dept++;
        else 
	 N_Int++;
       }
     }
    }

  /* int N_SlFE= FESpace->GetN_SlaveDegrees(), N_InnerBound=FESpace->GetInnerBound();
  printf("Rank %d N_U %d OwnDofs is %d N_ActiveSlave %d N_ActiveHalo %d N_InnerBound %d---------\n",rank,N_U,N_OwnDof,N_ActiveSlave,N_ActiveHalo,N_InnerBound);
  printf("Rank %d N_Master %d N_Dept %d N_Int %d N_Active %d N_SlFE %d---------\n",rank,N_Master,N_Dept,N_Int,N_Active,N_SlFE);
  MPI_Finalize();exit(0);*/
     //  printf("Rank %d      OwnDofs is %d---------\n",rank,N_OwnDof);
    sdisp = new int[size];
    rdisp = new int[size];
    sdisplH = new int[size];
    rdisplH  = new int[size];
    OwnDofs = new int[N_OwnDof];
    rdispl[0] = 0; rdisplH[0]=0; rdisp[0] = 0;

    MPI_Alltoall(N_DofRecv,1, MPI_INT,N_DofSend ,1, MPI_INT, Comm);
    MPI_Alltoall(N_DofRecvH,1, MPI_INT,N_DofSendH ,1, MPI_INT, Comm);

    for(i=1;i<size;i++)
    {
      rdispl[i] = rdispl[i-1] + N_DofRecv[i-1];
      rdisplH[i] = rdisplH[i-1] + N_DofRecvH[i-1];
    }
   
    sdispl[0] = 0; 
    for(i=0;i<size;i++)
    {
      SizeDofSend += N_DofSend[i] + N_DofSendH[i];
      N_SendDof += N_DofSend[i];
      N_SendDofH += N_DofSendH[i];
      sdispl[i] = sdispl[i-1] + N_DofSend[i-1];
      sdisplH[i] = sdisplH[i-1] + N_DofSendH[i-1];
    }
    
    DofSend  = new int[N_SendDof];
    DofSendH  = new int[N_SendDofH]; 
    DofRecv  = new int[N_Slave];
    DofRecvH  = new int[N_Halo];
   
    for(i=0;i<2;i++)
    {
      SlaveBuf[i]  = new int[N_Slave];
      MasterBuf[i] = new int[N_SendDof];
    }
     
    Reorder = new int[N_U];
    NewGN = new int[N_U];
    int xx=0,xy=N_Master,xz=N_Master+N_Int,yy=N_Master+N_Int+N_Dept,zz=N_Master+N_Int+N_Dept+N_ActiveSlave;
    memcpy (temp_arr,rdispl, size*SizeOfInt );
    m=0;
    for(i=0;i<N_U;i++)
    {
      if(Master[i] != rank)  
      {
	if(Verify[i]==-1)
	{
	  SlaveBuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
	  SlaveBuf[DOF_NO]   [temp_arr[Master[i]]] = MappingData[DOF_NO]   [i];
	  DofRecv[temp_arr[Master[i]]] = i;
          temp_arr[Master[i]]++;
	  if(i<N_Active) { Reorder[yy] = i ; NewGN[i]=yy; yy++;}
	}
	else
	{
	  if(i<N_Active) {Reorder[zz] = i ; NewGN[i]=zz; zz++;}
	}
      }
      else 
      {
	OwnDofs[m]=i;
	m++;
	if(i<N_Active) 
	{
	 if(Verify[i] == -1)
	 {
	  Reorder[xx] = i ; NewGN[i]=xx; xx++;
	 }
	 else if(Verify[i] == -2)
	 {
	  Reorder[xz] = i ; NewGN[i]=xz; xz++;
	 }
         else
         {
	  Reorder[xy] = i ; NewGN[i]=xy; xy++;
         }
	}
      }
    }

 MPI_Alltoallv(SlaveBuf[GLOBAL_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[GLOBAL_NO], N_DofSend, sdispl, MPI_INT, Comm);	
 MPI_Alltoallv(SlaveBuf[DOF_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[DOF_NO], N_DofSend, sdispl, MPI_INT, Comm);
 
    for(i=0;i<N_SendDof;i++)
    {
      temp_globalno = MasterBuf[GLOBAL_NO][i];
      temp_dofno    = MasterBuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
      DofSend[i] = temp; 
    }
 
    for(i=0;i<2;i++)
    {
      delete [] SlaveBuf[i];
      delete [] MasterBuf[i];
    }
    
    for(i=0;i<2;i++)
    {
      SlaveBuf[i]  = new int[N_Halo];
      MasterBuf[i] = new int[N_SendDofH];
    }
 memcpy (temp_arr,rdisplH, size*SizeOfInt );
    for(i=0;i<N_U;i++)
    {
      if(Master[i] != rank && Verify[i]!=-1)
      {
	SlaveBuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
	SlaveBuf[DOF_NO]   [temp_arr[Master[i]]] = MappingData[DOF_NO]   [i];
	  DofRecvH[temp_arr[Master[i]]] = i;
          temp_arr[Master[i]]++;
      }
    }
 
 MPI_Alltoallv(SlaveBuf[GLOBAL_NO],N_DofRecvH, rdisplH, MPI_INT,MasterBuf[GLOBAL_NO], N_DofSendH, sdisplH, MPI_INT, Comm);	
 MPI_Alltoallv(SlaveBuf[DOF_NO],N_DofRecvH, rdisplH, MPI_INT,MasterBuf[DOF_NO], N_DofSendH, sdisplH, MPI_INT, Comm);
 
    for(i=0;i<N_SendDofH;i++)
    {
      temp_globalno = MasterBuf[GLOBAL_NO][i];
      temp_dofno    = MasterBuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
      DofSendH[i] = temp; 
    }
   
    for(i=0;i<2;i++)
    {
      delete [] SlaveBuf[i];
      delete [] MasterBuf[i];
    }
    delete [] SlaveBuf;
    delete [] MasterBuf;

    //DofSendPos = new int[SizeDofSend];
    //DofRecvPos = new int[N_Slave];
   // for(i=0;i<SizeDofSend;i++) DofSendPos[i]=i;
   // for(i=0;i<N_Slave;i++) DofRecvPos[i]=i;
     //SortPos(DofSend,DofSendPos,SizeDofSend);
    //quickSort(DofSend,DofSendPos,0,SizeDofSend);
  if(N_SendDof>0)
   Send_Info = new double[N_SendDof];
  if(N_Slave>0)
   Recv_Info = new double[N_Slave];
    
    
 if(N_SendDofH>0)
   Send_InfoH = new double[N_SendDofH];
 if(N_Halo>0)
   Recv_InfoH = new double[N_Halo];
    
    if(TDatabase::ParamDB->SC_VERBOSE>2)
      if(rank==TDatabase::ParamDB->Par_P0)
	printf("ConstructDofMap done !!!\n");
       
}

void TParFECommunicator3D::SetNewGN()
{
 int rank, size, i,jj,ii, j, k, l, m;
 int *DOF, *GlobalNumbers, *BeginIndex,N_LocDof, *NewGlobalNumbers, *temp,N_Active;
 int N_U, N_OwnCells, N, N_Cells;
  TCollection *Coll;
  TBaseCell *cell, *Neib_cell;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
  
  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_U = FESpace->GetN_DegreesOfFreedom();
  N_Active = FESpace->GetN_ActiveDegrees();
  N_OwnCells = Coll->GetN_OwnCells();
  BeginIndex = FESpace->GetBeginIndex();		
  GlobalNumbers = FESpace->GetGlobalNumbers(); 
  N_LocDof = BeginIndex[1] - BeginIndex[0];

  temp = new int[N_U];
  memcpy (temp,Master, N_U*SizeOfInt );
 
  for(i=0; i<N_Cells; i++)
  {
     cell = Coll->GetCell(i);
     DOF = GlobalNumbers + BeginIndex[i];
      for(j=0; j<N_LocDof; j++)
      {
	N = DOF[j];
	if(N<N_Active) GlobalNumbers[BeginIndex[i]+j] = NewGN[N];
      }  // for(j=0; j<N_DOF; j++)
  } //  for(i=0; i<N_Cells

   // FESpace->SetGlobalNumbers(NewGlobalNumbers);
    
    for(i=0;i<N_SendDof;i++)  if(DofSend[i]<N_Active) DofSend[i] = NewGN[DofSend[i]];
    for(i=0;i<N_SendDofH;i++) if(DofSendH[i]<N_Active) DofSendH[i]= NewGN[DofSendH[i]];
    for(i=0;i<N_Slave;i++)    if(DofRecv[i]<N_Active) DofRecv[i] = NewGN[DofRecv[i]];
    for(i=0;i<N_Halo;i++)     if(DofRecvH[i]<N_Active) DofRecvH[i]= NewGN[DofRecvH[i]];
    
    for(i=0;i<N_Active;i++) Master[NewGN[i]]=temp[i];
    
    delete [] temp;

   temp = new int[N_OwnDof];
   memcpy (temp, OwnDofs, N_OwnDof*SizeOfInt );
       
   for(i=0;i<N_OwnDof;i++) if(OwnDofs[i]<N_Active) OwnDofs[i]=temp[NewGN[i]];
   delete [] temp;
   if(rank==TDatabase::ParamDB->Par_P0)  
       cout<<"New global numbers set"<<"\n";
}

void TParFECommunicator3D::ColorAndReorder(int start, int end, int &numColors, int *&ptrColors)
{
  int i,j,u,V,max,rank,nDof = (end-start);
  int *RowPtr,*KCol,*temp;
  
  RowPtr = sqstruct->GetRowPtr();
  KCol = sqstruct->GetKCol();
  V = nDof;
  MPI_Comm_rank(Comm, &rank);
  
 if(nDof)
 {
  int result[V];
  max = RowPtr[Reorder[start]+1]-RowPtr[Reorder[start]];
 
    // Assign the first color to first vertex
    result[0]  = 0;
 
    // Initialize remaining V-1 vertices as unassigned
    for (i = 1; i < V; i++)
    {
        result[i] = -1;  // no color is assigned to u
        u = Reorder[start+i];
        if(max<(RowPtr[u+1]-RowPtr[u])) 
             max=(RowPtr[u+1]-RowPtr[u]);
    }
  
    // A temporary array to store the available colors. True
    // value of available[cr] would mean that the color cr is
    // assigned to one of its adjacent vertices
    bool available[max];
    int ptr[max+1]; ptr[max]=0;
    for (int cr = 0; cr < max; cr++)
    {
        available[cr] = false; ptr[cr]=0;
    }
        ptr[1]=1;
    // Assign colors to remaining V-1 vertices
    for (j = 1; j < V; j++)
    {
        // Process all adjacent vertices and flag their colors as unavailable
        u = Reorder[start+j];
	
        for (i = RowPtr[u]; i < RowPtr[u+1]; ++i)
            if (NewGN[KCol[i]]>=start && NewGN[KCol[i]]<end && result[NewGN[KCol[i]]-start] != -1)
                available[result[NewGN[KCol[i]]-start]] = true;

        // Find the first available color
        int cr;
        for (cr = 0; cr < max; cr++)
            if (available[cr] == false)
                break;
 
        result[j] = cr; // Assign the found color
        ptr[cr+1]++; 
        // Reset the values back to false for the next iteration
        for (i = RowPtr[u]; i < RowPtr[u+1]; ++i)
            if (NewGN[KCol[i]]>=start && NewGN[KCol[i]]<end && result[NewGN[KCol[i]]-start] != -1)
                available[result[NewGN[KCol[i]]-start]] = false;
	  
    }
    
    for(i=1;i<max;i++) if(!ptr[i]) break;
    numColors = (i-1);
    
    ptrColors = new int[numColors+1];
   
      temp = new int[nDof];
      memcpy(temp,&Reorder[start],nDof*SizeOfInt);
     
      for(i=1;i<=max;i++) ptr[i] += ptr[i-1];
      memcpy(ptrColors,ptr,(numColors+1)*SizeOfInt);
      
      for(i=0;i<=numColors;i++) ptrColors[i]+=start;
      
//         if(rank==0) {
//     for(i=0;i<(numColors+1);i++) cout<<ptrColors[i]<<"\t";
//     cout<<"\n"; 
     //  cout<<numColors<<"nDof "<<nDof<<"\n";
// }
      
      for(i=0;i<nDof;i++)
      {
	Reorder[start+ptr[result[i]]] = temp[i];
	NewGN[temp[i]] = start+ptr[result[i]];
	ptr[result[i]]++;
      }
 }
 else
   numColors = 0;
  
 /*
    // print the result
    for (int u = 0; u < V; u++)
        cout << "Vertex " << u << " --->  Color "
             << result[u] << endl;*/
}

void TParFECommunicator3D::ConstructDofMap()
{
   
 int rank, size, i,jj,ii, j, k, l, m, m1, P, N_Cells, N_U,N_Active, N_LocDof, ID;
 int M, N, N_Vert, N_Joints, N_JointDOF, Neib_ID;
 int *DOF, *GlobalNumbers, *BeginIndex, *JointDof, Disp, N_EdgeDOF, *EdgeDof,*LocalIndex,*Verify;
 int N_CrossEdgeNeibs, *CrossEdgeNeibsRank, N_Edges, N_VertInCell,N_OwnCells;
 int N_VertCrossNeibs, *VertCrossNeibs, VertDof;
 int test, N_Dof;

 double x,y,z;

 bool UPDATE;

 TCollection *Coll;
 TBaseCell *cell, *Neib_cell;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  N_U = FESpace->GetN_DegreesOfFreedom();
  N_Active = FESpace->GetN_ActiveDegrees();
  N_OwnCells = Coll->GetN_OwnCells();
  
 // printf("N_U is %d----Rank %d\n",N_U,rank);
  
  
  for(i=0;i<N_OwnCells;i++)
  {
    cell = Coll->GetCell(i);
    if(cell->IsHaloCell()) printf("This shudnt happen---------------------------------------------\n");
  }
  for(i=N_OwnCells;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    if(!cell->IsHaloCell()) printf("This shudnt happen---------------------------------------------\n");
  }

  /** *************************************************************************************/
  /** ARRAY CONTAINING GLOBAL DOF numbers (GLOBAL means LOCALLY GLOBAL for the subdomain) */ 
  /** BeginIndex: START ID IN GlobalNumbers ARRAY Cellwise*********************************/
  BeginIndex = FESpace->GetBeginIndex();		
  GlobalNumbers = FESpace->GetGlobalNumbers();          
   
  /** Array containing number of ranks(process ID) associated (surrounding) with each DOF (including halo DOF) */ 
  N_DofRankIndex = new int[N_U];
 
  memset(N_DofRankIndex, 1, N_U*SizeOfInt);
  /** Array containing global number of all Local cells */
  LocalIndex = new int[N_Cells]; 
  
  int **MappingData;
  MappingData = new int*[2];
  for(i=0;i<2;i++)
  {
    MappingData[i] = new int[N_U];
    memset(MappingData[i], 0, N_U*SizeOfInt);
  }
  
  Master = new int[N_U];
  for(i=0;i<N_U;i++) Master[i]=rank;
    
  Verify = new int[N_U];
  memset(Verify, 0, N_U*SizeOfInt);
  
 /** START ---> [ DofRankIndex ------- N_DofRankIndex -------- N_DependentCells ] **/ 
 for(i=0; i<N_Cells; i++)
  {
     cell = Coll->GetCell(i);
     LocalIndex[i] = cell->GetGlobalCellNo();
     
     if(cell->IsDependentCell() && !cell->IsHaloCell())
     {
      DOF = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];

      for(j=0; j<N_LocDof; j++)
      {
	N = DOF[j];
	Verify[N] = -1;
      }  // for(j=0; j<N_DOF; j++)
      
     }
  } //  for(i=0; i<N_Cells
    

 for(i=N_OwnCells;i<N_Cells;i++)
 {
      cell = Coll->GetCell(i);
      ID = cell->GetSubDomainNo();
      DOF = GlobalNumbers + BeginIndex[i];
      N_LocDof = BeginIndex[i+1] - BeginIndex[i];

      for(j=0; j<N_LocDof; j++)
      {
	N = DOF[j];
	//if(N>=N_Active) continue;
	if(Verify[N] == 0) 
	{
	  Verify[N] = 1;
	  Master[N] = ID;
	}
	if(ID < Master[N]) Master[N] = ID;
	MappingData[GLOBAL_NO][N]    = cell->GetGlobalCellNo();
	MappingData[DOF_NO][N]       = j; 
      }
      
  }
 
   
    int temp,temp_globalno,temp_dofno;
  
  
   /** END ---> [ Master ------------ Verify] **/ 
  
    int *N_DOFtobeverified;	  /** Array containing how many DOF's need to be verified from other processors(ranks) */
    int *N_DOFverifiedbythisrank; /** Array containing how many DOF's need to be verified by this rank for other processors(ranks) */
    int VerifiedDOF = 0;	  /** Total no of DOF's that need to be verified from other ranks */
    int Verifiedbythisrank =0;    /** Total no of DOF verified by this rank */
    int **sendbuf,**recvbuf,*verrecvbuf,*versendbuf;
     
    /** PROTECTED VARIABLES  * int *sdispl,*rdispl;*/
    int *temp_arr;
	
    sendbuf 	= new int*[2];
    recvbuf 	= new int*[2];
    sdispl  	= new int[size];
    rdispl  	= new int[size];
    temp_arr	= new int[size];
	
    N_DOFtobeverified = new int[size];
    memset(N_DOFtobeverified,0,size*SizeOfInt);
    N_DOFverifiedbythisrank = new int[size];
    memset(N_DOFtobeverified,0,size*SizeOfInt);
	
    for(i=0;i<N_U;i++)
    {
      if(Verify[i] == 1)
      {
	N_DOFtobeverified[Master[i]]++;
	VerifiedDOF++;
      }
    }
		
    MPI_Alltoall(N_DOFtobeverified,1, MPI_INT,N_DOFverifiedbythisrank ,1, MPI_INT, Comm);
					
    for(i=0;i<size;i++)
      Verifiedbythisrank +=  N_DOFverifiedbythisrank[i];
		
    for(i=0;i<2;i++)
    {
      sendbuf[i] = new int[VerifiedDOF];
      recvbuf[i] = new int[Verifiedbythisrank];
    }
		
    //printf("Rank F %d----------%d\n",rank,Verifiedbythisrank);
    sdispl[0] = 0;
    for(i=1;i<size;i++)
    {
      sdispl[i] = N_DOFtobeverified[i-1] + sdispl[i-1];
      if(rank == 0);
      //printf(" %d \n",sdispl[i]);
    }
					
    rdispl[0] = 0;
    for(i=1;i<size;i++)
      rdispl[i] = N_DOFverifiedbythisrank[i-1] + rdispl[i-1];
					
					
    memcpy (temp_arr,sdispl, size*SizeOfInt );
					
    for(i=0;i<N_U;i++)
     {
      if(Verify[i] == 1)
      {
	//printf("%d ",temp_arr[Master[i]]);
	sendbuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
	sendbuf[DOF_NO][temp_arr[Master[i]]] = MappingData[DOF_NO][i];
	temp_arr[Master[i]]++;
      }
    }

    MPI_Alltoallv(sendbuf[GLOBAL_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[GLOBAL_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
    MPI_Alltoallv(sendbuf[DOF_NO],N_DOFtobeverified, sdispl, MPI_INT,recvbuf[DOF_NO], N_DOFverifiedbythisrank, rdispl, MPI_INT, Comm);
	      
    verrecvbuf  	= new int[Verifiedbythisrank];
    versendbuf  	= new int[VerifiedDOF];

    for(i=0;i<Verifiedbythisrank;i++)
    {
      temp = GetLocalIndex(N_Cells,LocalIndex,recvbuf[GLOBAL_NO][i]);
      verrecvbuf[i] = Master[GlobalNumbers[temp * N_LocDof + recvbuf[DOF_NO][i]]];
    }
	
    MPI_Alltoallv(verrecvbuf, N_DOFverifiedbythisrank, rdispl, MPI_INT,versendbuf,N_DOFtobeverified, sdispl, MPI_INT, Comm);
					
    for(i=0;i<VerifiedDOF;i++)
    {
      temp_globalno = sendbuf[GLOBAL_NO][i];
      temp_dofno    = sendbuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
      if(Verify[temp] != 1)
	printf("Error : This degree of Freedom (%d) didn't require verification\n",temp);
      Master[temp] = versendbuf[i];
    }
					
/**===================================== END : MASTER FOR ALL DEGREES OF FREEDOM HAS BEEN DECIDED ========================*/

    delete [] verrecvbuf;
    delete [] versendbuf;
    delete [] N_DOFtobeverified;
    delete [] N_DOFverifiedbythisrank;
    delete [] temp_arr;

		
    for(i=0;i<2;i++)
    {
      delete [] sendbuf[i];
      delete [] recvbuf[i];
    }

    delete [] sendbuf;
    delete [] recvbuf;
		
				
/**==================================START : NOW FIRST GATHERING INFORMATION REGARDING WHICH 
 *                                    DOF's have to be received from OTHER PROCESSORS==================*/
					
    int **SlaveBuf,**MasterBuf,SizeDofSend = 0;
    /** PROTECTED VARIABLES * int *N_DofSend,*N_DofRecv,*DofSend,*DofRecv;*/
    N_Slave = 0;
    N_OwnDof = 0;
					
    N_DofSend = new int[size];
    N_DofRecv = new int[size];
    temp_arr  = new int[size];
    SlaveBuf  = new int*[2];
    MasterBuf = new int*[2];
	
	//for(i=N_Active;i<N_U;i++) Master[i]=rank;
	
    memset(N_DofRecv,0,size*SizeOfInt);
			
    for(i=0;i<N_U;i++)
    {
     if(Master[i] != rank)
     {  
      N_Slave++;
      N_DofRecv[Master[i]]++;
     } 
     else 
       N_OwnDof++;
    }
  //  printf("Rank %d      OwnDofs is %d---------\n",rank,N_OwnDof);

    OwnDofs = new int[N_OwnDof];
    rdispl[0] = 0; 

    MPI_Alltoall(N_DofRecv,1, MPI_INT,N_DofSend ,1, MPI_INT, Comm);

    for(i=1;i<size;i++)
      rdispl[i] = rdispl[i-1] + N_DofRecv[i-1];
	
    sdispl[0] = 0;
    for(i=0;i<size;i++)
    {
      SizeDofSend += N_DofSend[i];
      sdispl[i] = sdispl[i-1] + N_DofSend[i-1];
    }


    DofSend  = new int[SizeDofSend];
    DofRecv  = new int[N_Slave];
    N_SendDof = SizeDofSend;
	
    for(i=0;i<2;i++)
    {
      SlaveBuf[i]  = new int[N_Slave];
      MasterBuf[i] = new int[SizeDofSend];
    }
	
    memcpy (temp_arr,rdispl, size*SizeOfInt );
	
    m=0;
    for(i=0;i<N_U;i++)
    {
      if(Master[i] != rank)
      {
	SlaveBuf[GLOBAL_NO][temp_arr[Master[i]]] = MappingData[GLOBAL_NO][i];
	SlaveBuf[DOF_NO]   [temp_arr[Master[i]]] = MappingData[DOF_NO]   [i];
	DofRecv[temp_arr[Master[i]]] = i;
	temp_arr[Master[i]]++;
      }
      else 
      {
	OwnDofs[m]=i;
	m++;
      }
    }
			
    MPI_Alltoallv(SlaveBuf[GLOBAL_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[GLOBAL_NO], N_DofSend, sdispl, MPI_INT, Comm);
    MPI_Alltoallv(SlaveBuf[DOF_NO],N_DofRecv, rdispl, MPI_INT,MasterBuf[DOF_NO], N_DofSend, sdispl, MPI_INT, Comm);
		
    for(i=0;i<SizeDofSend;i++)
    {
      temp_globalno = MasterBuf[GLOBAL_NO][i];
      temp_dofno    = MasterBuf[DOF_NO][i];
      temp = GetLocalIndex(N_Cells,LocalIndex,temp_globalno);
      temp = GlobalNumbers[temp*N_LocDof + temp_dofno];
      DofSend[i] = temp;  
    }
  
  //  if(TDatabase::ParamDB->SC_VERBOSE>1)   
	//printf(" Rank %d ------ NUMBER OF DOF's to be sent = %d -------- NUMBER OF DOF's to be recv = %d\n",rank,N_SendDof,N_Slave); 
    for(i=0;i<2;i++)
    {
      delete [] SlaveBuf[i];
      delete [] MasterBuf[i];
    }

    delete [] SlaveBuf;
    delete [] MasterBuf;
	
    DofSendPos = new int[SizeDofSend];
    DofRecvPos = new int[N_Slave];
    for(i=0;i<SizeDofSend;i++) DofSendPos[i]=i;
    for(i=0;i<N_Slave;i++) DofRecvPos[i]=i;
     //SortPos(DofSend,DofSendPos,SizeDofSend);
    //quickSort(DofSend,DofSendPos,0,SizeDofSend);
  if(N_SendDof>0)
   Send_Info = new double[N_SendDof];
  if(N_Slave>0)
   Recv_Info = new double[N_Slave];
    
    if(TDatabase::ParamDB->SC_VERBOSE>2)
      if(rank==TDatabase::ParamDB->Par_P0)
	printf("ConstructDofMap done !!!\n");
}


void TParFECommunicator3D::SetFENeibCommunicationSteps()
{
  int rank,size,i, j, k, l, m, N, P, Q;
  int *dofNeibIDs,N_DofNeibs_MaxAll;
  int *pos, *N_DofNeibs_Array, *DofNeibIDs_Array, N_Entries;
  bool UPDATE;
  
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);
 
  N_DofNeibs=0;
  dofNeibIDs= new int[size];
  m=0;
  for(i=0;i<size;i++)
  {
    if(N_DofRecv[i]>0)
    {
      N_DofNeibs++;
      dofNeibIDs[m]=i;
      m++;
    }
    else if(N_DofSend[i]>0)
    {
      N_DofNeibs++;
      dofNeibIDs[m]=i;
      m++;
    }
  }
  DofNeibIDs = new int[N_DofNeibs];
  memcpy(DofNeibIDs, dofNeibIDs, N_DofNeibs*SizeOfInt);
  delete [] dofNeibIDs;
  
  IndexOfNeibRank = new int[size];
  for(i=0; i<size; i++)
   IndexOfNeibRank[i] = -1;
  for(i=0; i<N_DofNeibs; i++)
   IndexOfNeibRank[DofNeibIDs[i]] = i;
  
  MPI_Allreduce(&N_DofNeibs, &N_DofNeibs_MaxAll, 1, MPI_INT, MPI_MAX, Comm); 
  if(rank==0)
  {
   pos = new int[size];
   for(i=0;i<size;i++)
     pos[i] = i*N_DofNeibs_MaxAll;
   N_DofNeibs_Array = new int[size];
   DofNeibIDs_Array= new int[size*N_DofNeibs_MaxAll];
  }

  MPI_Gather(&N_DofNeibs, 1, MPI_INT, N_DofNeibs_Array, 1, MPI_INT, 0, Comm);
  MPI_Gatherv(DofNeibIDs, N_DofNeibs, MPI_INT, DofNeibIDs_Array, N_DofNeibs_Array, pos,MPI_INT, 0, Comm);
 
  if(rank==0)
  {
   int step, *Entries, *ColInd, *RowPtr, *Collection, MaxPossibleSteps;
   N_Entries = 0;

   for(i=0;i<size;i++)
      N_Entries +=N_DofNeibs_Array[i];
   printf("No.of entries are %d\n",N_Entries);
   Entries = new int[N_Entries];         // adjacency matrix memory allocation
   ColInd     = new int[N_Entries];
   RowPtr = new int[size+1];

   RowPtr[0] = 0; m = 0;
   for(i=0;i<size;i++)
   {
     N =  N_DofNeibs_Array[i];
     RowPtr[i+1] = RowPtr[i] + N;
     for(j=0;j<N;j++)
     {
       ColInd[m] = DofNeibIDs_Array[i*N_DofNeibs_MaxAll + j];  // assmeble the matrix
       Entries[m] =-1;
       m++;
     }
   }
   
    // print adjacency matrix
   /*  for(i=0;i<size;i++)
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
        printf("a( %d , %d )=  %d \n",i, ColInd[j], Entries[j]);
    MPI_Finalize();exit(0);*/
    
   // find number of steps to communicate between neibs (edge coloring alogrithm)
   Collection = new int[size];
   Max_CommunicationSteps = -1;
   MaxPossibleSteps = (N_DofNeibs_MaxAll+1)*10;
   for(step=0;step<MaxPossibleSteps;step++)
   {
    m=0;
    for(k=0;k<size;k++)
    {
     UPDATE = TRUE;
     for(l=0;l<m;l++)
     if(Collection[l] == k)
     {
      UPDATE = FALSE;
       break;
     }
     if(UPDATE)
     {
       Collection[m] = k;
        m++;
     }
     else
     {
      continue;
     }
 
     for(l=RowPtr[k];l<RowPtr[k+1];l++)
     { 
        // only upper tirangular martix is enough (due to symmetry)
         if( ColInd[l]<k) continue;

         UPDATE = TRUE;
         for(P=0;P<m;P++)
         if(Collection[P] == ColInd[l])
         {
           // ColInd[l] is already in the collection
          UPDATE = FALSE;
          break;
         } 
       
         if(UPDATE)
         {
          if(Entries[l]<0)
          {
            Entries[l] = step;
            if(Max_CommunicationSteps<step+1) Max_CommunicationSteps=step+1;
            Collection[m] = ColInd[l];
            m++;
            break;
           }
         }
        } //  for(l=RowPtr[k];
      } // for(k=0;
    } //  for(step=0;step<

    printf("Total number of steps needed to communicate is %d\n",Max_CommunicationSteps);
   //  print adjacency matrix
//     for(i=0;i<size;i++)
//      for(j=RowPtr[i];j<RowPtr[i+1];j++)
//        printf("a( %d , %d )=  %d \n",i, ColInd[j], Entries[j]);
    

    N_CommunicationProcesses = new int[Max_CommunicationSteps];
    SendID_Index= new int[Max_CommunicationSteps*size];
    ReceiveID_Index = new int[Max_CommunicationSteps*size];

    for(j=0;j<Max_CommunicationSteps;j++)
     N_CommunicationProcesses[j] =0;

    for(i=0;i<Max_CommunicationSteps*size;i++)
     {
      SendID_Index[i] = -1;
      ReceiveID_Index[i] = -1;
     }

    for(i=0;i<size;i++)
     for(j=RowPtr[i];j<RowPtr[i+1];j++)
       if(i<ColInd[j])
        {
         if( Entries[j]<0)
          { 
           printf("Error in finding DOF communication steps %d\n",Entries[j]);
           MPI_Abort(MPI_COMM_WORLD,  0);
          }
         k =  Entries[j];
	// if(k==1) printf("i is %d j is %d\n",i,ColInd[j]);
         SendID_Index[k*size + N_CommunicationProcesses[k]] = i;
         ReceiveID_Index[k*size + N_CommunicationProcesses[k]] =  ColInd[j];
         N_CommunicationProcesses[k]++;
        }
   /*
     for(i=0;i<Max_CommunicationSteps;i++)
     {
       N= N_CommunicationProcesses[i];
       printf("Step %d: No.of processes is %d\n",i,N);
       for(j=0;j<N;j++)
        printf("Send from process %d to %d -------- and onto\n",
                   SendID_Index[i*size +j], ReceiveID_Index[i*size +j]);
       printf("\n");
     }   */

    delete [] Collection;
  }
   
   MPI_Bcast(&Max_CommunicationSteps, 1, MPI_INT, 0, Comm);
   MPI_Barrier(Comm);

   if(rank!=0)
    {
     N_CommunicationProcesses = new int[Max_CommunicationSteps];
     SendID_Index= new int[Max_CommunicationSteps*size];
     ReceiveID_Index = new int[Max_CommunicationSteps*size];
    }

   MPI_Bcast(N_CommunicationProcesses, Max_CommunicationSteps, MPI_INT, 0, Comm);
   MPI_Bcast(SendID_Index, Max_CommunicationSteps*size, MPI_INT, 0, Comm);
   MPI_Bcast(ReceiveID_Index, Max_CommunicationSteps*size, MPI_INT,  0, Comm);

 //  printf("Parallel SetIntraCommunicator ranks  %d, N_DofNeibs %d\n", rank, N_DofNeibs);
       //     MPI_Finalize();
       //     exit(0);
}


void TParFECommunicator3D::CommUpdate(double *sol, double *rhs)
{
 double timeC;
  int i, j, k, l, N,iter, rank, size;
 int sendID, recvID, *val_send, *val_recv;
 double t1,t2=0.0,t3,tSum=0.;
 
 MPI_Status status;
 MPI_Comm_rank(Comm,&rank);	
 MPI_Comm_size(Comm, &size);
  
timeC = 0;

//  if(N_SendDof>0)
//    Send_Info = new double[N_SendDof];
//  if(N_Slave>0)
//    Recv_Info = new double[N_Slave];
 
  t1=MPI_Wtime();
  for(i=0;i<N_SendDof;i++)
  {
   Send_Info[i]=sol[DofSend[i]];
   //Send_Info[DofSendPos[i]]=sol[DofSend[i]];
  }
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
 
  for(i=0; i<Max_CommunicationSteps; i++)
  {
    N= N_CommunicationProcesses[i];
    for(j=0;j<N;j++)
    {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if((rank!=sendID) && (rank!=recvID))  continue;
      
      if(rank==sendID && N_DofSend[recvID]>0)
        MPI_Send(&Send_Info[sdispl[recvID]], N_DofSend[recvID], MPI_DOUBLE,  recvID, 100, Comm);
      
      else if(rank==recvID && N_DofRecv[sendID]>0)
        MPI_Recv(&Recv_Info[rdispl[sendID]], N_DofRecv[sendID], MPI_DOUBLE,  sendID, 100, Comm, &status);


      if(rank==recvID && N_DofSend[sendID])
         MPI_Send(&Send_Info[sdispl[sendID]], N_DofSend[sendID], MPI_DOUBLE,  sendID, 200, Comm);
       
      else if(rank==sendID && N_DofRecv[recvID])
         MPI_Recv(&Recv_Info[rdispl[recvID]], N_DofRecv[recvID], MPI_DOUBLE,  recvID, 200, Comm, &status);
     } //for(j=0;j<N
   } //for(i=0; i<Max_CommunicationSteps

  for(i=0;i<N_Slave;i++)
  {  
      rhs[DofRecv[i]] = Recv_Info[i];
      sol[DofRecv[i]] = Recv_Info[i];
  }
  
//   if(N_SendDof>0)
//      delete [] Send_Info;
//   if(N_Slave>0)
//      delete [] Recv_Info;
  
   MPI_Barrier(MPI_COMM_WORLD);
 
  
  if(TDatabase::ParamDB->SC_SMOOTHER_SCALAR==5)
    CommUpdateH(sol,rhs);
  //CommUpdateAlltoAllv(sol, rhs);
}

void TParFECommunicator3D::CommUpdateM(double *sol, double *rhs)
{
  double timeC;
  int i, j, k, l, N,iter, rank, size;
 int sendID, recvID, *val_send, *val_recv;
 double t1,t2=0.0,t3,tSum=0.;
 
 MPI_Status status;
 MPI_Comm_rank(Comm,&rank);	
 MPI_Comm_size(Comm, &size);

//  if(N_SendDof>0)
//    Send_Info = new double[N_SendDof];
//  if(N_Slave>0)
//    Recv_Info = new double[N_Slave];
 
  t1=MPI_Wtime();
  for(i=0;i<N_SendDof;i++)
  {
   Send_Info[i]=sol[DofSend[i]];
   //Send_Info[DofSendPos[i]]=sol[DofSend[i]];
  }
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
 
  for(i=0; i<Max_CommunicationSteps; i++)
  {
    N= N_CommunicationProcesses[i];
    for(j=0;j<N;j++)
    {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if((rank!=sendID) && (rank!=recvID))  continue;
      
      if(rank==sendID && N_DofSend[recvID]>0)
        MPI_Send(&Send_Info[sdispl[recvID]], N_DofSend[recvID], MPI_DOUBLE,  recvID, 100, Comm);
      
      else if(rank==recvID && N_DofRecv[sendID]>0)
        MPI_Recv(&Recv_Info[rdispl[sendID]], N_DofRecv[sendID], MPI_DOUBLE,  sendID, 100, Comm, &status);


      if(rank==recvID && N_DofSend[sendID])
         MPI_Send(&Send_Info[sdispl[sendID]], N_DofSend[sendID], MPI_DOUBLE,  sendID, 200, Comm);
       
      else if(rank==sendID && N_DofRecv[recvID])
         MPI_Recv(&Recv_Info[rdispl[recvID]], N_DofRecv[recvID], MPI_DOUBLE,  recvID, 200, Comm, &status);
     } //for(j=0;j<N
   } //for(i=0; i<Max_CommunicationSteps

  for(i=0;i<N_Slave;i++)
  {  
      rhs[DofRecv[i]] = Recv_Info[i];
      sol[DofRecv[i]] = Recv_Info[i];
  }
  
//   if(N_SendDof>0)
//      delete [] Send_Info;
//   if(N_Slave>0)
//      delete [] Recv_Info;
  
   MPI_Barrier(MPI_COMM_WORLD);
 
   
  //CommUpdateAlltoAllv(sol, rhs);
}

void TParFECommunicator3D::CommUpdateH(double *sol, double *rhs)
{
    double timeC;
  int i, j, k, l, N,iter, rank, size;
 int sendID, recvID, *val_send, *val_recv;
 double t1,t2=0.0,t3,tSum=0.;
 
 MPI_Status status;
 MPI_Comm_rank(Comm,&rank);	
 MPI_Comm_size(Comm, &size);
 
//  if(N_SendDofH>0)
//    Send_Info = new double[N_SendDofH];
//  if(N_Halo>0)
//    Recv_Info = new double[N_Halo];
 
  t1=MPI_Wtime();
  for(i=0;i<N_SendDofH;i++)
  {
   Send_InfoH[i]=sol[DofSendH[i]];
   //Send_InfoH[DofSendPos[i]]=sol[DofSend[i]];
  }
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
 
  for(i=0; i<Max_CommunicationSteps; i++)
  {
    N= N_CommunicationProcesses[i];
    for(j=0;j<N;j++)
    {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if((rank!=sendID) && (rank!=recvID))  continue;
      
      if(rank==sendID && N_DofSendH[recvID]>0)
        MPI_Send(&Send_InfoH[sdisplH[recvID]], N_DofSendH[recvID], MPI_DOUBLE,  recvID, 100, Comm);
      
      else if(rank==recvID && N_DofRecvH[sendID]>0)
        MPI_Recv(&Recv_InfoH[rdisplH[sendID]], N_DofRecvH[sendID], MPI_DOUBLE,  sendID, 100, Comm, &status);


      if(rank==recvID && N_DofSendH[sendID])
         MPI_Send(&Send_InfoH[sdisplH[sendID]], N_DofSendH[sendID], MPI_DOUBLE,  sendID, 200, Comm);
       
      else if(rank==sendID && N_DofRecvH[recvID])
         MPI_Recv(&Recv_InfoH[rdisplH[recvID]], N_DofRecvH[recvID], MPI_DOUBLE,  recvID, 200, Comm, &status);
     } //for(j=0;j<N
   } //for(i=0; i<Max_CommunicationSteps

   
  for(i=0;i<N_Halo;i++)
  {  
      rhs[DofRecvH[i]] = Recv_InfoH[i];
      sol[DofRecvH[i]] = Recv_InfoH[i];
  }
   
    
//   if(N_SendDofH>0)
//      delete [] Send_Info;
//   if(N_Halo>0)
//      delete [] Recv_Info;
}

void TParFECommunicator3D::CommUpdate(double *sol)
{
 int i, j, k, l, N,iter, rank, size;
 int sendID, recvID, *val_send, *val_recv;
 double t1,t2=0.0,t3,tSum=0.;
  double timeC; 
 MPI_Status status;
 MPI_Comm_rank(Comm,&rank);	
 MPI_Comm_size(Comm, &size);
  

//  if(N_SendDof>0)
//    Send_Info = new double[N_SendDof];
//  if(N_Slave>0)
//    Recv_Info = new double[N_Slave];
 
  t1=MPI_Wtime();
  for(i=0;i<N_SendDof;i++)
  {
   Send_Info[i]=sol[DofSend[i]];
   //Send_Info[DofSendPos[i]]=sol[DofSend[i]];
  }
  t2=MPI_Wtime(); 
  timeC+=(t2-t1);
 
  for(i=0; i<Max_CommunicationSteps; i++)
  {
    N= N_CommunicationProcesses[i];
    for(j=0;j<N;j++)
    {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if((rank!=sendID) && (rank!=recvID))  continue;
      
      if(rank==sendID && N_DofSend[recvID]>0)
        MPI_Send(&Send_Info[sdispl[recvID]], N_DofSend[recvID], MPI_DOUBLE,  recvID, 100, Comm);
      
      else if(rank==recvID && N_DofRecv[sendID]>0)
        MPI_Recv(&Recv_Info[rdispl[sendID]], N_DofRecv[sendID], MPI_DOUBLE,  sendID, 100, Comm, &status);


      if(rank==recvID && N_DofSend[sendID])
         MPI_Send(&Send_Info[sdispl[sendID]], N_DofSend[sendID], MPI_DOUBLE,  sendID, 200, Comm);
       
      else if(rank==sendID && N_DofRecv[recvID])
         MPI_Recv(&Recv_Info[rdispl[recvID]], N_DofRecv[recvID], MPI_DOUBLE,  recvID, 200, Comm, &status);
     } //for(j=0;j<N
   } //for(i=0; i<Max_CommunicationSteps

  for(i=0;i<N_Slave;i++)
  {  
     // rhs[DofRecv[i]] = Recv_Info[i];
      sol[DofRecv[i]] = Recv_Info[i];
  }
  
//   if(N_SendDof>0)
//      delete [] Send_Info;
//   if(N_Slave>0)
//      delete [] Recv_Info;
  
   MPI_Barrier(MPI_COMM_WORLD);
 
  
//   if(TDatabase::ParamDB->SC_SMOOTHER_SCALAR==5)
//     CommUpdateH(sol,rhs);
  //CommUpdateAlltoAllv(sol, rhs);
}

void TParFECommunicator3D::CommUpdateAlltoAllv(double *sol, double *rhs)
{
   int i,rank;
   
  MPI_Comm_rank(Comm,&rank);		

//   Send_Info = new double[N_SendDof];
//   Recv_Info = new double[N_Slave];

  for(i=0;i<N_SendDof;i++)
  {
      Send_Info[i]=sol[DofSend[i]];
  }
  MPI_Alltoallv(Send_Info,N_DofSend,sdispl,MPI_DOUBLE,Recv_Info,N_DofRecv,rdispl,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_Slave;i++)
  {  
    rhs[DofRecv[i]] = Recv_Info[i];
    sol[DofRecv[i]] = Recv_Info[i];
  }
  
//   delete [] Send_Info;
//   delete [] Recv_Info;
  /*
  Send_Info = new double[N_SendDofH];
  Recv_Info = new double[N_Halo];

  for(i=0;i<N_SendDofH;i++)
  {
      Send_Info[i]=sol[DofSendH[i]];
  }
  MPI_Alltoallv(Send_Info,N_DofSendH,sdisplH,MPI_DOUBLE,Recv_Info,N_DofRecvH,rdisplH,MPI_DOUBLE,Comm);
  
  for(i=0;i<N_Halo;i++)
  {  
    rhs[DofRecvH[i]] = Recv_Info[i];
    sol[DofRecvH[i]] = Recv_Info[i];
  }
  
  delete [] Send_Info;
  delete [] Recv_Info;*/
}

void TParFECommunicator3D::SetSlaveDofRows(int *Row,int *KCol,double *Values,double *rhs)
{
  int i, j,begin, end;
  
  for(i=0;i<N_Slave;i++)
    {
      begin=Row[DofRecv[i]];
      end=Row[DofRecv[i]+1];
      for(j=begin;j<end;j++)
      {
         if(KCol[j] == DofRecv[i])
	   Values[j] = 1;
	 else
	   Values[j] = 0;
      }
     // rhs[DofRecv[i]] = 0;
    } 
  //  CommUpdate(rhs);
    /*
    for(i=0;i<N_Halo;i++)
    {
      begin=Row[DofRecvH[i]];
      end=Row[DofRecvH[i]+1];
      for(j=begin;j<end;j++)
      {
         if(KCol[j] == DofRecvH[i])
	   Values[j] = 1;
	 else
	   Values[j] = 0;
      }
      rhs[DofRecvH[i]] = 0;
    }*/
}

void TParFECommunicator3D::CommUpdateReduce(double *sol,double *rhs)
{
  int i,j,rank,size,sendID, recvID,N;
  MPI_Status status;
  MPI_Comm_rank(Comm,&rank);		
  MPI_Comm_size(Comm, &size);
//  Send_Info = new double[N_SendDof];
//  Recv_Info = new double[N_Slave];
  
  for(i=0;i<N_Slave;i++)
  {
      Recv_Info[i]=rhs[DofRecv[i]];
  }
 
//  MPI_Alltoallv(Recv_Info,N_DofRecv,rdispl,MPI_DOUBLE,Send_Info,N_DofSend,sdispl,MPI_DOUBLE,Comm);
  
    
  for(i=0; i<Max_CommunicationSteps; i++)
  {
    N= N_CommunicationProcesses[i];
    for(j=0;j<N;j++)
    {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];
     // if(rank==0) printf("rank %d step i is %d send ID %d recvID %d\n",rank,i,sendID,recvID);

      if((rank!=sendID) && (rank!=recvID))  continue;
     
      if(rank==sendID && N_DofRecv[recvID]>0)
      {
	//printf("rank is %d\n",rank);
        MPI_Send(&Recv_Info[rdispl[recvID]], N_DofRecv[recvID], MPI_DOUBLE,  recvID, 100, Comm);
      }
      
      else if(rank==recvID && N_DofSend[sendID]>0)
        MPI_Recv(&Send_Info[sdispl[sendID]], N_DofSend[sendID], MPI_DOUBLE,  sendID, 100, Comm, &status);


      if(rank==recvID && N_DofRecv[sendID])
         MPI_Send(&Recv_Info[rdispl[sendID]], N_DofRecv[sendID], MPI_DOUBLE,  sendID, 200, Comm);
       
      else if(rank==sendID && N_DofSend[recvID])
         MPI_Recv(&Send_Info[sdispl[recvID]], N_DofSend[recvID], MPI_DOUBLE,  recvID, 200, Comm, &status);
      
     } //for(j=0;j<N
   } //for(i=0; i<Max_CommunicationSteps*/
 
  for(i=0;i<N_SendDof;i++)
  {  
      rhs[DofSend[i]] += Send_Info[i];
      sol[DofSend[i]] += Send_Info[i];
  }

/** MASTER HAS FINAL VALUE **/  
  for(i=0;i<N_SendDof;i++)
  {
      Send_Info[i]=rhs[DofSend[i]];
  }
  
 // MPI_Alltoallv(Send_Info,N_DofSend,sdispl,MPI_DOUBLE,Recv_Info,N_DofRecv,rdispl,MPI_DOUBLE,Comm);
  
    
  for(i=0; i<Max_CommunicationSteps; i++)
  {
    N= N_CommunicationProcesses[i];
    for(j=0;j<N;j++)
    {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];
     // if(rank==0) printf("rank %d step i is %d send ID %d recvID %d\n",rank,i,sendID,recvID);

      if((rank!=sendID) && (rank!=recvID))  continue;
     
      if(rank==sendID && N_DofSend[recvID]>0)
      {
	//printf("rank is %d\n",rank);
        MPI_Send(&Send_Info[sdispl[recvID]], N_DofSend[recvID], MPI_DOUBLE,  recvID, 100, Comm);
      }
      
      else if(rank==recvID && N_DofRecv[sendID]>0)
        MPI_Recv(&Recv_Info[rdispl[sendID]], N_DofRecv[sendID], MPI_DOUBLE,  sendID, 100, Comm, &status);


      if(rank==recvID && N_DofSend[sendID])
         MPI_Send(&Send_Info[sdispl[sendID]], N_DofSend[sendID], MPI_DOUBLE,  sendID, 200, Comm);
       
      else if(rank==sendID && N_DofRecv[recvID])
         MPI_Recv(&Recv_Info[rdispl[recvID]], N_DofRecv[recvID], MPI_DOUBLE,  recvID, 200, Comm, &status);
      
     } //for(j=0;j<N
   } //for(i=0; i<Max_CommunicationSteps*/
  
  for(i=0;i<N_Slave;i++)
  {  
      rhs[DofRecv[i]] = Recv_Info[i];
      sol[DofRecv[i]] = Recv_Info[i];
  }
  
//   delete [] Send_Info;
//   delete [] Recv_Info;  
}


/** Scheduling for communication between processor based on used FESpace */
void TParFECommunicator3D::ScheduleParFEComm3D()
{
 int rank, size;
 int i, j, k, l, m, S_ID, R_ID, N_U, M, N, P, MaxNeibs_All, N_Entries;
 int *pos, *N_DofNeibs_Array, *DofNeibs_Array, N_Update;

 bool UPDATE;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  // first find the N_Neibs of each procesor and its neib ranks
  N_U = FESpace->GetN_DegreesOfFreedom();

  if(rank==0)
   {
    pos = new int[size]; 
    N_DofNeibs_Array = new int[size];
   }

  MPI_Allreduce(&N_Neibs, &MaxNeibs_All, 1, MPI_INT, MPI_MAX, Comm); 

  if(rank==0)
   {
    for(i=0;i<size;i++)
     pos[i] = i*MaxNeibs_All;

    DofNeibs_Array= new int[size*MaxNeibs_All];
   }

   MPI_Gather(&N_Neibs, 1, MPI_INT, N_DofNeibs_Array, 1, MPI_INT, 0, Comm);
   MPI_Gatherv(NeibsRank, N_Neibs, MPI_INT, DofNeibs_Array, N_DofNeibs_Array, pos, MPI_INT, 0, Comm);

   if(rank==0)
    {
     int step, *Entries, *ColInd, *RowPtr, *Collection, MaxPossibleSteps;

     N_Entries = 0;
     for(i=0;i<size;i++)
      N_Entries +=N_DofNeibs_Array[i];

     // adjacency matrix memory allocation
     Entries = new int[N_Entries];
     ColInd  = new int[N_Entries];
     RowPtr = new int[size+1];

     RowPtr[0] = 0;
     m = 0;

     for(i=0;i<size;i++)
      {
       N = N_DofNeibs_Array[i];
       RowPtr[i+1] = RowPtr[i] + N;
       for(j=0;j<N;j++)
        {
         ColInd[m] = DofNeibs_Array[i*MaxNeibs_All + j];
         // assmeble the matrix
         Entries[m] =-1;
         m++;
        }
      }// for(i=0;i<size;

//     // print adjacency matrix
//     for(i=0;i<size;i++)
//      for(j=RowPtr[i];j<RowPtr[i+1];j++)
//        printf("a( %d , %d )=  %d \n",i, ColInd[j], Entries[j]);

    // find number of steps to communicate between neibs (edge coloring alogrithm)
     Collection = new int[size];
     Max_CommunicationSteps = -1;

     N_Update = 0;
     for(i=0;i<size;i++)
      for(j=RowPtr[i];j<RowPtr[i+1];j++)
       if(i<ColInd[j])
        N_Update++;

//         printf("a( %d , %d )=  %d \n",i, ColInd[j], Entries[j]); 

     //if max. order of any vertex is N, the max. no. diff color is N+1
     MaxPossibleSteps = 2*MaxNeibs_All; // excluding own rank

     // printf("Possible steps %d N_Update %d\n", MaxPossibleSteps, N_Update);

     for(step=0;step<MaxPossibleSteps;step++)
      {
       m=0;
       for(k=0;k<size;k++)
        {
         UPDATE = TRUE;
         for(l=0;l<m;l++)
           if(Collection[l] == k)
            {
             UPDATE = FALSE;
             break;
            }

         if(UPDATE)
          {
           Collection[m] = k;
           m++;
          }
         else
          { continue; }
 
         for(l=RowPtr[k];l<RowPtr[k+1];l++)
          {
           // only upper tirangular martix is enough (due to symmetry)
           if( ColInd[l]<k) continue;

           UPDATE = TRUE;
           for(P=0;P<m;P++)
            if(Collection[P] == ColInd[l])
             {
              // ColInd[l] is already in the collection
              UPDATE = FALSE;
              break;
             } 

            if(UPDATE)
             {
              if(Entries[l]<0)
               {
                Entries[l] = step;
                if(Max_CommunicationSteps<step+1) Max_CommunicationSteps=step+1;
                Collection[m] = ColInd[l];
                m++;

                N_Update--;
//                 printf("a( %d , %d )=  %d \n",k, ColInd[l], Entries[l]);
                break;
               }
             }
          } //  for(l=RowPtr[k];
        } // for(k=0;

        if(N_Update==0) break;
       } //  for(step=0;step<

  if(Max_CommunicationSteps<0)
     {Max_CommunicationSteps=0;}

//   printf("Possible steps %d N_Update %d\n", step, N_Update);
  printf("Total number of steps needed to communicate is %d\n",Max_CommunicationSteps);

    // print adjacency matrix
    for(i=0;i<size;i++)
     for(j=RowPtr[i];j<RowPtr[i+1];j++)
      if(i<ColInd[j] && i==18)
       printf("a( %d , %d )=  %d \n",i, ColInd[j], Entries[j]);


  if(Max_CommunicationSteps>0)
   {
    N_ActiveCommuPross = new int[Max_CommunicationSteps];
    SendID_Index= new int[Max_CommunicationSteps*size];
    ReceiveID_Index = new int[Max_CommunicationSteps*size];

    for(j=0;j<Max_CommunicationSteps;j++)
     N_ActiveCommuPross[j] =0;

    for(i=0;i<Max_CommunicationSteps*size;i++)
     {
      SendID_Index[i] = -1;
      ReceiveID_Index[i] = -1;
     }

    for(i=0;i<size;i++)
     for(j=RowPtr[i];j<RowPtr[i+1];j++)
      {
       if(i<ColInd[j])
        {
         if( Entries[j]<0)
          {
           printf("Error in finding DOF communication steps (%d, %d) : %d\n", i, ColInd[j], Entries[j]);
           MPI_Abort(MPI_COMM_WORLD,  0);
          }
         k =  Entries[j];
         SendID_Index[k*size + N_ActiveCommuPross[k]] = i;
         ReceiveID_Index[k*size + N_ActiveCommuPross[k]] =  ColInd[j];
         N_ActiveCommuPross[k]++;
        }
       }

     for(i=0;i<Max_CommunicationSteps;i++)
      {
       printf("Step %d:\n",i);
       N= N_ActiveCommuPross[i];
       for(j=0;j<N;j++)
        printf("Send from process %d to %d -------- and onto\n",
                   SendID_Index[i*size +j], ReceiveID_Index[i*size +j]);
       printf("\n");
      }
     } // if(Max_CommunicationSteps>0)

     delete [] Entries;
     delete [] ColInd;
     delete [] RowPtr; 
     delete [] Collection;
     delete [] pos;
     delete [] N_DofNeibs_Array;
     delete [] DofNeibs_Array;
    }//if(rank==0)


  MPI_Bcast(&Max_CommunicationSteps, 1, MPI_INT, 0, Comm);
     
   if(rank!=0 && Max_CommunicationSteps>0)
    {
     N_ActiveCommuPross = new int[Max_CommunicationSteps];
     SendID_Index= new int[Max_CommunicationSteps*size];
     ReceiveID_Index = new int[Max_CommunicationSteps*size];
    }

  MPI_Bcast(N_ActiveCommuPross, Max_CommunicationSteps, MPI_INT, 0, Comm);
  MPI_Bcast(SendID_Index, Max_CommunicationSteps*size, MPI_INT, 0, Comm);
  MPI_Bcast(ReceiveID_Index, Max_CommunicationSteps*size, MPI_INT,  0, Comm);

  N=0;
  if(Max_CommunicationSteps>0 && N_Neibs>0)
   {
    delete [] NeibsRank;
    NeibsRank = new int[N_Neibs];

    for(i=0;i<Max_CommunicationSteps;i++)
     {
      M = N_ActiveCommuPross[i];

      for(j=0;j<M;j++)
       {
        R_ID = ReceiveID_Index[i*size +j];
        S_ID = SendID_Index[i*size +j];

        UPDATE = FALSE;
         if( (rank==S_ID || rank==R_ID) && S_ID!=R_ID  )
           UPDATE = TRUE;

        if(UPDATE)
         {

          if(S_ID!=rank)
           { NeibsRank[N]=S_ID; }
          else
           { NeibsRank[N]=R_ID; }
          N++;
         } // if(UPDATE)
       } //  for(j=0;j<M;j++)
      } // for(i=0;i<Max_Communicati
   } // if(rank!=0)


  if(N_Neibs!=0 && N_Neibs!=N+1)
   { 
    printf("Rank %d Error in finding N_Neibs  %d\n",rank, N_Neibs);
    MPI_Abort(MPI_COMM_WORLD,  0);
   }

  N_Neibs = N;

  if(N_Neibs)
   {
    IndexOfNeibRank = new int[size];

    for(i=0; i<size; i++)
     IndexOfNeibRank[i] = -1;

    for(i=0; i<N_Neibs; i++)
     IndexOfNeibRank[NeibsRank[i]] = i;
   }

 
     
     
   // only for testing
   // check the N_DepDofmapping
   MPI_Status status;
   for(i=0;i<Max_CommunicationSteps;i++)
   {
    N = N_ActiveCommuPross[i];

//          if(rank==size-1)
//         for(j=0;j<N;j++)
//            printf("Send from process %d to %d -------- and onto\n",
//                    SendID_Index[i*size +j], ReceiveID_Index[i*size +j]); 
//     
//           MPI_Finalize();
//      exit(0);  
    
    for(j=0;j<N;j++)
     {
      S_ID = SendID_Index[i*size +j];
      R_ID = ReceiveID_Index[i*size +j];

      M=0;
      P=0;
      if(rank==S_ID)
       {
        for(k=0; k<N_U; k++)
         for(l=0; l<N_DofRankIndex[k]; l++)
          if(DofRankIndex[k*MaxSubDomainPerDof + l]==R_ID)
           M++;
        MPI_Send(&M, 1, MPI_INT,  R_ID, 100, Comm);
       }
     else if(rank==R_ID)
       {
        for(k=0; k<N_U; k++)
         for(l=0; l<N_DofRankIndex[k]; l++)
          if(DofRankIndex[k*MaxSubDomainPerDof + l]==S_ID)
           M++;

        MPI_Recv(&P, 1, MPI_INT,  S_ID, 100, Comm, &status);

        if(M !=P)
         {
          printf("Error in dept. Dof Mapping, check DofRankIndex in ScheduleParFEComm3D !!!\n");
          printf("Rank %d NeibRank %d P %d M %d!!!\n", rank, S_ID, P, M);
          MPI_Abort(Comm, 0);
         }
       }
     } // for(j=0;j<N;j++)
    //if(rank==TDatabase::ParamDB->Par_P5)
      //printf("\n");
    } //for(i=0;i<Max_CommunicationSteps;i++)

  if(TDatabase::ParamDB->SC_VERBOSE>4)
  if(rank==TDatabase::ParamDB->Par_P0)
   printf("ScheduleParFEComm3D done !!!\n");
   //   MPI_Finalize();
   // exit(0);
} // ScheduleParFEComm3D


/** Mapping of DOFs between subdomains */
void TParFECommunicator3D::MapDofFromNeib3D()
{
 int i, ii, j, jj, k, kk, l, m, n, P, rank, size, N_Cells, M, N;
 int *GlobalNumbers, *BeginIndex, *DOF, ID, N_Joints;
 int *JointDof, Neib_ID, N_JointDOF, *N_SubDJointsOfNeib;
 int Max_N_SubDJointsOfNeib, MaxN_JointDofs, MapType, joint_No;
 int **SendBuf, **RecvBuf, **NeibLocalDof, *MapPair, N1, N2;
 int N_Edges, N_EdgeDOF;
 int N_CrossEdgeNeibs, *CrossEdgeNeibsRank, *CrossEdgeNeibsGlobalNo;
 int *CrossEdgeNeibsLocalEdgeNo, *CrossEdgeNeibsMaptype, *EdgeDof;
 int *N_SubDEdgesOfNeib, Max_N_SubDEdgesOfNeib, **EdgeSendBuf, **NeibEdgeLocalDof;
 int **EdgeRecvBuf, **RecvEdgeLocalDof, *N_CrossVertNeibs, *CrossVertNeibs;
 int N_VertInCell, Max_N_CrossVertNeibs, N_CrossVertNeibCells, *CrossVertNeibCellRank;
 int *CrossVertNeibCellGlobalNo, *CrossVertNeibCellLocVertNo;
 int VertDof;

 bool UPDATE;

 TCollection *Coll;
 TBaseCell *cell, *NeibCell;
 TVertex *Vert;
 TJoint *Joint;
 FE3D FeId;
 TFEDesc3D *FeDesc;
 TEdge *edge;
 TFE3D *FE;
 FEDesc3D FEDesc0, FEDesc1;
 Refinements Ref0, Ref1;
 TFE3DMapper *Mapper;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  BeginIndex = FESpace->GetBeginIndex();
  GlobalNumbers = FESpace->GetGlobalNumbers();
  MaxN_JointDofs = 0;
  Max_N_SubDJointsOfNeib=0;
  Max_N_SubDEdgesOfNeib=0;
  Max_N_CrossVertNeibs=0;

//   MasterDofCounter = new int[size];

//====================================================================
/** mapping of SubDomain Joint Dof start */
//====================================================================

  if(N_Neibs)
   {
    N_SubDJointsOfNeib = new int[N_Neibs];
    memset(N_SubDJointsOfNeib, 0, N_Neibs*SizeOfInt);

   for(ii=0; ii<N_DependentCells; ii++)
    {
     i = DependentCellIndex[ii];
     cell = Coll->GetCell(i);
     ID = cell->GetSubDomainNo();

     N_Joints = cell->GetN_Joints();
     FeId = FESpace->GetFE3D(i, cell);
     FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
     N_JointDOF = FeDesc->GetN_JointDOF();

     // Discspace space
     if(N_JointDOF<=0)
      break; 

      for(j=0; j<N_Joints; j++)
       {
        Joint = cell->GetJoint(j);
        if(MaxN_JointDofs < FeDesc->GetN_JointDOF())
         MaxN_JointDofs = FeDesc->GetN_JointDOF();

         // put neib joint subdomain ID in all DOF in the joint
        if(Joint->GetType() == SubDomainJoint)
         {
          Neib_ID = ((TSubDomainJoint *)Joint)->GetNeibRank();

          k=0;
          while(Neib_ID != NeibsRank[k]) k++;

          N_SubDJointsOfNeib[k]++;
         }// if(Joint->GetType(
       } // for(j=0; j<N_Joint
     } // for(ii=0; ii<N_DependentCells; ii++)

    // find max. of all N_SubDJointsOfNeib
    for(k=0; k<N_Neibs; k++)
     if(Max_N_SubDJointsOfNeib < N_SubDJointsOfNeib[k])
       Max_N_SubDJointsOfNeib = N_SubDJointsOfNeib[k];

 if(Max_N_SubDJointsOfNeib)
  {
   SendBuf = new int*[4];
   RecvBuf = new int*[4];

   for(j=0; j<4; j++)
    {
     SendBuf[j] = new int[N_Neibs*Max_N_SubDJointsOfNeib];
     RecvBuf[j] = new int[N_Neibs*Max_N_SubDJointsOfNeib];
    }

   memset(N_SubDJointsOfNeib, 0, N_Neibs*SizeOfInt);

   // first find SubDomain Joint neibs RefDesc and FeDesc
   for(ii=0; ii<N_DependentCells; ii++)
    {
     i = DependentCellIndex[ii];
     cell = Coll->GetCell(i);
     ID = cell->GetSubDomainNo();

      DOF = GlobalNumbers + BeginIndex[i];
      N_Joints = cell->GetN_Joints();
      FeId = FESpace->GetFE3D(i, cell);
      FE = TFEDatabase3D::GetFE3D(FeId);
      FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
      N_JointDOF = FeDesc->GetN_JointDOF();

      // Discspace space
       if(N_JointDOF<=0)
        break;  

      for(j=0; j<N_Joints; j++)
       {
        Joint = cell->GetJoint(j);

         // put neib joint subdomain ID in all DOF in the joint
        if(Joint->GetType() ==  SubDomainJoint)
         {
          Neib_ID = ((TSubDomainJoint *)Joint)->GetNeibRank();
          N =((TSubDomainJoint *)Joint)->GetNeibLocalJointNo();
          NeibCell = Joint->GetNeighbour(cell);

           k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
           M = k*Max_N_SubDJointsOfNeib + N_SubDJointsOfNeib[k];

           SendBuf[0][M] = NeibCell->GetGlobalCellNo(); // SubDJntNeibGlobalNo
           SendBuf[1][M] = N;                            //  SubDJntNeibLocalEdgeNo
           SendBuf[2][M] = (int)(FE->GetFEDesc3D_ID()); //  SubDJntOwnFEDesc
           SendBuf[3][M] = (int)(cell->GetRefDesc()->GetType()); //  SubDJntOwnRefDesc

           N_SubDJointsOfNeib[k]++;

         }// if(Joint->GetType(
       } // for(j=0; j<N_Joint
     } // for(ii=0; ii<N_DependentCells; ii++)
    } // if(Max_N_SubDJointsOfNeib)
   } // if(N_Neibs)

    // communicate between neibs 
    FECommunicateNeib(SendBuf, Max_N_SubDJointsOfNeib, N_SubDJointsOfNeib, RecvBuf,  N_SubDJointsOfNeib, 4);

//     check the neib subdomain cells are using same FEDisc
   if(Max_N_SubDJointsOfNeib)
    for(i=0; i<N_Neibs; i++)
     {
      N = N_SubDJointsOfNeib[i];
      for(j=0; j<N; j++)
       {
        M = i*Max_N_SubDJointsOfNeib + j;
        ii= RecvBuf[0][M];

        k = GetLocalCellIndex(N_DependentCells, DeptCellGlobalNo, DependentCellIndex, ii);
        cell =  Coll->GetCell(k);

        if(ii != cell->GetGlobalCellNo())
         {
           printf("Error  check MapDofFromNeib3D   \n" );
           printf("Neib Rank %d  CellIndex %d ii, Global No %d   %d \n", NeibsRank[i], k, ii,  Coll->GetCell(k)->GetGlobalCellNo());
           MPI_Finalize();
           exit(0);
         }

        DOF = GlobalNumbers + BeginIndex[k];
        N_Joints = cell->GetN_Joints();
        FeId = FESpace->GetFE3D(k, cell);
        FE = TFEDatabase3D::GetFE3D(FeId);
        FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
        Joint = cell->GetJoint(RecvBuf[1][M]);

        if(Joint->GetType() !=  SubDomainJoint)
         {
          printf("Error  check MapDofFromNeib3D   \n" );
          printf("Neib Rank %d  CellIndex %d SubDomainJoint  %d \n", NeibsRank[i], k, SubDomainJoint);
          MPI_Finalize();
          exit(0);
         }

        FEDesc0 = FE->GetFEDesc3D_ID();
        FEDesc1 = (FEDesc3D)RecvBuf[2][M];
        if(FEDesc0 !=  FEDesc1)
         {
          printf("Error  heterogeneous FE discretisation not yet allowed btween subdomains \n" );
          printf("Rank %d FEDesc0 %d  FEDesc1 %d \n", rank, FEDesc0, FEDesc1);
          MPI_Finalize();
          exit(0);
         }

        Ref0 = cell->GetRefDesc()->GetType();
        Ref1 = (Refinements)RecvBuf[3][M];

        if(Ref0 !=  Ref1)
         {
          printf("Error  NonMatching refinments are not yet allowed btween subdomains \n" );
          printf("Rank %d Ref0 %d  Ref1 %d \n", rank, Ref0, Ref1);
          MPI_Finalize();
          exit(0);
         }
       } // for(j=0; j<N;
     } // for(i=0; i<N_Nei
 
//   if(N_Neibs>0)
//    MasterOfSubDomainDofs = new int[N_JointDOF];

   //mapping starts, first mapping between SubDomain Joints
  if(N_Neibs>0 &&  Max_N_SubDJointsOfNeib>0)
   {
    for(j=0; j<4; j++)
     delete []  SendBuf[j];

    delete [] SendBuf;

   SendBuf = new int*[1];
   NeibLocalDof = new int*[1];
   SendBuf[0] = new int[N_Neibs*Max_N_SubDJointsOfNeib*N_JointDOF];
   NeibLocalDof[0] = new int[N_Neibs*Max_N_SubDJointsOfNeib*N_JointDOF];

   memset(N_SubDJointsOfNeib, 0, N_Neibs*SizeOfInt);

   for(ii=0; ii<N_DependentCells; ii++)
    {
     i = DependentCellIndex[ii];
     cell = Coll->GetCell(i);
     ID = cell->GetSubDomainNo();

      DOF = GlobalNumbers + BeginIndex[i];
      N_Joints = cell->GetN_Joints();
      FeId = FESpace->GetFE3D(i, cell);
      FE = TFEDatabase3D::GetFE3D(FeId);
      FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
      FEDesc0 = FE->GetFEDesc3D_ID();
      Mapper=TFEDatabase3D::GetFE3DMapper(FEDesc0, FEDesc0);

      if(FeDesc->GetN_JointDOF() !=  N_JointDOF)
       {
        printf("Error  N_JointDOF should be equal in all cells   \n" );
        printf("Rank %d  N_JointDOF %d N_JointDOF  %d \n", rank, FeDesc->GetN_JointDOF(), N_JointDOF);
        MPI_Finalize();
        exit(0);
       }

      for(j=0; j<N_Joints; j++)
       {
        Joint = cell->GetJoint(j);

        if(Joint->GetType() ==  SubDomainJoint)
         {
          Neib_ID = ((TSubDomainJoint *)Joint)->GetNeibRank();
          k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
          M = k*Max_N_SubDJointsOfNeib + N_SubDJointsOfNeib[k];

          JointDof = FeDesc->GetJointDOF(j);
          MapType = Joint->GetMapType();
          MapPair = Mapper->GetPairs(MapType);

          if(Mapper->GetN_Pairs() !=  N_JointDOF)
           {
            printf("Error  N_JointDOF should be equal N_Pairs   \n" );
            printf("Rank %d  N_JointDOF %d N_Pairs  %d \n", rank, N_JointDOF, Mapper->GetN_Pairs());
            MPI_Finalize();
            exit(0);
           }

          N1 = M*N_JointDOF;
          for(jj=0; jj<N_JointDOF; jj++)
           {
            N2 = MapPair[2*jj + 1]-N_JointDOF;
            SendBuf[0][N1 + N2] =  DOF[ JointDof[jj] ];
           } // for(jj=0; jj<N_JointDOF; jj

          N_SubDJointsOfNeib[k]++;
         }// if(Joint->GetType(
       } // for(j=0; j<N_Joint
     } // for(ii=0; ii<N_DependentCells; ii++)
   } //if(N_Neibs)
 

//    N2=0;
//   if(rank==TDatabase::ParamDB->Par_P6)
//     for(i=0; i<N_Neibs; i++)
//      {
//       N = N_SubDJointsOfNeib[i];
//       if(NeibsRank[i]==TDatabase::ParamDB->Par_P5)
//        for(j=0; j<N; j++)
//        {
//         M = i*Max_N_SubDJointsOfNeib + j;
//         N1 = M*N_JointDOF;
//         for(jj=0; jj<N_JointDOF; jj++)
//          printf("M %d, Neib Rank %d OwnDof %d \n", N2++, NeibsRank[i],  SendBuf[0][N1 + jj]  );
//        } // for(j=0; j<N_Joints; 
//      }
// 
//  MPI_Finalize();
//  exit(0);

    // communicate JointDof between neibs
    for(i=0; i<N_Neibs; i++)
     N_SubDJointsOfNeib[i] *=N_JointDOF;
    FECommunicateNeib(SendBuf, Max_N_SubDJointsOfNeib*N_JointDOF, N_SubDJointsOfNeib, NeibLocalDof, N_SubDJointsOfNeib, 1);

    for(i=0; i<N_Neibs; i++)
     N_SubDJointsOfNeib[i] /=N_JointDOF;

    if(N_Neibs>0 && Max_N_SubDJointsOfNeib>0)
      memset(N_DeptDofNeibs, 0, N_DeptDofs*SizeOfInt);

   if(Max_N_SubDJointsOfNeib>0)
    for(i=0; i<N_Neibs; i++)
     {
      ID = NeibsRank[i];
      N = N_SubDJointsOfNeib[i];
      for(j=0; j<N; j++)
       {
        M = i*Max_N_SubDJointsOfNeib + j;

        ii= RecvBuf[0][M]; // global cell no
        joint_No = RecvBuf[1][M]; // local joint no of global cell

        k = GetLocalCellIndex(N_DependentCells, DeptCellGlobalNo, DependentCellIndex, ii);

        DOF = GlobalNumbers + BeginIndex[k];
        cell =  Coll->GetCell(k);
        FeId = FESpace->GetFE3D(k, cell);
        FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
        JointDof = FeDesc->GetJointDOF(joint_No);
        N1 = M*N_JointDOF;

//         //identify the master of all dofs on this joint
//         ShapeDesc = cell->GetShapeDesc();
//         ShapeDesc->GetFaceVertex(TmpFV, TmpLen, MaxLen);
//         N_FV = TmpLen[joint_No];
// 
// 
//         for(jj=0; jj<N_JointDOF; jj++)
//          MasterOfSubDomainDofs[jj] = -1;
// 
//         for(jj=0; jj<size; jj++)
//          MasterDofCounter[jj] = 0;
// 
//         for(jj=0; jj<N_FV; jj++)
//          {
//           m = TmpFV[joint_No*MaxLen+jj];  // local vert No
//           MasterID = (cell->GetVertex(m))->GetRank_ID();
//           n = FeDesc->GetVertDOF(m);
//           kk=0;
//           while(JointDof[kk] != n ) kk++;
//           MasterOfSubDomainDofs[kk] = MasterID;
//           MasterDofCounter[MasterID]++;
// 
//           if(MasterID!=rank || MasterID!=ID)
//             printf("Rank %d  Neib %d MasterID %d \n", rank, ID, MasterID);
//          }
// 
//         // vert rank with min no. of indices will take all unassined (inner) dofs
//         m=-1;
//         for(jj=0; jj<size; jj++)
//          if(m < MasterDofCounter[jj])
//           {
//            m = MasterDofCounter[jj];
//            MasterID = jj;
//           }

        for(jj=0; jj<N_JointDOF; jj++)
         {
          m = DOF[JointDof[jj]];
//           n = GetDeptIndex(N_DeptDofs, DeptDofs, m);
          n = DepDofIndexOfLocDof[m];
          P = N_DeptDofNeibs[n];

          UPDATE=TRUE;
          for(l=0; l<P; l++)
           if(DeptDofNeibRanks[n*MaxSubDomainPerDof + l] == ID)
           {
            UPDATE=FALSE;
            break;
           }

          if(UPDATE)
           {
            DeptDofNeibRanksLocalDOF[n*MaxSubDomainPerDof + P] = NeibLocalDof[0][N1 + jj];
            DeptDofNeibRanks[n*MaxSubDomainPerDof + P] = ID;
            N_DeptDofNeibs[n]++;

//             if(MasterOfSubDomainDofs[jj]==-1)
//              {
// //               DeptDofMaster[n] = MasterID;
//               }
//             else
//              { DeptDofMaster[n] = MasterOfSubDomainDofs[jj]; }

           }
         } //  for(jj=0; jj<N_JointDOF;
       } // for(j=0; j<N;
     } // for(i=0; i<N_Nei

  if(N_Neibs>0 && Max_N_SubDJointsOfNeib>0)
   {
    delete []  SendBuf[0];

    for(j=0; j<4; j++)
     delete []  RecvBuf[j];

    delete [] SendBuf;
    delete [] RecvBuf;
   }


//====================================================================
/** mapping of SubDomain Joint Dof --- end */

/** now map the cross edge DOFs --- start */
// =====================================================================

  if(N_Neibs)
   {
    N_SubDEdgesOfNeib = new int[N_Neibs];
    memset(N_SubDEdgesOfNeib, 0, N_Neibs*SizeOfInt);

    /** for cross edge cells set clioboard (if any) */
    for(ii=0; ii<N_DependentCells; ii++)
     {
      i = DependentCellIndex[ii];
      cell = Coll->GetCell(i);
      if(cell->IsCrossEdgeCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(ii=0; ii<N_DependentCells; ii++)


   for(ii=0; ii<N_DependentCells; ii++)
    {
     i = DependentCellIndex[ii];
     cell = Coll->GetCell(i);
     ID = cell->GetSubDomainNo();

//      if(rank!=ID)
//       continue;

     if(cell->IsCrossEdgeCell())
      {
       DOF = GlobalNumbers + BeginIndex[i];
       N_Edges=cell->GetN_Edges();
       FeId = FESpace->GetFE3D(i, cell);
       FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
       N_EdgeDOF = FeDesc->GetN_EdgeDOF();

       // Discspace and non-conforming space
       if(N_EdgeDOF<=0)
        break;
       
       for(j=0;j<N_Edges;j++)
        {
         edge = cell->GetEdge(j);

          if( edge->IsSubDomainCrossEdge() &&  (edge->GetClipBoard()==-1) )
           {
            edge->SetClipBoard(5);

            edge->GetCrossEdgeNeibs(N_CrossEdgeNeibs, CrossEdgeNeibsRank, CrossEdgeNeibsGlobalNo,
                                    CrossEdgeNeibsLocalEdgeNo, CrossEdgeNeibsMaptype);

            for(jj=0;jj<N_CrossEdgeNeibs;jj++)
             {
              Neib_ID = CrossEdgeNeibsRank[jj];
              k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
              N_SubDEdgesOfNeib[k]++;
             } // for(jj=0;jj<N_CrossEdgeNe
            } // if( edge->IsSubDomainCrossEdge() &&  (edge
         } // for(j=0;j<N_Edges
      } //  if(cell->IsCrossEdgeCell())
    } // for(ii=0; ii<N_DependentCells; ii

    // find max. of all N_SubDEdgesOfNeib
    for(k=0; k<N_Neibs; k++)
     if(Max_N_SubDEdgesOfNeib < N_SubDEdgesOfNeib[k])
       Max_N_SubDEdgesOfNeib = N_SubDEdgesOfNeib[k];

   memset(N_SubDEdgesOfNeib, 0, N_Neibs*SizeOfInt);
    
//    printf("  rank %d  N_SubDEdgesOfNeib  %d\n", rank, Max_N_SubDEdgesOfNeib); 

   if(Max_N_SubDEdgesOfNeib>0)
    {
     EdgeSendBuf = new int*[2];
     EdgeRecvBuf = new int*[2];
     NeibEdgeLocalDof = new int*[1];
     RecvEdgeLocalDof = new int*[1];

     EdgeSendBuf[0] = new int[N_Neibs*Max_N_SubDEdgesOfNeib];
     EdgeSendBuf[1] = new int[N_Neibs*Max_N_SubDEdgesOfNeib];
     EdgeRecvBuf[0] = new int[N_Neibs*Max_N_SubDEdgesOfNeib];
     EdgeRecvBuf[1] = new int[N_Neibs*Max_N_SubDEdgesOfNeib];

     NeibEdgeLocalDof[0] = new int[N_Neibs*Max_N_SubDEdgesOfNeib*N_EdgeDOF];
     RecvEdgeLocalDof[0] = new int[N_Neibs*Max_N_SubDEdgesOfNeib*N_EdgeDOF];

    /** for cross edge cells set clioboard (if any) */
    for(ii=0; ii<N_DependentCells; ii++)
     {
      i = DependentCellIndex[ii];
      cell = Coll->GetCell(i);
      if(cell->IsCrossEdgeCell())
       {
        N_Edges=cell->GetN_Edges();
        for(j=0;j<N_Edges;j++)
         (cell->GetEdge(j))->SetClipBoard(-1);
       } // if(cell->IsDependentCell())
     }// for(i=0;i<N_OwnCel

   for(ii=0; ii<N_DependentCells; ii++)
    {
     i = DependentCellIndex[ii];
     cell = Coll->GetCell(i);
     ID = cell->GetSubDomainNo();

//      if(rank!=ID)
//       continue;

     if(cell->IsCrossEdgeCell())
      {
       DOF = GlobalNumbers + BeginIndex[i];
       N_Edges=cell->GetN_Edges();
       FeId = FESpace->GetFE3D(i, cell);
       FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
       N_EdgeDOF = FeDesc->GetN_EdgeDOF();

       for(j=0;j<N_Edges;j++)
        {
          edge = cell->GetEdge(j);

          if( edge->IsSubDomainCrossEdge() &&  (edge->GetClipBoard()==-1) )
           {
            edge->SetClipBoard(5);
            EdgeDof = FeDesc->GetEdgeDOF(j);

            edge->GetCrossEdgeNeibs(N_CrossEdgeNeibs, CrossEdgeNeibsRank, CrossEdgeNeibsGlobalNo,
                                    CrossEdgeNeibsLocalEdgeNo, CrossEdgeNeibsMaptype);

            for(jj=0;jj<N_CrossEdgeNeibs;jj++)
             {
              Neib_ID = CrossEdgeNeibsRank[jj];
              k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
              M = k*Max_N_SubDEdgesOfNeib + N_SubDEdgesOfNeib[k];

              EdgeSendBuf[0][M] = CrossEdgeNeibsGlobalNo[jj];
              EdgeSendBuf[1][M] = CrossEdgeNeibsLocalEdgeNo[jj];
 
              MapType = CrossEdgeNeibsMaptype[jj];
              N1 = M*N_EdgeDOF;

              for(kk=0; kk<N_EdgeDOF; kk++)
               if(MapType==1)
                { NeibEdgeLocalDof[0][N1 + kk] =  DOF[ EdgeDof[kk] ]; }
               else if(MapType==-1)
                { NeibEdgeLocalDof[0][N1 + N_EdgeDOF-1 - kk] =  DOF[ EdgeDof[kk] ]; }
               else
                {
                 printf("Error  CrossEdgeNeibsMaptype rank %d\n", rank );
                 MPI_Finalize();
                 exit(0);
                }

              N_SubDEdgesOfNeib[k]++;

             } // for(jj=0;jj<N_CrossEdgeNe
            } // if( edge->IsSubDomainCrossEdge() &&  (edge
         } // for(j=0;j<N_Edges
      } //  if(cell->IsCrossEdgeCell())
     } // for(ii=0; ii<N_DependentCells; ii

//   N2=0;
//   if(rank==TDatabase::ParamDB->Par_P6)
//     for(i=0; i<N_Neibs; i++)
//      {
//       N = N_SubDEdgesOfNeib[i];
//       if(NeibsRank[i]==TDatabase::ParamDB->Par_P5)
//        for(j=0; j<N; j++)
//        {
//         M = i*Max_N_SubDEdgesOfNeib + j;
//         printf("Index %d, Rank %d, NeibRank %d NeibGlobal No %d Edge %d \n", 
// 	       N2++, rank, NeibsRank[i],  EdgeSendBuf[0][M], EdgeSendBuf[1][M] );
//        } // for(j=0; j<N_Joints; 
//      } 
     
    } // if(Max_N_SubDEdgesOfNeib>0)    
   } //  if(N_Neibs)
   
  
   // communicate between neibs 
   FECommunicateNeib(EdgeSendBuf, Max_N_SubDEdgesOfNeib, N_SubDEdgesOfNeib, EdgeRecvBuf, N_SubDEdgesOfNeib, 2);

    //check the cross edge neib 
   if(Max_N_SubDEdgesOfNeib>0)    
    for(i=0; i<N_Neibs; i++)
     {
      N = N_SubDEdgesOfNeib[i];
      for(j=0; j<N; j++)
       {
        M = i*Max_N_SubDEdgesOfNeib + j;
        ii= EdgeRecvBuf[0][M];

        k = GetLocalCellIndex(N_DependentCells, DeptCellGlobalNo, DependentCellIndex, ii);
        cell =  Coll->GetCell(k);

        if(ii != cell->GetGlobalCellNo())
         {
           printf("Error  check MapDofFromNeib3D SubDEdgesOfNeib  \n" );
           printf("Neib Rank %d  CellIndex %d ii, Global No %d   %d \n", NeibsRank[i], k, ii,  Coll->GetCell(k)->GetGlobalCellNo());
           MPI_Finalize();
           exit(0);
         }

         edge = cell->GetEdge(EdgeRecvBuf[1][M]);

        if(!(edge->IsSubDomainCrossEdge()))
         {
          printf("Error  check MapDofFromNeib3D SubDEdgesOfNeib  \n" );
          printf("Neib Rank %d  edge must be SubDomainCrossEdge \n", NeibsRank[i]);
          MPI_Finalize();
          exit(0);
         }
       } // for(j=0; j<N;
      } // for(i=0; i<N_Nei   

//    if(rank==TDatabase::ParamDB->Par_P5)
//    for(i=0; i<N_DeptDofs; i++)
//     printf("Rankk %d i %d  DeptDofs %d \n",rank, i, DeptDofs[i] ); 

    // communicate JointDof between neibs
   if(Max_N_SubDEdgesOfNeib>0)
    for(i=0; i<N_Neibs; i++)
     N_SubDEdgesOfNeib[i] *=N_EdgeDOF;
    FECommunicateNeib(NeibEdgeLocalDof, Max_N_SubDEdgesOfNeib*N_EdgeDOF, N_SubDEdgesOfNeib, RecvEdgeLocalDof, N_SubDEdgesOfNeib, 1);

   if(Max_N_SubDEdgesOfNeib>0)
    for(i=0; i<N_Neibs; i++)
     N_SubDEdgesOfNeib[i] /=N_EdgeDOF;

    //put the cross edge neibs local DOF
    if(Max_N_SubDEdgesOfNeib>0)
    for(i=0; i<N_Neibs; i++)
     {
      ID = NeibsRank[i];
      N = N_SubDEdgesOfNeib[i];
      for(j=0; j<N; j++)
       {
        M = i*Max_N_SubDEdgesOfNeib + j;
        ii= EdgeRecvBuf[0][M];

        k = GetLocalCellIndex(N_DependentCells, DeptCellGlobalNo, DependentCellIndex, ii);
        cell =  Coll->GetCell(k);

        FeId = FESpace->GetFE3D(k, cell);
        FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
        N_EdgeDOF = FeDesc->GetN_EdgeDOF();
        EdgeDof = FeDesc->GetEdgeDOF(EdgeRecvBuf[1][M]);

        DOF = GlobalNumbers + BeginIndex[k];
        N1 = M*N_EdgeDOF;

//         //start vert rank will be master for all but end dofs in this edge
//         ShapeDesc = cell->GetShapeDesc();
//         ShapeDesc->GetEdgeVertex(TmpFV);
//         m = TmpFV[2*EdgeRecvBuf[1][M]];
//         MasterID = (cell->GetVertex(m))->GetRank_ID();
// 
//         m = TmpFV[2*EdgeRecvBuf[1][M] +1];
//         MasterID_1 = (cell->GetVertex(m))->GetRank_ID();
// 
//        //max rank will take all dofs in this edge
//        if(MasterID < MasterID_1)
//          MasterID = MasterID_1;

        for(kk=0; kk<N_EdgeDOF; kk++)
         {
          m = DOF[EdgeDof[kk]];
//           n = GetDeptIndex(N_DeptDofs, DeptDofs, m);
          n = DepDofIndexOfLocDof[m];  
          P = N_DeptDofNeibs[n];

          UPDATE=TRUE;
          for(l=0; l<P; l++)
           if(DeptDofNeibRanks[n*MaxSubDomainPerDof + l] == ID)
           {
            UPDATE=FALSE;
            break;
           }

          if(UPDATE)
           {
            DeptDofNeibRanksLocalDOF[n*MaxSubDomainPerDof + P] = RecvEdgeLocalDof[0][N1 + kk] ;
            DeptDofNeibRanks[n*MaxSubDomainPerDof + P] = ID;
            N_DeptDofNeibs[n]++;
//             DeptDofMaster[n] = MasterID;
           }
         } //  for(kk=0; kk<N_EdgeDOF; kk+
       } // for(j=0; j<N;
      } // for(i=0; i<N_Nei

  if(N_Neibs)
   delete [] N_SubDEdgesOfNeib;

  if(Max_N_SubDEdgesOfNeib>0)
   {

    delete [] EdgeSendBuf[0];
    delete [] EdgeSendBuf[1];

    delete [] EdgeRecvBuf[0];
    delete [] EdgeRecvBuf[1];

    delete [] NeibEdgeLocalDof[0];
    delete [] RecvEdgeLocalDof[0];

    delete [] EdgeSendBuf;
    delete [] EdgeRecvBuf;
    delete [] NeibEdgeLocalDof;
    delete [] RecvEdgeLocalDof;
   }

// ====================================================================
/** now map the cross edge DOFs --- end */
/** now map the cross vertex DOFs  --- start  */
// ====================================================================
  if(N_Neibs)
   {
    N_CrossVertNeibs = new int[N_Neibs];
    memset(N_CrossVertNeibs, 0, N_Neibs*SizeOfInt);
   }

//     printf("MapDofFromNeib3D done %d\n", N_DependentCells);

   for(ii=0; ii<N_DependentCells; ii++)
    {
     i = DependentCellIndex[ii];
     cell = Coll->GetCell(i);

     if(cell->IsCrossVertexCell())
      {
       N_VertInCell = cell->GetN_Vertices();
       for(j=0; j<N_VertInCell; j++)
         (cell->GetVertex(j))->SetClipBoard(-1);
      } // if(cell->IsCrossVertexCell())
     } // for(ii=0; ii<N_DependentCells; ii++)

//    if(rank==TDatabase::ParamDB->Par_P0)
//     printf("rank %d MapDofFromNeib3D done %d\n", rank, N_DependentCells);
// 
//  MPI_Finalize();
//  exit(0);


   for(ii=0; ii<N_DependentCells; ii++)
    {
     i = DependentCellIndex[ii];
     cell = Coll->GetCell(i);
     ID = cell->GetSubDomainNo();

//      if(rank!=ID)
//       continue;

     if(cell->IsCrossVertexCell())
       {
        N_VertInCell = cell->GetN_Vertices();
        FeId = FESpace->GetFE3D(i, cell);
        FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
        M = FeDesc->GetN_VertDOF();

       // Discspace and non-conforming space
       if(M<=0)
        break;

        for(j=0;j<N_VertInCell;j++)
         {
          Vert=cell->GetVertex(j);

          if( Vert->IsCrossVert() &&  (Vert->GetClipBoard()==-1) )
           {
            Vert->SetClipBoard(5);
            Vert->GetCrossNeibsInfo(N_CrossVertNeibCells, CrossVertNeibCellRank, 
                                    CrossVertNeibCellGlobalNo, CrossVertNeibCellLocVertNo);

             for(jj=0; jj<N_CrossVertNeibCells; jj++)
              {
               Neib_ID = CrossVertNeibCellRank[jj];
               k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
               N_CrossVertNeibs[k]++;
               // if(rank==6)
               // printf("Rankk %d N_Neibs %d Neib_ID %d\n",rank,  N_CrossVertNeibs[k], Neib_ID);
              } // for(k=0; k<N_CrossVertNeibCells
           } //  if( Vert->IsCrossVert() &&  (V
         } //  for(j=0;j<N_VertIn
       } // if(cell->IsCrossVertexCell())
   } //  for(i=0; i<N_Cells




   // find max. of all Max_N_CrossVertNeibs
   for(k=0; k<N_Neibs; k++)
    if(Max_N_CrossVertNeibs < N_CrossVertNeibs[k])
      Max_N_CrossVertNeibs = N_CrossVertNeibs[k];

  if(N_Neibs)
   memset(N_CrossVertNeibs, 0, N_Neibs*SizeOfInt);

   if(Max_N_CrossVertNeibs)
    {
     SendBuf = new int*[3];
     RecvBuf = new int*[3];

     for(j=0; j<3; j++)
      {
       SendBuf[j] = new int[N_Neibs*Max_N_CrossVertNeibs];
       RecvBuf[j] = new int[N_Neibs*Max_N_CrossVertNeibs];
      }

   for(ii=0; ii<N_DependentCells; ii++)
    {
     i = DependentCellIndex[ii];
     cell = Coll->GetCell(i);

     if(cell->IsCrossVertexCell())
      {
       N_VertInCell = cell->GetN_Vertices();
       for(j=0; j<N_VertInCell; j++)
        (cell->GetVertex(j))->SetClipBoard(-1);
      } // for(j=0; j<N_Vert
     } // for(ii=0; ii<N_DependentCells; ii++)

   for(ii=0; ii<N_DependentCells; ii++)
    {
     i = DependentCellIndex[ii];
     cell = Coll->GetCell(i);
     ID = cell->GetSubDomainNo();
/*
     if(rank!=ID)
      continue;*/

     if(cell->IsCrossVertexCell())
       {
        N_VertInCell = cell->GetN_Vertices();
        FeId = FESpace->GetFE3D(i, cell);
        FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);
        M = FeDesc->GetN_VertDOF();
        DOF = GlobalNumbers + BeginIndex[i];

        for(j=0;j<N_VertInCell;j++)
         {
          Vert=cell->GetVertex(j);

          if( Vert->IsCrossVert() &&  (Vert->GetClipBoard()==-1.) )
           {
            Vert->SetClipBoard(5);
            VertDof = FeDesc->GetVertDOF(j);
            Vert->GetCrossNeibsInfo(N_CrossVertNeibCells, CrossVertNeibCellRank, 
                                    CrossVertNeibCellGlobalNo, CrossVertNeibCellLocVertNo);
             for(jj=0; jj<N_CrossVertNeibCells; jj++)
              {
               Neib_ID = CrossVertNeibCellRank[jj];
               k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID

               M = k*Max_N_CrossVertNeibs + N_CrossVertNeibs[k];

               SendBuf[0][M] = CrossVertNeibCellGlobalNo[jj];
               SendBuf[1][M] = CrossVertNeibCellLocVertNo[jj];
               SendBuf[2][M] = DOF[VertDof];

               N_CrossVertNeibs[k]++;
              } // for(k=0; k<N_CrossVertNeibCells
           } //  if( Vert->IsCrossVert() &&  (V
         } //  for(j=0;j<N_VertIn
       } // if(cell->IsCrossVertexCell())
     } //  for(i=0; i<N_Cells
    } // if(Max_N_CrossVertNeibs)

//   N2=0;
//   if(rank==TDatabase::ParamDB->Par_P6)
//     for(i=0; i<N_Neibs; i++)
//      {
//       N = N_CrossVertNeibs[i];
//       
//            printf("rank %d  N_CrossVertNeibs  %d \n", rank, N);     
//       
//       
//       if(NeibsRank[i]==TDatabase::ParamDB->Par_P5)
//        for(j=0; j<N; j++)
//        {
//         M = i*Max_N_CrossVertNeibs + j;
//         printf("M %d, Neib Rank %d Global No %d vert %d \n", N2++, NeibsRank[i],  SendBuf[0][M], SendBuf[1][M] );
//        } // for(j=0; j<N_Joints; 
//      }

//  MPI_Finalize();
//  exit(0);

    // communicate between neibs 
    FECommunicateNeib(SendBuf, Max_N_CrossVertNeibs, N_CrossVertNeibs, RecvBuf, N_CrossVertNeibs, 3);

   if(Max_N_CrossVertNeibs>0)
    {
     for(i=0; i<N_Neibs; i++)
     {
      ID = NeibsRank[i];
      N = N_CrossVertNeibs[i];
      for(j=0; j<N; j++)
       {
        M = i*Max_N_CrossVertNeibs + j;
        ii= RecvBuf[0][M];

        k = GetLocalCellIndex(N_DependentCells, DeptCellGlobalNo, DependentCellIndex, ii);
        cell =  Coll->GetCell(k);

        if(!( cell->IsCrossVertexCell() ))
         {
          printf("Error  check cell must be Cross Vertex Cell  \n" );
//           printf("Rank: %d Neib Rank %d  own global cell %d neib global cell %d local vert no %d \n", rank, NeibsRank[i], cell->GetGlobalCellNo(), ii, RecvBuf[1][M]);
           
          printf("Rank: %d localcellNo %d GlobalcellNo %d  Neib Rank %d  \n",rank, k, cell->GetGlobalCellNo(),  NeibsRank[i] );
          MPI_Finalize();
          exit(0);
         }

        FeId = FESpace->GetFE3D(k, cell);
        FeDesc = TFEDatabase3D::GetFEDesc3DFromFE3D(FeId);

//         MasterID = (cell->GetVertex(RecvBuf[1][M]))->GetRank_ID();

        // should be only one DOF on each vert
        if(FeDesc->GetN_VertDOF())
         {
          VertDof = FeDesc->GetVertDOF(RecvBuf[1][M]);

          DOF = GlobalNumbers + BeginIndex[k];
          m = DOF[VertDof];
//           n = GetDeptIndex(N_DeptDofs, DeptDofs, m);
          n = DepDofIndexOfLocDof[m];	  
          P = N_DeptDofNeibs[n];

          UPDATE=TRUE;
          for(l=0; l<P; l++)
           if(DeptDofNeibRanks[n*MaxSubDomainPerDof + l] == ID)
           {
            UPDATE=FALSE;
            break;
           }

          if(UPDATE)
           {
            DeptDofNeibRanksLocalDOF[n*MaxSubDomainPerDof + P] = RecvBuf[2][M];
            DeptDofNeibRanks[n*MaxSubDomainPerDof + P] = ID;
            N_DeptDofNeibs[n]++;
//             DeptDofMaster[n] = MasterID;
           }
          } //  if(FeDesc->GetN_VertDOF()) 
       } // for(j=0; j<N;
      } // for(i=0; i<N_Nei

     for(j=0; j<3; j++)
      {
       delete [] SendBuf[j];
       delete [] RecvBuf[j];
      }

     delete [] SendBuf;
     delete [] RecvBuf;
    } // if(Max_N_CrossVertNeibs)   

  if(N_Neibs)
   {
    delete [] N_CrossVertNeibs;
    delete [] NeibLocalDof[0];

    delete [] NeibLocalDof;
    delete [] N_SubDJointsOfNeib;
//     delete [] MasterOfSubDomainDofs;
   }

//  delete [] MasterDofCounter;

//    if(rank==TDatabase::ParamDB->Par_P5)
//    for(i=0; i<N_DeptDofs; i++)
//     if(N_DeptDofNeibs[i]- N_DeptDofNeibs[i] )
//     printf("Rankk %d i %d  N_DeptDofNeibs %d test %d \n",rank, i, N_DeptDofNeibs[i], N_DeptDofNeibs[i] ); 

//  M=0;
//   if(rank==TDatabase::ParamDB->Par_P6)
//    for(i=0; i<N_DeptDofs; i++)
//       for(j=0; j<N_DeptDofNeibs[i]; j++)
//     if(DeptDofNeibRanks[i*MaxSubDomainPerDof + j]==TDatabase::ParamDB->Par_P5)
//         printf("Rankk %d index %d dof %d Neib %d NeibDof %d  \n",
// 	   rank,  M++, DeptDofs[i], DeptDofNeibRanks[i*MaxSubDomainPerDof + j], DeptDofNeibRanksLocalDOF[i*MaxSubDomainPerDof + j] );

//  MPI_Finalize();
//  exit(0);

  if(TDatabase::ParamDB->SC_VERBOSE>4)
   if(rank==TDatabase::ParamDB->Par_P0)
    printf("MapDofFromNeib3D done !!!\n");

//  MPI_Finalize();
//  exit(0);

} //  TParFECommunicator3D::MapDofFromNeib3D()

/** Construct Global Dof From Neib dofs, min. rank of dept. dof will be their master*/
void TParFECommunicator3D::ConstructGlobalDofFromNeib3D()
{
 int i, j, k, m1, m2, M, N, rank, size, ID, Neib_ID, N_Cells;
 int N_U, GblDOfStart;
 int N_SlaveDof, N_MasterDof, N_SelfDof;
 int *N_SendDofs, *N_RecevDofs, **sendbuf, **recevbuf;
 int *BeginIndex, N_LocDof, *GlobalNumbers, *DofMarker, *pos, *DipDeoDof;

 bool Found;
 TCollection *Coll;

   MPI_Comm_rank(Comm, &rank);
   MPI_Comm_size(Comm, &size);

   N_U = FESpace->GetN_DegreesOfFreedom();
   Coll = FESpace->GetCollection();
   N_Cells = Coll->GetN_Cells();
   BeginIndex = FESpace->GetBeginIndex();
   N_LocDof = BeginIndex[1] - BeginIndex[0];  // assume that  all cells in fespace have same FE
   GlobalNumbers = FESpace->GetGlobalNumbers();

   N_SelfDof = 0;  // dof not connected with neib subdomains
   N_MasterDof = 0;// dept. dof, but this process is the master
   N_SlaveDof = 0; // dept. dof, but master is neib process
   N_OwnDof = 0; // N_SelfDof + N_MasterDof

    for(i=0;i<N_DeptDofs;i++)
     {
      DipDeoDof = DeptDofNeibRanks + (i*MaxSubDomainPerDof);

      DeptDofMaster[i] = rank;
      ID = rank;

      // find the small rank index among this dept. dofs
      for(j=0;j<N_DeptDofNeibs[i];j++)
       if(ID < DipDeoDof[j])
        ID = DipDeoDof[j];

      if(ID!=rank)
       {
        DeptDofMaster[i] = ID;
        N_SlaveDof++;
       }

//        if(rank==TDatabase::ParamDB->Par_P5 && ID==TDatabase::ParamDB->Par_P6)
//          printf("Error: rank %d Dof %d masterrank %d\n", rank, DeptDofs[i],  DeptDofMaster[i]);
     } //for(i=0;i<N_DeptDofs;

//   MPI_Finalize();
//   exit(0);  
//      if(rank==TDatabase::ParamDB->Par_P5)
//       {
//        i = GetDeptIndex(N_DeptDofs, DeptDofs, 1);
// 
//        if(rank==TDatabase::ParamDB->Par_P5)
//          printf("Error: rank %d Dof %d masterrank %d\n", rank, 12,  DeptDofMaster[i]);
//      }

   N_SelfDof = N_U - N_DeptDofs;
   N_MasterDof = N_DeptDofs - N_SlaveDof;
   N_OwnDof = N_SelfDof + N_MasterDof;

// printf("Rank %d, N_U %d N_OwnDof %d\n", rank, N_U, N_OwnDof+N_SlaveDof);


//   if(rank==0)
//    {
//     printf("Rank %d, N_OwnDof %d\n", rank, N_U);
//     N_OwnDof = 0;
//    }

  MPI_Allreduce(&N_OwnDof, &N_GlobalDegreesOfFreedom, 1, MPI_INT, MPI_SUM, Comm);

//   FESpace->SetN_GlobalDegreesOfFreedom(N_OwnDofSum);

//   N_OwnDofSum = FESpace->GetN_GlobalDegreesOfFreedom();
// 
//   if(rank==0)
    printf("Rank %d, N_OwnDof %d\n", rank, N_GlobalDegreesOfFreedom);
// 
//   MPI_Finalize();
//   exit(0);  

  N_DistDofAll = new int[size];

  MPI_Allgather(&N_OwnDof, 1, MPI_INT, N_DistDofAll, 1, MPI_INT, Comm);

   GblDOfStart = 0;
   for(i=0;i<rank;i++)
    GblDOfStart +=N_DistDofAll[i];

   GlobalDofOfLocalDof = new int[N_U];
   for(i=0;i<N_U;i++)
    GlobalDofOfLocalDof[i] = -1;


   DofMarker = new int[N_U];
   memset(DofMarker, 0, N_U*SizeOfInt);

   // first fill the dept. dof
   for(i=0;i<N_DeptDofs;i++)
    if(DeptDofMaster[i] == rank)
      DofMarker[DeptDofs[i]] = 2;

   for(i=0;i<N_U;i++)
    if(N_DofRankIndex[i]==1)
      DofMarker[i] = 1;


   if(N_OwnDof)
    {
     OwnDofs = new int[N_OwnDof];

//      // first fill the dept. dof
//      m1=GblDOfStart+N_SelfDof;
//      for(i=0;i<N_DeptDofs;i++)
//       if(DeptDofMaster[i] == rank)
//         GlobalDofOfLocalDof[DeptDofs[i]] = m1++;
// 
// 
//      m1=GblDOfStart;
//      for(i=0;i<N_U;i++)
//       if(N_DofRankIndex[i]==1)
//         GlobalDofOfLocalDof[i] = m1++;

     m1=GblDOfStart;
     for(i=0;i<N_U;i++)
      if(DofMarker[i]>0)
       {
//  if(rank==TDatabase::ParamDB->Par_P5)
//   printf("Row %d:  GlobalDofOfLocalDof  %d  \n", i,  m1 );
        GlobalDofOfLocalDof[i] = m1++;
       }

     if(m1-GblDOfStart !=  N_OwnDof)
      {
       printf("Error: in GlobalDOf construction Rank %d i %d \n",rank, i);
       MPI_Finalize();
       exit(0);
      }

     m1 = 0;
     for(i=0;i<N_U;i++)
      if(DofMarker[i]>0) // self or master
        OwnDofs[m1++] = i;

     if(m1 !=  N_OwnDof)
      {
       printf("Error: in GlobalDOf construction Rank %d i %d \n",rank, i);
       MPI_Finalize();
       exit(0);
      }


    } // if(N_OwnDof)

    delete [] DofMarker;

  MPI_Allreduce(&N_MasterDof, &MaxN_MasterDofs_All, 1, MPI_INT, MPI_MAX, Comm);


  MPI_Allreduce(&N_SlaveDof, &MaxN_SlaveDofs_All, 1, MPI_INT, MPI_MAX, Comm);

//printf("Rank %d MaxN_SlaveDofs_All %d  N_SlaveDof %d \n", rank, MaxN_SlaveDofs_All, N_SlaveDof);

  if(N_Neibs)
   {
    N_SendDofs = new int[N_Neibs];
    N_RecevDofs = new int[N_Neibs];

    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);

    sendbuf = new int*[2];
    recevbuf = new int*[2];
    sendbuf[0] = new int[N_Neibs*MaxN_DeptDofs_All];
    sendbuf[1] = new int[N_Neibs*MaxN_DeptDofs_All];
    recevbuf[0] = new int[N_Neibs*MaxN_DeptDofs_All];
    recevbuf[1] = new int[N_Neibs*MaxN_DeptDofs_All];

    // first built DeptDofLocalIndexOfNeibRanks
    for(i=0;i<N_DeptDofs;i++)
     {
      N = N_DeptDofNeibs[i];
      for(j=0;j<N;j++)
       {
        Neib_ID = DeptDofNeibRanks[i*MaxSubDomainPerDof + j];
        k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
        M = k*MaxN_DeptDofs_All + N_SendDofs[k];
        N_SendDofs[k]++;

        sendbuf[0][M] = DeptDofNeibRanksLocalDOF[i*MaxSubDomainPerDof + j];
       } // for(j=0;j<N;j++)
     } //for(i=0;i<N_DeptDofs;
   } // if(N_Neibs)

   FECommunicateNeib(sendbuf, MaxN_DeptDofs_All, N_SendDofs, recevbuf, N_SendDofs, 1);

  if(N_Neibs)
   {
    memset(N_RecevDofs, 0, N_Neibs*SizeOfInt);

    for(i=0;i<N_DeptDofs;i++)
     {
      N = N_DeptDofNeibs[i];
      for(j=0;j<N;j++)
       {
        m1 = i*MaxSubDomainPerDof + j;
        Neib_ID = DeptDofNeibRanks[m1];
        k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
        M = k*MaxN_DeptDofs_All + N_RecevDofs[k];
        N_RecevDofs[k]++;
        DeptDofLocalIndexOfNeibRanks[m1] = recevbuf[0][M];
       } // for(j=0;j<N;j++)
     } //for(i=0;i<N_DeptDofs;
   } // if(N_Neibs)


  if(N_Neibs)
   {
    memset(N_SendDofs, 0, N_Neibs*SizeOfInt);
    memset(N_RecevDofs, 0, N_Neibs*SizeOfInt);

    for(i=0;i<N_DeptDofs;i++)
     {
      ID = DeptDofMaster[i];

      if(ID==rank) // master of this dof
       {
        N = N_DeptDofNeibs[i];
        for(j=0;j<N;j++)
         {
          Neib_ID = DeptDofNeibRanks[i*MaxSubDomainPerDof + j];
          k = IndexOfNeibRank[Neib_ID]; // LocindexOfNeib_ID
          M = k*MaxN_MasterDofs_All + N_SendDofs[k];
          N_SendDofs[k]++;

          sendbuf[0][M] = DeptDofNeibRanksLocalDOF[i*MaxSubDomainPerDof + j];
          sendbuf[1][M] = GlobalDofOfLocalDof[DeptDofs[i]];
         } // for(j=0;j<N;j++)
       } // if(ID==rank) 
      else // slave of this dof, info have to be obtained from neibs
       {
        k = IndexOfNeibRank[ID]; // LocindexOfNeib_ID
        N_RecevDofs[k]++;
       } // elsei(ID==rank) 
     } //for(i=0;i<N_DeptDofs;
   } // if(N_Neibs)

  // type =1, i.e., higher to lower ranks
  FECommunicateOneWay(sendbuf, MaxN_MasterDofs_All, N_SendDofs, recevbuf, N_RecevDofs, 2, 1);

  for(i=0; i<N_Neibs; i++)
   {
    N = N_RecevDofs[i];
    for(j=0;j<N;j++)
     {
      M = i*MaxN_MasterDofs_All + j;

      if(GlobalDofOfLocalDof[recevbuf[0][M]]!=-1)
       printf("Error: rank %d NeibsRank %d Dof %d RecevGlobalDof %d GlobalDof %d\n",
              rank, NeibsRank[i], recevbuf[0][M], recevbuf[1][M], GlobalDofOfLocalDof[recevbuf[0][M]]);


      GlobalDofOfLocalDof[recevbuf[0][M]] = recevbuf[1][M];
     } // for(j=0;j<N;j++)
    } // for(i=0; i<N_Nei*/


//    if(rank!=0)
   for(i=0;i<N_U;i++)
    if(GlobalDofOfLocalDof[i] == -1)
     {
      printf("Error: in GlobalDOf construction Rank %d i %d \n",rank, i);
      MPI_Finalize();
      exit(0);
    }

//    if(rank==1)
//    for(i=0;i<N_U;i++)
//       printf("%d Rank %d  GlobalDofOfLocalDof   %d \n",i, rank, GlobalDofOfLocalDof[i]);

  // Mapping of global dof from all subdomains to root
//  int *SubDomainParentCellNo, *SubDomainDofGlobalNumbers;
//  int *SubDomainGlobalNumbers, *SubDomainBeginIndex;
//  int *IntArray = new int[2];

   pos = new int[size];
   for(i=0;i<size;i++)
    pos[i] = i*MaxN_LocalDofAllRank;

   GlobalDofOFLocalDofAllRank = new int[size*MaxN_LocalDofAllRank];
   N_LocalDofAllRank = new int[size];
   MPI_Allgather(&N_U, 1, MPI_INT, N_LocalDofAllRank, 1, MPI_INT, Comm);

   MPI_Allgatherv(GlobalDofOfLocalDof, N_U, MPI_INT, GlobalDofOFLocalDofAllRank, N_LocalDofAllRank,
                  pos, MPI_INT, Comm);

   delete [] pos;

//  MPI_Status status;
//   if(rank!=0)
//    {
// //      IntArray[0] = N_Cells;
// //      IntArray[1] = N_U;
//      MPI_Send(&N_U, 1, MPI_INT, 0, 100, Comm);
// 
// //      SubDomainParentCellNo = Coll->GetGlobalIndex();
// //      MPI_Isend(SubDomainParentCellNo, N_Cells, MPI_INT, 0, 200, Comm, &request001);
//      MPI_Send(GlobalDofOfLocalDof, N_U, MPI_INT, 0, 300, Comm);
// //      MPI_Isend(GlobalNumbers, N_Cells*N_LocDof, MPI_INT, 0, 400, Comm, &request003);
// //      MPI_Isend(BeginIndex, N_Cells, MPI_INT, 0, 500, Comm, &request004);
//    }
//   else
//    {
// //     int l, P, N_SubDomainCells, N_SubDomainDOF, *GlobNo, *DOF, *SubDomainDOF;
// //     int *NewGlobalNumbers, *NewDOF;
//     int N_SubDomainDOF, *GlobNo;
// 
// 
// //      NewGlobalNumbers = new int[N_Cells*N_LocDof];
// 
//      N_LocalDofAllRank = new int[size];
//      GlobalDofOFLocalDofAllRank = new int[size*MaxN_LocalDofAllRank];
//      N_LocalDofAllRank[rank] =  N_U;
// 
//      //copy own GlobalDofOfLocalDof values first
//      memcpy(GlobalDofOFLocalDofAllRank, GlobalDofOfLocalDof,  N_U*SizeOfInt);
// 
//      for(j=1;j<size;j++)
//       {
//        MPI_Recv(&N_SubDomainDOF, 1, MPI_INT, j, 100, Comm, &status);
// 
// //        N_SubDomainCells = IntArray[0]; // including halo cells if any
// //        N_SubDomainDOF = IntArray[1];   // including halo cells DOF if any 
//        N_LocalDofAllRank[j] = N_SubDomainDOF;
// 
// //        if(j>1)
// //         {
// //          delete [] SubDomainParentCellNo;
// //          delete [] SubDomainGlobalNumbers;
// //          delete [] SubDomainBeginIndex;
// //         }
// 
//        // /** get the subdomain cells' global cell numbers */
// //        SubDomainParentCellNo = new int[N_SubDomainCells];
// //        MPI_Recv(SubDomainParentCellNo, N_SubDomainCells, MPI_INT, j, 200, Comm, &status);
// 
//        GlobNo =  GlobalDofOFLocalDofAllRank +j*MaxN_LocalDofAllRank;
//        MPI_Recv(GlobNo, N_SubDomainDOF, MPI_INT, j, 300, Comm, &status);
// 
// //        SubDomainGlobalNumbers = new int[N_SubDomainCells*N_LocDof];
// //     MPI_Recv(SubDomainGlobalNumbers, N_SubDomainCells*N_LocDof, MPI_INT, j, 400, Comm, &status);
// // 
// //        SubDomainBeginIndex = new int[N_SubDomainCells];
// //        MPI_Recv(SubDomainBeginIndex, N_SubDomainCells, MPI_INT, j, 500, MPI_COMM_WORLD, &status);
// 
// //        for(k=0;k<N_SubDomainCells;k++)
// //         {
// //          M = SubDomainParentCellNo[k];
// //          DOF = GlobalNumbers + BeginIndex[M];
// // //          NewDOF = NewGlobalNumbers + BeginIndex[M];
// // 
// //          SubDomainDOF = SubDomainGlobalNumbers + SubDomainBeginIndex[k];
// // 
// //          for(l=0; l<N_LocDof; l++)
// //           {
// //            N=DOF[l];
// //            P=GlobNo[SubDomainDOF[l]];
// // 
// //            if(GlobalDofOfLocalDof[N]==-1)
// //             {
// //              GlobalDofOfLocalDof[N] = P;
// //             }
// //            else if(GlobalDofOfLocalDof[N] != P)
// //             {
// //              printf("Root: error in GlobalDOf mapping Rank %d   old %d new %d\n",j,  GlobalDofOfLocalDof[N], P);
// //              MPI_Abort(Comm,  0);
// //             }
// // 
// // //            NewDOF[l]=P;
// //           } // for(l=0; l<N_LocDof; l++)
// //         } //
//       }// for(j=1;j<size;j++)
// 
// //        if(j>1)
// //         {
// //          delete [] SubDomainParentCellNo;
// //          delete [] SubDomainGlobalNumbers;
// //          delete [] SubDomainBeginIndex;
// //         }
// 
// //     // new global numbers assigned by subdomains will be put in the root
// //     memcpy(GlobalNumbers, NewGlobalNumbers, N_Cells*N_LocDof*SizeOfInt);
// //     delete [] NewGlobalNumbers;
// 
// //    for(i=0;i<N_U;i++)
// //     if(GlobalDofOfLocalDof[i] == -1)
// //      {
// //       printf("Error: in GlobalDOf construction Rank %d i %d \n",rank, i);
// //       MPI_Finalize();
// //       exit(0);
// //     }
// 
// // //     // no need in root, since GlobalNumbers and GlobalDofOfLocalDof are same
// // // //     delete [] GlobalDofOfLocalDof;
//    } // else if(rank==0)



  if(N_Neibs)
   {
    delete [] N_SendDofs;
    delete [] N_RecevDofs;


    for(i=0;i<2;i++)
     {
      delete [] sendbuf[i];
      delete [] recevbuf[i];
     }

    delete [] sendbuf;
    delete [] recevbuf;
   }

//    delete [] N_OwnDofAll;
//    delete [] IntArray;

//  if(rank==TDatabase::ParamDB->Par_P5)
//   printf("Rank %d, N_MasterDof_All %d\n", rank, MaxN_MasterDofs_All);

  if(TDatabase::ParamDB->SC_VERBOSE>4)
  if(rank==TDatabase::ParamDB->Par_P0)
   printf("ConstructGlobalDofFromNeib3D done !!!\n");

//      MPI_Finalize();
//      exit(0);

}


/** Mapping of global dof from root to all subdomains */
void TParFECommunicator3D::MapDofFromRoot3D()
{
 int i, j, k, l, M, rank, size, N_Cells;
 int N_U, N_LocDof;
 int *GlobalNumbers, *BeginIndex, *DOF, *SubDomainDOF, P, N;
 int *IntArray = new int[2];
 int *SubDomainParentCellNo, *SubDomainDofGlobalNumbers;
 int *SubDomainGlobalNumbers, *SubDomainBeginIndex;

 TCollection *Coll;
 TBaseCell *cell;

 MPI_Status status;

  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  N_U = FESpace->GetN_DegreesOfFreedom();
  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = FESpace->GetBeginIndex();
  GlobalNumbers = FESpace->GetGlobalNumbers();
  N_LocDof = BeginIndex[1] - BeginIndex[0];  // assume that  all cells in fespace have same FE

   if(rank==0)
    {
     int N_SubDomainCells, N_SubDomainDOF;

     N_LocalDofAllRank = new int[size];
     N_LocalDofAllRank[rank] =  N_U;
     GlobalDofOFLocalDofAllRank = new int[(size-1)*MaxN_LocalDofAllRank];

     for(j=1;j<size;j++)
      {
       MPI_Recv(IntArray, 2, MPI_INT, j, 100, Comm, &status);

       N_SubDomainCells = IntArray[0]; // including halo cells if any
       N_SubDomainDOF = IntArray[1];   // including halo cells DOF if any   
       N_LocalDofAllRank[j] = N_SubDomainDOF;

       if(j>1)
        {
         delete [] SubDomainParentCellNo;
         delete [] SubDomainGlobalNumbers;
         delete [] SubDomainBeginIndex;
        }

       /** get the subdomain cells' global cell numbers */
       SubDomainParentCellNo = new int[N_SubDomainCells];
       MPI_Recv(SubDomainParentCellNo, N_SubDomainCells, MPI_INT, j, 200, Comm, &status);

       SubDomainGlobalNumbers = new int[N_SubDomainCells*N_LocDof];
       MPI_Recv(SubDomainGlobalNumbers, N_SubDomainCells*N_LocDof, MPI_INT, j, 300, Comm, &status);

       SubDomainBeginIndex = new int[N_SubDomainCells];
       MPI_Recv(SubDomainBeginIndex, N_SubDomainCells, MPI_INT, j, 400, MPI_COMM_WORLD, &status);

       /** mapping begin*/
       for(k=0; k<N_SubDomainCells; k++)
        {
         M = SubDomainParentCellNo[k];
         DOF = GlobalNumbers + BeginIndex[M];
         SubDomainDOF = SubDomainGlobalNumbers + SubDomainBeginIndex[k];

         for(l=0; l<N_LocDof; l++)
          {
           N=DOF[l];
           P = SubDomainDOF[l];

           GlobalDofOFLocalDofAllRank[(j-1)*MaxN_LocalDofAllRank +  P] = N;
          } // for(l=0; l<N_LocDof; l++)
        } // for(i=0; i<N_SubDomainCells; i++)

       if(j>1)
        MPI_Wait(&request006, MPI_STATUS_IGNORE);

//         MPI_Isend(GlobalDofOFLocalDofAllRank+(j-1)*MaxN_LocalDofAllRank, N_SubDomainDOF, MPI_INT,
//                   j, 500, Comm, &request006);
        MPI_Send(GlobalDofOFLocalDofAllRank+(j-1)*MaxN_LocalDofAllRank, N_SubDomainDOF, MPI_INT,
                  j, 500, Comm);
      }//for(j=1;j<size

     delete [] SubDomainParentCellNo;
     delete [] SubDomainGlobalNumbers;
     delete [] SubDomainBeginIndex;
    }
   else
    {
     IntArray[0] = N_Cells;
     IntArray[1] = N_U;
     MPI_Send(IntArray, 2, MPI_INT, 0, 100, Comm);

     SubDomainParentCellNo = Coll->GetGlobalIndex();
     MPI_Send(SubDomainParentCellNo, N_Cells, MPI_INT, 0, 200, Comm);

     MPI_Send(GlobalNumbers, N_Cells*N_LocDof, MPI_INT, 0, 300, Comm);
     MPI_Send(BeginIndex, N_Cells, MPI_INT, 0, 400, Comm);

     GlobalDofOfLocalDof = new int[N_U];
//      MPI_Irecv(GlobalDofOfLocalDof, N_U, MPI_INT, 0, 500, Comm, &request006);
     MPI_Recv(GlobalDofOfLocalDof, N_U, MPI_INT, 0, 500, Comm, MPI_STATUS_IGNORE);
// // //       WaitForMapDofFromRoot3D();
//      MPI_Wait(&request006, MPI_STATUS_IGNORE);
//       if(rank==1)
//        {
//         for(i=0; i<N_U; i++)
//          printf("In rank LocalDof %d:  GlobalDofOfLocalDof  : %d \n",i,  GlobalDofOfLocalDof[i]);
//        }
    } // else if(rank==0)

  delete [] IntArray;

//  printf("Rankk %d MapDofFromRoot3D  %d  \n",rank, N_Neibs );
//  MPI_Finalize();
// exit(0);
}


/** Communication between processor based on defined schedule */
/** arrays (send & recev) with same (senddisp) size */
int TParFECommunicator3D::FECommunicateNeib(int **sendbuf, int senddisp, int *sendlen,
                                              int **recvbuf, int *recevlen, int N_array)
{
 int i, j, k, l, m, M, N, rank, size, sendID, recvID;
 int totaldisp, *buffer, *tmp, pos, *Sendbuf, *Recvbuf, max;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  max = N_array*senddisp;
  MPI_Allreduce(&max, &totaldisp, 1, MPI_INT, MPI_MAX, Comm);

//  if(TDatabase::ParamDB->Par_P6==5)
//    printf(" rank %d FECommunicateNeib    totaldisp %d\n", rank, totaldisp);

//   if(max==0)  
//    totaldisp=0;

  if(totaldisp!=0)
  {
   buffer = new int[totaldisp];
   tmp = new int[totaldisp];

   for(i=0;i<Max_CommunicationSteps;i++)
   {
    N = N_ActiveCommuPross[i];
    for(j=0;j<N;j++)
     {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank
        pos = 0;

        if(sendlen[M]>0)
         {
         // pack all arrays and send 

         for(l=0;l<N_array; l++)
          {
           Sendbuf = sendbuf[l];
           MPI_Pack(Sendbuf+M*senddisp, sendlen[M], MPI_INT, buffer, totaldisp*SizeOfInt,  &pos, Comm);
          }

          MPI_Send(buffer, pos, MPI_PACKED,  recvID, 100, Comm);
         } // // if(sendlen[M]>0)
       }
      else if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; 

        if(recevlen[M]>0)
         {
          MPI_Recv(buffer, totaldisp*SizeOfInt, MPI_PACKED,  sendID, 100, Comm, &status);

          pos = 0;
          for(l=0;l<N_array; l++)
           {
            MPI_Unpack(buffer, totaldisp*SizeOfInt, &pos, tmp, recevlen[M], MPI_INT, Comm);
            Recvbuf = recvbuf[l];
            for(m=0;m<recevlen[M]; m++)
              Recvbuf[M*senddisp + m] = tmp[m];
           } //  for(l=0;l<N_ar
         } // if(recevlen[M]>0)
       } //  else if(rank==recvID

      if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank
        pos = 0;

        if(sendlen[M]>0)
         {
          // pack all arrays and send 
          for(l=0;l<N_array; l++)
           {
            Sendbuf = sendbuf[l];
            MPI_Pack(Sendbuf+M*senddisp, sendlen[M], MPI_INT, buffer, totaldisp*SizeOfInt,  &pos, Comm);
           }

          MPI_Send(buffer, pos, MPI_PACKED,  sendID, 200, Comm);
         } // if(sendlen[M]>0)
       }
      else if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(recevlen[M]>0)
         {
          MPI_Recv(buffer, totaldisp*SizeOfInt, MPI_PACKED,  recvID, 200, Comm, &status);

          pos = 0;
          for(l=0;l<N_array; l++)
           {
            MPI_Unpack(buffer, totaldisp*SizeOfInt, &pos, tmp, recevlen[M], MPI_INT, Comm);
            Recvbuf = recvbuf[l];
            for(m=0;m<recevlen[M]; m++)
              Recvbuf[M*senddisp + m] = tmp[m];
           } //  for(l=0;l<N_ar
          } // if(recevlen[M]>0
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max

   delete [] buffer;
   delete [] tmp;
  }

//   printf(" rank %d FECommunicateNeib senddisp %d totaldisp %d\n", rank, senddisp, totaldisp);
//   MPI_Finalize();
//   exit(0);  

  return 0;
} // FECommunicateNeib



/** Communication between processor based on defined schedule */
/** send and recev arrays with same size */
int TParFECommunicator3D::FECommunicateNeib(double **sendbuf, int senddisp, int *sendlen,
                                            double **recvbuf, int *recevlen, int N_array)
{
 int i, j, k, l, m, M, N, rank, size, sendID, recvID;
 int totaldisp,  pos,  max;

 double *buffer, *tmp, *Sendbuf, *Recvbuf;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  max = N_array*senddisp;
  MPI_Allreduce(&max, &totaldisp, 1, MPI_INT, MPI_MAX, Comm);

  if(totaldisp!=0)
  {
   buffer = new double[totaldisp];
   tmp = new double[totaldisp];

   for(i=0;i<Max_CommunicationSteps;i++)
   {
    N = N_ActiveCommuPross[i];
    for(j=0;j<N;j++)
     {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank
        pos = 0;

        if(sendlen[M]>0)
         {
         // pack all arrays and send 

         for(l=0;l<N_array; l++)
          {
           Sendbuf = sendbuf[l];
           MPI_Pack(Sendbuf+M*senddisp, sendlen[M], MPI_DOUBLE, buffer, totaldisp*SizeOfDouble,  &pos, Comm);
          }

          MPI_Send(buffer, pos, MPI_PACKED,  recvID, 100, Comm);
         } // // if(sendlen[M]>0)
       }
      else if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; 

        if(recevlen[M]>0)
         {
          MPI_Recv(buffer, totaldisp*SizeOfDouble, MPI_PACKED,  sendID, 100, Comm, &status);

          pos = 0;
          for(l=0;l<N_array; l++)
           {
            MPI_Unpack(buffer, totaldisp*SizeOfDouble, &pos, tmp, recevlen[M], MPI_DOUBLE, Comm);
            Recvbuf = recvbuf[l];
            for(m=0;m<recevlen[M]; m++)
              Recvbuf[M*senddisp + m] = tmp[m];
           } //  for(l=0;l<N_ar
         } // if(recevlen[M]>0)
       } //  else if(rank==recvID

      if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank
        pos = 0;

        if(sendlen[M]>0)
         {
          // pack all arrays and send 
          for(l=0;l<N_array; l++)
           {
            Sendbuf = sendbuf[l];
            MPI_Pack(Sendbuf+M*senddisp, sendlen[M], MPI_DOUBLE, buffer, totaldisp*SizeOfDouble,  &pos, Comm);
           }

          MPI_Send(buffer, pos, MPI_PACKED,  sendID, 200, Comm);
         } // if(sendlen[M]>0)
       }
      else if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(recevlen[M]>0)
         {
          MPI_Recv(buffer, totaldisp*SizeOfDouble, MPI_PACKED,  recvID, 200, Comm, &status);

          pos = 0;
          for(l=0;l<N_array; l++)
           {
            MPI_Unpack(buffer, totaldisp*SizeOfDouble, &pos, tmp, recevlen[M], MPI_DOUBLE, Comm);
            Recvbuf = recvbuf[l];
            for(m=0;m<recevlen[M]; m++)
              Recvbuf[M*senddisp + m] = tmp[m];
           } //  for(l=0;l<N_ar
          } // if(recevlen[M]>0
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max

   delete [] buffer;
   delete [] tmp;
  }

//   printf(" rank %d FECommunicateNeib senddisp %d totaldisp %d\n", rank, senddisp, totaldisp);
//   MPI_Finalize();
//   exit(0);  

  return 0;
} // FECommunicateNeib


int TParFECommunicator3D::FECommunicateNeib(double *sendbuf, int senddisp, int *sendlen,
                                            double *recvbuf, int N_dim)
{
 int i, j, M, N, rank, size, sendID, recvID;
 double *Sendbuf, *Recvbuf;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  for(i=0;i<Max_CommunicationSteps;i++)
   {
    N = N_ActiveCommuPross[i];
    for(j=0;j<N;j++)
     {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Sendbuf = sendbuf+N_dim*M*senddisp;
//           MPI_Send(Sendbuf, N_dim*sendlen[M], MPI_DOUBLE,  recvID, 100, Comm);
          MPI_Send(Sendbuf, N_dim*senddisp, MPI_DOUBLE,  recvID, 100, Comm);
         } // // if(sendlen[M]>0)
       }
      else if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Recvbuf = recvbuf+N_dim*M*senddisp;
//           MPI_Recv(Recvbuf, N_dim*sendlen[M], MPI_DOUBLE,  sendID, 100, Comm, &status);
          MPI_Recv(Recvbuf, N_dim*senddisp, MPI_DOUBLE,  sendID, 100, Comm, &status);
         } // if(sendlen[M]>0)
       } //  else if(rank==recvID

      if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Sendbuf = sendbuf+N_dim*M*senddisp;
//           MPI_Send(Sendbuf, N_dim*sendlen[M], MPI_DOUBLE,  sendID, 100, Comm);
          MPI_Send(Sendbuf, N_dim*senddisp, MPI_DOUBLE,  sendID, 100, Comm);
         } // if(sendlen[M]>0)
       }
      else if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Recvbuf = recvbuf+N_dim*M*senddisp;
//           MPI_Recv(Recvbuf, N_dim*sendlen[M], MPI_DOUBLE,  recvID, 100, Comm, &status);
          MPI_Recv(Recvbuf, N_dim*senddisp, MPI_DOUBLE,  recvID, 100, Comm, &status);
          } // if(sendlen[M]>0
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max

  return 0;
} // FECommunicateNeib

int TParFECommunicator3D::FECommunicateNeib(double *sendbuf, int senddisp, int *sendlen,
                                            double *recvbuf, int *recevlen, int N_dim)
{
 int i, j, M, N, rank, size, sendID, recvID;
 double *Sendbuf, *Recvbuf;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  for(i=0;i<Max_CommunicationSteps;i++)
   {
    N = N_ActiveCommuPross[i];
    for(j=0;j<N;j++)
     {
      sendID = SendID_Index[i*size +j];
      recvID = ReceiveID_Index[i*size +j];

      if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Sendbuf = sendbuf+N_dim*M*senddisp;
          MPI_Send(Sendbuf, N_dim*sendlen[M], MPI_DOUBLE,  recvID, 100, Comm);
//           MPI_Send(Sendbuf, N_dim*senddisp, MPI_DOUBLE,  recvID, 100, Comm);
         } // // if(sendlen[M]>0)
       }
      else if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank

        if(recevlen[M]>0)
         {
          Recvbuf = recvbuf+N_dim*M*senddisp;
          MPI_Recv(Recvbuf, N_dim*recevlen[M], MPI_DOUBLE,  sendID, 100, Comm, &status);
//           MPI_Recv(Recvbuf, N_dim*senddisp, MPI_DOUBLE,  sendID, 100, Comm, &status);
         } // if(recevlen[M]>0)
       } //  else if(rank==recvID

      if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank

        if(sendlen[M]>0)
         {
          Sendbuf = sendbuf+N_dim*M*senddisp;
          MPI_Send(Sendbuf, N_dim*sendlen[M], MPI_DOUBLE,  sendID, 100, Comm);
//           MPI_Send(Sendbuf, N_dim*senddisp, MPI_DOUBLE,  sendID, 100, Comm);
         } // if(sendlen[M]>0)
       }
      else if(rank==sendID)
       {
        M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank

        if(recevlen[M]>0)
         {
          Recvbuf = recvbuf+N_dim*M*senddisp;
          MPI_Recv(Recvbuf, N_dim*recevlen[M], MPI_DOUBLE,  recvID, 100, Comm, &status);
//           MPI_Recv(Recvbuf, N_dim*senddisp, MPI_DOUBLE,  recvID, 100, Comm, &status);
          } // if(recevlen[M]>0
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max

  return 0;
} // FECommunicateNeib

int TParFECommunicator3D::FECommunicateOneWay(int **sendbuf, int disp, int *sendlen,
                                            int **recvbuf, int *recevlen, int N_array, int type)
{
 int i, j, k, l, m, M, N, rank, size, sendID, recvID;
 int totaldisp, *buffer, *tmp, pos, *Sendbuf, *Recvbuf;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

  totaldisp = N_array*disp;

// printf(" rank %d FECommunicateNeib    totaldisp %d\n", rank, totaldisp);

 if(totaldisp!=0)
  {
   buffer = new int[totaldisp];
   tmp = new int[totaldisp];

   for(i=0;i<Max_CommunicationSteps;i++)
    {
     N = N_ActiveCommuPross[i];
     for(j=0;j<N;j++)
      {
       if(type==0)
        {
         sendID = SendID_Index[i*size +j];
         recvID = ReceiveID_Index[i*size +j];
        }
       else
        {
         recvID = SendID_Index[i*size +j];
         sendID = ReceiveID_Index[i*size +j];
        }

       if(rank==sendID)
        {
         M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank
         pos = 0;

        if(sendlen[M]>0)
         {
         // pack all arrays and send 
         for(l=0;l<N_array; l++)
          {
           Sendbuf = sendbuf[l];
           MPI_Pack(Sendbuf+M*disp, sendlen[M], MPI_INT, buffer, totaldisp*SizeOfInt,  &pos, Comm);
          }

          MPI_Send(buffer, pos, MPI_PACKED,  recvID, 100, Comm);
         } // // if(sendlen[M]>0)
       }
      else if(rank==recvID)
       {
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank

        if(recevlen[M]>0)
         {
          MPI_Recv(buffer, totaldisp*SizeOfInt, MPI_PACKED,  sendID, 100, Comm, &status);

          pos = 0;
          for(l=0;l<N_array; l++)
           {
            MPI_Unpack(buffer, totaldisp*SizeOfInt, &pos, tmp, recevlen[M], MPI_INT, Comm);
            Recvbuf = recvbuf[l];
            for(m=0;m<recevlen[M]; m++)
              Recvbuf[M*disp + m] = tmp[m];
           } //  for(l=0;l<N_ar
         } // if(sendlen[M]>0)
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max

   delete [] buffer;
   delete [] tmp;
  } // if(tot

// //   printf(" rank %d FECommunicateNeib senddisp %d totaldisp %d\n", rank, senddisp, totaldisp);
// //   MPI_Finalize();
// //   exit(0);  
// 
//   return 0;
} // FECommunicateNeib

/** send array of disp size between neibs */
int TParFECommunicator3D::FECommunicateOneWay(double *sendbuf, int disp, double *recvbuf,  int type)
{
 int i, j, k, l, m, M, N, rank, size, sendID, recvID;
 double *Sendbuf, *Recvbuf;

  MPI_Status status;
  MPI_Comm_rank(Comm, &rank);
  MPI_Comm_size(Comm, &size);

 if(disp!=0)
  {
   for(i=0;i<Max_CommunicationSteps;i++)
    {
     N = N_ActiveCommuPross[i];
     for(j=0;j<N;j++)
      {
       if(type==0)
        {
         sendID = SendID_Index[i*size +j];
         recvID = ReceiveID_Index[i*size +j];
        }
       else
        {
         recvID = SendID_Index[i*size +j];
         sendID = ReceiveID_Index[i*size +j];
        }
       if(rank==sendID)
        {
         M = IndexOfNeibRank[recvID]; // LocindexOfRecvRank         
         Sendbuf = sendbuf+M*disp;
//printf("rank %d i %d j %d  recvID %d\n", rank, i, j,  recvID);
         MPI_Send(Sendbuf, disp, MPI_DOUBLE,  recvID , 100, Comm);
       }
      else if(rank==recvID)
       { 
        M = IndexOfNeibRank[sendID]; // LocindexOfRecvRank
        Recvbuf = recvbuf+M*disp;
//printf("rank %d  i %d j %d   sendID %d \n", rank, i, j, sendID);
        MPI_Recv(Recvbuf, disp, MPI_DOUBLE,  sendID, 100, Comm, &status);
       } //  else if(rank==recvID
     } //  for(j=0;j<N;j++)
   } //  for(i=0; i<Max
  } // if(tot

// //   printf(" rank %d FECommunicateNeib senddisp %d totaldisp %d\n", rank, senddisp, totaldisp);
// //   MPI_Finalize();
// //   exit(0);  
// 
//   return 0;
} // FECommunicateNeib


// root should not be called this function
int TParFECommunicator3D::WaitForGlobalDofFromNeib3D()
{
//    MPI_Wait(&request001, MPI_STATUS_IGNORE);
//    MPI_Wait(&request002, MPI_STATUS_IGNORE);
//    MPI_Wait(&request003, MPI_STATUS_IGNORE);
//    MPI_Wait(&request004, MPI_STATUS_IGNORE);

  return 0;
}


int TParFECommunicator3D::WaitForMapDofFromRoot3D()
{
   MPI_Wait(&request006, MPI_STATUS_IGNORE);

  return 0;
}


#endif

