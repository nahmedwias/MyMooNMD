// =======================================================================
// @(#)MainUtilities.C        1.11 07/03/00
//
// Purpose:     main utilities
//
// Authors:     Gunar Matthies  18.10.99
//
// =======================================================================

#include <MainUtilities.h>
#include <Constants.h>
#include <Database.h>
#include <Enumerations.h>
#include <Convolution.h>
#include <LinAlg.h>
#include <ConvDiff.h>

#include <AllRefTrans.h>

#ifdef __2D__
#include <DiscreteForm2D.h>
#include <FEDatabase2D.h>
#include <FEVectFunct2D.h>
#include <TNSE2D_Routines.h>
#include <NSE2DSUPG.h>
#endif

#ifdef __3D__
#include <FEDatabase3D.h>
#include <FEVectFunct2D.h>
#include <FEVectFunct3D.h>
#include <CommonRoutineTNSE3D.h>
#include <BoundFace.h>
#endif

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
// #include <malloc.h>
#include <MooNMD_Io.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <stdio.h>

#ifdef __MAC64__
#include <malloc/malloc.h> // defines _MALLOC_MALLOC_H_
#else
#include <malloc.h>
#endif

double GetTime()
{
  struct rusage usage;
  double ret;
    
  if(getrusage(RUSAGE_SELF, &usage) == -1) 
    cerr << "Error in GetTime!" << endl;

  ret = ((double) usage.ru_utime.tv_usec)/1000000;

  ret += usage.ru_utime.tv_sec;

  return ret;
}

int GetMemory()
{
#ifdef _MALLOC_MALLOC_H_
 struct mstats info;
 info = mstats();
 return info.bytes_free+info.chunks_free;
#else  
  struct mallinfo MALLINFO;
  MALLINFO = mallinfo();
  return MALLINFO.usmblks+MALLINFO.uordblks;
#endif
}

//Print all memory information available through mallinfo.
void display_mallinfo(const std::string& program_part)
{
#ifdef _MALLOC_MALLOC_H_
  Output::print("Memory usage info called in program part: ", program_part);
  Output::print(GetMemory());  
#else
   
  struct mallinfo mi;
  
  mi = mallinfo();
  Output::print("Memory usage info called in program part: ", program_part);
  Output::print("Total non-mmapped bytes (arena):       ",mi.arena);
  Output::print("# of free chunks (ordblks):            ",mi.ordblks);
  Output::print("# of free fastbin blocks (smblks):     ",mi.smblks);
  Output::print("# of mapped regions (hblks):           ",mi.hblks);
  Output::print("Bytes in mapped regions (hblkhd):      ",mi.hblkhd);
  Output::print("Max. total allocated space (usmblks):  ",mi.usmblks);
  Output::print("Free bytes held in fastbins (fsmblks): ",mi.fsmblks);
  Output::print("Total allocated space (uordblks):      ",mi.uordblks);
  Output::print("Total free space (fordblks):           ",mi.fordblks);
  Output::print("Topmost releasable block (keepcost):   ",mi.keepcost);
#endif
}

#ifdef __2D__
int VertexNumber(int IEH, int NVE)
{
  int ret;

  if(NVE==4)
  {
    switch(IEH)
    {
      case 2: ret=3; break;
      case 3: ret=2; break;
      default: ret=IEH;
    }
  }
  else
  {
    ret = IEH;
  }

  return ret;
}

void StreamFunction(const TFESpace2D *velo, double *u1, double *u2,
                    const TFESpace2D *stream, double *psi)
{
  int IEH,IEL,k,l, N_Cells, NVE, NVE2, NVT, VertNum, m, l1;
  TCollection *Coll;
  TBaseCell *cell, *neigh, **CellList;
  int ListPointer, ListInput;
  int *KVIND;
  FE2D ele;
  TFE2D *element;
  TFEDesc2D *desc;
  int *VeloGlobalNumbers, *VeloBeginIndex;
  int *StreamGlobalNumbers, *StreamBeginIndex;
  int IH, JVE, IHV, IND, INDH, IVTH, IVT, IELH;
  double px1, px2, py1, py2, dn1, dn2;
  BaseFunct2D basefunct;
  double v1, v2;
  int *JointDOFs;

  VeloGlobalNumbers = velo->GetGlobalNumbers();
  VeloBeginIndex = velo->GetBeginIndex();

  StreamGlobalNumbers = stream->GetGlobalNumbers();
  StreamBeginIndex = stream->GetBeginIndex();

  Coll = velo->GetCollection();
  N_Cells = Coll->GetN_Cells();
  NVT = stream->GetN_DegreesOfFreedom();
  KVIND = new int [NVT];
  memset(KVIND, 0, SizeOfInt*NVT);
  memset(psi, 0, SizeOfDouble*NVT);
  CellList = new TBaseCell* [N_Cells];
  memset(CellList, 0, sizeof(TBaseCell*)*N_Cells);

  for(IEH=0;IEH<N_Cells;IEH++)
    Coll->GetCell(IEH)->SetClipBoard(IEH);

  l = StreamGlobalNumbers[0];
  psi[l] = 0.0;
  KVIND[l] = 1;

  cell = Coll->GetCell(0);
  CellList[0] = cell;
  l = -(10+0);
  cell->SetClipBoard(l);
  ListPointer = 0; 
  ListInput = 1;

  while(ListPointer<ListInput)
  {
    cell = CellList[ListPointer];

    l = cell->GetClipBoard();
    IEL = IEH = -(10+l);
    
    // cout << "cell: " << IEL << endl;
    IH = 0;
    ele = velo->GetFE2D(IEL, cell);
    NVE2 = cell->GetN_Edges();

    for(k=0;k<NVE2;k++)
    {
      VertNum = StreamBeginIndex[IEL]+VertexNumber(k, NVE2);
      JVE = StreamGlobalNumbers[VertNum];
      IHV = KVIND[JVE];
      IH = IH+IHV;
      if(IHV >= 1) IND = k;
    } // endfor k

    if((IH >= NVE2) || (IH == 0))
    {
      // already all value computed in this cell
      cell->SetClipBoard(-1);

      for(k=0;k<NVE2;k++)
      {
        neigh = cell->GetJoint(k)->GetNeighbour(cell);
   
        if(neigh)
        {
          if(neigh->GetClipBoard()>=0)
          { 
            // neigh was not handled
            IELH = neigh->GetClipBoard();
            IH = 0;
            NVE = neigh->GetN_Edges();
  
            for(l=0;l<NVE;l++)
            {
              VertNum = StreamBeginIndex[IELH]+VertexNumber(l, NVE);
              JVE = StreamGlobalNumbers[VertNum];
              IHV = KVIND[JVE];
              IH = IH+IHV;
              if(IHV>=1) INDH = l;
            } // endfor l
    
            if((IH<NVE) && (IH>0))
            {
              // add cell to list since not all values were calculated
              CellList[ListInput] = neigh;
              l = -(10+neigh->GetClipBoard());
              if(l<0)
                neigh->SetClipBoard(l);
              ListInput++;
            }
          } // >=0
        } // if neigh
      } // endfor k;
    }
    else
    {
      // at least one value must be calculated
      // cout << "cell number: " << IEL << endl;
      for(k=0;k<NVE2-1;k++)
      {
        INDH = (IND+1) % NVE2; 
        VertNum = StreamBeginIndex[IEL]+VertexNumber(INDH, NVE2);
        IVTH = StreamGlobalNumbers[VertNum];

        element = TFEDatabase2D::GetFE2D(ele);
        basefunct = element->GetBaseFunct2D_ID();
        desc = element->GetFEDesc2D();
        // int N_JointDOFs = desc->GetN_JointDOF();
        if(KVIND[IVTH]<1)
        {
          KVIND[IVTH]=1;

          l = VertexNumber(IND, NVE2);
          IVT = StreamGlobalNumbers[StreamBeginIndex[IEL]+l];
  
          // do calculation
          cell->GetVertex(IND)->GetCoords(px1, py1);
          cell->GetVertex(INDH)->GetCoords(px2, py2);
          dn1 = py2-py1;
          dn2 = px1-px2;
  
          l1 = 0;
          while(StreamGlobalNumbers[StreamBeginIndex[IEL]+
                  VertexNumber(l1,NVE2)]!=IVT) l1++;
          // cout << "l1= " << l1 << " IVT: " << IVT;
          // cout << " IVTH: " << IVTH << endl;

          switch(basefunct)
          {
            case BF_N_T_P1_2D:
            case BF_N_Q_Q1_2D:
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+IND];
                psi[IVTH] = psi[IVT] + u1[m]*dn1 + u2[m]*dn2;
              break;
            case BF_C_T_P0_2D:
            case BF_C_Q_Q0_2D:
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]];
                psi[IVTH] = psi[IVT] + u1[m]*dn1 + u2[m]*dn2;
              break;
            case BF_C_T_P1_2D:
            case BF_C_Q_Q1_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[0]];
                v1 = u1[m];
                v2 = u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[1]];
                v1 += u1[m];
                v2 += u2[m];
                v1 /= 2;
                v2 /= 2;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            case BF_C_T_P2_2D:
            case BF_C_Q_Q2_2D:
            case BF_C_T_B2_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[0]];
                v1 = u1[m];
                v2 = u2[m];
                // cout << m << "  " << u1[m] << endl;
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[1]];
                v1 += 4*u1[m];
                v2 += 4*u2[m];
                // cout << m << "  " << u1[m] << endl;
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[2]];
                v1 += u1[m];
                v2 += u2[m];
                // cout << m << "  " << u1[m] << endl;
                v1 /= 6;
                v2 /= 6;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
                // cout << "IVT: " << IVT << endl;
                // cout << "IVTH: " << IVTH << endl;
                // cout << "v1: " << v1 << "     " << "v2: " << v2 << endl;
                // cout << "delta psi: " << v1*dn1 + v2*dn2 << endl;
                // cout << endl;
              break;
            case BF_C_T_P3_2D:
            case BF_C_Q_Q3_2D:
            case BF_C_T_B3_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[0]];
                v1 = u1[m];
                v2 = u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[1]];
                v1 += 3*u1[m];
                v2 += 3*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[2]];
                v1 += 3*u1[m];
                v2 += 3*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[3]];
                v1 += u1[m];
                v2 += u2[m];
                v1 /= 8;
                v2 /= 8;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            case BF_C_T_P4_2D:
            case BF_C_Q_Q4_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[0]];
                v1 = 7*u1[m];
                v2 = 7*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[1]];
                v1 += 32*u1[m];
                v2 += 32*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[2]];
                v1 += 12*u1[m];
                v2 += 12*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[3]];
                v1 += 32*u1[m];
                v2 += 32*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[4]];
                v1 += 7*u1[m];
                v2 += 7*u2[m];
                v1 /= 45;
                v2 /= 45;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            case BF_C_T_P5_2D:
            case BF_C_Q_Q5_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[0]];
                v1 = 19*u1[m];
                v2 = 19*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[1]];
                v1 += 75*u1[m];
                v2 += 75*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[2]];
                v1 += 50*u1[m];
                v2 += 50*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[3]];
                v1 += 50*u1[m];
                v2 += 50*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[4]];
                v1 += 75*u1[m];
                v2 += 75*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[5]];
                v1 += 19*u1[m];
                v2 += 19*u2[m];
                v1 /= 144;
                v2 /= 144;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            case BF_C_T_P6_2D:
            case BF_C_Q_Q6_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[0]];
                v1 = 41*u1[m];
                v2 = 41*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[1]];
                v1 += 216*u1[m];
                v2 += 216*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[2]];
                v1 += 27*u1[m];
                v2 += 27*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[3]];
                v1 += 272*u1[m];
                v2 += 272*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[4]];
                v1 += 27*u1[m];
                v2 += 27*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[5]];
                v1 += 216*u1[m];
                v2 += 216*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[6]];
                v1 += 41*u1[m];
                v2 += 41*u2[m];
                v1 /= 420;
                v2 /= 420;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            case BF_C_T_P7_2D:
            case BF_C_Q_Q7_2D:
                JointDOFs = desc->GetJointDOF(l1);
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[0]];
                v1 = 751*u1[m];
                v2 = 751*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[1]];
                v1 += 3577*u1[m];
                v2 += 3577*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[2]];
                v1 += 1323*u1[m];
                v2 += 1323*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[3]];
                v1 += 2989*u1[m];
                v2 += 2989*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[4]];
                v1 += 2989*u1[m];
                v2 += 2989*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[5]];
                v1 += 1323*u1[m];
                v2 += 1323*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[6]];
                v1 += 3577*u1[m];
                v2 += 3577*u2[m];
                m=VeloGlobalNumbers[VeloBeginIndex[IEL]+JointDOFs[7]];
                v1 += 751*u1[m];
                v2 += 751*u2[m];
                v1 /= 8640;
                v2 /= 8640;
  
                psi[IVTH] = psi[IVT] + v1*dn1 + v2*dn2;
              break;
            default:
                cerr << "no stream function calculation possible for";
                cerr << " the base function set: " << basefunct << endl;
          }
        } // endif

        IND = INDH;
      } // endfor k

      // cout << "NVE2: " << NVE2 << endl;
      for(k=0;k<NVE2;k++)
      {
        neigh = cell->GetJoint(k)->GetNeighbour(cell);
   
        if(neigh)
        {
          if(neigh->GetClipBoard()>=0)
          { 
            // neigh was not handled
            IELH = neigh->GetClipBoard();
            IH = 0;
            NVE = neigh->GetN_Edges();
  
            for(l=0;l<NVE;l++)
            {
              VertNum = StreamBeginIndex[IELH]+VertexNumber(l, NVE);
              JVE = StreamGlobalNumbers[VertNum];
              IHV = KVIND[JVE];
              IH = IH+IHV;
              if(IHV>=1) INDH = l;
            } // endfor l
    
            if((IH<NVE) && (IH>0))
            {
              // add cell to list since not all values were calculated
              CellList[ListInput] = neigh;
              l = -(10+neigh->GetClipBoard());
              if(l<0)
                neigh->SetClipBoard(l);
              ListInput++;
            }
          } // !=-1
        } // if neigh
      } // endfor k;
      cell->SetClipBoard(-2);
    }
    ListPointer++;
  } // endwhile

/*
  for(IEH=0;IEH<NVT;IEH++)
    if(KVIND[IEH]==0)
      cerr << IEH << "   " << KVIND[IEH] << endl;
 */

  delete KVIND;

  delete CellList;
}

void ComputeVorticityDivergence(const TFESpace2D *, TFEFunction2D *u1,
                                TFEFunction2D *u2,
                                const TFESpace2D *vorticity_space,
                                double *vort,  double *div)
{
  int i, j, N_Cells, index, *DOF, N_loc_dofVort, N_Vort;
  int *GlobalNumbersVort, *BeginIndexVort, *N_Found;
  TCollection *Coll;
  TBaseCell *cell;
  FE2D CurrentElementVort;
  TFE2D *Element;
  TNodalFunctional2D *nf;
  RefTrans2D RefTrans;
  const double *xi_ref, *eta_ref;
  double AbsDetjkVort[MaxN_QuadPoints_2D];
  double X_orig[MaxN_PointsForNodal2D], Y_orig[MaxN_PointsForNodal2D];
  double val[3];

  int N_LocalDOF;
  double PointValuesDiv[MaxN_PointsForNodal2D];
  double FunctionalValuesDiv[MaxN_PointsForNodal2D];
  double PointValuesVort[MaxN_PointsForNodal2D];
  double FunctionalValuesVort[MaxN_PointsForNodal2D];

  // get information of the numbering of the degrees of freedom
  // of the vorticity
  GlobalNumbersVort = vorticity_space->GetGlobalNumbers();
  BeginIndexVort = vorticity_space->GetBeginIndex();
  N_Vort = vorticity_space->GetN_DegreesOfFreedom();
  memset(vort,0,N_Vort*SizeOfDouble);
  memset(div,0,N_Vort*SizeOfDouble);
  //doubel *div1 =  new double[N_Vort];
  //memset(div1,0,N_Vort*SizeOfDouble);
  N_Found = new int[N_Vort];
  memset(N_Found,0,N_Vort*SizeOfInt);

  // get pointer to set of mesh cells which define the fe space
  Coll = vorticity_space->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  // loop over all mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElementVort = vorticity_space->GetFE2D(i, cell);
    // get fe from its id
    Element = TFEDatabase2D::GetFE2D(CurrentElementVort);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    TFEDatabase2D::SetCellForRefTrans(cell,RefTrans);
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional2D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofVort, xi_ref, eta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase2D::GetOrigFromRef(RefTrans,N_loc_dofVort, xi_ref, 
                                  eta_ref, X_orig, Y_orig, AbsDetjkVort);
   
    DOF = GlobalNumbersVort + BeginIndexVort[i];
    
    memset(PointValuesDiv, 0, MaxN_PointsForNodal2D*SizeOfDouble);
    memset(PointValuesVort, 0, MaxN_PointsForNodal2D*SizeOfDouble);
    
    for (j=0;j<N_loc_dofVort;j++)
    {
      // values computed with FindGradient have to averaged on a periodic boundary !!! 

      u1->FindGradientLocal(cell,i,X_orig[j],Y_orig[j],val);
      PointValuesVort[j] -= val[2];
      PointValuesDiv[j] += val[1];
      // vort[index] -= val[2];
      // div[index] += val[1];

      u2->FindGradientLocal(cell,i,X_orig[j],Y_orig[j],val); 
      PointValuesVort[j] += val[1];
      PointValuesDiv[j] += val[2];

      //OutPut(values[2] << " " << vort[index] << endl)
      // OutPut(vort[index] << endl);
    }    

    nf->GetAllFunctionals(Coll, cell, PointValuesVort,
                          FunctionalValuesVort);
    nf->GetAllFunctionals(Coll, cell, PointValuesDiv,
                          FunctionalValuesDiv);

    N_LocalDOF = Element->GetN_DOF();
    for(j=0;j<N_LocalDOF;j++)
    {
      index = DOF[j];
      vort[index] += FunctionalValuesVort[j];
      div[index] += FunctionalValuesDiv[j];
      N_Found[index]++;
    }
  } 
 
  for (i=0;i<N_Vort;i++)
  {
    vort[i]/=(double)N_Found[i];
    div[i]/=(double)N_Found[i];
  }
  delete [] N_Found;

}

// determine L2 and H1 error, 2D
void L2H1Errors(int N_Points, std::array<double*, 2>, double *AbsDetjk, 
                const double *Weights, double,
                double **Der, double **Exact,
                double **, double *LocError)
{
  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(int i=0;i<N_Points;i++)
  {
    double *deriv = Der[i];
    double *exactval = Exact[i];
    double w = Weights[i]*AbsDetjk[i];

    double t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;

    t = deriv[1]-exactval[1];
    LocError[1] += w*t*t;
      
    t = deriv[2]-exactval[2];
    LocError[1] += w*t*t;
  } // endfor i
//   cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// determine L2-error, divergence error and H1 error, 2D
void L2DivH1Errors(int N_Points, std::array<double*, 2>, double *AbsDetjk, 
                  const double *Weights, double,
                  double **Der, double **Exact,
                  double **, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;
  LocError[2] = 0.0;
 // cout << endl;
  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];
    
    t = deriv[0]-exactval[0]; // x-component
    LocError[0] += w*t*t;
    t = deriv[3]-exactval[4]; // y-component
    LocError[0] += w*t*t;

    // L2-error of divergence
    t  = deriv[1]-exactval[1]; // x-derivative of x-component
    t += deriv[5]-exactval[6]; // y-derivative of y-component
    LocError[1] += w*t*t;
    
    // H1 semi norm
    t = deriv[1]-exactval[1]; // x-derivative of x-component
    LocError[2] += w*t*t;
    t = deriv[2]-exactval[2]; // y-derivative of x-component
    LocError[2] += w*t*t;
    t = deriv[4]-exactval[5]; // x-derivative of y-component
    LocError[2] += w*t*t;
    t = deriv[5]-exactval[6]; // y-derivative of y-component
    LocError[2] += w*t*t;
  } // endfor i
  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}



void SDFEMErrors(int N_Points, std::array<double*, 2>, double *AbsDetjk, 
                 const double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w;
  double *coeff, c0, c1, c2, c3, c5;
  double e0, e1, e2, e3;
  double loc0, loc1, loc2, loc3, loc4;

  //static double delta0 = TDatabase::ParamDB->DELTA0;
  //static double delta1 = TDatabase::ParamDB->DELTA1;
  //double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double delta;

  loc0 = 0.0;
  loc1 = 0.0;
  loc2 = 0.0;
  loc3 = 0.0;
  loc4 = 0.0;

  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
    c0 = coeff[0];
    c1 = coeff[1];
    c2 = coeff[2];
    c3 = coeff[3];
    c5 = MAX(fabs(c1),fabs(c2));

    // if (X[i]>TDatabase::ParamDB->P6) continue;
    // if (Y[i]>TDatabase::ParamDB->P6) continue;

    // SUPG parameter
    delta = Compute_SDFEM_delta<2>(hK, c0, {{c1, c2}}, c3, c5);
    
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    // error in solution
    e0 = deriv[0]-exactval[0];
    if (fabs(e0)>loc4)
      loc4 = fabs(e0);
    loc0 += w*e0*e0;
    if (fabs(e0) > loc3)
	     loc3 = fabs(e0);
    
    // error in derivative of solution
    e1 = deriv[1]-exactval[1];
    loc1 += w*e1*e1;
    e2 = deriv[2]-exactval[2];
    loc1 += w*e2*e2;
    
    // sd error
    // THIS IS ONLY CORRECT IF DIV(b) = 0
    e3 = c1*e1+c2*e2;
    loc2 += w*(c0*(e1*e1+e2*e2) + c3*e0*e0 + delta*e3*e3);
    
  } // endfor i
  LocError[0] = loc0;
  LocError[1] = loc1;
  LocError[2] = loc2;
  LocError[3] = loc3;

  //cout << "LocError[3]: " << LocError[3] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// determine L2, H1 and SDFEM error, in (0,P6)^2
void SDFEMErrorsSmooth(int N_Points, double *X, double *Y, double *AbsDetjk, 
                 const double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError)
{
  int i, sd_type = TDatabase::ParamDB->SDFEM_TYPE;
  double *deriv, *exactval, w;
  double *coeff, c0, c1, c2, c5, alpha;
  double e0, e1, e2, e3;
  double loc0, loc1, loc2;

  static double delta0 = TDatabase::ParamDB->DELTA0;
  static double delta1 = TDatabase::ParamDB->DELTA1;
  double delta;

  loc0 = 0.0;
  loc1 = 0.0;
  loc2 = 0.0;

  for(i=0;i<N_Points;i++)
  {
    if (X[i]>TDatabase::ParamDB->P6) continue;
    if (Y[i]>TDatabase::ParamDB->P6) continue;

    coeff = coeffs[i];
    c0 = coeff[0];
    c1 = coeff[1];
    c2 = coeff[2];
    c5 = MAX(fabs(c1),fabs(c2));

    // if (X[i]>TDatabase::ParamDB->P6) continue;
    // if (Y[i]>TDatabase::ParamDB->P6) continue;

    if (sd_type==0)
    {
       if(c0 < hK*c5)
          delta = delta0 * hK/c5;
       else
          delta = delta1 *hK*hK/c0 ;
    }
    else
    {
       if (c5>0)
       {
          alpha = c5*hK/(2*c0);
          delta = hK*(1/tanh(alpha) - 1/alpha)/(2*c5);
       }
       else
          delta = 0;
    }
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    e0 = deriv[0]-exactval[0];
    loc0 += w*e0*e0;

    e1 = deriv[1]-exactval[1];
    loc1 += w*e1*e1;
    e2 = deriv[2]-exactval[2];
    loc1 += w*e2*e2;

    e3 = c1*e1+c2*e2;
    loc2 += w*(c0*(e1*e1+e2*e2) + e0*e0 + delta*e3*e3);

  } // endfor i

  LocError[0] = loc0;
  LocError[1] = loc1;
  LocError[2] = loc2;

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// determine L2, H1 and SDFEM error for smooth region in the
// example JohnMaubachTobiska1997 (x-0.5)^2+(y-0.5)^2 > r^2 
void SDFEMErrorsSmooth_JohnMaubachTobiska1997
(int N_Points, double *X, double *Y, double *AbsDetjk, 
 const double *Weights, double hK, double **Der, double **Exact,
 double **coeffs, double *LocError)
{
    int i, sd_type = TDatabase::ParamDB->SDFEM_TYPE;
    double *deriv, *exactval, w;
    double *coeff, c0, c1, c2, c3, c5, alpha;
    double e0, e1, e2, e3;
    double loc0, loc1, loc2;
    
    static double delta0 = TDatabase::ParamDB->DELTA0;
    static double delta1 = TDatabase::ParamDB->DELTA1;
    double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
    double delta, norm_b, x, y, r2=0.3*0.3;
    
    loc0 = 0.0;
    loc1 = 0.0;
    loc2 = 0.0;
    
    for(i=0;i<N_Points;i++)
    {
	x = X[i]-0.5;
	y = Y[i]-0.5;   
	if (x*x+y*y<=r2)
	  continue;
	coeff = coeffs[i];
	c0 = coeff[0];
	c1 = coeff[1];
	c2 = coeff[2];
	c3 = coeff[3];
	c5 = MAX(fabs(c1),fabs(c2));
		
	// SUPG parameter
	switch(sd_type)
	{
	    case 0:
		if(c0 < hK*c5)
		    delta = delta0 * hK/c5;
		else
		    delta = delta1 *hK*hK/c0;
		break;
	    case 1:
		norm_b = sqrt(c1*c1+c2*c2);
		//norm_b = linfb;
		if (norm_b > 0)
	    {
		alpha = norm_b*hK/(2*c0);
		delta = hK*(1/tanh(alpha) - 1/alpha)/(2*norm_b);
	    }
		else
		    delta = 0;
		break;
	    case 9:
		if (c0 <= hK)
		{
		    delta = delta0*time_step;	
                    //delta = delta0*time_step*time_step;
		}
		else
		{
		    delta = delta0*time_step;
                    //delta = delta0*time_step*time_step;
		}
		break;
	    case 10: 
		// second estimate in paper with Julia Novo
		// get the unscaled parameters
		if(c0 <= hK)
		    delta = delta0 * hK;
		else
		    delta = delta0 *hK*hK/c0 ;
		break;
	    case 11:
		// for estimate in paper with Julia Novo
		norm_b = sqrt(c1*c1+c2*c2);
		delta = delta0 * hK * sqrt(time_step)/norm_b;
		break;
	    default:
	      //OutPut("CHECK IF CORRECT DELTA IN ERROR COMPUTATION !!!"<<endl);
		if(c0 < hK*c5)
		    delta = delta0 * hK/c5;
		else
		    delta = delta1 *hK*hK/c0;
		break;
	}
	deriv = Der[i];
	exactval = Exact[i];
	w = Weights[i]*AbsDetjk[i];
	
	// error in solution
	e0 = deriv[0]-exactval[0];
	loc0 += w*e0*e0;
	
	// error in derivative of solution
	e1 = deriv[1]-exactval[1];
	loc1 += w*e1*e1;
	e2 = deriv[2]-exactval[2];
	loc1 += w*e2*e2;
	
	// sd error
	// THIS IS ONLY CORRECT IF DIV(b) = 0
	e3 = c1*e1+c2*e2;
	loc2 += w*(c0*(e1*e1+e2*e2) + c3*e0*e0 + delta*e3*e3);
    } // endfor i
    LocError[0] = loc0;
    LocError[1] = loc1;
    LocError[2] = loc2;
}

// determine errors to interpolant
// paper with Julia Novo
void SDFEMErrorsInterpolant(int N_Points, double *, double *, double *AbsDetjk, 
                 const double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w;
  double *coeff, c0, c1, c2, c5;
  double e0, e1, e2, e3, e4;
  double loc0, loc1, loc2, loc3;
  double delta, c_inv;

  loc0 = 0.0;
  loc1 = 0.0;
  loc2 = 0.0;
  loc3 = 0.0;

  switch (TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE)
  {
    case 1:  c_inv = 1.0;
    break;
    case 2:
        // triangle
        c_inv = sqrt(48.0);
        // quad
        c_inv = sqrt(24.0);
        break;
    case 3:
        // triangle
        c_inv = sqrt((435+sqrt(26025.0))/4.0);
        // quad
        c_inv = sqrt((244+sqrt(9136.0))/3.0);
        break;
    default:
      OutPut("c_inv not defined " << endl);
      exit(4711);
  }

  for(i=0;i<N_Points;i++)
  {
    // compute SUPG parameter
    coeff = coeffs[i];
    c0 = coeff[0];
    c1 = coeff[1];
    c2 = coeff[2];
    // double c3 = coeff[3];
    c5 = MAX(fabs(c1),fabs(c2));
    delta = Compute_SDFEM_delta<2>(hK, coeff[0], {{coeff[1], coeff[2]}}, coeff[3], c5);
    if (TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 120814)
      delta = TDatabase::ParamDB->INTERNAL_P1_Array[TDatabase::ParamDB->INTERNAL_LEVEL];

    // NOTE: order in derivatives in Derivatives_SD
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    // error in solution
    e0 = deriv[2]-exactval[0];
    loc0 += w*e0*e0/delta;
    
    // error in streamline derivative of solution
    e1 = c1*(deriv[0]-exactval[1]);
    loc1 += w*delta*e1*e1;
    e2 = c2*(deriv[1]-exactval[2]);
    loc1 += w*delta*e2*e2;

    // first additional term, with gradient
    e3 = c0 * c_inv * (deriv[0]-exactval[1])/ hK;
    loc2 += w * 16.0 * delta * e3 * e3;
    e3 = c0 * c_inv * (deriv[1]-exactval[2])/ hK;
    loc2 += w * 16.0 * delta * e3 * e3;

    // second additional term, with Laplacian
     e4 = c0 * (deriv[3] + deriv[4] - exactval[3]);
    loc3 += w * 8.0 * delta * e4 * e4;
  } // endfor i
  LocError[0] = loc0;
  LocError[1] = loc1;
  LocError[2] = loc2;
  LocError[3] = loc3;
 //cout << "LocError[3]: " << LocError[3] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// determine L2, H1 and SDFEM error for Oseen
void SPGErrorsOseen(int N_Points, double *, double *, double *AbsDetjk, 
                 const double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w;
  double *coeff, eps, u1, u2, c;
  double e0, e1, e2, e3;
  double loc0, loc1, loc2, loc3;

  //double delta0 = TDatabase::ParamDB->DELTA0;
  //double delta1 = TDatabase::ParamDB->DELTA1;
  //double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double delta;

  loc0 = 0.0;
  loc1 = 0.0;
  loc2 = 0.0;
  loc3 = 0.0;

  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
    eps = coeff[0]; 
    u1 = coeff[3];  
    u2 = coeff[4];
    c = coeff[5];

    // stabilization parameter
    delta =  SUPG_Parameter(hK, eps, u1, u2, c);
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    // error in solution
    e0 = deriv[0]-exactval[0];
    loc0 += w*e0*e0;
    if (fabs(e0) > loc3)
       loc3 = fabs(e0);
    
    // error in derivative of solution
    e1 = deriv[1]-exactval[1];
    loc1 += w*e1*e1;
    e2 = deriv[2]-exactval[2];
    loc1 += w*e2*e2;
    // sd error
    e3 = u1*e1+u2*e2;
    loc2 += w*(eps*(e1*e1+e2*e2) + c*e0*e0 + delta*e3*e3);
   } // endfor i
  LocError[0] = loc0;
  LocError[1] = loc1;
  LocError[2] = loc2;
  LocError[3] = loc3;
  //cout << "LocError[3]: " << LocError[3] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// determine L2, H1 and pressure part of SUPG error for Oseen
void SPGErrorsOseenPressure(int N_Points, double *, double *, double *AbsDetjk, 
                 const double *Weights, double hK, double **Der, double **Exact,
                 double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w;
  double *coeff, eps, u1, u2, c;
  double e0, e1, e2, e3;
  double loc0, loc1, loc2;

  //double delta0 = TDatabase::ParamDB->DELTA0;
  //double delta1 = TDatabase::ParamDB->DELTA1;
  //double time_step = TDatabase::TimeDB->CURRENTTIMESTEPLENGTH;
  double delta;

  loc0 = 0.0;
  loc1 = 0.0;
  loc2 = 0.0;
 
  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
    eps = coeff[0]; 
    u1 = coeff[3];  
    u2 = coeff[4];
    c = coeff[5];

    // stabilization parameter
    delta =  SUPG_Parameter(hK, eps, u1, u2, c);
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    // error in solution
    e0 = deriv[0]-exactval[0];
    loc0 += w*e0*e0;
    
    // error in derivative of solution
    e1 = deriv[1]-exactval[1];
    loc1 += w*e1*e1;
    e2 = deriv[2]-exactval[2];
    loc1 += w*e2*e2;
    // sd error
    e3 = u1*e1+u2*e2;
    loc2 += w*delta*e3*e3;
    
  } // endfor i
  LocError[0] = loc0;
  LocError[1] = loc1;
  LocError[2] = loc2;
  //cout << "LocError[3]: " << LocError[3] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}


// determine L1 error, 2D
void L1Error(int N_Points, std::array<double*, 2> xy, double *AbsDetjk, 
             const double *Weights, double hK, 
             double **Der, double **Exact,
             double **, double *LocError)
{
  int i, index;
  double *deriv, *exactval, w, t, area, v[3], va[3], xa[3], ya[3];
  double a[2], b[2], fac, areaa;
  LocError[0] = 0.0;

  // conforming linear finite elemens
  if (hK==-4711)
  {
      // area of the mesh cell
      area = AbsDetjk[0]/2;
      // compute values of the fe function in the vertices of the mesh cell
      // giveen are the values in the midpoints of the edges
      // Der[0][0] - edge between P0 and P1
      // Der[1][0] - edge between P1 and P2
      // Der[2][0] - edge between P2 and P0
      v[0] = Der[0][0] - Der[1][0] +   Der[2][0];
      v[1] = Der[1][0] - Der[2][0] +   Der[0][0];
      v[2] = Der[2][0] - Der[0][0] +   Der[1][0];
      //for (i=0;i<3;i++)
      //          OutPut(X[i] << " " << Y[i] << " " << v[i] << endl);
      // four different cases
      // case one and two: no change of sign
      if ((v[0] >= 0)&& (v[1] >= 0) && (v[2] >= 0))
      {
          LocError[0] = area*(v[0]+v[1]+v[2])/3;
          return;
      }
      if ((v[0] <= 0)&& (v[1] <= 0) && (v[2] <= 0))
      {
          LocError[0] = -area*(v[0]+v[1]+v[2])/3;
          return;
      }
      // change of sign, one value has different sign than the other two ones
      // this value gets index 2
      t = v[0]*v[1]*v[2];
      // two negative, one positive values
      if (t>0)
      {
          for (i=0;i<3;i++)
          {
              if (v[i] > 0)
              {
                  index = i;
                  break;
              }
          }
      }
      else
      // two positive, one negative values
      {
          for (i=0;i<3;i++)
          {
              if (v[i] < 0)
              {
                  index = i;
                  break;
              }
          }
      }
      for (i=0;i<3;i++)
      {
          va[i] = v[(index+1+i)%3]; 
          xa[i] = xy[0][(index+1+i)%3]; 
          ya[i] = xy[1][(index+1+i)%3]; 
      }     
      //for (i=0;i<3;i++)
      //          OutPut(xa[i] << " " << ya[i] << " " << va[i] << endl);
      // compute area of triangle with different sign of the fe function
      a[0] = xa[0] - xa[2];
      a[1] = ya[0] - ya[2];
      fac = -va[2]/(va[0]-va[2]);
      a[0] *= fac;
      a[1] *= fac;

      b[0] = xa[1] - xa[2];
      b[1] = ya[1] - ya[2];
      fac = -va[2]/(va[1]-va[2]);
      b[0] *= fac;
      b[1] *= fac;

      areaa = fabs(a[0]*b[1]-b[0]*a[1])/2;

      //OutPut(a[0] << " " << a[1] << " " << b[0] << " " << b[1] << " " << areaa<< endl);
      LocError[0] = area*(v[0]+v[1]+v[2])/3 - 2*areaa*va[2]/3;
      LocError[0] = fabs(LocError[0]);
  }
  else
  // use quadrature rule for other types of finite elements
  {
      for(i=0;i<N_Points;i++)
      {
          deriv = Der[i];
          exactval = Exact[i];
          w = Weights[i]*AbsDetjk[i];
          
          t = deriv[0]-exactval[0];
          LocError[0] += w*fabs(t);
      } // endfor i
  }
  //OutPut("LocError(L1): " << LocError[0] << endl);
}

// determine deformation tensor error
void DeformationTensorError(int N_Points, double *, double *,
                double *AbsDetjk, const double *Weights, double,
                double **Der, double **Exact,
                double **, double *LocError)
{
  int i;
  double *deriv_x, *exactval, *deriv_y, *exactval1, w, t;

  LocError[0] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    // first component
    deriv_x = Der[i];
    exactval = Exact[i];
    // second component
    deriv_y = Der[i+N_Points];
    exactval1 = Exact[i+N_Points];
    // weight
    w = Weights[i]*AbsDetjk[i];

    // left upper term
    t = deriv_x[1]-exactval[1];
    LocError[0] += w*t*t;
    // right upper and left lower term
    t = ((deriv_x[2]-exactval[2])+(deriv_y[1]-exactval1[1]));
    LocError[0] += w*t*t/2.0;
    // right lower term
    t = deriv_y[2]-exactval1[2];
    LocError[0] += w*t*t;

  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// determine L2 and H1 error, 2D
void H1Norm(int N_Points, double *, double *, double *AbsDetjk, 
            const double *Weights, double,
            double **Der, double **Exact,
            double **, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t;

  LocError[0] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;

    t = deriv[1]-exactval[1];
    LocError[0] += w*t*t;
    t = deriv[2]-exactval[2];
    LocError[0] += w*t*t;
  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}
// compute the error in the divergence
void DivergenceError(int N_Points, double *, double *,
		     double *AbsDetjk, const double *Weights, double,
		     double **Der, double **,
		     double **, double *LocError)
{
  int i;
  double *deriv_x, *deriv_y, w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;
 
  for(i=0;i<N_Points;i++)
  {
    // first component
    deriv_x = Der[i];
    // second component
    deriv_y = Der[i+N_Points];
    // weight
    w = Weights[i]*AbsDetjk[i];

     // u_x + v_y
    t = fabs(deriv_x[1] + deriv_y[2]);
    LocError[0] += w*t;
    LocError[1] += w*t*t;    
  } // endfor i

  //cout << "LocError[0]: " << LocError[0] << endl;
}

// compute the error in the grad-div term for Oseen
void DivergenceErrorGradDivOseen(int N_Points, double *, double *,
         double *AbsDetjk, const double *Weights, double hK, 
         double **Der, double **,
         double **coeffs, double *LocError)
{
  int i;
  double *deriv_x, *deriv_y, w, t, nu, b1, b2, mu, *coeff;

  LocError[0] = 0.0;
 
  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
    nu = coeff[0];
    b1 = coeff[3];
    b2 = coeff[4];

   // get stabilization parameters
    mu = graddiv_parameterOseen(hK, nu, b1, b2);
   // first component
    deriv_x = Der[i];
    // second component
    deriv_y = Der[i+N_Points];
    // weight
    w = Weights[i]*AbsDetjk[i];
    
     // u_x + v_y
    t = fabs(deriv_x[1] + deriv_y[2]);
    LocError[0] += w*mu*t*t;    
  } // endfor i

  //cout << "LocError[0]: " << LocError[0] << endl;
}


// mesh cell parameters for shock capturing scheme DC_CD
void Parameters_DC_CD(int N_Points, double *, double *, double *AbsDetjk, 
           const double *Weights, double,
           double **Der, double **Exact,
           double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t, *coeff;
  double eps, b1, b2, c, f;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
    eps = coeff[0];
    b1 = coeff[1];
    b2 = coeff[2];
    c = coeff[3];
    f = coeff[4];

    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[2]-exactval[2];
   // LocError[0] += w*t*t;

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;
    t = deriv[1]-exactval[1];
    LocError[0] += w*t*t;
   
    //   t= -eps*deriv[3]+ b1*deriv[1]+b2*deriv[2] + c*deriv[0] - f;
    t= -eps*(deriv[3] + deriv[4]) + b1*deriv[0]+b2*deriv[1] + c*deriv[2] - f;
    LocError[1] += w*t*t;
    

  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// mesh cell values for gradient and residual 
void Parameters_Gradient_Residual(int N_Points, double *, double *, double *AbsDetjk,
           const double *Weights, double,
           double **Der, double **Exact,
           double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t, *coeff;
  double eps, b1, b2, c, f;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
    eps = coeff[0];
    b1 = coeff[1];
    b2 = coeff[2];
    c = coeff[3];
    f = coeff[4];

    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    // L^2 norm of gradient gradient
    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;
    t = deriv[1]-exactval[1];
    LocError[0] += w*t*t;

    // L^2 norm of residual
    t= -eps*(deriv[3] + deriv[4]) + b1*deriv[0]+b2*deriv[1] + c*deriv[2] - f;
    LocError[1] += w*t*t;
  } // endfor i
}

#endif // 2D

#ifdef __3D__
void ComputeVorticityDivergence(TFESpace3D *, TFEFunction3D *u1, 
                                TFEFunction3D *u2, TFEFunction3D *u3,
                                TFESpace3D *vorticity_space, 
                                TFEFunction3D *vort1, 
                                TFEFunction3D *vort2, TFEFunction3D *vort3,
                                double *div)
                                
{
  int i, j, N_Cells, index, *DOF, N_loc_dofVort, N_Vort;
  int *GlobalNumbersVort, *BeginIndexVort, *N_Found;
  TCollection *Coll;
  TBaseCell *cell;
  FE3D CurrentElementVort;
  TFE3D *Element;
  TNodalFunctional3D *nf;
  RefTrans3D RefTrans;
  const double *xi_ref, *eta_ref, *zeta_ref;  //values[4];
  double AbsDetjkVort[MaxN_QuadPoints_3D];
  double X_orig[MaxN_PointsForNodal3D], Y_orig[MaxN_PointsForNodal3D];
  double Z_orig[MaxN_PointsForNodal3D];
  double val[4], *v1, *v2, *v3;  // *div1;

  // get information of the numbering of the degrees of freedom
  // of the vorticity
  GlobalNumbersVort = vorticity_space->GetGlobalNumbers();
  BeginIndexVort = vorticity_space->GetBeginIndex();
  N_Vort = vorticity_space->GetN_DegreesOfFreedom();
  v1 =  vort1->GetValues();
  memset(v1,0,N_Vort*SizeOfDouble);
  v2 =  vort2->GetValues();
  memset(v2,0,N_Vort*SizeOfDouble);
  v3 =  vort3->GetValues();
  memset(v3,0,N_Vort*SizeOfDouble);
  memset(div,0,N_Vort*SizeOfDouble);

  N_Found = new int[N_Vort];
  memset(N_Found,0,N_Vort*SizeOfInt);

  // get pointer to set of mesh cells which define the fe space
  Coll = vorticity_space->GetCollection();
  // get number of mesh cells
  N_Cells = Coll->GetN_Cells();

  // loop over all mesh cells
  for(i=0;i<N_Cells;i++)
  {
    // get current mesh cell
    cell = Coll->GetCell(i);        
    
    // compute geometric positions of the fe nodes
    // get id of finite element in current mesh cell
    CurrentElementVort = vorticity_space->GetFE3D(i, cell);
    // get fe from its id
    Element = TFEDatabase3D::GetFE3D(CurrentElementVort);
    // get reference transformation 
    RefTrans = Element->GetRefTransID();
    TFEDatabase3D::SetCellForRefTrans(cell,RefTrans);
    // get pointer to the nodal functionals (fe nodes) of the fe 
    // (in ref. cell)
    nf = Element->GetNodalFunctional3D();
    // get number and coordinates of local dof in ref cell
    // xi_ref, eta_ref are pointers
    nf->GetPointsForAll(N_loc_dofVort, xi_ref, eta_ref, zeta_ref);
      
    // get coordinates of fe nodes in original cell
    // input: RefTrans, N_loc_dof, xi_ref, eta_ref
    // output : X_orig, Y_orig - pointers to arrays with coordinates
    //          AbsDetjk - same as above    
    TFEDatabase3D::GetOrigFromRef(RefTrans,N_loc_dofVort, 
                                  xi_ref, eta_ref, zeta_ref,
                                  X_orig, Y_orig, Z_orig,
                                  AbsDetjkVort);
   
    DOF = GlobalNumbersVort + BeginIndexVort[i];
    for (j=0;j<N_loc_dofVort;j++)
    {
      // values computed with FindGradient have to averaged on a periodic boundary !!! 
      // use FindGradientLocal
      index = DOF[j];
      u1->FindGradientLocal(cell,i,X_orig[j],Y_orig[j],Z_orig[j],val);
      v2[index] += val[3]; // u1_z
      v3[index] -= val[2]; // u1_y
      div[index] += val[1]; // u1_x
      //OutPut("val1 " << val[1] << " " << val[2] << " " << val[3]<< endl);
      u2->FindGradientLocal(cell,i,X_orig[j],Y_orig[j],Z_orig[j],val); 
      v1[index] -= val[3];
      v3[index] += val[1];
      div[index] +=val[2]; // u2_y
      //OutPut("val2 " << val[1] << " " << val[2] << " " << val[3]<< endl);

      u3->FindGradientLocal(cell,i,X_orig[j],Y_orig[j],Z_orig[j],val); 
      v1[index]  += val[2];
      v2[index]  -= val[1];
      div[index] += val[3]; // u3_z
      //OutPut("val3 " << val[1] << " " << val[2] << " " << val[3]<< endl);

      N_Found[index]++;
    }    
  } 
 
  for (i=0;i<N_Vort;i++)
  {
    v1[i]/=N_Found[i];
    v2[i]/=N_Found[i];
    v3[i]/=N_Found[i];
    //OutPut( v1[i] << " " << v2[i] << " " << v3[i] << endl);
    div[i]/=N_Found[i];
  }
  delete N_Found;

}



// determine L2 and H1 error, 3D
void L2H1Errors(int N_Points, std::array<double*, 3>,
                double *AbsDetjk, 
                const double *Weights, double,
                double **Der, double **Exact,
                double **, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;

    t = deriv[1]-exactval[1];
    LocError[1] += w*t*t;
    t = deriv[2]-exactval[2];
    LocError[1] += w*t*t;
    t = deriv[3]-exactval[3];
    LocError[1] += w*t*t;
  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}
void L2H1ErrorsSmooth(int N_Points, double *X, double *Y, double *Z, 
                double *AbsDetjk, 
                const double *Weights, double,
                double **Der, double **Exact,
                double **, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    if (X[i]>TDatabase::ParamDB->P6) continue;
    if (Y[i]>TDatabase::ParamDB->P6) continue;
    if (Z[i]>TDatabase::ParamDB->P6) continue;
    
    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;

    t = deriv[1]-exactval[1];
    LocError[1] += w*t*t;
    t = deriv[2]-exactval[2];
    LocError[1] += w*t*t;
    t = deriv[3]-exactval[3];
    LocError[1] += w*t*t;
  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

void L2DivH1Errors(int N_Points, std::array<double*, 3>,
                   double *AbsDetjk, const double *Weights, double,
                   double **Der, double **Exact, double **,
                   double *LocError)
{
  LocError[0] = 0.0;
  LocError[1] = 0.0;
  LocError[2] = 0.0;
  
  for(int i = 0; i < N_Points; i++)
  {
    double w = Weights[i]*AbsDetjk[i];
    
    // L2-error
    double t = Der[i][0] - Exact[i][0]; // x-component
    LocError[0] += w*t*t;
    t = Der[i][4] - Exact[i][5];        // y-component
    LocError[0] += w*t*t;
    t = Der[i][8] - Exact[i][10];       // z-component
    LocError[0] += w*t*t;
    
    // L2-error of divergence
    t  = Der[i][1] - Exact[i][1];  // x-derivative of x-component
    t += Der[i][6] - Exact[i][7];  // y-derivative of y-component
    t += Der[i][11]- Exact[i][13]; // z-derivative of z-component
    LocError[1] += w*t*t;
    
    // H1 semi norm
    t = Der[i][1] - Exact[i][1];  // x-derivative of x-component
    LocError[2] += w*t*t;
    t = Der[i][2] - Exact[i][2];  // y-derivative of x-component
    LocError[2] += w*t*t;
    t = Der[i][3] - Exact[i][3];  // z-derivative of x-component
    LocError[2] += w*t*t;
    t = Der[i][5] - Exact[i][6];  // x-derivative of y-component
    LocError[2] += w*t*t;
    t = Der[i][6] - Exact[i][7];  // y-derivative of y-component
    LocError[2] += w*t*t;
    t = Der[i][7] - Exact[i][8];  // z-derivative of y-component
    LocError[2] += w*t*t;
    t = Der[i][9] - Exact[i][11]; // x-derivative of z-component
    LocError[2] += w*t*t;
    t = Der[i][10]- Exact[i][12]; // y-derivative of z-component
    LocError[2] += w*t*t;
    t = Der[i][11]- Exact[i][13]; // z-derivative of z-component
    LocError[2] += w*t*t;
  }
}

// determine L1 error, 3D
void L1Error(int, std::array<double*, 3>, double *, const double *, double,
             double **, double **, double **, double *LocError)
{
    OutPut("computation of L1-error not implemented !!!" <<endl);
    LocError[0] = 0;
}

// determine deformation tensor error
void DeformationTensorError(int N_Points, double *, double *, double *,
                double *AbsDetjk, const double *Weights, double, 
                double **Der, double **Exact,
                double **, double *LocError)
{
  int i;
  double *deriv_x, *exactval, *deriv_y, *exactval1, *deriv_z, *exactval2, w, t;

  LocError[0] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    // first component
    deriv_x = Der[i];
    exactval = Exact[i];
    // second component
    deriv_y = Der[i+N_Points];
    exactval1 = Exact[i+N_Points];
    // third component
    deriv_z = Der[i+2*N_Points];
    exactval2 = Exact[i+2*N_Points];
    // weight
    w = Weights[i]*AbsDetjk[i];

    // left upper term
    t = deriv_x[1]-exactval[1];
    LocError[0] += w*t*t;
    // 12 and 21 term 
    t = ((deriv_x[2]-exactval[2])+(deriv_y[1]-exactval1[1]));
    LocError[0] += w*t*t/2.0;
    // 22 term
    t = deriv_y[2]-exactval1[2];
    LocError[0] += w*t*t;
    // 13 and 31 term
    t = ((deriv_x[3]-exactval[3])+(deriv_z[1]-exactval2[1]));
    LocError[0] += w*t*t/2.0;
    // 23 and 32 term
    t = ((deriv_y[3]-exactval1[3])+(deriv_z[2]-exactval2[2]));
    LocError[0] += w*t*t/2.0;
    // 33 term
    t = deriv_z[3]-exactval2[3];
    LocError[0] += w*t*t;       
  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

// compute the error in the divergence
void DivergenceError(int N_Points, double *, double *, double *,
		     double *AbsDetjk, const double *Weights, double,
		     double **Der, double **,
		     double **, double *LocError)
{
  int i;
  double *deriv_x, *deriv_y, *deriv_z, w, t;

  LocError[0] = 0.0;
  LocError[1] = 0.0;
 
  for(i=0;i<N_Points;i++)
  {
    // first component
    deriv_x = Der[i];
    // second component
    deriv_y = Der[i+N_Points];
    // third component
    deriv_z = Der[i+2*N_Points];
    // weight
    w = Weights[i]*AbsDetjk[i];

     // u_x + v_y + w_z
    t = fabs(deriv_x[1] + deriv_y[2] + deriv_z[3]);
    LocError[0] += w*t;
    LocError[1] += w*t*t;    
  } // endfor i

  //cout << "LocError[0]: " << LocError[0] << endl;
}
// mesh cell parameters for shock capturing scheme DC_CD
void Parameters_DC_CD(int N_Points, double *, double *, double *,
                      double *AbsDetjk, 
                      const double *Weights, double,
                      double **Der, double **Exact,
                      double **coeffs, double *LocError)
{
  int i;
  double *deriv, *exactval, w, t, *coeff;
  double b1, b2, b3, c, f;  //eps;

  LocError[0] = 0.0;
  LocError[1] = 0.0;

  for(i=0;i<N_Points;i++)
  {
    coeff = coeffs[i];
//    eps = coeff[0];
    b1 = coeff[1];
    b2 = coeff[2];
    b3 = coeff[3];
    c = coeff[4];
    f = coeff[5];

    deriv = Der[i];
    exactval = Exact[i];
    w = Weights[i]*AbsDetjk[i];

    t = deriv[3]-exactval[3];
    LocError[0] += w*t*t;

    t = deriv[0]-exactval[0];
    LocError[0] += w*t*t;
    t = deriv[1]-exactval[1];
    LocError[0] += w*t*t;
    t = deriv[2]-exactval[2];
    LocError[0] += w*t*t;
    
    //t= -eps*(deriv[4] + deriv[5] + deriv[6]) + b1*deriv[1]+ b2*deriv[2] + b3*deriv[3] + c*deriv[0] - f;
    t= b1*deriv[0]+ b2*deriv[1] + b3*deriv[2] + c*deriv[3] - f;
    LocError[1] += w*t*t;
    

  } // endfor i

  // cout << "LocError[0]: " << LocError[0] << endl;
  // cout << "LocError[1]: " << LocError[1] << endl;
}

/*******************************************************************************/
//
// compute the Q criterion for a flow field
//
/*******************************************************************************/
void Q_criterion(TCollection *Coll,
TFEFunction3D *velocity1, TFEFunction3D *velocity2,
TFEFunction3D *velocity3, double *Qcrit)
{
  int i, j, N_V, N_Cells;
  double x[8],y[8],z[8],values[4], eps = 1e-6;
  double grad_velo_xx, grad_velo_xy, grad_velo_xz,  grad_velo_yx, grad_velo_yy,  grad_velo_yz;
  double grad_velo_zx, grad_velo_zy, grad_velo_zz; 
  double x_sp, y_sp, z_sp, val;
 
  TBaseCell *cell;

   // number of cells
  N_Cells = Coll->GetN_Cells();
  //loop over the mesh cells of the global grid
  for(i=0;i<N_Cells;i++)
  {
    // get cell
    cell = Coll->GetCell(i);
    // number of vertices per cell
    N_V = cell->GetN_Vertices();
    // compute barycenter
    x_sp=0.;
    y_sp=0.;
    z_sp=0.;
    for (j=0;j<N_V;j++)
    {
      // read coordinates of the mesh cell
      cell->GetVertex(j)->GetCoords(x[j], y[j], z[j]);
      x_sp += x[j];
      y_sp += y[j];
      z_sp += z[j];
     }
     x_sp /= N_V;
     y_sp /= N_V;
     z_sp /= N_V;     
      // compute gradient in barycenter
      velocity1->FindGradientLocal(cell,i,x_sp,y_sp,z_sp,values);
      grad_velo_xx = values[1];
      grad_velo_xy = values[2];
      grad_velo_xz = values[3];
      velocity2->FindGradientLocal(cell,i,x_sp,y_sp,z_sp,values);
      grad_velo_yx = values[1];
      grad_velo_yy = values[2];
      grad_velo_yz = values[3];
      velocity3->FindGradientLocal(cell,i,x_sp,y_sp,z_sp,values);
      grad_velo_zx = values[1];
      grad_velo_zy = values[2];
      grad_velo_zz = values[3];
      
     // Q criterion : 0.5 * (Q:Q-S:S)
     val=-0.5*((grad_velo_xx*grad_velo_xx)
                 +(grad_velo_yy*grad_velo_yy)
                 +(grad_velo_zz*grad_velo_zz)
                 +(grad_velo_xy*grad_velo_yx)
                 +(grad_velo_xz*grad_velo_zx)
                 +(grad_velo_zy*grad_velo_yz)
                 +(grad_velo_yx*grad_velo_xy)
                 +(grad_velo_zx*grad_velo_xz)
                 +(grad_velo_yz*grad_velo_zy));

     OutPut(val << " ");
     if (fabs(val) < eps)
     {
  Qcrit[i]=0;
     }
     else
     {
       if (val > 0)
  Qcrit[i] = 1; 
  else 
    Qcrit[i]=-1;
     }
  }
}

#endif // 3D

#ifdef __2D__
// ========================================================================
// put l infinity norm of u in coeff5
// ========================================================================
void LInfU(int N_Points, double **, double **Params, TBaseCell *)
{
  int i;
  double max, u1, u2, u, *param;

  max = -1;

  for(i=0;i<N_Points;i++)
  {
    param = Params[i];
    u1 = param[0];
    u2 = param[1];

    u = MAX(fabs(u1), fabs(u2));

    if(u>max) max = u;
  }

  for(i=0;i<N_Points;i++)
    Params[i][2] = max;
}

void linfb(int N_Points, double **Coeffs, double **, TBaseCell *)
{
  int i;
  double max, *coeff, b1, b2, b;

  max = -1;

  for(i=0;i<N_Points;i++)
  {
    coeff = Coeffs[i];
    b1 = coeff[1];
    b2 = coeff[2];

    b = MAX(fabs(b1), fabs(b2));
    if(b>max) max = b;
  }

  for(i=0;i<N_Points;i++)
    Coeffs[i][5] = max;
}
void ave_l2b_quad_points(int N_Points, double **Coeffs, double **, TBaseCell *)
{
  int i;
  double max, *coeff, b1, b2, b;

  max = 0;

  for(i=0;i<N_Points;i++)
  {
    coeff = Coeffs[i];
    b1 = coeff[1];
    b2 = coeff[2];

    b= sqrt(b1*b1+b2*b2); 
    max += b;
  }

  max /= N_Points;

  for(i=0;i<N_Points;i++)
    Coeffs[i][5] = max;
}

#endif // 2D

#ifdef __3D__
void linfb(int N_Points, double **Coeffs, double **, TBaseCell *)
{
  int i;
  double max, *coeff, b1, b2, b3, b;

  max = -1;

  for(i=0;i<N_Points;i++)
  {
    coeff = Coeffs[i];
    b1 = coeff[1];
    b2 = coeff[2];
    b3 = coeff[3];

    b = MAX(fabs(b1), fabs(b2));
    b = MAX(fabs(b3), b);

    if(b>max) max = b;
  }

  for(i=0;i<N_Points;i++)
    Coeffs[i][6] = max;
}

void ave_l2b_quad_points(int N_Points, double **Coeffs, double **, TBaseCell *)
{
  int i;
  double max, *coeff, b1, b2, b3, b;

  max = 0;

  for(i=0;i<N_Points;i++)
  {
    coeff = Coeffs[i];
    b1 = coeff[1];
    b2 = coeff[2];
    b3 = coeff[3];

    b =  sqrt(b1*b1+b2*b2+b3*b3);

    max += b;
  }

  max /= N_Points;

  for(i=0;i<N_Points;i++)
    Coeffs[i][6] = max;
}
#endif // 3D


void ExactNull(double, double, double, double *values)
{
  values[0] =0;
  values[1] =0;
  values[2] =0;
  values[3] =0;
  values[4] =0;
}

void ExactNull(double, double, double *values)
{
  values[0] =0;
  values[1] =0;
  values[2] =0;
  values[3] =0;
  //values[0] = x*(1-x)*y*(1-y);
}

void BoundConditionVMM(int, double, BoundCond &cond)
{
   cond = NEUMANN;
}

void BoundConditionNoBoundCondition(int, double, BoundCond &cond)
{
   cond = NEUMANN;
}
void BoundConditionNoBoundCondition(double, double, double, BoundCond &cond)
{
   cond = NEUMANN;
}
void BoundaryValueHomogenous(int, double, double &value)
{
  value = 0;
}
void BoundaryValueHomogenous(double, double, double, double &value)
{
  value = 0;
}
void BoundaryValueNoBoundaryValue(int, double, double &value)
{
  value = 0;
}

void BoundConditionNSE(int, double, BoundCond &cond)
{
   cond = DIRICHLET;
}

void BoundaryConditionPressSep(int, double, BoundCond &cond)
{
   cond = NEUMANN;
}

void BoundaryValuePressSep(int, double, double &value)
{
  value = 0;
}
void BoundaryConditionPressSep3D(double, double, double, BoundCond &cond)
{
  cond = NEUMANN;
}

void BoundaryValuePressSep3D(double, double, double, double &value)
{
  value = 0;
}

void BoundaryConditionNewton(double, double, double, BoundCond &cond)
{
  cond = DIRICHLET;
}

void BoundaryValueNewton(double, double, double, double &value)
{
  value = 0;
}

// boundary condition for assembling in the FEM-FCT scheme
void BoundCondition_FEM_FCT(int, double, BoundCond &cond)
{
    cond = NEUMANN;
}

void BoundValue_FEM_FCT(int, double, double &value)
{
    value = 0;
}
void BoundCondition_FEM_FCT(double, double, double, BoundCond &cond)
{
    cond = NEUMANN;
}
void BoundValue_FEM_FCT(double, double, double, double &value)
{
  value = 0;
}

extern "C" {
  void Cpp_Print(char *s)
  {
    OutFile << s;
  }
}

// ========================================================================
// boundary values for auxiliary problem in Galdi/Layton model
// ========================================================================

// void BoundConditionAuxProblem(int i, double t, BoundCond &cond)
// {
//   cond = NEUMANN;
//   //cond = DIRICHLET;
// }
// 
// void BoundValueAuxProblem(int BdComp, double Param, double &value)
// {
//   value = 0;
// }

// ========================================================================
// boundary values for higher order fe in VMS
// ========================================================================

// void ho_BoundCondition(int i, double t, BoundCond &cond)
// {
//   cond = DIRICHLET;
// }
// 
// void ho_BoundValue(int BdComp, double Param, double &value)
// {
//   value = 0;
// }


void SetPolynomialDegree()
{
    if ((TDatabase::ParamDB->ANSATZ_ORDER >0)
        && (TDatabase::ParamDB->ANSATZ_ORDER <10))
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE =
            TDatabase::ParamDB->ANSATZ_ORDER;
        return;
    }
    
    if ((TDatabase::ParamDB->ANSATZ_ORDER <0)
        && (TDatabase::ParamDB->ANSATZ_ORDER > -10))
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE =
            -TDatabase::ParamDB->ANSATZ_ORDER;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == -101)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 1;
        return;
    }

    // Discontious Galerkin Methods
    if ((TDatabase::ParamDB->ANSATZ_ORDER <-10)
        && (TDatabase::ParamDB->ANSATZ_ORDER > -15))
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE =
            -TDatabase::ParamDB->ANSATZ_ORDER - 10;
        return;
    }
    if ((TDatabase::ParamDB->ANSATZ_ORDER <-100)
        && (TDatabase::ParamDB->ANSATZ_ORDER > -150))
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE =
            (-TDatabase::ParamDB->ANSATZ_ORDER - 100)/10;
        return;
    }

    //==========LOCALPROJECTION=================
    if (TDatabase::ParamDB->ANSATZ_ORDER == 100)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 1;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 201
      || TDatabase::ParamDB->ANSATZ_ORDER == 200
      || TDatabase::ParamDB->ANSATZ_ORDER == 211
      || TDatabase::ParamDB->ANSATZ_ORDER == 221
      || TDatabase::ParamDB->ANSATZ_ORDER == 222)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 2;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 302
      || TDatabase::ParamDB->ANSATZ_ORDER == 301
      || TDatabase::ParamDB->ANSATZ_ORDER == 312
      || TDatabase::ParamDB->ANSATZ_ORDER == 322)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 3;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 403
      || TDatabase::ParamDB->ANSATZ_ORDER == 402
      || TDatabase::ParamDB->ANSATZ_ORDER == 413
      || TDatabase::ParamDB->ANSATZ_ORDER == 423)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 4;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 504
      || TDatabase::ParamDB->ANSATZ_ORDER == 503
      || TDatabase::ParamDB->ANSATZ_ORDER == 514
      || TDatabase::ParamDB->ANSATZ_ORDER == 524)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 5;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 605
      || TDatabase::ParamDB->ANSATZ_ORDER == 604
      || TDatabase::ParamDB->ANSATZ_ORDER == 615
      || TDatabase::ParamDB->ANSATZ_ORDER == 625)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 6;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 706
      || TDatabase::ParamDB->ANSATZ_ORDER == 705
      || TDatabase::ParamDB->ANSATZ_ORDER == 716
      || TDatabase::ParamDB->ANSATZ_ORDER == 726)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 7;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 807
      || TDatabase::ParamDB->ANSATZ_ORDER == 806
      || TDatabase::ParamDB->ANSATZ_ORDER == 817
      || TDatabase::ParamDB->ANSATZ_ORDER == 827)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 8;
        return;
    }
    if (TDatabase::ParamDB->ANSATZ_ORDER == 908
      || TDatabase::ParamDB->ANSATZ_ORDER == 907
      || TDatabase::ParamDB->ANSATZ_ORDER == 918
      || TDatabase::ParamDB->ANSATZ_ORDER == 928)
    {
        TDatabase::ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 9;
        return;
    }
    //===========================================



    OutPut("INTERNAL_POLYNOMIAL_DEGREE not defined for ANSATZ_ORDER "
           << TDatabase::ParamDB->ANSATZ_ORDER << endl);
    exit(1);
}

void CheckMaximumPrinciple(TSquareMatrix *A, double *sol, int N_Active,
                           double *errors)
{
  int i,j,k,index,ii,jj,kk,found;
  double tol=1e-8, minimum, maximum, val, maxerror;
  const int *RowPtr, *KCol;
  int N_DOF, maxprinciple = 0; 

  RowPtr = A->GetRowPtr();
  KCol = A->GetKCol();
  //double *Entries = A->GetEntries();
  N_DOF = A->GetN_Rows();
  maxerror = 0;

  j = RowPtr[0];
  for(i=0;i<N_DOF;i++)
  {
      k = RowPtr[i+1];
      val = sol[i];
      maximum = -1e10;
      minimum = 1e10;
      if (i<N_Active)
      {
          for(;j<k;j++)
          {
              index = KCol[j];
              if (index==i)
                  continue;
              if (sol[index]>maximum)
                  maximum = sol[index];
              if (sol[index]<minimum)
                  minimum = sol[index];
          }
      }
      else
          // Dirichlet nodes
      {
          jj = RowPtr[0];
          for (ii=0;ii<N_Active;ii++)
          {
              found = 0;
              kk = RowPtr[ii+1];
              for (;jj<kk;jj++)
              {
                  index = KCol[jj];
                  // off diagonal entry in row ii and column i
                  if (index==i)
                  {
                      found++;
                      if (sol[ii]>maximum)
                          maximum = sol[ii];
                      if (sol[ii]<minimum)
                          minimum = sol[ii];
                  }
                  // check for neighbour Dirichlet nodes
                  if (found)
                  {
                      jj = RowPtr[ii];
                      for (;jj<kk;jj++)
                      {
                          index = KCol[jj];
                          // off diagonal entry in row ii and column >= N_Active
                          if ((index>=N_Active)&&(index!=i))
                          {
                              if (sol[index]>maximum)
                                  maximum = sol[index];
                              if (sol[index]<minimum)
                                  minimum = sol[index];
                          }
                          
                      }
                  }
              }
          }
      }   
      // computation of maximum and minimum of neighbours done
      minimum = minimum - tol;
      maximum = maximum + tol;
      if ((minimum<=val) && (maximum>=val))
      {
          maxprinciple++;
          continue;
      }
      if ((maximum>-2) &&  (val-maximum > maxerror))
          maxerror = val-maximum;
      if ((minimum<2) && (minimum-val > maxerror))
          maxerror = minimum-val;
  } // endfor i
  OutPut("# dof " << N_DOF << " max principle fulfilled " << maxprinciple <<
         " " << maxprinciple*100.0/N_DOF << "% maxerror " << maxerror << endl);
  errors[0] = maxprinciple*100.0/N_DOF ;
  errors[1] = maxerror;
} // end CheckMaximumPrinciple

//save array of double pointers in a file
void SaveData(char *name, int N_Array, double **sol, int *N_Unknowns)
{
    int i;
    char OldString[] = "_old", MvString[] = "mv ", SpaceString[] = " ", tmp[200];
    
    // move the previous data set
    strcpy(tmp,MvString);
    strcat(tmp,name);
    strcat(tmp,SpaceString);
    strcat(tmp,name);
    strcat(tmp,OldString);
    OutPut(tmp<<endl);
    //system(tmp);

  std::ofstream dat(name);

  if(!dat)
  {
    OutPut("CANNOT OPEN FILE '" << name << "' FOR SAVING DATA!" << endl);
    return;
  }
  for (i=0;i<N_Array;i++)
  {
      OutPut(N_Unknowns[i] << " ");
   	  dat.write((char *)sol[i],sizeof(double)*N_Unknowns[i]);
  }  
  dat.close();
  
  OutPut("wrote output into file: " << name << endl);
}
void ReadData(char *name, int N_Array, double **sol, int *N_Unknowns)
{
    int i;
  std::ifstream dat(name);

  if(!dat)
  {
    OutPut("CANNOT OPEN FILE '" << name << "' FOR READING DATA!" << endl);
    exit(4711);
  }

  for (i=0;i<N_Array;i++)
  {
   	  dat.read((char *)sol[i],sizeof(double)*N_Unknowns[i]);
	  //SwapDoubleArray(sol[i], N_Unknowns[i]); 
  }  
  dat.close();
  
  OutPut("read input from file: " << name << endl);
}

// save sol into a file
void SaveData(std::string basename, double *sol, int nDOF)
{
  std::string filename = basename + "init";
  std::ofstream ofile;
  ofile.open( filename.c_str() , std::ios::out | std::ios::trunc );
  if( !ofile.is_open() )
  {
    OutPut("CANNOT OPEN FILE '" << filename << "' FOR SAVING DATA!" << endl);
    exit( 4711 );
  }
  ofile << setprecision( 12 );
  ofile.write((char *)sol,sizeof(double)*nDOF);
  ofile.close();
  OutPut("saving data into file: " << filename << endl);
}

void ReadData(std::string filename, double *sol, int nDOF)
{
  std::ifstream ifile(filename.c_str());
  if(!ifile)
  {
    OutPut("CANNOT OPEN FILE '" << filename << "' FOR READING DATA!" << endl);
    exit(4711);
  }
  ifile.read((char *)sol,sizeof(double)*nDOF);  
  ifile.close();
  OutPut("reading data from file: " << filename << endl);
}


/******************************************************************************/
//
// sets the nodes with global dof neum_to_diri to Dirichlet nodes 
// by setting the matrix and the rhs
//
/******************************************************************************/

#ifdef __2D__
void SetDirichletNodesFromNeumannNodes(TSquareMatrix2D **SQMATRICES, 
				       double *rhs, double *sol,
				       int N_neum_to_diri,
				       int *neum_to_diri,
				       int *neum_to_diri_bdry,
				       double *neum_to_diri_param,
				       BoundValueFunct2D *BoundaryValue)
#endif    
#ifdef __3D__
void SetDirichletNodesFromNeumannNodes(TSquareMatrix3D **SQMATRICES, 
				       double *rhs, double *sol,
				       int N_neum_to_diri,
				       int *neum_to_diri,
				       double *neum_to_diri_x,
				       double *neum_to_diri_y,
				       double *neum_to_diri_z,
				       BoundValueFunct3D *BoundaryValue)
#endif    

{
#ifdef __2D__
    TSquareMatrix2D *MatrixA;
#endif    
#ifdef __3D__
    TSquareMatrix3D *MatrixA;
#endif    
    double *Entries_A;
    int i, l, l0, l1, index;
    
    //OutPut(N_neum_to_diri << endl);
    if (N_neum_to_diri == 0)
	return;

    MatrixA = SQMATRICES[0];
    const int * RowPtr_A      = MatrixA->GetRowPtr();
    const int * KCol_A        = MatrixA->GetKCol();
    Entries_A     = MatrixA->GetEntries();
    // loop over dof to change
    for (i=0;i<N_neum_to_diri;i++)
    {
	index = neum_to_diri[i];
	l0 = RowPtr_A[index];
	l1 = RowPtr_A[index+1];
	for (l=l0;l<l1;l++)
	{
	    // diagonal entry
	    if (KCol_A[l]==index)
		Entries_A[l] = 1.;  
	    else
		Entries_A[l] = 0.;
	}
	// set boundary condition
#ifdef __2D__
	BoundaryValue(neum_to_diri_bdry[i], neum_to_diri_param[i], rhs[index]);
#endif    
#ifdef __3D__
	BoundaryValue(neum_to_diri_x[i], neum_to_diri_y[i],  neum_to_diri_z[i], rhs[index]);
	//OutPut("b "<< neum_to_diri_bdry[i]  << " " << neum_to_diri_param[i] << " " << rhs[index] << endl);
#endif    
	

        sol[index] = rhs[index];
    }
}


double graddiv_parameterOseen(double hK, double, double, double)
{
  double tau;
  
  switch(TDatabase::ParamDB->DIV_DIV_STAB_TYPE)
  {
    // constant
    case 0:
      tau = TDatabase::ParamDB->DIV_DIV_STAB_C1;
      break;
      // depending on mesh width 
    case 1:
      //  if (nu < hK)
        tau = TDatabase::ParamDB->DIV_DIV_STAB_C1 * pow(hK,TDatabase::ParamDB->DIV_DIV_STAB_C2);
      // else
      //  tau = TDatabase::ParamDB->DIV_DIV_STAB_C1;
      break;
    default:
      OutPut("DIV_DIV_STAB_TYPE " << TDatabase::ParamDB->DIV_DIV_STAB_TYPE << " not implemented" << endl);
      exit(4711);
  }
  return(tau);
}
