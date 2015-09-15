#ifdef __2D__

#include <FreeSurface2D.h>
#include <IsoBoundEdge.h>
#include <IsoInterfaceJoint.h>
#include <FEVectFunct2D.h>
#include <Database.h>
#include <FEDatabase2D.h>
#include <LinAlg.h>
#include <MooNMD_Io.h>
#include <QuadAffin.h>
#include <QuadBilinear.h>
#include <QuadIsoparametric.h>
#include <TriaAffin.h>
#include <TriaIsoparametric.h>
#include <InterfaceJoint.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>

void SolveGridEquation(double **Entries, double *sol, double *rhs,
                       int *KCol, int *RowPtr, int N_DOF)
{
  int i,j,k, col, Diognal;
  double *Entries11, *Entries12, *Entries21, *Entries22;
  double sum1, sum2, max_sum1, max_sum2;
  int start, end;

  double max_error, error=1.e-12;
  int iter;

  Entries11 = Entries[0];
  Entries12 = Entries[1];
  Entries21 = Entries[2];
  Entries22 = Entries[3];
  
  max_error = 1.; iter = 0;
  while(max_error>error)
  {
    max_error = 0.0; iter++;
    for(i=0;i<N_DOF;i++)
    {
      start = RowPtr[i];
      end = RowPtr[i+1];
      sum1 = rhs[i];
      sum2 = rhs[i+N_DOF];
      for(k=start;k<end;k++)
      {
        col = KCol[k];
        if (col==i) Diognal = k;
        sum1 -= Entries11[k] * sol[col]
              +Entries12[k] * sol[col+N_DOF];
        sum2 -= Entries21[k] * sol[col]
              +Entries22[k] * sol[col+N_DOF];
      } // endfor k
     // sol[i] += sum1/Entries11[start];
     // sol[i+N_DOF] += sum2/Entries22[start];
        sol[i] += sum1/Entries11[Diognal];
        sol[i+N_DOF] += sum2/Entries22[Diognal];
      if(max_error<fabs(sum1/Entries11[Diognal])) max_error = fabs(sum1/Entries11[Diognal]);
      if(max_error<fabs(sum2/Entries22[Diognal])) max_error = fabs(sum2/Entries22[Diognal]);
    } // endfor i
    if(iter == 1000) break;
  } // end while
// OutPut("Grid Solver: Number iteration "<<iter<<endl);
// exit(0);

}


void Solver_3dia(int N_Splines, double *a, double *b, double *c, double *rhs, double *sol)
{
  double *alpha, *beta, *y;
  int i, N;

  N = N_Splines+1;
  alpha = new double[N]; beta = new double[N]; y = new double[N];

  alpha[0] = a[0]; y[0] = rhs[0];
  for(i=1;i<N;i++)
  {
    beta[i] = b[i]/alpha[i-1];
    alpha[i] = a[i]-beta[i]*c[i-1];
    y[i] = rhs[i]-beta[i]*y[i-1];
  }

  sol[N-1] = y[N-1]/alpha[N-1];
  for(i=N-2;i>=0;i--)
    sol[i] = (y[i]-c[i]*sol[i+1])/alpha[i];

  delete [] alpha; delete [] beta; delete [] y;
}



// ====================================================================
// modify matrices due to integrals on free surface
// ====================================================================
void FreeSurfInt(TSquareMatrix2D *A11, TSquareMatrix2D *A12,
                 TSquareMatrix2D *A21, TSquareMatrix2D *A22,
                 double *rhs1, double *rhs2,
                 BoundCondFunct2D *BoundaryCondition,
                 double dt, double factor)
{
  int i,j,k,l,m;
  TBaseCell *cell;
  TCollection *Coll;
  int N_Cells, N_Vertices, N_Edges;
  TJoint *joint;
  TIsoBoundEdge *isoboundedge;
  TBoundComp *BoundComp;
  int comp;
  double t0, t1, n0, n1, normn;
  BoundCond Cond0, Cond1;
  int JointNumbers[MAXN_JOINTS], IJoint, N_IsoJoints;
  FE2D FEId;
  TFE2D *ele;
  RefTrans2D RefTrans;
  TRefTrans2D *F_K;
  const TFESpace2D *fespace;
  BF2DRefElements RefElement;
  QuadFormula1D LineQuadFormula;
  TQuadFormula1D *qf1;
  int N_LinePoints;
  double *LineWeights, *zeta;
  double x0, y0, x1, y1;
  int N_BaseFunct, *N_BaseFuncts;
  double **uref, **uxiref, **uetaref;
  double uorig[MaxN_BaseFunctions2D], uxorig[MaxN_BaseFunctions2D];
  double uyorig[MaxN_BaseFunctions2D];
  BaseFunct2D *BaseFuncts;
  double r2, r;
  int *KCol, *RowPtr;
  double *ValuesA11, *ValuesA12, *ValuesA21, *ValuesA22;
  int *BeginIndex, *GlobalNumbers, *DOF, TestDOF, AnsatzDOF;
  int index1, index2;
  double val;

  double Ca = TDatabase::ParamDB->P9;

  BaseFuncts = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D();
  N_BaseFuncts = TFEDatabase2D::GetN_BaseFunctFromFE2D();

  fespace = A11->GetFESpace();
  Coll = fespace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  BeginIndex = fespace->GetBeginIndex();
  GlobalNumbers = fespace->GetGlobalNumbers();

  RowPtr = A11->GetRowPtr();
  KCol = A11->GetKCol();

  ValuesA11 = A11->GetEntries();
  ValuesA12 = A12->GetEntries();
  ValuesA21 = A21->GetEntries();
  ValuesA22 = A22->GetEntries();

  for(i=0;i<N_Cells;i++)
  {
    // cout << endl << "CELL number: " << i << endl;
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    IJoint = 0;
    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge) 
      {
        isoboundedge = (TIsoBoundEdge *)joint;
        BoundComp = isoboundedge->GetBoundComp();
        isoboundedge->GetParameters(t0, t1);
        comp=BoundComp->GetID();
        BoundaryCondition(comp, t0, Cond0);
        BoundaryCondition(comp, t1, Cond1);

        if(Cond0 == FREESURF)
        {
          JointNumbers[IJoint] = j;
          IJoint++;
        }
      } // endif
    } // endfor j

    N_IsoJoints = IJoint;
    if(N_IsoJoints > 0)
    {
      // cout << "Cell " << i << " has free surface." << endl;
      FEId = fespace->GetFE2D(i, cell);
      ele = TFEDatabase2D::GetFE2D(FEId);
      RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TQuadIsoparametric *)F_K)->SetCell(cell);
        break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
          ((TTriaIsoparametric *)F_K)->SetCell(cell);
        break;
      } // endswitch

      l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
      LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
      qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
      qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
      TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                  ->MakeRefElementData(LineQuadFormula);

      DOF = GlobalNumbers + BeginIndex[i];
      N_BaseFunct = N_BaseFuncts[FEId];

      for(j=0;j<N_IsoJoints;j++)
      {
        IJoint = JointNumbers[j];
        // cout << "joint number: " << IJoint << endl;
        cell->GetVertex(IJoint)->GetCoords(x0, y0);
        cell->GetVertex((IJoint+1) % N_Edges)->GetCoords(x1, y1);

        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetTangent(IJoint, zeta[k], t0, t1);
          normn = sqrt(t0*t0+t1*t1);
          n0 =  t1/normn;
          n1 = -t0/normn;
          // cout << "zeta: " << zeta[k] << endl;
          // cout << "k= " << k << "  tangent: " << t0 << " " << t1 << endl;
          // cout << "length: " << sqrt(t0*t0+t1*t1) << endl;
          uref = TFEDatabase2D::GetJointValues2D(BaseFuncts[FEId], 
                        LineQuadFormula, IJoint);
          uxiref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], 
                        LineQuadFormula, IJoint, D10);
          uetaref = TFEDatabase2D::GetJointDerivatives2D(BaseFuncts[FEId], 
                        LineQuadFormula, IJoint, D01);
          switch(RefElement)
          {
            case BFUnitSquare:
              ((TQuadIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
            break;

            case BFUnitTriangle:
              ((TTriaIsoparametric *)F_K)->GetOrigValues(IJoint,  zeta[k],
                        N_BaseFunct, uref[k], uxiref[k], uetaref[k],
                        uorig, uxorig, uyorig);
            break;
          } // endswitch

          // modify matrices
          r2 = 1/(t0*t0+t1*t1);
          r = dt * sqrt(t0*t0+t1*t1)/Ca;
          for(l=0;l<N_BaseFunct;l++)
          {
            TestDOF = DOF[l];

            // updating rhs
            val = -r2* t0 * (uxorig[l]*t0+uyorig[l]*t1);
            val *= LineWeights[k]*r;
            // cout << "Rhs1: " << TestDOF << " " << val;
            rhs1[TestDOF] += val;
            // cout << " = " << rhs1[TestDOF] << endl;

            val = -r2* t1 * (uxorig[l]*t0+uyorig[l]*t1);
            val *= LineWeights[k]*r;
            // cout << "Rhs2: " << TestDOF << " " << val;
            rhs2[TestDOF] += val;
            // cout << " = " << rhs2[TestDOF] << endl;

            val = sqrt(t0*t0+t1*t1) * LineWeights[k] *
                        n1*n1 * factor * uorig[l]*dt;
            // cout << "l= " << l << " bf: " << uorig[l] << endl;
            rhs1[TestDOF] += val*n0;
            rhs2[TestDOF] += val*n1;

            index2 = RowPtr[TestDOF+1];
            for(m=0;m<N_BaseFunct;m++)
            {
              AnsatzDOF = DOF[m];
              // cout << AnsatzDOF << " -- " << TestDOF << endl;
              index1 = RowPtr[TestDOF];
              if(index1+1 == index2) continue;
              while(KCol[index1] != AnsatzDOF) index1++;

              val = r2*dt*( 
                        (uxorig[m]*t0 + uyorig[m]*t1)  *
                        (uxorig[l]*t0 + uyorig[l]*t1) );
              val *= LineWeights[k]*r;
              // cout << "A11: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA11[index1] += val;

              val = r2*dt*( 
                        (uxorig[m]*t0 + uyorig[m]*t1)  *
                        (uxorig[l]*t0 + uyorig[l]*t1) );
              val *= LineWeights[k]*r;
              // cout << "A22: " << TestDOF << " ";
              // cout << AnsatzDOF << " " << val << endl;
              ValuesA22[index1] += val;
            } // endfor m
          } // endfor l
        } // endfor k
      } // endfor j
    } // end (N_IsoJoints > 0)
    else
    {
      // cout << "Cell " << i << " has NO free surface." << endl;
    }
  } // endfor i
}

// ====================================================================
// determine grid velocity in whole domain
// ====================================================================
void GetGridVelocity(TMultiGrid2D *GridMG, TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *AuxGridPos,
                     double *Nx, double *Ny,
                     TFEVectFunct2D *Velocity, double dt,
                     TFEVectFunct2D *GridVelocity)
{
  int i,j,k,l,m;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels;
  TMGLevel2D *Level;
  int *DOF, *JointDOF;
  const TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  bool OnBoundary;
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y;
  double x0, x1, y0, y1;
  int GridLength;
  double *Rhs, *Sol;
  int N_BoundaryNodes;
  double res, oldres;
  TJoint *joint;
  TVertex **Vertices;
  double *gridvelo;
  double *LineWeights, *zeta;
  int N_LinePoints;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  double normalx, normaly, nx, ny, tx, ty, hE;
  int IIso;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  int N_Inner, N_;
  double un;

  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
  // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

  NewValuesX = AuxGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // determine outer normal vectors
  IIso = N_BoundaryNodes;
  for(i=0;i<N_Cells;i++)
  {
    // cout << "cell: " << i << endl;
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        // cout << "joint: " << j << endl;
        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t*LineWeights[k]*nx;
          normaly += t*LineWeights[k]*ny;
          hE += t*LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k 

        DOF = GridGlobalNumbers + GridBeginIndex[i];

/*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
        
              case 3:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
*/
// /*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
        
              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
// */

/* 
        // not needed
        if(cell->GetJoint(j)->GetType() == IsoBoundEdge)
        {
          FEId = VelocitySpace->GetFE2D(i, cell);
          FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
          N_LocalDOFs = FEDesc->GetN_JointDOF();
          for(k=1;k<N_LocalDOFs-1;k++)
          {
            Nx[IIso] += normalx;
            Ny[IIso] += normaly;
            IIso++;
          } // endfor
        } // endif
*/
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;
  }

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = false;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
        OnBoundary = true;
    } // endfor j

    if(OnBoundary)
    {
      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);
  
      DOF = VeloGlobalNumbers + VeloBeginIndex[i];
  
      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];
        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();      
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          X[0] = X[1] = X[2] = 0;
          Y[0] = Y[1] = Y[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOF !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          X[0] = X[1] = X[2] = X[3] = 0;
          Y[0] = Y[1] = Y[2] = Y[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      DOF = GridGlobalNumbers + GridBeginIndex[i];
  
      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if(TDatabase::ParamDB->P5 > 0)
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
            NewValuesX[l] = ValuesX[l] + dt*VX[j];
            NewValuesY[l] = ValuesY[l] + dt*VY[j];
          }
        }
      } // endfor j
    } // endif
  } // endfor i

  N_Levels = GridMG->GetN_Levels();
  Level = GridMG->GetLevel(N_Levels-1);
  Rhs = Level->GetRhs();
  Sol = Level->GetSolution();
  
  memset(Rhs, 0, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Rhs + (GridLength-N_BoundaryNodes), 
         NewValuesX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  memcpy(Sol, NewValuesX, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), 
         NewValuesX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  res = 1;
  j=0;
  while(res>TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR && j<100)
  {
    GridMG->Cycle(N_Levels-1, res);
    if(TDatabase::ParamDB->SC_VERBOSE > 1)
    {
      OutPut("residual after iteration " << j << ": " << res);
      if(j>0)
      {
        OutPut(" rate: " << res/oldres << endl);
      }
      else
      {
        OutPut(endl);
      }
    } // SC_VERBOSE

    oldres = res;
    j++;
  }

  memcpy(NewValuesX, Sol, GridLength*SizeOfDouble);

  memset(Rhs, 0, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Rhs + (GridLength-N_BoundaryNodes), 
         NewValuesY+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  memcpy(Sol, NewValuesY, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), 
         NewValuesY+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  res = 1;
  j=0;
  while(res>TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR && j<100)
  {
    GridMG->Cycle(N_Levels-1, res);
    if(TDatabase::ParamDB->SC_VERBOSE > 1)
    {
      OutPut("residual after iteration " << j << ": " << res);
      if(j>0)
      {
        OutPut(" rate: " << res/oldres << endl);
      }
      else
      {
        OutPut(endl);
      }
    } // SC_VERBOSE

    oldres = res;
    j++;
  }

  memcpy(NewValuesY, Sol, GridLength*SizeOfDouble);

  // put solution into grid position
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    DOF = GridGlobalNumbers + GridBeginIndex[i];
    N_Vertices = cell->GetN_Vertices();

    switch(N_Vertices)
    {
      case 3:
        for(j=0;j<N_Vertices;j++)
        {
          k = DOF[j];
          cell->GetVertex(j)->SetCoords(NewValuesX[k], NewValuesY[k]);
        }
      break;
      
      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewValuesX[k], NewValuesY[k]);
      break;
    } // endswitch
  } // endfor i

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, NewValuesX, 2*GridLength*SizeOfDouble);
  Daxpy(2*GridLength, -1, ValuesX, gridvelo);
  Dscal(2*GridLength, 1/dt, gridvelo);

  // for(i=0;i<GridLength;i++)
  //   cout << gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;

} // GetGridVelocity

// ====================================================================
// determine new grid position
// ====================================================================
void MoveGrid(TMultiGrid2D *GridMG, TFEVectFunct2D *GridPos,
              double *IsoX, double *IsoY,
              double *Nx, double *Ny,
              TFEVectFunct2D *Velocity, double dt,
              TFEVectFunct2D *NewGridPos, 
              double *NewIsoX, double *NewIsoY,
              TFEVectFunct2D *GridVelocity)
{
  int i,j,k,l,m;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels;
  TMGLevel2D *Level;
  int *DOF, *JointDOF;
  const TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll;
  TBaseCell *cell;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  bool OnBoundary;
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y;
  double x0, x1, y0, y1;
  int GridLength;
  double *Rhs, *Sol;
  int N_BoundaryNodes;
  double res, oldres;
  TJoint *joint;
  TIsoBoundEdge *isojoint;
  TVertex **Vertices;
  int IIso;
  double *gridvelo;
  double *LineWeights, *zeta;
  int N_LinePoints;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  double normalx, normaly, nx, ny, tx, ty, hE;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  int N_Inner, N_;
  double un;
  int polydegree;
  QuadFormula2D QuadFormula;

  if(TDatabase::ParamDB->P5>0)
  {
    OutPut("Boundary update using (u.n)n" << endl);
  }
  else
  {
    OutPut("Boundary update using u" << endl);
  }

  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();

  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - N_Inner;
  // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

  NewValuesX = NewGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

  Coll = VelocitySpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // determine outer normal vectors
  IIso = N_BoundaryNodes;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();

    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
      {
        FEId = VelocitySpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(2*polydegree-1);
            ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
            ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t * LineWeights[k] * nx;
          normaly += t * LineWeights[k] * ny;
          hE += t * LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k 

        DOF = GridGlobalNumbers + GridBeginIndex[i];
/*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
        
              case 3:
                l = DOF[2] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
                l = DOF[0] - N_Inner;
                if(l>=0)
                {
                  Nx[l] += normalx;
                  Ny[l] += normaly;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
*/

// /*
        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
        
              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
// */

        if(cell->GetJoint(j)->GetType() == IsoBoundEdge)
        {
          FEId = VelocitySpace->GetFE2D(i, cell);
          FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
          N_LocalDOFs = FEDesc->GetN_JointDOF();
          t = 2.0/(N_LocalDOFs-1);
          for(k=1;k<N_LocalDOFs-1;k++)
          {
            /*
            Nx[IIso] += normalx;
            Ny[IIso] += normaly;
            */
            // /*
            s = -1.0 + k*t;
            F_K->GetOuterNormal(j, s, nx, ny);
            Nx[IIso] += nx;
            Ny[IIso] += ny;
            // */
            IIso++;
          } // endfor
        } // endif
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;

    // cout << setw(5) << i << "n = (" << Nx[i] << ", " << Ny[i] << ")" << endl;
  }

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = false;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) )
        OnBoundary = true;
    } // endfor j

    if(OnBoundary)
    {
      FEId = VelocitySpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);
  
      DOF = VeloGlobalNumbers + VeloBeginIndex[i];
  
      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];
        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();      
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      DOF = GridGlobalNumbers + GridBeginIndex[i];
  
      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if(TDatabase::ParamDB->P5 > 0)
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
            NewValuesX[l] = ValuesX[l] + dt*VX[j];
            NewValuesY[l] = ValuesY[l] + dt*VY[j];
          }
        }
      } // endfor j
    } // endif
  } // endfor i

  N_Levels = GridMG->GetN_Levels();
  Level = GridMG->GetLevel(N_Levels-1);
  Rhs = Level->GetRhs();
  Sol = Level->GetSolution();
  
  memset(Rhs, 0, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Rhs + (GridLength-N_BoundaryNodes), 
         NewValuesX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  memcpy(Sol, NewValuesX, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), 
         NewValuesX+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  res = 1;
  j=0;
  while(res>TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR && j<100)
  {
    GridMG->Cycle(N_Levels-1, res);
    if(TDatabase::ParamDB->SC_VERBOSE > 1)
    {
      OutPut("residual after iteration " << j << ": " << res);
      if(j>0)
      {
        OutPut(" rate: " << res/oldres << endl);
      }
      else
      {
        OutPut(endl);
      }
    } // SC_VERBOSE

    oldres = res;
    j++;
  }

  memcpy(NewValuesX, Sol, GridLength*SizeOfDouble);

  memset(Rhs, 0, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Rhs + (GridLength-N_BoundaryNodes), 
         NewValuesY+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  memcpy(Sol, NewValuesY, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes), 
         NewValuesY+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  res = 1;
  j=0;
  while(res>TDatabase::ParamDB->SC_LIN_RES_NORM_MIN_SCALAR && j<100)
  {
    GridMG->Cycle(N_Levels-1, res);
    if(TDatabase::ParamDB->SC_VERBOSE > 1)
    {
      OutPut("residual after iteration " << j << ": " << res);
      if(j>0)
      {
        OutPut(" rate: " << res/oldres << endl);
      }
      else
      {
        OutPut(endl);
      }
    } // SC_VERBOSE
    oldres = res;
    j++;
  }

  memcpy(NewValuesY, Sol, GridLength*SizeOfDouble);

  // put solution into grid position
  IIso = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    DOF = GridGlobalNumbers + GridBeginIndex[i];
    N_Vertices = cell->GetN_Vertices();

    switch(N_Vertices)
    {
      case 3:
        for(j=0;j<N_Vertices;j++)
        {
          k = DOF[j];
          cell->GetVertex(j)->SetCoords(NewValuesX[k], NewValuesY[k]);
        }
      break;
      
      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewValuesX[k], NewValuesY[k]);
      break;
    } // endswitch

    N_Edges = cell->GetN_Edges();
    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        isojoint = (TIsoBoundEdge *)joint;
        k = isojoint->GetN_Vertices();
        Vertices = isojoint->GetVertices();
        FEId = VelocitySpace->GetFE2D(i, cell);
        FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        m = FEDesc->GetN_JointDOF();
        if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[i];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX[IIso], IsoY[IIso]);
            if(TDatabase::ParamDB->P5 > 0)
            {
              un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes] 
                  + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
              IsoX[IIso] += dt*un*Nx[IIso+N_BoundaryNodes];
              IsoY[IIso] += dt*un*Ny[IIso+N_BoundaryNodes];
              // cout << "U:   " << ValuesVX[m] << " " << ValuesVY[m] << endl;
              // cout << "N:   " << Nx[IIso+N_BoundaryNodes] << " "
              //                 << Ny[IIso+N_BoundaryNodes] << endl;
              // cout << "UNN: " << un*Nx[IIso+N_BoundaryNodes] << " " 
              //                 << un*Ny[IIso+N_BoundaryNodes] << endl;
            }
            else
            {
              IsoX[IIso] += dt*ValuesVX[m];
              IsoY[IIso] += dt*ValuesVY[m];
            }
            Vertices[l]->SetCoords(IsoX[IIso], IsoY[IIso]);
            IIso++;
          } // endfor l
        }
        else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        }
      } // endif
    } // endfor j
  } // endfor i

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, NewValuesX, 2*GridLength*SizeOfDouble);
  Daxpy(2*GridLength, -1, ValuesX, gridvelo);
  Dscal(2*GridLength, 1/dt, gridvelo);

} // MoveGrid


void GetGridVelocity(double **Entries, double *Sol, double *Rhs,
                     int *KCol, int *RowPtr,
                     TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *AuxGridPos,
                     TFEVectFunct2D *Velocity, double dt,
                     TFEVectFunct2D *GridVelocity, int *Velo_CellNo)
{
  int i,j,k,l,m;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels;
  TMGLevel2D *Level;
  int *DOF, *JointDOF;
  const TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll, *Velo_Coll;
  TBaseCell *cell;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  bool OnBoundary;
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY, *d;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y;
  double x0, x1, y0, y1;
  int GridLength;
  int N_BoundaryNodes;
  double res, oldres;
  TJoint *joint;
  TVertex **Vertices;
  double *gridvelo, *Nx, *Ny;
  double *LineWeights, *zeta;
  int N_LinePoints, Phase_No;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  double normalx, normaly, nx, ny, tx, ty;
  int IIso, Velo_N_Cells;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  int N_Inner, N_, Velo_i;
  double un, hE;

  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();
  Velo_Coll = VelocitySpace->GetCollection();
  Velo_N_Cells = Velo_Coll->GetN_Cells();


  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
//   cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

//   cout << GridLength << " N_BoundaryNodes:" << endl;
// 
// exit(0);

  d = new double[2*GridLength];
  Nx = new double[N_BoundaryNodes];
  Ny = new double[N_BoundaryNodes];

  memset(Nx, 0, N_BoundaryNodes*SizeOfDouble);
  memset(Ny, 0, N_BoundaryNodes*SizeOfDouble);

//   cout << GridLength << " N_BoundaryNodes:" << endl;
  NewValuesX = AuxGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

//   Coll = VelocitySpace->GetCollection();

  Coll = GridSpace->GetCollection();
  N_Cells = Coll->GetN_Cells();


//   for(i=0;i<N_Cells;i++)
//     {
// cout << "cell : " << i << " number " << Velo_CellNo[i] <<endl;
// }

  // determine outer normal vectors
  IIso = N_BoundaryNodes;
 if(TDatabase::ParamDB->P5 > 0)
  {
  // determine outer normal vectors
   for(i=0;i<N_Cells;i++)
    {
//     cout << "cell: " << i << endl;
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
    N_Edges = cell->GetN_Edges();


    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) 
            || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) 
             || (cell->GetJoint(j)->GetType() == InterfaceJoint) )
      {
        // cout << "joint: " << j << endl;
        Phase_No = cell->GetPhase_ID();
        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;

//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
//        
        Velo_i = cell->GetGlobalCellNo();

        FEId = VelocitySpace->GetFE2D(Velo_i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t*LineWeights[k]*nx;
          normaly += t*LineWeights[k]*ny;
          hE += t*LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];

        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
  if(Phase_No==1)  // outer phase has inward normal
   {
    Nx[i] /= -t;
    Ny[i] /= -t;
    }
   else
    {
    Nx[i] /= t;
    Ny[i] /= t;
    }

  }
 }

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = false;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) 
         || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) )
        OnBoundary = true;
    } // endfor j

    if(OnBoundary)
    {
//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
      Velo_i = cell->GetGlobalCellNo();

      FEId = VelocitySpace->GetFE2D(Velo_i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);
  
      DOF = VeloGlobalNumbers + VeloBeginIndex[Velo_i];
  
      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];
        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          X[0] = X[1] = X[2] = 0;
          Y[0] = Y[1] = Y[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOF !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          X[0] = X[1] = X[2] = X[3] = 0;
          Y[0] = Y[1] = Y[2] = Y[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      DOF = GridGlobalNumbers + GridBeginIndex[i];
  
      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
         if(TDatabase::ParamDB->P5 > 0)
         {
          un = VX[j]*Nx[k] + VY[j]*Ny[k];
          cout << l <<"  ---  "<< VX[j] << "  ---  " << VY[j] << endl;
          cout << "update normal direction is not good in rising bubble" << endl;
          NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
          NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
         }
         else
         {
//           if(fabs(ValuesX[l])<1e-8 && fabs(VX[j])>1e-16)   
//            cout << " x " << ValuesX[l] <<  " y " << ValuesY[l] << " VX " << VX[j] << endl;

          if(ValuesX[l]!=0)
           NewValuesX[l] = ValuesX[l] + dt*VX[j]; //axial bd no need x update
          NewValuesY[l] = ValuesY[l] + dt*VY[j];
         }
        }
      } // endfor j
    } // endif
  } // endfor i


  memset(Rhs, 0, 2*GridLength*SizeOfDouble);

   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes),
          d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes),
          d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes),
         d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes),
         d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  SolveGridEquation(Entries, Sol, Rhs, KCol, RowPtr, GridLength);

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, Sol, 2*GridLength*SizeOfDouble);
  Dscal(2*GridLength, 1./dt, gridvelo);


//    for(i=0;i<GridLength;i++)
//      cout << gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;
// exit(0);

  delete [] d;
  delete [] Nx;
  delete [] Ny;

//    for(i=0;i<GridLength;i++)
//      cout << gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;
// exit(0);

} // GetGridVelocity

void GetGridVelo_outer(double **Entries, double *Sol, double *Rhs,
                       int *KCol, int *RowPtr,
                       TFEVectFunct2D *GridPos,
                       TFEVectFunct2D *AuxGridPos,
                       TFEVectFunct2D *Velocity, double dt,
                       TFEVectFunct2D *GridVelocity, int *Velo_CellNo)
{
  int i,j,k,l,m;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int N_Levels;
  TMGLevel2D *Level;
  int *DOF, *JointDOF;
  const TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll, *Velo_Coll;
  TBaseCell *cell;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  bool OnBoundary;
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY, *d;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y;
  double x0, x1, y0, y1;
  int GridLength;
  int N_BoundaryNodes;
  double res, oldres;
  TJoint *joint;
  TVertex **Vertices;
  double *gridvelo, *Nx, *Ny;
  double *LineWeights, *zeta;
  int N_LinePoints, Phase_No;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  double normalx, normaly, nx, ny, tx, ty;
  int IIso, Velo_N_Cells;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  int N_Inner, N_, Velo_i;
  double un, hE;

  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();
  Velo_Coll = VelocitySpace->GetCollection();
  Velo_N_Cells = Velo_Coll->GetN_Cells();


  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - GridSpace->GetN_Inner();
//   cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;

//   cout << GridLength << " N_BoundaryNodes:" << endl;


  d = new double[2*GridLength];
  Nx = new double[N_BoundaryNodes];
  Ny = new double[N_BoundaryNodes];

  memset(Nx, 0, N_BoundaryNodes*SizeOfDouble);
  memset(Ny, 0, N_BoundaryNodes*SizeOfDouble);

//   cout << GridLength << " N_BoundaryNodes:" << endl;
  NewValuesX = AuxGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);

//   Coll = VelocitySpace->GetCollection();

  Coll = GridSpace->GetCollection();
  N_Cells = Coll->GetN_Cells();


//   for(i=0;i<N_Cells;i++)
//     {
// cout << "cell : " << i << " number " << Velo_CellNo[i] <<endl;
// }

  // determine outer normal vectors
  IIso = N_BoundaryNodes;
 if(TDatabase::ParamDB->P5 > 0)
  {
  // determine outer normal vectors
   for(i=0;i<N_Cells;i++)
    {
//     cout << "cell: " << i << endl;
    cell = Coll->GetCell(i);
    cell->SetClipBoard(i);
    N_Edges = cell->GetN_Edges();


    for(j=0;j<N_Edges;j++)
    {
     if( !(cell->GetJoint(j)->InnerJoint()) 
        || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) 
        || (cell->GetJoint(j)->GetType() == InterfaceJoint) )
      {
        // cout << "joint: " << j << endl;
        Phase_No = cell->GetPhase_ID();
        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;

//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
        Velo_i = cell->GetGlobalCellNo();
        FEId = VelocitySpace->GetFE2D(Velo_i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t*LineWeights[k]*nx;
          normaly += t*LineWeights[k]*ny;
          hE += t*LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];

        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
  if(Phase_No==1)  // outer phase has inward normal
   {
    Nx[i] /= -t;
    Ny[i] /= -t;
    }
   else
    {
    Nx[i] /= t;
    Ny[i] /= t;
    }

  }
 }

  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell  = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = false;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) 
        || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) 
        || (cell->GetJoint(j)->GetType() == InterfaceJoint) )
        OnBoundary = true;
    } // endfor j

    if(OnBoundary)
    {
//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
      Velo_i = cell->GetGlobalCellNo();
      FEId = VelocitySpace->GetFE2D(Velo_i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);
  
      DOF = VeloGlobalNumbers + VeloBeginIndex[Velo_i];
  
      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];
        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          X[0] = X[1] = X[2] = 0;
          Y[0] = Y[1] = Y[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOF !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          X[0] = X[1] = X[2] = X[3] = 0;
          Y[0] = Y[1] = Y[2] = Y[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      DOF = GridGlobalNumbers + GridBeginIndex[i];
  
      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
         if(TDatabase::ParamDB->P5 > 0)
         {
          un = VX[j]*Nx[k] + VY[j]*Ny[k];
          cout << l <<"  ---  "<< VX[j] << "  ---  " << VY[j] << endl;
          cout << "update normal direction is not good in raising bubble" << endl;
          NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
          NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
         }
         else
         {
/*          if(fabs(ValuesX[l])<1e-8 && fabs(VX[j])>1e-16)   
           cout << "GetGridVelo_outer x " << ValuesX[l] <<  " y " << ValuesY[l] << " VX " << VX[j] << endl;*/  
          if(ValuesX[l]!=0) 
           NewValuesX[l] = ValuesX[l] + dt*VX[j];
          NewValuesY[l] = ValuesY[l] + dt*VY[j];
         }
        }
      } // endfor j
    } // endif
  } // endfor i


  memset(Rhs, 0, 2*GridLength*SizeOfDouble);

   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes),
          d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes),
          d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes),
         d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes),
         d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  SolveGridEquation(Entries, Sol, Rhs,
                    KCol, RowPtr, GridLength);

  gridvelo = GridVelocity->GetValues();
  memcpy(gridvelo, Sol, 2*GridLength*SizeOfDouble);
  Dscal(2*GridLength, 1/dt, gridvelo);


//    for(i=0;i<GridLength;i++)
//      cout << gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;
// exit(0);

  delete [] d;
  delete [] Nx;
  delete [] Ny;

//    for(i=0;i<GridLength;i++)
//      cout << gridvelo[i] << "  ---  " << gridvelo[i+GridLength] << endl;
// exit(0);

} // GetGridVelocity




void MoveGrid_2Phase(double **Entries, double *Sol, double *Rhs,
                     int *KCol, int *RowPtr,
                     TFEVectFunct2D *GridPos,
                     TFEVectFunct2D *NewGridPos,
                     TFEVectFunct2D *Velocity, double dt,  
                     int *Velo_CellNo, int isoupdate)
{
  int i,j,k,l,m;
  int *VeloBeginIndex, *VeloGlobalNumbers;
  int *GridBeginIndex, *GridGlobalNumbers;
  int N_Cells, N_Vertices, N_Edges, N_LocalDOFs;
  int *DOF, *JointDOF;
  const TFESpace2D *VelocitySpace, *GridSpace;
  TCollection *Coll, *Velo_Coll;
  TBaseCell *cell;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  bool OnBoundary;
  double xi[4], eta[4], X[4], Y[4], VX[4], VY[4];
  double FunctValues[4][MaxN_BaseFunctions2D];
  double FEValuesX[MaxN_BaseFunctions2D];
  double FEValuesY[MaxN_BaseFunctions2D];
  double *ValuesX, *ValuesY, *d;
  double *ValuesVX, *ValuesVY;
  double *NewValuesX, *NewValuesY;
  double s, t, x, y;
  double x0, x1, y0, y1;
  int GridLength;
  int N_BoundaryNodes;
  double res, oldres;
  TJoint *joint;
  TIsoInterfaceJoint *isojoint;
  TVertex **Vertices;
  int IIso;
  double *gridvelo, *Nx, *Ny;
  double *LineWeights, *zeta;
  int N_LinePoints;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  double normalx, normaly, nx, ny, tx, ty, hE;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  int N_Inner, N_, Velo_N_Cells;
  double un;
  int polydegree, Velo_i;
  QuadFormula2D QuadFormula;
  double IsoX, IsoY;


  VelocitySpace = Velocity->GetFESpace2D();
  VeloBeginIndex = VelocitySpace->GetBeginIndex();
  VeloGlobalNumbers = VelocitySpace->GetGlobalNumbers();
  ValuesVX = Velocity->GetValues();
  ValuesVY = ValuesVX + Velocity->GetLength();
  Velo_Coll = VelocitySpace->GetCollection();
  Velo_N_Cells =  Velo_Coll->GetN_Cells();

  GridSpace = GridPos->GetFESpace2D();
  GridBeginIndex = GridSpace->GetBeginIndex();
  GridGlobalNumbers = GridSpace->GetGlobalNumbers();
  GridLength = GridPos->GetLength();
  ValuesX = GridPos->GetValues();
  ValuesY = ValuesX + GridLength;

  N_Inner = GridSpace->GetN_Inner();
  N_BoundaryNodes = GridLength - N_Inner;
  // cout << "N_BoundaryNodes: " << N_BoundaryNodes << endl;
  d = new double[ 2*GridLength];

  NewValuesX = NewGridPos->GetValues();
  NewValuesY = NewValuesX + GridLength;

  memcpy(NewValuesX, ValuesX, 2*GridLength*SizeOfDouble);
  Nx = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
  Ny = new double[2*N_BoundaryNodes]; // additional values for edge midpoints
//   Coll = VelocitySpace->GetCollection();

  Coll = GridSpace->GetCollection();
  N_Cells = Coll->GetN_Cells();

  // determine outer normal vectors
  IIso = N_BoundaryNodes;
 if(TDatabase::ParamDB->P5>0)
  {
   for(i=0;i<N_Cells;i++)
    {
     cell = Coll->GetCell(i);
     N_Edges = cell->GetN_Edges();

     for(j=0;j<N_Edges;j++)
      {
       if( !(cell->GetJoint(j)->InnerJoint())
            || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) )
       {

//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell == Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
        Velo_i = cell->GetGlobalCellNo();
        FEId = VelocitySpace->GetFE2D(Velo_i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(2*polydegree-1);
            ((TTriaIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TTriaIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
            QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
            ((TQuadIsoparametric *)F_K)->SetApproximationOrder(polydegree);
            ((TQuadIsoparametric *)F_K)->SetQuadFormula(QuadFormula);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch

        normalx = 0;
        normaly = 0;
        hE = 0;
        for(k=0;k<N_LinePoints;k++)
        {
          F_K->GetOuterNormal(j, zeta[k], nx, ny);
          F_K->GetTangent(j, zeta[k], tx, ty);
          t = sqrt(tx*tx+ty*ty);
          normalx += t * LineWeights[k] * nx;
          normaly += t * LineWeights[k] * ny;
          hE += t * LineWeights[k];
          // cout << "k= " << k << " " << nx << " " << ny << endl;
        } // endfor k

        DOF = GridGlobalNumbers + GridBeginIndex[i];

        switch(N_Edges)
        {
          case 3:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;

          case 4:
            switch(j)
            {
              case 0:
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 1:
                l = DOF[1] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 2:
                l = DOF[3] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;

              case 3:
                l = DOF[2] - N_Inner;
                F_K->GetOuterNormal(j, -1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
                l = DOF[0] - N_Inner;
                F_K->GetOuterNormal(j, 1, nx, ny);
                if(l>=0)
                {
                  Nx[l] += nx*hE;
                  Ny[l] += ny*hE;
                } // endif
              break;
            } // endswitch j
          break;
        } // endswitch N_Edges

        if(cell->GetJoint(j)->GetType() == IsoBoundEdge ||
           (cell->GetJoint(j)->GetType() == IsoInterfaceJoint))
        {
          FEId = VelocitySpace->GetFE2D(Velo_CellNo[i], cell);
          FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
          N_LocalDOFs = FEDesc->GetN_JointDOF();
          t = 2.0/(N_LocalDOFs-1);
          for(k=1;k<N_LocalDOFs-1;k++)
          {
            /*
            Nx[IIso] += normalx;
            Ny[IIso] += normaly;
            */
            // /*
            s = -1.0 + k*t;
            F_K->GetOuterNormal(j, s, nx, ny);
            Nx[IIso] += nx;
            Ny[IIso] += ny;
            // */
            IIso++;
          } // endfor
        } // endif
      } // !InnerJoint
    } // endfor j
  } // endfor i

  N_ = IIso;
  // normalize normal vector
  for(i=0;i<N_;i++)
  {
    x = Nx[i];
    y = Ny[i];
    t = sqrt(x*x+y*y);
    Nx[i] /= t;
    Ny[i] /= t;

    // cout << setw(5) << i << "n = (" << Nx[i] << ", " << Ny[i] << ")" << endl;
  }
} // if(TDatabase::ParamDB->P5 > 0)


//          if(TY[0]>=0 && TY[1]>=0 )
//             cout<< " x0: "<< setw(10)<< TX[0]<< setw(10)<<" x1: " << TX[1]<< setw(10)<<
// 	                 "y1 - y0: "<< setw(10) << TX[1] - TX[0] <<endl;



  // determine new position of boundary vertices
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    OnBoundary = false;
    for(j=0;j<N_Edges;j++)
    {
      if( !(cell->GetJoint(j)->InnerJoint()) 
         || (cell->GetJoint(j)->GetType() == IsoInterfaceJoint) )
        OnBoundary = true;
    } // endfor j

    if(OnBoundary)
    {
//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
      Velo_i = cell->GetGlobalCellNo();
      FEId = VelocitySpace->GetFE2D(Velo_i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      switch(N_Edges)
      {
        case 3:
          xi[0]  = 0; xi[1]  = 1; xi[2]  = 0;
          eta[0] = 0; eta[1] = 0; eta[2] = 1;
          VX[0] = VX[1] = VX[2] = 0;
          VY[0] = VY[1] = VY[2] = 0;
          N_Vertices = N_Edges;
        break;

        case 4:
          // NOTE: sorted due to number of DOFs !!!
          xi[0]  = -1; xi[1]  =  1; xi[2]  = -1; xi[3]  =  1;
          eta[0] = -1; eta[1] = -1; eta[2] =  1; eta[3] =  1;
          VX[0] = VX[1] = VX[2] = VX[3] = 0;
          VY[0] = VY[1] = VY[2] = VY[3] = 0;
          N_Vertices = N_Edges;
        break;

        default:
          Error("only triangles and quadrilaterals are allowed!" << endl);
          exit(-1);
      } // endswitch

      for(j=0;j<N_Vertices;j++)
        bf->GetDerivatives(D00, xi[j], eta[j], FunctValues[j]);

      DOF = VeloGlobalNumbers + VeloBeginIndex[Velo_i];
  
      for(j=0;j<N_LocalDOFs;j++)
      {
        k = DOF[j];
        s = ValuesVX[k];
        t = ValuesVY[k];
        for(l=0;l<N_Vertices;l++)
        {
          VX[l] += FunctValues[l][j]*s;
          VY[l] += FunctValues[l][j]*t;
        } // endfor l
      } // endfor j

      FEId = GridSpace->GetFE2D(i, cell);
      Element = TFEDatabase2D::GetFE2D(FEId);
      BF = Element->GetBaseFunct2D_ID();      
      if( (BF != BF_C_T_P1_2D) && (BF != BF_C_Q_Q1_2D) )
      {
        Error("Grid Space must be conforming and of first order!" << endl);
        exit(-1);
      }  // endif
      bf = Element->GetBaseFunct2D();
      N_LocalDOFs = Element->GetN_DOF();

      DOF = GridGlobalNumbers + GridBeginIndex[i];
  
      for(j=0;j<N_Vertices;j++)
      {
        l = DOF[j];
        k = l - N_Inner;
        if(k>=0)
        {
          if(TDatabase::ParamDB->P5 > 0)
          {
            un = VX[j]*Nx[k] + VY[j]*Ny[k];
            NewValuesX[l] = ValuesX[l] + dt*un*Nx[k];
            NewValuesY[l] = ValuesY[l] + dt*un*Ny[k];
          }
          else
          {
           if(ValuesX[l]!=0)
             NewValuesX[l] = ValuesX[l] + dt*VX[j];
           NewValuesY[l] = ValuesY[l] + dt*VY[j]; 
          }
        }
      } // endfor j
    } // endif
  } // endfor i

  memset(Rhs, 0, 2*GridLength*SizeOfDouble);

   memcpy(d, NewValuesX, 2*GridLength*SizeOfDouble);
   Daxpy(2*GridLength, -1, ValuesX, d);
   memcpy(Rhs + (GridLength-N_BoundaryNodes),
          d+(GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);
   memcpy(Rhs + (2*GridLength-N_BoundaryNodes),
          d+(2*GridLength-N_BoundaryNodes),
          N_BoundaryNodes*SizeOfDouble);

  memset(Sol, 0 , 2*GridLength*SizeOfDouble);
  memcpy(Sol + (GridLength-N_BoundaryNodes),
         d+(GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);
  memcpy(Sol + (2*GridLength-N_BoundaryNodes),
         d+(2*GridLength-N_BoundaryNodes),
         N_BoundaryNodes*SizeOfDouble);

  SolveGridEquation(Entries, Sol, Rhs,
                    KCol, RowPtr, GridLength);

  memcpy(d, ValuesX, 2*GridLength*SizeOfDouble);
  Daxpy(2*GridLength, 1, Sol, d);
  memcpy(NewValuesX, d, (GridLength-N_BoundaryNodes)*SizeOfDouble);
  memcpy(NewValuesY, d+GridLength, (GridLength-N_BoundaryNodes)*SizeOfDouble);


/*
  for(i=0;i<GridLength;i++)
   cout << i << "  ---  "<<Sol[i] << "  ---  " << Sol[i+GridLength] << endl;
*/
  // put solution into grid position
  IIso = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    DOF = GridGlobalNumbers + GridBeginIndex[i];
    N_Vertices = cell->GetN_Vertices();

    switch(N_Vertices)
    {
      case 3:
        for(j=0;j<N_Vertices;j++)
        {
          k = DOF[j];
          cell->GetVertex(j)->SetCoords(NewValuesX[k], NewValuesY[k]);
        }
      break;

      case 4:
        k = DOF[0];
        cell->GetVertex(0)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[1];
        cell->GetVertex(1)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[3];
        cell->GetVertex(2)->SetCoords(NewValuesX[k], NewValuesY[k]);
        k = DOF[2];
        cell->GetVertex(3)->SetCoords(NewValuesX[k], NewValuesY[k]);
      break;
    } // endswitch


  if(isoupdate)
   {
    N_Edges = cell->GetN_Edges();
    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge || joint->GetType() == IsoInterfaceJoint)
      {
        isojoint = (TIsoInterfaceJoint *)joint;
        k = isojoint->GetN_Vertices();
//         cout <<" k " << k << endl;
        Vertices = isojoint->GetVertices();

//      for(Velo_i=0;Velo_i<Velo_N_Cells;Velo_i++)
//       {
//        if(cell==Velo_Coll->GetCell(Velo_i))
//          break;
//        
//        if(Velo_i==Velo_N_Cells-1)
//         {
//         cout<< "error in getvelocity function grid cell not found in velo cells" <<endl;
//         exit(-1);
//         }
//        }
        Velo_i = cell->GetGlobalCellNo();
        
        FEId = VelocitySpace->GetFE2D(Velo_i, cell);
        FEDesc = TFEDatabase2D::GetFEDesc2DFromFE2D(FEId);
        m = FEDesc->GetN_JointDOF();
        if(m == k+2)
        {
          JointDOF = FEDesc->GetJointDOF(j);
          DOF =  VeloGlobalNumbers+VeloBeginIndex[Velo_i];
          for(l=0;l<k;l++)
          {
            m = DOF[JointDOF[l+1]];
            Vertices[l]->GetCoords(IsoX, IsoY);
    
           if(TDatabase::ParamDB->P5>0)
            {
             un = ValuesVX[m]*Nx[IIso+N_BoundaryNodes]
                 + ValuesVY[m]*Ny[IIso+N_BoundaryNodes];
             IsoX += dt*un*Nx[IIso+N_BoundaryNodes];
             IsoY += dt*un*Ny[IIso+N_BoundaryNodes];
            }
           else
            {
             IsoX += dt*ValuesVX[m];
             IsoY += dt*ValuesVY[m];
            }
       
//            if(IsoX<=0.01 && IsoY<0.5)
//             {
// //              cout << IsoY << "  IsoX less than zero !!!!!!!!!!"  << IsoX << endl;
//               cell->GetVertex(j)->GetCoords(x0, y0);
//               cell->GetVertex((j+1) % N_Edges)->GetCoords(x1, y1);           
//  
//              cout << "X: " << x0  << " " <<  IsoX  << " " <<x1 << " "  << (x0+x1)/2.<< " "  << IsoX - (x0+x1)/2.  << endl; 
//              cout << "Y: " << y0  << " " <<  IsoY  << " " <<y1 << " "  <<  (y0+y1)/2. << " "  <<  IsoY - (y0+y1)/2. << endl; 
// 	     
//              cout << "U1: " << ValuesVX[DOF[JointDOF[l]]]  << " " <<  ValuesVX[DOF[JointDOF[l+1]]]  << " " << ValuesVX[DOF[JointDOF[l+2]]]  << " " << endl; 
//              cout << "U2: " << ValuesVY[DOF[JointDOF[l]]]  << " " <<  ValuesVY[DOF[JointDOF[l+1]]]  << " " << ValuesVY[DOF[JointDOF[l+2]]]  << endl;    
//             }
       
            Vertices[l]->SetCoords(IsoX, IsoY);
            IIso++;

          } // endfor l
        }
        else
        {
          // approximation order of isoparametric boundary and velocity
          // element must be the same
          Error("No match in isoparametric case" << endl);
          exit(-1);
        }
      } // endif
    } // endfor j
  }
  } // endfor i

 // gridvelo = GridVelocity->GetValues();
//  memcpy(gridvelo, NewValuesX, 2*GridLength*SizeOfDouble);
//  Daxpy(2*GridLength, -1, ValuesX, gridvelo);
//  Dscal(2*GridLength, 1/dt, gridvelo);
  delete [] d;
  delete [] Nx;
  delete [] Ny;
} // MoveGrid

// ====================================================================
// Get the inner angles of the cells in whole domain
// ====================================================================
void Getcellangle(TFESpace2D *Space, double *MinMaxAngle)
{
 int i,j,k,l, N_Cells, N_Edges;
 int found,  N_LinePoints;

 double TX[4], TY[4], hE[4], Theta, tx, ty, Test, MQI=0.;
 TBaseCell *cell;
 FE2D FEId;
 BF2DRefElements RefElement;
 TRefTrans2D *F_K;
 RefTrans2D RefTrans;
 TCollection *Cells;

  MinMaxAngle[0] = 180;  // Min_Angel = 180
  MinMaxAngle[1] = 0;  // Max_Angel = 0
  Cells = Space->GetCollection();
  N_Cells = Cells->GetN_Cells();
     
//      TX      = new double[4];  // Max no edges in 2d
//      TY      = new double[4];  // Max no edges in 2d
//      hE      = new double[4];  // Max no edges in 2d

  for(i=0;i<N_Cells;i++)
   {
     cell    = Cells->GetCell(i);
     N_Edges = cell->GetN_Edges();

     FEId = Space->GetFE2D(i, cell);
     RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);

     switch(RefElement)
        {
         case BFUnitTriangle:

            RefTrans = TriaAffin;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaAffin *)F_K)->SetCell(cell);

          break;

          case BFUnitSquare:

            RefTrans = QuadAffin;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadAffin *)F_K)->SetCell(cell);

          break;

          default:
            Error("only triangles and quadrilaterals are allowed" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
          } // endswitch

     for(j=0;j<N_Edges;j++)
      {
        F_K->GetTangent(j, 0, tx, ty);
        TX[j] = tx;
        TY[j] = ty;
        hE[j] = sqrt(tx*tx+ty*ty);

   // cout <<"cell : " <<i << "  j= " << j << ": " <<TX[j]<< "------ " << TY[j] << endl;
       } // endfor j

//      Test = 0;
      k = N_Edges -1;
      for(j=0;j<N_Edges;j++)
      {
       if(hE[j]==0.0 || hE[k]== 0.0 )
        Theta = 0.0;
       else
        Theta = acos(-(TX[j]*TX[k]+TY[j]*TY[k])/(hE[j]*hE[k]))*(180/3.141592654);

       k = j;
//        Test +=Theta;
       if(MinMaxAngle[0]>Theta) MinMaxAngle[0] = Theta;
       if(MinMaxAngle[1]<Theta) MinMaxAngle[1] = Theta;
//        cout <<"cell : " <<i << "  j= " << j << ": " << " Theta : " << Theta << endl;
//  *****************************************************
//  Grid test

      MQI += (60. - Theta)*(60. - Theta);
//  *****************************************************

     }
//       cout <<"cell : " <<i <<  " sum of 3 angels : " << Test << endl;
     //  cout<<endl;

   } // endfor i

   MQI /=double(3*N_Cells);
   MQI = sqrt(MQI);

// OutPut("Mesh Quality Indicator: "<< MQI<< endl);
//    delete [] TX;
//    delete [] TY;
//    delete [] hE;
 //cout<< " Min_Angel: "<< MinMaxAngle[0]<< "  Max_Angel : "<<MinMaxAngle[1]<< endl;
// exit(0);
}



double Volume(TFESpace2D *FESpace)
{
  TCollection *Coll;
  double vol, locvol;
  TBaseCell *cell;
  int i, j, N_Cells, N_Edges;
  FE2D FEId;
  RefTrans2D RefTrans;
  TRefTrans2D *rt;
  QuadFormula2D QuadFormula;
  TQuadFormula2D *qf2;
  int polydegree;
  bool IsIsoparametric;
  TJoint *joint;
  JointType jointtype;
  BoundTypes bdtype;

  Coll = FESpace->GetCollection();
  N_Cells = Coll->GetN_Cells();
  vol = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = Coll->GetCell(i);
    FEId = FESpace->GetFE2D(i, cell);

    RefTrans = TFEDatabase2D::GetRefTrans2D_IDFromFE2D(FEId);
    N_Edges = cell->GetN_Edges(); 

    IsIsoparametric = false;
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(j=0;j<N_Edges;j++)
      {
        joint = cell->GetJoint(j);
        jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == InterfaceJoint)
        {
          bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == IsoInterfaceJoint ||
           jointtype == IsoBoundEdge)
          IsIsoparametric = true;
      }
    } // endif 

    if(IsIsoparametric)
    {
      switch(N_Edges)
      {
        case 4:
          RefTrans = QuadIsoparametric;
        break;

        case 3:
          RefTrans = TriaIsoparametric;
        break;
      }
    } // endif IsIsoparametric

    rt = TFEDatabase2D::GetRefTrans2D(RefTrans);
    switch(RefTrans)
    {
      case TriaAffin:
        ((TTriaAffin *)rt)->SetCell(cell);
        locvol = ((TTriaAffin *)rt)->GetVolume();
      break;

      case TriaIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        QuadFormula = TFEDatabase2D::GetQFTriaFromDegree(2*polydegree-1);
        ((TTriaIsoparametric *)rt)->SetApproximationOrder(polydegree);
        ((TTriaIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TTriaIsoparametric *)rt)->SetCell(cell);
        locvol = ((TTriaIsoparametric *)rt)->GetVolume();
      break;

      case QuadAffin:
        ((TQuadAffin *)rt)->SetCell(cell);
        locvol = ((TQuadAffin *)rt)->GetVolume();
      break;

      case QuadBilinear:
        ((TQuadBilinear *)rt)->SetCell(cell);
        locvol = ((TQuadBilinear *)rt)->GetVolume();
      break;

      case QuadIsoparametric:
        polydegree = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        QuadFormula = TFEDatabase2D::GetQFQuadFromDegree(3*polydegree);
        ((TQuadIsoparametric *)rt)->SetApproximationOrder(polydegree);
        ((TQuadIsoparametric *)rt)->SetQuadFormula(QuadFormula);
        ((TQuadIsoparametric *)rt)->SetCell(cell);
        locvol = ((TQuadIsoparametric *)rt)->GetVolume();
      break;
    }

    vol += locvol;
  } // endfor i

  return vol;
}


void ReParametrize_pts(int &N_Edges, TBaseCell **cell, int *EdgeNo, double h_min, double **FreePts)
{
  int i, j, ISpline, N_Splines, N_V, ORDER, VSP, m, k, i3, N_E = N_Edges;
  double *h, *t;
  double *a, *b, *c, *x, *y, teta;
  double *rhs, *Mx, *My, *Params, *Param9;
  double phi1, phi2, phi3, phi4, X, Y, T;
  double dx0, dy0, dx1, dy1, Isoteta;
  TIsoBoundEdge *isojoint;
  TIsoInterfaceJoint *isoIntjoint;
  TVertex **IsoVertices;
  TJoint *Joint;
  TBaseCell *Me;

  ORDER = 0;
  VSP = TDatabase::ParamDB->VELOCITY_SPACE;

    if (abs(VSP) > 20)
   {ORDER = abs(VSP) - 20;}
  else if ( abs(VSP) > 10)
   {ORDER = abs(VSP) - 10;}
  else ORDER = abs(VSP);

  N_V = N_E+ N_E*(ORDER-1) + 1; // add end point for spline

  N_Splines = N_V-1;
  h = new double[N_Splines+1];
  t = new double[N_Splines+1];
  a = new double[N_Splines+1];
  b = new double[N_Splines+1];
  c = new double[N_Splines+1];
  rhs = new double[N_Splines+1];
  Mx = new double[N_Splines+1];
  My = new double[N_Splines+1];
  Params = new double [10*N_Splines];
  Param9 = new double [N_Splines+1];

  x = new double[N_V];
  y = new double[N_V];

   m = 0;
   for(i=0;i<N_E;i++)
   {
    Me = cell[i];
    Me->GetVertex(EdgeNo[i])->GetCoords(x[m], y[m]);
    m++;

    Joint = cell[i]->GetJoint(EdgeNo[i]);
    isojoint = (TIsoBoundEdge *)Joint;
    k = isojoint->GetN_Vertices();
    if(k==ORDER-1)
     {
      IsoVertices = isojoint->GetVertices();
      for(i3=0;i3<k;i3++)
       {   
        IsoVertices[i3]->GetCoords(x[m], y[m]);
//          cout<< i<<" FreeGaus " << (180/Pi)*atan2(y[m], x[m]) <<endl;
        m++;
       } 
     }
    else
     {
      // only second order conforming elements implimented
      cout<< " No match in isopoints per free edge in Remesc2D.C "<<endl;
      exit(0);
     }
   }

//   end vertex of the freeboundary
  k = cell[N_E-1]->GetN_Edges();
  cell[N_E-1]->GetVertex((EdgeNo[N_E-1]+1) % k)->GetCoords(x[m], y[m]);

 h[0] = 0.0; t[0] = 0.0;

 for(i=1;i<=N_Splines;i++)
  {
    h[i] = sqrt((x[i]-x[i-1])*(x[i]-x[i-1])+(y[i]-y[i-1])*(y[i]-y[i-1]));
    t[i] = t[i-1] + h[i];
  }

  dx0 = (x[1]-x[0])/h[1];
  dy0 = (y[1]-y[0])/h[1];

  dx1 = (x[N_Splines]-x[N_Splines-1])/h[N_Splines];
  dy1 = (y[N_Splines]-y[N_Splines-1])/h[N_Splines];


  a[0] = 2.; c[0] = 1.; rhs[0] = -6./h[1]*(dx0 - (x[1]-x[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    a[i] = 2.;
    b[i] = h[i]/(h[i]+h[i+1]);
    c[i] = h[i+1]/(h[i]+h[i+1]);
    rhs[i] = 6./(h[i]+h[i+1])*((x[i+1]-x[i])/h[i+1]-(x[i]-x[i-1])/h[i]);
  }
  b[N_Splines] = 1.; a[N_Splines] = 2.;
  rhs[N_Splines] = 6./h[N_Splines]*(dx1 - (x[N_Splines]-x[N_Splines-1])/h[N_Splines]);
  
  Solver_3dia(N_Splines, a, b, c, rhs, Mx);


  rhs[0] = -6./h[1]*(dy0 - (y[1]-y[0])/h[1]);
  for(i=1;i<N_Splines;i++)
  {
    rhs[i] = 6./(h[i]+h[i+1])*((y[i+1]-y[i])/h[i+1]-(y[i]-y[i-1])/h[i]);
  }
  rhs[N_Splines] = 6./h[N_Splines]*(dy1 - (y[N_Splines]-y[N_Splines-1])/h[N_Splines]);
  
  Solver_3dia(N_Splines, a, b, c, rhs, My);

  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;
    Params[ISpline    ] = x[i]; 
    Params[ISpline + 1] = y[i];
    Params[ISpline + 2] = x[i+1]; 
    Params[ISpline + 3] = y[i+1];
    Params[ISpline + 4] = -Mx[i]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

    Params[ISpline + 5] = -My[i]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 6] = Mx[i+1]*h[i+1]*h[i+1]/2. +
                          ((x[i+1]-x[i])/h[i+1]-h[i+1]/6.*(Mx[i+1]-Mx[i]))*h[i+1];

    Params[ISpline + 7] = My[i+1]*h[i+1]*h[i+1]/2. +
                          ((y[i+1]-y[i])/h[i+1]-h[i+1]/6.*(My[i+1]-My[i]))*h[i+1];
    Params[ISpline + 8] = t[i+1]/t[N_Splines];
    Params[ISpline + 9] = 0.;
  
   //cout<<"  "<<Params[ISpline + 8]<<'\t'<<Params[ISpline + 9]<<endl;
  }

  /*
  for(i=0;i<N_Splines;i++)
  {
    ISpline = i*10;
    //cout<<i+1<< " subspline"<<endl;
    cout<<setprecision(15);
    cout<<"  "<<Params[ISpline    ]<<'\t'<<Params[ISpline + 1]<<endl;
    cout<<"  "<<Params[ISpline + 2]<<'\t'<<Params[ISpline + 3]<<endl;
    cout<<"  "<<Params[ISpline + 4]<<'\t'<<Params[ISpline + 5]<<endl;
    cout<<"  "<<Params[ISpline + 6]<<'\t'<<Params[ISpline + 7]<<endl;
    cout<<"  "<<Params[ISpline + 8]<<'\t'<<Params[ISpline + 9]<<endl;
  }
  */

   delete [] FreePts[0];
   delete [] FreePts[1];

   N_E = int(t[N_Splines]/h_min);
  // minimum points will not be less than the initial number of points (used in stbubble problem)
  //    if(N_E<int(TDatabase::ParamDB->P6))
 //         N_E = int(TDatabase::ParamDB->P6);

   N_Edges = N_E;
   FreePts[0] = new double[N_E+1];
   FreePts[1] = new double[N_E+1];

   teta = 1.0/N_E;

   T = 0.;
   Param9[0] = 0.;
   for(i=1;i<=N_Splines;i++) 
    Param9[i] = Params[(i-1)*10+8];

   m = 0;
   for(j=0;j<=N_E;j++)
    {
     T = double(j)*teta;
     for(i=1;i<=N_Splines;i++)
      {
       ISpline = (i-1)*10;
       if((T>=Param9[i-1]) && (T<=Param9[i]))
        {
      // further T must be from [0;1] on a subspline
//           cout<< ISpline << ' ' << T;
         T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
//           cout<< ' ' << T;
         break;
        }
      }
  
    phi1 = (2.*T*T - 3.*T)*T + 1.;
    phi2 = (-2.*T + 3.)*T*T;
    phi3 = (T*T - 2.*T + 1.)*T;
    phi4 = (T - 1)*T*T;

    X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
        Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
    Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
        Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

    if(fabs(X)<1.e-12)
      X = 0.;
    
    FreePts[0][j] = X;
    FreePts[1][j] = Y;
//        cout<< j  << " X " << X << " Y " << Y<<endl; 
    if(j<N_E)
     {
     // need to set it for surface interpolation maping
     Joint = cell[j]->GetJoint(EdgeNo[j]);
     isoIntjoint = (TIsoInterfaceJoint *)Joint;
     k = isojoint->GetN_Vertices();
     IsoVertices = isojoint->GetVertices();
     Isoteta = 1./(double)(k+1);

     for(i3=0;i3<k;i3++)
      {   
       T = ( (double)(j) + ((double)(i3+1.))*Isoteta )*teta;
       for(i=1;i<=N_Splines;i++)
        {
         ISpline = (i-1)*10;
         if((T>=Param9[i-1]) && (T<=Param9[i]))
          {
      // further T must be from [0;1] on a subspline
//          cout<< ISpline << ' ' << T;
           T = (T-Param9[i-1])/(Param9[i]-Param9[i-1]);
//          cout<< ' ' << T <<endl;
          break;
         }
	}
      phi1 = (2.*T*T - 3.*T)*T + 1.;
      phi2 = (-2.*T + 3.)*T*T;
      phi3 = (T*T - 2.*T + 1.)*T;
      phi4 = (T - 1)*T*T;

      X = Params[ISpline    ]*phi1 + Params[ISpline + 2]*phi2 +
	  Params[ISpline + 4]*phi3 + Params[ISpline + 6]*phi4;
      Y = Params[ISpline + 1]*phi1 + Params[ISpline + 3]*phi2 +
	  Params[ISpline + 5]*phi3 + Params[ISpline + 7]*phi4;

      IsoVertices[i3]->SetCoords(X, Y);
      //     cout<< j << "ISO X " << X << " ISO Y " << Y<<endl;
     } //   for(i3=0;i3<k;i3
    } // if(j<N_E)

   }  // for j
//   end vertex of the freeboundary
 

   delete [] h; delete [] t; delete [] a; delete [] b;
   delete [] c; delete [] rhs; delete [] Mx; delete [] My;
   delete [] Params; delete [] Param9;  delete [] x; delete [] y;
}




void IntUn(TFEVectFunct2D *u, double *Nx, double *Ny)
{
  int i,j,k,l,m,n;
  TBaseCell *cell;
  TCollection *coll;
  int N_Cells, N_Edges, N_DOF;
  int N_LinePoints;
  double *zeta, *LineWeights;
  double *u1, *u2;
  const TFESpace2D *USpace;
  TJoint *joint;
  FE2D FEId;
  TFE2D *Element;
  TFEDesc2D *FEDesc;
  BaseFunct2D BF;
  TBaseFunct2D *bf;
  TQuadFormula1D *qf1;
  QuadFormula1D LineQuadFormula;
  BF2DRefElements RefElement;
  TRefTrans2D *F_K;
  RefTrans2D RefTrans;
  double **JointValues, *Values;
  double x0, x1, y0, y1, s, t;
  double un, int_un, loc_un;
  double n1, n2, t1, t2;
  double u1loc, u2loc;
  int *DOF;
  int *BeginIndex, *GlobalNumbers;

  USpace = u->GetFESpace2D();
  coll = USpace->GetCollection();
  N_Cells = coll->GetN_Cells();
  N_DOF = u->GetLength();
  u1 = u->GetValues();
  u2 = u1+N_DOF;
  BeginIndex = USpace->GetBeginIndex();
  GlobalNumbers = USpace->GetGlobalNumbers();

  int_un = 0;
  for(i=0;i<N_Cells;i++)
  {
    cell = coll->GetCell(i);
    N_Edges = cell->GetN_Edges();
    DOF = GlobalNumbers + BeginIndex[i];

    for(j=0;j<N_Edges;j++)
    {
      joint = cell->GetJoint(j);
      if(joint->GetType() == IsoBoundEdge)
      {
        loc_un = 0;

        FEId = USpace->GetFE2D(i, cell);
        l = TFEDatabase2D::GetPolynomialDegreeFromFE2D(FEId);
        LineQuadFormula =  TFEDatabase2D::GetQFLineFromDegree(2*l);
        qf1 = TFEDatabase2D::GetQuadFormula1D(LineQuadFormula);
        qf1->GetFormulaData(N_LinePoints, LineWeights, zeta);
        TFEDatabase2D::GetBaseFunct2DFromFE2D(FEId)
                    ->MakeRefElementData(LineQuadFormula);
        BF = TFEDatabase2D::GetBaseFunct2D_IDFromFE2D(FEId);
        N_DOF = TFEDatabase2D::GetN_BaseFunctFromFE2D(FEId);

        RefElement = TFEDatabase2D::GetRefElementFromFE2D(FEId);
        switch(RefElement)
        {
          case BFUnitTriangle:
            RefTrans = TriaIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TTriaIsoparametric *)F_K)->SetCell(cell);
          break;

          case BFUnitSquare:
            RefTrans = QuadIsoparametric;
            F_K = TFEDatabase2D::GetRefTrans2D(RefTrans);
            ((TQuadIsoparametric *)F_K)->SetCell(cell);
          break;

          default:
            Error("only triangles and quadrilaterals are allowes" << endl);
            Error("file: " << __FILE__ << " line " << __LINE__ << endl);
            exit(-1);
        } // endswitch RefElement

        cell->GetVertex(j)->GetCoords(x0, y0);
        cell->GetVertex((j+1)%N_Edges)->GetCoords(x1, y1);
        t = x1-x0;
        s = y1-y0;

        JointValues = TFEDatabase2D::GetJointValues2D(BF, LineQuadFormula, j);
        for(k=0;k<N_LinePoints;k++)
        {
          un = 0;
          Values = JointValues[k];
          F_K->GetOuterNormal(j, zeta[k], n1, n2);
          u1loc = u2loc = 0;
          for(l=0;l<N_DOF;l++)
          {
            // cout << "l= " << l << ": " << Values[l];
            // cout << "zeta[k]: " << zeta[k] << endl;
            m = DOF[l];
            u1loc += u1[m]*Values[l];
            u2loc += u2[m]*Values[l];
          } // endfor l
          un += u1loc*n1 + u2loc*n2;
          F_K->GetTangent(j, zeta[k], t1, t2);
          loc_un += sqrt(t1*t1+t2*t2)*LineWeights[k]*un;
        } // endfor k
        int_un += loc_un;
      } // endif
    } // endfor j
  } // endfor i

  OutPut("Int(u.n) = " << int_un << endl);
  // exit(-1);
}

#endif
