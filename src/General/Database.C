// =======================================================================
// @(#)Database.C        1.37 06/27/00
// 
// Class:       TDatabase
// Purpose:     database of needed refinement, mapping and
//              shape descriptors
//              parameter database
//
// Author:      Volker Behns  29.07.97
//
// =======================================================================
#ifdef _MPI
#  include "mpi.h"
#endif

#include <Database.h>
#include <Constants.h>
#include <Mapper.h>
#include <Line.h>
#include <ParameterDatabase.h>
#include <MainUtilities.h>
// more explicit path due to name conflict with the 'triangle' library on mac
#include <../Geometry/Triangle.h> 
#include <Quadrangle.h>
#include <Parallelogram.h>
#include <Rectangle.h>
#include <RefNoRef.h>
#include <RefLineDesc.h>
#include <RefQuadRegDesc.h>
#include <RefQuadBis0Desc.h>
#include <RefQuadBis1Desc.h>
#include <RefQuad1Conf0Desc.h>
#include <RefQuad1Conf1Desc.h>
#include <RefQuad1Conf2Desc.h>
#include <RefQuad1Conf3Desc.h>
#include <RefQuad2Conf0Desc.h>
#include <RefQuad2Conf1Desc.h>
#include <RefQuad2Conf2Desc.h>
#include <RefQuad2Conf3Desc.h>
#include <RefQuadToTri0Desc.h>
#include <RefQuadToTri1Desc.h>
#include <RefTriRegDesc.h>
#include <RefTriBaryDesc.h>
#include <RefTriBis0Desc.h>
#include <RefTriBis1Desc.h>
#include <RefTriBis2Desc.h>
#include <RefTriBis01Desc.h>
#include <RefTriBis02Desc.h>
#include <RefTriBis10Desc.h>
#include <RefTriBis12Desc.h>
#include <RefTriBis20Desc.h>
#include <RefTriBis21Desc.h>
#include <It_Between.h>
#include <It_EQ.h>
#include <It_Finest.h>
#include <It_EQLevel.h>
#include <It_LELevel.h>
#include <It_OCAF.h>
#include <Utilities.h>
#include <MooNMD_Io.h>
#include <string.h>

#ifdef __MORTAR__
  #include <It_Mortar.h>
  #include <RefMortar0Desc.h>
  #include <RefMortar1Desc.h>
  #include <RefMortarLineDesc.h>
#endif

#ifdef __3D__
  #include <Tetrahedron.h>
  #include <Hexahedron.h>
  #include <Brick.h>
  #include <RefTetraRegDesc.h>
  #include <RefTetraBaryDesc.h>
  #include <RefTetraReg0Desc.h>
  #include <RefTetraReg1Desc.h>
  #include <RefTetraReg2Desc.h>
  #include <RefTetraBis0Desc.h>
  #include <RefTetraBis1Desc.h>
  #include <RefTetraBis2Desc.h>
  #include <RefTetraBis3Desc.h>
  #include <RefTetraBis4Desc.h>
  #include <RefTetraBis5Desc.h>
  #include <RefTetraBis01Desc.h>
  #include <RefTetraBis02Desc.h>
  #include <RefTetraBis03Desc.h>
  #include <RefTetraBis04Desc.h>
  #include <RefTetraBis05Desc.h>
  #include <RefTetraBis10Desc.h>
  #include <RefTetraBis12Desc.h>
  #include <RefTetraBis13Desc.h>
  #include <RefTetraBis14Desc.h>
  #include <RefTetraBis15Desc.h>
  #include <RefTetraBis20Desc.h>
  #include <RefTetraBis21Desc.h>
  #include <RefTetraBis23Desc.h>
  #include <RefTetraBis24Desc.h>
  #include <RefTetraBis25Desc.h>
  #include <RefTetraBis30Desc.h>
  #include <RefTetraBis32Desc.h>
  #include <RefTetraBis34Desc.h>
  #include <RefTetraBis35Desc.h>
  #include <RefTetraBis40Desc.h>
  #include <RefTetraBis41Desc.h>
  #include <RefTetraBis43Desc.h>
  #include <RefTetraBis45Desc.h>
  #include <RefTetraBis51Desc.h>
  #include <RefTetraBis52Desc.h>
  #include <RefTetraBis53Desc.h>
  #include <RefTetraBis54Desc.h>
  #include <RefTetraQuad0Desc.h>
  #include <RefTetraQuad1Desc.h>
  #include <RefTetraQuad2Desc.h>
  #include <RefTetraQuad3Desc.h>
  #include <RefHexaRegDesc.h>
#endif

// Constructors
TDatabase::TDatabase()
{
  // allocate databases
  ShapeDB = new TShapeDesc*[N_SHAPES];
  RefDescDB = new TRefDesc*[N_SHAPES + N_REFDESC + 2*N_MORTARDESC];
  MapperDB = new TMapper*[N_MAPPER];
  IteratorDB = new TIterator*[N_ITERATORS];
  ParamDB = new TParamDB;
  TimeDB = new TTimeDB;

  // initialize shape descriptors
  ShapeDB[S_Line] = new TLine();
  RefDescDB[S_Line] = new TRefNoRef(ShapeDB[S_Line]);

  ShapeDB[Triangle] = new TTriangle();
  RefDescDB[Triangle] = new TRefNoRef(ShapeDB[Triangle]);

  ShapeDB[Quadrangle] = new TQuadrangle();
  RefDescDB[Quadrangle] = new TRefNoRef(ShapeDB[Quadrangle]);

  ShapeDB[Parallelogram] = new TParallelogram();
  RefDescDB[Parallelogram] = new TRefNoRef(ShapeDB[Parallelogram]);

  ShapeDB[Rectangle] = new TRectangle();
  RefDescDB[Rectangle] = new TRefNoRef(ShapeDB[Rectangle]);

  #ifdef __3D__
    ShapeDB[Tetrahedron] = new TTetrahedron();
    RefDescDB[Tetrahedron] = new TRefNoRef(ShapeDB[Tetrahedron]);

    ShapeDB[Hexahedron] = new THexahedron();
    RefDescDB[Hexahedron] = new TRefNoRef(ShapeDB[Hexahedron]);

    ShapeDB[Brick] = new TBrick();
    RefDescDB[Brick] = new TRefNoRef(ShapeDB[Brick]);
  #endif

  // initialize refinement descriptors
  RefDescDB[N_SHAPES + LineReg] = new TRefLineDesc(ShapeDB[S_Line]);
  RefDescDB[N_SHAPES + TriReg]  = new TRefTriRegDesc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBary]  = new TRefTriBaryDesc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis0] = new TRefTriBis0Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis1] = new TRefTriBis1Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis2] = new TRefTriBis2Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis01]= new TRefTriBis01Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis02]= new TRefTriBis02Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis10]= new TRefTriBis10Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis12]= new TRefTriBis12Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis20]= new TRefTriBis20Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + TriBis21]= new TRefTriBis21Desc(ShapeDB[Triangle]);
  RefDescDB[N_SHAPES + QuadReg] = new TRefQuadRegDesc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + ParallReg] = new TRefQuadRegDesc(ShapeDB[Parallelogram]);
  RefDescDB[N_SHAPES + RectReg] = new TRefQuadRegDesc(ShapeDB[Rectangle]);
  RefDescDB[N_SHAPES + QuadBis0] = new TRefQuadBis0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + QuadBis1] = new TRefQuadBis1Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf0] = new TRefQuad1Conf0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf1] = new TRefQuad1Conf1Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf2] = new TRefQuad1Conf2Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad1Conf3] = new TRefQuad1Conf3Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf0] = new TRefQuad2Conf0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf1] = new TRefQuad2Conf1Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf2] = new TRefQuad2Conf2Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES+Quad2Conf3] = new TRefQuad2Conf3Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + QuadToTri0] = new
      TRefQuadToTri0Desc(ShapeDB[Quadrangle]);
  RefDescDB[N_SHAPES + QuadToTri1] = new
      TRefQuadToTri1Desc(ShapeDB[Quadrangle]);

  #ifdef __3D__
    RefDescDB[N_SHAPES + TetraReg] =
         new TRefTetraRegDesc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBary] =
         new TRefTetraBaryDesc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraReg0] =
         new TRefTetraReg0Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraReg1] =
         new TRefTetraReg1Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraReg2] =
         new TRefTetraReg2Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis0] = new TRefTetraBis0Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis1] = new TRefTetraBis1Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis2] = new TRefTetraBis2Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis3] = new TRefTetraBis3Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis4] = new TRefTetraBis4Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis5] = new TRefTetraBis5Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis01] = new TRefTetraBis01Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis02] = new TRefTetraBis02Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis03] = new TRefTetraBis03Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis04] = new TRefTetraBis04Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis05] = new TRefTetraBis05Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis10] = new TRefTetraBis10Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis12] = new TRefTetraBis12Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis13] = new TRefTetraBis13Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis14] = new TRefTetraBis14Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis15] = new TRefTetraBis15Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis20] = new TRefTetraBis20Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis21] = new TRefTetraBis21Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis23] = new TRefTetraBis23Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis24] = new TRefTetraBis24Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis25] = new TRefTetraBis25Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis30] = new TRefTetraBis30Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis32] = new TRefTetraBis32Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis34] = new TRefTetraBis34Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis35] = new TRefTetraBis35Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis40] = new TRefTetraBis40Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis41] = new TRefTetraBis41Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis43] = new TRefTetraBis43Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis45] = new TRefTetraBis45Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis51] = new TRefTetraBis51Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis52] = new TRefTetraBis52Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis53] = new TRefTetraBis53Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraBis54] = new TRefTetraBis54Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad0] = new TRefTetraQuad0Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad1] = new TRefTetraQuad1Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad2] = new TRefTetraQuad2Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + TetraQuad3] = new TRefTetraQuad3Desc(ShapeDB[Tetrahedron]);
    RefDescDB[N_SHAPES + HexaReg]  =
         new TRefHexaRegDesc(ShapeDB[Hexahedron]);
    RefDescDB[N_SHAPES + BrickReg]  =
         new TRefHexaRegDesc(ShapeDB[Brick]);
  #endif

  #ifdef __3D__
    //initialize mapper
    MapperDB[MapTriReg0] = new TMapper(MapTriReg0);
    MapperDB[MapTriReg1] = new TMapper(MapTriReg1);
    MapperDB[MapTriReg2] = new TMapper(MapTriReg2);

    MapperDB[MapTriBis00] = new TMapper(MapTriBis00);
    MapperDB[MapTriBis01] = new TMapper(MapTriBis01);
    MapperDB[MapTriBis02] = new TMapper(MapTriBis02);
    MapperDB[MapTriBis10] = new TMapper(MapTriBis10);
    MapperDB[MapTriBis11] = new TMapper(MapTriBis11);
    MapperDB[MapTriBis12] = new TMapper(MapTriBis12);
    MapperDB[MapTriBis20] = new TMapper(MapTriBis20);
    MapperDB[MapTriBis21] = new TMapper(MapTriBis21);
    MapperDB[MapTriBis22] = new TMapper(MapTriBis22);
    MapperDB[MapTriBis010] = new TMapper(MapTriBis010);
    MapperDB[MapTriBis011] = new TMapper(MapTriBis011);
    MapperDB[MapTriBis012] = new TMapper(MapTriBis012);
    MapperDB[MapTriBis020] = new TMapper(MapTriBis020);
    MapperDB[MapTriBis021] = new TMapper(MapTriBis021);
    MapperDB[MapTriBis022] = new TMapper(MapTriBis022);
    MapperDB[MapTriBis100] = new TMapper(MapTriBis100);
    MapperDB[MapTriBis101] = new TMapper(MapTriBis101);
    MapperDB[MapTriBis102] = new TMapper(MapTriBis102);
    MapperDB[MapTriBis120] = new TMapper(MapTriBis120);
    MapperDB[MapTriBis121] = new TMapper(MapTriBis121);
    MapperDB[MapTriBis122] = new TMapper(MapTriBis122);
    MapperDB[MapTriBis200] = new TMapper(MapTriBis200);
    MapperDB[MapTriBis201] = new TMapper(MapTriBis201);
    MapperDB[MapTriBis202] = new TMapper(MapTriBis202);
    MapperDB[MapTriBis210] = new TMapper(MapTriBis210);
    MapperDB[MapTriBis211] = new TMapper(MapTriBis211);
    MapperDB[MapTriBis212] = new TMapper(MapTriBis212);

    MapperDB[MapQuadReg0] = new TMapper(MapQuadReg0);
    MapperDB[MapQuadReg1] = new TMapper(MapQuadReg1);
    MapperDB[MapQuadReg2] = new TMapper(MapQuadReg2);
    MapperDB[MapQuadReg3] = new TMapper(MapQuadReg3);
  #endif

  // initialize iterators
  IteratorDB[It_EQ] = new TIt_EQ();
  IteratorDB[It_LE] = new TIt_LE();
  IteratorDB[It_Finest] = new TIt_Finest();
  IteratorDB[It_EQLevel] = new TIt_EQLevel();
  IteratorDB[It_LELevel] = new TIt_LELevel();
  IteratorDB[It_Between] = new TIt_Between();
  IteratorDB[It_OCAF] = new TIt_OCAF();

  #ifdef __MORTAR__
    IteratorDB[It_Mortar1] = new TIt_Mortar();
    IteratorDB[It_Mortar2] = new TIt_Mortar();
  #endif  
  // Initialization of the default parameters
  SetDefaultParameters();
}

TShapeDesc **TDatabase::ShapeDB = nullptr;
TRefDesc   **TDatabase::RefDescDB = nullptr;
TMapper    **TDatabase::MapperDB = nullptr;
TIterator  **TDatabase::IteratorDB = nullptr;
TParamDB   *TDatabase::ParamDB = nullptr;
TTimeDB    *TDatabase::TimeDB = nullptr;

// Methods

#ifdef __MORTAR__

void TDatabase::AddMortar0(int Mortar_Ni, int N)
{
  RefDescDB[N_SHAPES + Mortar + Mortar_Ni] = new
               TRefMortar0Desc(ShapeDB[Quadrangle], Mortar_Ni, N);
  RefDescDB[N_SHAPES + MortarLine + Mortar_Ni] = new
               TRefMortarLineDesc(ShapeDB[S_Line], N);
}

void TDatabase::AddMortar1(int Mortar_Ni, int N)
{
  RefDescDB[N_SHAPES + Mortar + Mortar_Ni] = new
               TRefMortar1Desc(ShapeDB[Quadrangle], Mortar_Ni, N);
  RefDescDB[N_SHAPES + MortarLine + Mortar_Ni] = new
               TRefMortarLineDesc(ShapeDB[S_Line], N);
}
#endif

void TDatabase::SetDefaultParameters()
{
  char *tmp;
  ParamDB->VERSION = 1;

  tmp = new char[12];
  strcpy(tmp,"NO_MAP_FILE");
  ParamDB->MAPFILE=tmp;
  tmp = new char[25];
  strcpy(tmp,"MooN_MD_default_outfile");
  ParamDB->OUTFILE=tmp;
  
  ParamDB->timeprofiling = 0; //time profiling
  ParamDB->MapperType = 1;
  ParamDB->DSType = 1;		//Parallel Direct Solver Type

  ParamDB->MESHGEN_ALLOW_EDGE_REF=0;
  ParamDB->MESHGEN_REF_QUALITY=30;
 
  ParamDB->RE_NR=1.0;
  ParamDB->RA_NR=1.0;
  ParamDB->ROSSBY_NR=0.0;
  ParamDB->START_RE_NR= -4711;
  ParamDB->RE_NR_INCREMENT=1.0;
  ParamDB->FLOW_PROBLEM_TYPE = 0;
  ParamDB->OSEEN_ZERO_ORDER_COEFF = 0.0;

  ParamDB->FR_NR=1.0;
  ParamDB->WB_NR= 1.0;
  ParamDB->PR_NR=1.0;
  ParamDB->PE_NR=1.0;  
  ParamDB->BI_NR = 0;
  ParamDB->Axial3D = 0;
  ParamDB->Axial3DAxis = 0;  

  ParamDB->ANSATZ_ORDER = 2;
  ParamDB->TEST_ORDER = 2;

  ParamDB->VELOCITY_SPACE = 22;
  ParamDB->PRESSURE_SPACE = -4711;
  ParamDB->PRESSURE_SEPARATION = 0;
  ParamDB->OMPNUMTHREADS=1;

  
  ParamDB->REFINEMENT = 0;
  ParamDB->GRID_TYPE = 0;
  ParamDB->GRID_TYPE_1 = 0;
  ParamDB->GRID_TYPE_2 = 0;
  ParamDB->CHANNEL_GRID_STRETCH = 2.75;
  ParamDB->ADAPTIVE_REFINEMENT_CRITERION = 3;
  ParamDB->ERROR_CONTROL = 1;
  ParamDB->REFINE_STRATEGY = 0;
  ParamDB->MAX_CELL_LEVEL = 1000;
  ParamDB->REFTOL = 0.5;
  ParamDB->COARSETOL = 0.0;  
  ParamDB->MIN_FRACTION_TO_CHANGE = 0.1;
  ParamDB->DECREASE_REFTOL_FACTOR = 0.8;
  ParamDB->INCREASE_COARSETOL_FACTOR = 1.1;
  ParamDB->FRACTION_OF_ERROR = 0.25;
  ParamDB->CONVERT_QUAD_TO_TRI = 0;
  ParamDB->N_CELL_LAYERS = 1;

  ParamDB->INTL_DISCTYPE = 1;
  ParamDB->UPWIND_ORDER = 0;
  ParamDB->UPWIND_FLUX_DAMP = 1;
  ParamDB->UPWIND_APPLICATION = 0;
  ParamDB->SHISHKIN_MESH = 0;
  ParamDB->SHISHKIN_DIAM = 1.0;
  ParamDB->NSTYPE = 1;
  ParamDB->DARCYTYPE = 1;
    
  ParamDB->SIGMA_PERM = 1;

  ParamDB->equal_order_stab_weight_PkPk = 0;
  ParamDB->grad_div_stab_weight = 0;
  ParamDB->SIGN_MATRIX_BI = 1;
  ParamDB->l_T = 1;
  ParamDB->L_0 = 1;
  ParamDB->SOURCE_SINK_FUNCTION = false;

  
  ParamDB->LAPLACETYPE = 0;
  ParamDB->USE_ISOPARAMETRIC = 1;
  ParamDB->VMM_COARSE_LEVEL = 4711;
  ParamDB->VMM_COARSE_SPACE_ORDER = 1;
  ParamDB->RFB_SUBMESH_LAYERS = 3;
  ParamDB->DEFECT_CORRECTION_TYPE = 0;
  ParamDB->CELL_MEASURE = 0;

  ParamDB->FACE_SIGMA = 1;
  ParamDB->WEAK_BC_SIGMA = 1;
  ParamDB->WEAK_BC = 0;
  ParamDB->TAU=1.0;
  ParamDB->TAU2=1.0;
  ParamDB->TAU3=1.0;
 
  ParamDB->DELTA0=1.0;
  ParamDB->SDFEM_POWER0=1.0;
  ParamDB->DELTA1=1.0;
  ParamDB->SDFEM_TYPE=2;
  ParamDB->SDFEM_NORM_B=0; // l_infty
  ParamDB->CIP_TYPE=0;
  ParamDB->ADJOINT_FACTOR_4_OMEGA_EQ_0 = 10.0;

  ParamDB->FILTER_WIDTH_CONSTANT = 2;
  ParamDB->FILTER_WIDTH_POWER = 1;
  ParamDB->GAUSSIAN_GAMMA = 6;
  ParamDB->CONVOLUTE_SOLUTION = 0;
  
  ParamDB->TURBULENT_VISCOSITY_TYPE = 1;
  ParamDB->TURBULENT_VISCOSITY_TENSOR = 0;
  ParamDB->TURBULENT_VISCOSITY_CONSTANT = 0.01;
  ParamDB->TURBULENT_VISCOSITY_POWER = 1;
  ParamDB->TURBULENT_VISCOSITY_SIGMA = 6;

  ParamDB->ARTIFICIAL_VISCOSITY_CONSTANT = 1; // parameters for VMS
  ParamDB->ARTIFICIAL_VISCOSITY_POWER = 1;

  ParamDB->FRICTION_CONSTANT = 0.0;      // free slip 
  ParamDB->FRICTION_POWER = 0.0;         // free slip
  ParamDB->FRICTION_TYPE = 0;            // friction type
  ParamDB->FRICTION_U0 = 1.0;            // U_0
  ParamDB->PENETRATION_CONSTANT = 1e12;  // no penetration
  ParamDB->PENETRATION_POWER = -2;        // no penetration

  ParamDB->DIV_DIV_STAB_TYPE = 0;        // stabilization for div-div term 
  ParamDB->DIV_DIV_STAB_C1 = 2;
  ParamDB->DIV_DIV_STAB_C2 = 1;

  ParamDB->LP_FULL_GRADIENT = 1;
  ParamDB->LP_STREAMLINE = 0;
  ParamDB->LP_DIVERGENCE = 0;
  ParamDB->LP_PRESSURE = 0;
  ParamDB->LP_COEFF_TYPE = 0;

  ParamDB->LP_FULL_GRADIENT_COEFF = 1.0;
  ParamDB->LP_STREAMLINE_COEFF= 1.0;
  ParamDB->LP_DIVERGENCE_COEFF = 1.0;
  ParamDB->LP_PRESSURE_COEFF = 1.0;

  ParamDB->LP_FULL_GRADIENT_EXPONENT = 1.0;
  ParamDB->LP_STREAMLINE_EXPONENT = 1.0;
  ParamDB->LP_DIVERGENCE_EXPONENT = 1.0;
  ParamDB->LP_PRESSURE_EXPONENT = 1.0;

  ParamDB->LP_ORDER_DIFFERENCE = 1;
  ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = -123;
  ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = -123;
  ParamDB->LP_DIVERGENCE_ORDER_DIFFERENCE = -123;
  ParamDB->LP_PRESSURE_ORDER_DIFFERENCE = -123;
  
  ParamDB->LP_CROSSWIND = 0;
  ParamDB->LP_CROSSWIND_COEFF_TYPE = 1;
  ParamDB->LP_CROSSWIND_COEFF = 1.0;
  ParamDB->LP_CROSSWIND_EXPONENT = 1.0;

  //======================================================================
  /** parameter for a posteriori parameter computation with adjoint problem */
  //======================================================================
  ParamDB->SOLVE_ADJOINT_PROBLEM = false; 
  ParamDB->SOLD_ADJOINT = 0;
  ParamDB->N_STAGES_ADJOINT = 1;
  ParamDB->SC_NONLIN_ITE_ADJOINT = 1000;
  ParamDB->OPTIMIZATION_ITE_TYPE_ADJOINT = 0;
  ParamDB->BFGS_VECTORS_ADJOINT = 25;
  ParamDB->RELATIVE_DECREASE_ADJOINT = 1e-4;
  ParamDB->PENALTY_ADJOINT = 0;
  ParamDB->PENALTY_VALUE_AT_ZERO_ADJOINT = 1e6;
  ParamDB->PENALTY_SMALLEST_PARAM_FAC_ADJOINT = 1e-2;
  ParamDB->PENALTY_LARGEST_PARAM_FAC_ADJOINT = 1e1;
  ParamDB->WEIGHT_RESIDUAL_L1_ADJOINT = 0;
  ParamDB->WEIGHT_RESIDUAL_L2_ADJOINT = 1;
  ParamDB->WEIGHT_GRADIENT_L1_ADJOINT = 0;
  ParamDB->WEIGHT_GRADIENT_L2_ADJOINT = 0;
  ParamDB->WEIGHT_STREAM_DER_L1_ADJOINT = 0;
  ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_ADJOINT = 0;
  ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT = 0;
  ParamDB->REG_POINT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT = 1.0;
  ParamDB->MIN_VAL_ADJOINT = -1e20;
  ParamDB->MAX_VAL_ADJOINT = 1e20;
  ParamDB->MIN_MAX_EXPONENT_ONE_ADJOINT = 1;
  ParamDB->MIN_MAX_EXPONENT_TWO_ADJOINT = 1;
  ParamDB->MIN_MAX_FACTOR_ONE_ADJOINT = 1;
  ParamDB->MIN_MAX_FACTOR_TWO_ADJOINT = 1;
  ParamDB->WEIGHT_RESIDUAL_LP_ADJOINT = 0;
  ParamDB->WEIGHT_RESIDUAL_EXP_LP_ADJOINT = 0;
  ParamDB->WEIGHT_RESIDUAL_CW_ADJOINT = 0;
  ParamDB->RESIDUAL_LP_ADJOINT = 2.0;

  ParamDB->MIN_MAX_ADJOINT = 0;
  ParamDB->INITIAL_STEEPEST_DESCENT_ADJOINT = 0;

  tmp = new char[40];
  strcpy(tmp,"MooN_MD_default_save_data_filename");
  ParamDB->SAVE_DATA_FILENAME=tmp;
  tmp = new char[40];
  strcpy(tmp,"MooN_MD_default_read_data_filename");
  ParamDB->READ_DATA_FILENAME=tmp;
  tmp = new char[40];
  strcpy(tmp, "NO_SMESH_FILE");
  ParamDB->SMESHFILE = tmp;
  tmp = new char[40];
  strcpy(tmp,"MooNMD_default_pod_filename");
  ParamDB->POD_FILENAME=tmp;
  //file for storing snapshots (ROM, reduced order modeling)
  tmp = new char[40];
  strcpy(tmp,"MooNMD_default_snap_filename");
  ParamDB->SNAP_FILENAME=tmp;


   /** parameters for SOLD schemes */
  ParamDB->SOLD_TYPE = 0;
  ParamDB->SOLD_PARAMETER_TYPE = 11;
  ParamDB->SOLD_CONST = 1.0;
  ParamDB->SOLD_POWER = 1.0;
  ParamDB->SOLD_S = 1.0;
  ParamDB->SOLD_U0 = 1.0;
  ParamDB->SOLD_PARAMETER_SCALING = 0;
  ParamDB->SOLD_PARAMETER_SCALING_FACTOR = 1.0;

  /** parameters for controling the program */
  ParamDB->WRITE_GRAPE = false; 
  ParamDB->WRITE_GMV = false; 
  ParamDB->WRITE_AMIRA = false; 
  ParamDB->WRITE_GNU = false; 
  ParamDB->WRITE_CASE = false; 
  ParamDB->SAVE_DATA = false; 
  ParamDB->READ_DATA = false; 
  ParamDB->READ_GRAPE_FILE = false; 
  ParamDB->ESTIMATE_ERRORS = false; 
  ParamDB->SOLVE_ADJOINT_PROBLEM = false; 
  ParamDB->COMPUTE_VORTICITY_DIVERGENCE = false;

  /** the following parameters are for individual use */
  ParamDB->P2 = 1.0;

  ParamDB->SC_NONLIN_ITE_ADJOINT = 1000;


  // ******** parameters for saddle point system *********//

  // parameters for nonlinear iteration
  ParamDB->SC_NONLIN_ITE_TYPE_SADDLE = 0;


    
  /** parameters for weakly imposing boundary/interface conditions */
  ParamDB->n_neumann_boundary = 0.;
  ParamDB->neumann_boundary_id.clear();
  ParamDB->neumann_boundary_value.clear();
    
  ParamDB-> n_g_v_boundary = 0.;
  ParamDB-> g_v_boundary_id.clear();
  ParamDB->g_v_boundary_value.clear();
    
  ParamDB-> n_unvn_boundary = 0.;
  ParamDB-> unvn_boundary_id.clear();
  ParamDB->unvn_boundary_value.clear();

  ParamDB-> n_gradunv_boundary = 0.;
  ParamDB-> gradunv_boundary_id.clear();
  ParamDB->gradunv_boundary_value.clear();
    
  ParamDB-> n_u_v_boundary = 0.;
  ParamDB-> u_v_boundary_id.clear();
  ParamDB->u_v_boundary_value.clear();
    
  ParamDB-> n_p_v_n_boundary = 0.;
  ParamDB-> p_v_n_boundary_id.clear();
  ParamDB-> p_v_n_boundary_value.clear();

    
  // Nitsche Combination - Weak Dirichlet Boundary Conditions 
  ParamDB-> n_nitsche_boundary = 0.;
  ParamDB-> nitsche_boundary_id.clear();
  ParamDB-> nitsche_penalty.clear();
  ParamDB-> s1 = 0;
  ParamDB-> s2 = 0;

  
  ParamDB->TETGEN_QUALITY = 0.0;
  ParamDB->TETGEN_VOLUMEN = 0.0;
  ParamDB->TETGEN_STEINER = 0;

  ParamDB->CHAR_L0=1.;
  ParamDB->D_VISCOSITY=1.0;
  ParamDB->SURF_TENSION=0.;
  ParamDB->IMPACT_ANGLE=90.;
  ParamDB->Area = 1.;

  // initialize TimeDB
  TimeDB->CURRENTTIME = 0;
  TimeDB->CURRENTTIMESTEPLENGTH = 1;
  TimeDB->TIMESTEPLENGTH = 1;
  TimeDB->MIN_TIMESTEPLENGTH = 1E-4;
  TimeDB->MAX_TIMESTEPLENGTH = 0.5;
  TimeDB->TIMESTEPLENGTH_TOL = 1e-3;
  TimeDB->TIMESTEPLENGTH_CONTROL = 0;
  TimeDB->TIMESTEPLENGTH_CONTROLLER = 0;  // mlh
  TimeDB->TIMESTEPLENGTH_PARA_KK_I = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_KK_P = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_KK_E = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_KK_R = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_KK_D = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_FAC = 0.8;
  TimeDB->TIMESTEPLENGTH_PARA_FAC_MAX = 5.0;
  TimeDB->TIMESTEPLENGTH_PARA_FAC_MIN = 0.2;
  TimeDB->TIMESTEPLENGTH_PARA_TOL = 1.0;
  TimeDB->TIMESTEPLENGTH_PARA_ATOL = 0.001;
  TimeDB->TIMESTEPLENGTH_PARA_RTOL = 0.001;
  TimeDB->RESET_CURRENTTIME = 0;
  TimeDB->STEADY_STATE_TOL = 1e-3;
  TimeDB->SCALE_DIVERGENCE_CONSTRAINT = -1.0;

  TimeDB->CONTROL=0;
  TimeDB->CONTROL_ALPHA=1;
  TimeDB->CONTROL_BETA=1;
  TimeDB->CONTROL_GAMMA=1;
  TimeDB->CONTROL_SAFTY=0.9;
  TimeDB->CONTROL_MINSCALE=0.1;
  TimeDB->CONTROL_MAXSCALE=5.0;
  
  // parameters for implicit Euler method
  TimeDB->THETA1 = 1;
  TimeDB->THETA2 = 0;
  TimeDB->THETA3 = 0;
  TimeDB->THETA4 = 1;
  TimeDB->TIME_DISC = 2;
  TimeDB->TIME_DISC2 = -1;

  TimeDB->EXTRAPOLATE_WEIGHT = 1;
  TimeDB->EXTRAPOLATE_STEPS = 0;
  TimeDB->EXTRAPOLATE_PRESSURE = 0;
  TimeDB->EXTRAPOLATE_VELOCITY = 0;

  TimeDB->T0 = 0;
  TimeDB->T1 = 0;
  TimeDB->T2 = 0;
  TimeDB->T3 = 0;
  TimeDB->T4 = 0;
  TimeDB->T5 = 0;
  TimeDB->T6 = 0;
  TimeDB->T7 = 0;
  TimeDB->T8 = 0;
  TimeDB->T9 = 0;

  TimeDB->STEPS_PER_IMAGE = 1;
  TimeDB->STEPS_PER_SNAP = 1;

  TimeDB->RB_TYPE = 3;
  TimeDB->RB_TYPE2 = -1;
  TimeDB->RB_SSC = 0;
  TimeDB->RB_SSC_TOL = 1;
  TimeDB->RB_SSC_ALPHA = 0.9;
  TimeDB->RB_SSC_ALPHA_MIN = 1.0;
  TimeDB->RB_SSC_ALPHA_MAX = 1.0;
  TimeDB->RB_SSC_MAX_ERROR = 1.0;

  TimeDB->RB_APPROX_J = 0;
  TimeDB->RB_APPROX_C = 0;
  
  // parameters for higher order Galerkin-type methods
  TimeDB->INTERNAL_SYSTEMSIZE = 0;
  TimeDB->INTERNAL_ALPHA = nullptr;
  TimeDB->INTERNAL_BETA = nullptr;
  TimeDB->VALUE_AT_ONE = nullptr;
  TimeDB->VAL_AT_QUAD_POINTS = nullptr;
  TimeDB->DER_AT_QUAD_POINTS = nullptr;
  TimeDB->CORR_AT_QUAD_POINTS = nullptr;
  TimeDB->DER_CORR_AT_QUAD_POINTS = nullptr;
  TimeDB->DER_AT_START = nullptr;
  TimeDB->DER_AT_ONE = nullptr;
  TimeDB->DER_COR_AT_ONE = nullptr;
  TimeDB->NORMW = 0;
  TimeDB->N_QUADPOINTS = 0;
  TimeDB->N_DEGREES = 0;
  TimeDB->ZETA = nullptr;
  TimeDB->WEIGHTS = nullptr;
  TimeDB->ALPHA0 = nullptr;
  TimeDB->BETA0 = nullptr;
  TimeDB->GAMMA0 = nullptr;
  TimeDB->CORRECTION = nullptr;
  TimeDB->POINTS = nullptr;

  ParamDB->INPUT_QUAD_RULE = 0;
  ParamDB->INTERNAL_PROBLEM_LINEAR = 0;
  ParamDB->INTERNAL_PRESSURE_SPACE = 0;
  ParamDB->INTERNAL_SLIP_WITH_FRICTION = 0;
  ParamDB->INTERNAL_SLIP_WITH_FRICTION_IDENTITY = 0;
  ParamDB->INTERNAL_QUAD_HEXA = 0;
  ParamDB->INTERNAL_QUAD_TETRA = 0;
  ParamDB->INTERNAL_QUAD_QUAD = 0;
  ParamDB->INTERNAL_QUAD_TRIA = 0;
  ParamDB->INTERNAL_QUAD_RULE = 0;
  ParamDB->INTERNAL_LOCAL_DOF = 0;
  ParamDB->INTERNAL_PERIODIC_IDENTITY = 0;
  ParamDB->INTERNAL_PROBLEM_IDENTITY = 0;
  ParamDB->INTERNAL_LEVEL = 0;
  ParamDB->INTERNAL_CONVECTION_EQ_VELOFIELD = 0;
  ParamDB->INTERNAL_STEADY_STATE_MATRICES_OR_RHS = 0;
  ParamDB->INTERNAL_AMG_SOLVES = 0;
  ParamDB->INTERNAL_AMG_PREPARE_TIME = 0.0;
  ParamDB->INTERNAL_GMRES_INFO = 1;
  ParamDB->INTERNAL_POLYNOMIAL_DEGREE = 1;
  ParamDB->INTERNAL_LINEAR_SCHEME = 1;
  ParamDB->INTERNAL_SOLD_ACTIVE = 0;
  ParamDB->INTERNAL_SORT_AMG = 1;
  ParamDB->INTERNAL_COERCIVITY = -4711.0;
  ParamDB->INTERNAL_FACE_INTEGRALS = 0;
  ParamDB->INTERNAL_NO_ESTIMATE_DIRICHLET_CELLS = 0;
  ParamDB->INTERNAL_WRONG_NEUMANN_CHECKED = 0;
  ParamDB->INTERNAL_BFGS_RESTART_ADJOINT = 0;
  ParamDB->INTERNAL_NEW_MATRICES_B = 1;
  ParamDB->INTERNAL_FULL_MATRIX_STRUCTURE = 0;
  ParamDB->INTERNAL_DISC_FLAG = 0;

  ParamDB->INTERNAL_FESPACE_CONSTRUCT = 0;
  ParamDB->INTERNAL_DO_NOT_RESPECT_DIRICHLET_BC = 0;
  ParamDB->INTERNAL_P1_Array = nullptr;
  ParamDB->INTERNAL_START_PARAM = 0;

  /** parameters for free surface calculation */
  ParamDB->FS_MAGNETLAW = 0;

  ParamDB->FS_L = 1;
  ParamDB->FS_U = 1;
  ParamDB->FS_T = 1;

  ParamDB->FS_ETA = 1;
  ParamDB->FS_RHO = 1;
  ParamDB->FS_ALPHA = 1;
  ParamDB->FS_G = 1;

  ParamDB->FS_MS = 1;
  ParamDB->FS_CHI0 = 1;

  ParamDB->FS_HM = 1;
  ParamDB->FS_DELTA_H = 0;
  ParamDB->FS_F = 0;
  
  ParamDB->FS_H0 = 1;
  ParamDB->FS_H1 = 1;
  ParamDB->FS_H2 = 1;
  
  ParamDB->FS_LH = 1;
  ParamDB->FS_GAMMA = 1;
  ParamDB->FS_HT = 1;

  ParamDB->FS_WE = 1;

  ParamDB->FS_WRITE = 0;
  ParamDB->FS_READ = 0;

  tmp = new char[12];
  strcpy(tmp,"FS_INNAME");
  ParamDB->FS_INNAME = tmp;

  tmp = new char[12];
  strcpy(tmp,"FS_OUTNAME");
  ParamDB->FS_OUTNAME = tmp;

  ParamDB->HEAT_TANGENTIAL_STRESS_FACTOR = 0;  
  ParamDB->HEAT_SOLID_SURFACE_FACTOR = 1.;    
  ParamDB->EQ_CONTACT_ANGLE = 0;   
  ParamDB->AD_CONTACT_ANGLE = 0;   
  ParamDB->RE_CONTACT_ANGLE = 0;  
  ParamDB->DY_CONTACT_ANGLE = 0;     
  ParamDB->CONTACT_ANGLE_TYPE = 0;  
   
  // ******** parameters for VMS *********//
  ParamDB->VMS_LARGE_VELOCITY_SPACE = 0;
  ParamDB->VMS_COARSE_MG_SMAGO = 1;
  // constants in AdaptProjectionSpace 
  ParamDB->VMS_ADAPT_LOWER = 0.25;
  ParamDB->VMS_ADAPT_MIDDLE = 1.0;
  ParamDB->VMS_ADAPT_UPPER = 3.0;
  ParamDB->VMS_ADAPT_STEPS = 1;
  ParamDB->VMS_ADAPT_COMP = 1;

  ParamDB->SUPERCONVERGENCE_ORDER = 0;
  ParamDB->WENO_TYPE = 0;

  /* the following parameters are for membrane REACTOR */
  ParamDB->REACTOR_P0 = 0.0;
  ParamDB->REACTOR_P1 = 0.0;
  ParamDB->REACTOR_P2 = 0.0;
  ParamDB->REACTOR_P3 = 0.0;
  ParamDB->REACTOR_P4 = 0.0;
  ParamDB->REACTOR_P5 = 0.0;
  ParamDB->REACTOR_P6 = 0.0;
  ParamDB->REACTOR_P7 = 0.0;
  ParamDB->REACTOR_P8 = 0.0;
  ParamDB->REACTOR_P9 = 0.0;
  ParamDB->REACTOR_P10 = 0.0;
  ParamDB->REACTOR_P11 = 0.0;
  ParamDB->REACTOR_P12 = 0.0;
  ParamDB->REACTOR_P13 = 0.0;
  ParamDB->REACTOR_P14 = 0.0;
  ParamDB->REACTOR_P15 = 0.0;
  ParamDB->REACTOR_P16 = 0.0;
  ParamDB->REACTOR_P17 = 0.0;
  ParamDB->REACTOR_P18 = 0.0;
  ParamDB->REACTOR_P19 = 0.0;
  ParamDB->REACTOR_P20 = 0.0;
  ParamDB->REACTOR_P21 = 0.0;
  ParamDB->REACTOR_P22 = 0.0;
  ParamDB->REACTOR_P23 = 0.0;
  ParamDB->REACTOR_P24 = 0.0;
  ParamDB->REACTOR_P25 = 0.0;
  ParamDB->REACTOR_P26 = 0.0;
  ParamDB->REACTOR_P27 = 0.0;
  ParamDB->REACTOR_P28 = 0.0;
  ParamDB->REACTOR_P29 = 0.0;
  ParamDB->REACTOR_P30 = 0.0;

// parameters for turbulent channel flow
  ParamDB->CHANNEL_STATISTICS2_WITH_MODEL = 0;

// parameters for turbulent flow around a squared cylinder
  ParamDB->CYLINDER_22000_YPLUS_SIDES = 1500;
  ParamDB->CYLINDER_22000_YPLUS_FRONT = 2300;
  ParamDB->CYLINDER_22000_YPLUS_BACK  = 1000;
  
// parameters for BULK computations
  ParamDB->BULK_REACTION_DISC = 0;
  ParamDB->BULK_PB_DISC = 0;
  ParamDB->BULK_PB_DISC_STAB = 0;
  ParamDB->BULK_PB_DISC_FCT_GROUP = 0;
  ParamDB->BULK_COUPLING = 0;
  ParamDB->BULK_GROWTH_RATE = 0;
  ParamDB->BULK_REACTION_MASS_LUMPING = 0;
  ParamDB->BULK_METHODS_OF_MOMENTS = 0;
  ParamDB->BULK_MOM_DISC = 0;
  ParamDB->BULK_SOLD_PARAMETER_TYPE = 0;
  ParamDB->N_CELL_LAYERS_PSD = 0;
  ParamDB->N_CELL_LAYERS_PSD_2 = 0;
  ParamDB->OUTPUT_NODE_LAYER_PSD = 0;
  ParamDB->BULK_REACTION_C_CUT = 0.0;

  ParamDB->BULK_l_infty = 0.0;
  ParamDB->BULK_u_infty = 0.0;
  ParamDB->BULK_c_infty = 0.0;
  ParamDB->BULK_c_C_infty_sat = 0.0;
  ParamDB->BULK_c_C_infty = 0.0;
  ParamDB->BULK_f_infty = 0.0;

  ParamDB->BULK_density = 0.0;
  ParamDB->BULK_dynamic_viscosity = 0.0;

  ParamDB->BULK_C_g = 0.0;
  ParamDB->BULK_C_nuc = 0.0;
  ParamDB->BULK_C_sat = 0.0;
  ParamDB->BULK_C_2 = 0.0;
  ParamDB->BULK_D_A = 0.0;
  ParamDB->BULK_D_P_0 = 0.0;
  ParamDB->BULK_D_P_MAX = 0.0;
  ParamDB->BULK_k_g = 0.0;
  ParamDB->BULK_k_r = 0.0;
  ParamDB->BULK_k_nuc = 0.0;
  ParamDB->BULK_D_P_MIN = 0.0;

// parameters for shear slip mesh update method computations
  ParamDB->SSMUM_MP_X = 0.5;
  ParamDB->SSMUM_MP_Y = 0.5;
  ParamDB->SSMUM_INNER_RADIUS = 0.25;
  ParamDB->SSMUM_OUTER_RADIUS = 0.35;
  ParamDB->SSMUM_ROT_PER_SECOND = 0.5;
  ParamDB->SSMUM_ANGLE = 0.0;
  ParamDB->SSMUM_MAX_CELLS_LAYERS = 1024;
  ParamDB->SSMUM_INTERPOLATION = 0;


  /** general parameters for population balances */
  ParamDB->PB_DISC_TYPE = 3;
  ParamDB->PB_TIME_DISC = 100;

  
  /** parameters for matlab output */
  tmp = new char[40];
  strcpy(tmp,"MooN_MD_default_matlab_matrix");
  ParamDB->MATLAB_MATRIX=tmp;
  
  ParamDB->WRITE_MATLAB = false;
  ParamDB->WRITE_MATLAB_MATRIX = false;

  /** parameter for non-conforming elements
      number corresponds to Apel-Matthies paper 2006 */
  ParamDB->NC_TYPE = 3;
  
  /** parameters for ROM*/
  ParamDB->WRITE_SNAPSHOTS = false;
  ParamDB->DO_ROM = false;
  ParamDB->DO_ROM_P = false;
  ParamDB->RANK_OF_BASIS = 0;
  ParamDB->RANK_OF_BASIS_P = 0;
  ParamDB->POD_INNER_PRODUCT = 1;
  ParamDB->POD_INNER_PRODUCT_P = 1;
  
  ParamDB->BUILD_PODFILE = false;
  ParamDB->POD_FLUCT_FIELD = false;
  ParamDB->POD_FLUCT_FIELD_P = false;
  ParamDB->P_ROM_METHOD = 1;

  tmp = new char[12];
  strcpy(tmp,"NO_PODFILE");
  ParamDB->PODFILE= tmp;
  
  /** parameter for type of projection method (NSE)**/
  ParamDB->PROJECTION_METHOD = 1;

  /** parameters for individual use */
  ParamDB->P0 = 0.;
  ParamDB->P1 = 0.;
  ParamDB->P2 = 0.;
  ParamDB->P3 = 0.;
  ParamDB->P4 = 0.;
  ParamDB->P5 = 0.;
  ParamDB->P6 = 0.;
  ParamDB->P7 = 0.;
  ParamDB->P8 = 0.;
  ParamDB->P9 = 0.;
  ParamDB->P10 = 0.;
  ParamDB->P11 = 0.;
  ParamDB->P12 = 0.;
  ParamDB->P13 = 0.;
  ParamDB->P14 = 0.;
  ParamDB->P15 = 0.;
  
  
  /** parameters for individual use in parallel computations */
  ParamDB->Par_P0 = 0;
  ParamDB->Par_P1 = 0;
  ParamDB->Par_P2 = 0;
  ParamDB->Par_P3 = 1;
  ParamDB->Par_P4 = 0;
  ParamDB->Par_P5 = 10;
  ParamDB->Par_P6 = 1;
  ParamDB->Par_P7 = 0;
  ParamDB->Par_P8 = 0;
  ParamDB->Par_P9 = 0;
  ParamDB->Par_P10 = -1.0;
  ParamDB->Par_P11 = 0;
  ParamDB->Par_P12 = 0;
  ParamDB->Par_P13 = 0;
  ParamDB->Par_P14 = 0;
  ParamDB->Par_P15 = 0;
  ParamDB->Par_P16 = 0;
  ParamDB->Par_P17 = 0;
  ParamDB->Par_P18 = 0;
  ParamDB->Par_P19 = 0;
  ParamDB->Par_P20 = 0;
  
  /** parameters for population balance computations */
  ParamDB->PBE_P0 = 0;
  ParamDB->PBE_P1 = 0;
  ParamDB->PBE_P2 = 0;
  ParamDB->PBE_P3 = 1;
  ParamDB->PBE_P4 = 0;
  ParamDB->PBE_P5 = 0;
  ParamDB->PBE_P6 = 0;
  ParamDB->PBE_P7 = 0;
  ParamDB->PBE_P8 = 0;
  ParamDB->PBE_P9 = 0;
  
   /** parameters for population balance computations */
  ParamDB->DG_P0 = 0.;
  ParamDB->DG_P1 = 0.;
  ParamDB->DG_P2 = 0.;
  ParamDB->DG_P3 = 0.;
  ParamDB->DG_P4 = 0.;
  ParamDB->DG_P5 = 0.;
  ParamDB->DG_P6 = 0.;
  ParamDB->DG_P7 = 0.;
  ParamDB->DG_P8 = 0.;
  ParamDB->DG_P9 = 0.;
  
  /** parameters for moving domains */
  ParamDB->MOVING_BOUNDARY = 0;

  
  /** parameters for moving domains */
  ParamDB->DEPENDENT_BASIS = 0;
  ParamDB->DEPENDENT_BASIS_Q1 = 0;
  ParamDB->DEPENDENT_BASIS_Q2 = 0;
 
    
  #ifdef _MPI
  ParamDB->Comm = MPI_COMM_WORLD;    
  #endif
  return;
}

void TDatabase::WriteParamDB(char *ExecutedFile)
{
  using namespace Output; // printToFile
  
  printToFile(">>> Start printing old database (ParamDB) <<<");
  printToFile("HOSTNAME: ", utilities::get_host_name(), " started on ", 
              utilities::get_date_and_time());
  printToFile("EXECUTED FILE: ", ExecutedFile);
  printToFile("VERSION: ", ParamDB->VERSION);
  printToFile("MAPFILE: ", ParamDB->MAPFILE);
  printToFile("OUTFILE: ", ParamDB->OUTFILE);
  printToFile("profiling: ", ParamDB->timeprofiling);
  printToFile("MapperType: ", ParamDB->MapperType);
  printToFile("DSType: ", ParamDB->DSType);
  
  printToFile("MESHGEN_ALLOW_EDGE_REF: ", ParamDB->MESHGEN_ALLOW_EDGE_REF);
  printToFile("MESHGEN_REF_QUALITY: ", ParamDB->MESHGEN_REF_QUALITY);

  printToFile("ANSATZ_ORDER: ", ParamDB->ANSATZ_ORDER);
  printToFile("TEST_ORDER: ", ParamDB->TEST_ORDER);
  
  printToFile("VELOCITY_SPACE: ", ParamDB->VELOCITY_SPACE);
  printToFile("PRESSURE_SPACE: ", ParamDB->PRESSURE_SPACE);
  printToFile("PRESSURE_SEPARATION: ", ParamDB->PRESSURE_SEPARATION);

  printToFile("OMPNUMTHREADS: ", ParamDB->OMPNUMTHREADS);

  printToFile("REFINEMENT: ", ParamDB->REFINEMENT);
  printToFile("GRID_TYPE: ", ParamDB->GRID_TYPE);
  printToFile("GRID_TYPE_1: ", ParamDB->GRID_TYPE_1);
  printToFile("GRID_TYPE_2: ", ParamDB->GRID_TYPE_2);
  printToFile("CHANNEL_GRID_STRETCH: ", ParamDB->CHANNEL_GRID_STRETCH);

  printToFile("ADAPTIVE_REFINEMENT_CRITERION: ", ParamDB->ADAPTIVE_REFINEMENT_CRITERION);
  printToFile("ERROR_CONTROL: ", ParamDB->ERROR_CONTROL);
  printToFile("REFINE_STRATEGY: ", ParamDB->REFINE_STRATEGY);
  printToFile("MAX_CELL_LEVEL: ", ParamDB->MAX_CELL_LEVEL);
  printToFile("REFTOL: ", ParamDB->REFTOL);
  printToFile("COARSETOL: ", ParamDB->COARSETOL);  
  printToFile("MIN_FRACTION_TO_CHANGE: ", ParamDB->MIN_FRACTION_TO_CHANGE);
  printToFile("DECREASE_REFTOL_FACTOR: ", ParamDB->DECREASE_REFTOL_FACTOR);
  printToFile("INCREASE_COARSETOL_FACTOR: ", ParamDB->INCREASE_COARSETOL_FACTOR);
  printToFile("FRACTION_OF_ERROR: ", ParamDB->FRACTION_OF_ERROR);
  printToFile("CONVERT_QUAD_TO_TRI: ", ParamDB->CONVERT_QUAD_TO_TRI);
  
  printToFile("INTL_DISCTYPE: ", ParamDB->INTL_DISCTYPE);
  printToFile("UPWIND_ORDER: ", ParamDB->UPWIND_ORDER);
  printToFile("UPWIND_FLUX_DAMP: ", ParamDB->UPWIND_FLUX_DAMP);
  printToFile("UPWIND_APPLICATION: ", ParamDB->UPWIND_APPLICATION);
  printToFile("SHISHKIN_MESH: ", ParamDB->SHISHKIN_MESH);
  printToFile("SHISHKIN_DIAM: ", ParamDB->SHISHKIN_DIAM);
  printToFile("NSTYPE: ", ParamDB->NSTYPE);
  printToFile("DARCYTYPE: ", ParamDB->DARCYTYPE);
  printToFile("SIGMA_PERM: ", ParamDB->SIGMA_PERM);
  printToFile("LAPLACETYPE: ", ParamDB->LAPLACETYPE);
  printToFile("USE_ISOPARAMETRIC: ", ParamDB->USE_ISOPARAMETRIC);
  printToFile("VMM_COARSE_LEVEL: ", ParamDB->VMM_COARSE_LEVEL);
  printToFile("VMM_COARSE_SPACE_ORDER: ", ParamDB->VMM_COARSE_SPACE_ORDER);
  printToFile("RFB_SUBMESH_LAYERS: ", ParamDB->RFB_SUBMESH_LAYERS);
  printToFile("DEFECT_CORRECTION_TYPE: ", ParamDB->DEFECT_CORRECTION_TYPE);
  printToFile("CELL_MEASURE: ", ParamDB->CELL_MEASURE);
  printToFile("FACE_SIGMA: ", ParamDB->FACE_SIGMA);
  printToFile("WEAK_BC_SIGMA: ", ParamDB->WEAK_BC_SIGMA);
  printToFile("WEAK_BC: ", ParamDB->WEAK_BC);

  printToFile("EQUAL_ORDER_STAB_WEIGHT_PkPk: ", ParamDB->equal_order_stab_weight_PkPk); 
  printToFile("GRAD_DIV_STAB_WEIGHT: ", ParamDB->grad_div_stab_weight); 
  printToFile("SIGN_MATRIX_BI: ", ParamDB->SIGN_MATRIX_BI);
  printToFile("l_T: ", ParamDB->l_T);
  printToFile("L_0: ", ParamDB->L_0); 
  printToFile("SOURCE_SINK_FUNCTION: ", ParamDB->SOURCE_SINK_FUNCTION);

  printToFile("RE_NR: ", ParamDB->RE_NR);
  printToFile("RA_NR: ", ParamDB->RA_NR);
  printToFile("ROSSBY_NR: ", ParamDB->ROSSBY_NR);
  printToFile("START_RE_NR: ", ParamDB->START_RE_NR);
  printToFile("RE_NR_INCREMENT: ", ParamDB->RE_NR_INCREMENT);
  printToFile("FLOW_PROBLEM_TYPE: ", ParamDB->FLOW_PROBLEM_TYPE);
  printToFile("FR_NR: ", ParamDB->FR_NR);
  printToFile("WB_NR: ", ParamDB->WB_NR);
  printToFile("PR_NR: ", ParamDB->PR_NR);
  printToFile("PE_NR: ", ParamDB->PE_NR);  
  printToFile("BI_NR: ", ParamDB->BI_NR);
  printToFile("Axial3D: ", ParamDB->Axial3D);
  printToFile("Axial3DAxis: ", ParamDB->Axial3DAxis);  
  printToFile("DELTA0: ", ParamDB->DELTA0);
  printToFile("SDFEM_POWER0: ", ParamDB->SDFEM_POWER0);
  printToFile("SDFEM_TYPE: ", ParamDB->SDFEM_TYPE);
  printToFile("SDFEM_NORM_B: ", ParamDB->SDFEM_NORM_B);
  printToFile("DELTA1: ", ParamDB->DELTA1);
  printToFile("CIP_TYPE: ", ParamDB->CIP_TYPE);
  printToFile("ADJOINT_FACTOR_4_OMEGA_EQ_0: ", ParamDB->ADJOINT_FACTOR_4_OMEGA_EQ_0);

  printToFile("SOLD_TYPE: ", ParamDB->SOLD_TYPE);
  printToFile("SOLD_PARAMETER_TYPE: ", ParamDB->SOLD_PARAMETER_TYPE);
  printToFile("SOLD_CONST: ", ParamDB->SOLD_CONST);
  printToFile("SOLD_POWER: ", ParamDB->SOLD_POWER);
  printToFile("SOLD_S: ", ParamDB->SOLD_S);
  printToFile("SOLD_U0: ", ParamDB->SOLD_U0);
  printToFile("SOLD_PARAMETER_SCALING: ", ParamDB->SOLD_PARAMETER_SCALING);
  printToFile("SOLD_PARAMETER_SCALING_FACTOR: ", ParamDB->SOLD_PARAMETER_SCALING_FACTOR);
  

  printToFile("FILTER_WIDTH_CONSTANT: ", ParamDB->FILTER_WIDTH_CONSTANT);
  printToFile("FILTER_WIDTH_POWER: ", ParamDB->FILTER_WIDTH_POWER);
  printToFile("GAUSSIAN_GAMMA: ", ParamDB->GAUSSIAN_GAMMA );
  printToFile("CONVOLUTE_SOLUTION: ", ParamDB->CONVOLUTE_SOLUTION);

  printToFile("TURBULENT_VISCOSITY_TYPE: ", ParamDB->TURBULENT_VISCOSITY_TYPE);
  printToFile("TURBULENT_VISCOSITY_TENSOR: ", ParamDB->TURBULENT_VISCOSITY_TENSOR);
  printToFile("TURBULENT_VISCOSITY_CONSTANT: ", ParamDB->TURBULENT_VISCOSITY_CONSTANT);
  printToFile("TURBULENT_VISCOSITY_POWER: ", ParamDB->TURBULENT_VISCOSITY_POWER);
  printToFile("TURBULENT_VISCOSITY_SIGMA: ", ParamDB->TURBULENT_VISCOSITY_SIGMA);

  printToFile("ARTIFICIAL_VISCOSITY_CONSTANT: ", ParamDB->ARTIFICIAL_VISCOSITY_CONSTANT);
  printToFile("ARTIFICIAL_VISCOSITY_POWER: ", ParamDB->ARTIFICIAL_VISCOSITY_POWER);

  printToFile("FRICTION_CONSTANT: ", ParamDB->FRICTION_CONSTANT);
  printToFile("FRICTION_POWER: ", ParamDB->FRICTION_POWER);
  printToFile("FRICTION_TYPE: ", ParamDB->FRICTION_TYPE);
  printToFile("FRICTION_U0: ", ParamDB->FRICTION_U0);
  printToFile("PENETRATION_CONSTANT: ", ParamDB->PENETRATION_CONSTANT);
  printToFile("PENETRATION_POWER: ", ParamDB->PENETRATION_POWER);
  printToFile("DIV_DIV_STAB_TYPE: ", ParamDB->DIV_DIV_STAB_TYPE); 
  printToFile("DIV_DIV_STAB_C1: ", ParamDB->DIV_DIV_STAB_C1); 
  printToFile("DIV_DIV_STAB_C2: ", ParamDB->DIV_DIV_STAB_C2); 
  printToFile("OSEEN_ZERO_ORDER_COEFF: ", ParamDB->OSEEN_ZERO_ORDER_COEFF);

  printToFile("LP_FULL_GRADIENT: ", ParamDB->LP_FULL_GRADIENT);
  printToFile("LP_STREAMLINE: ", ParamDB->LP_STREAMLINE);
  printToFile("LP_DIVERGENCE: ", ParamDB->LP_DIVERGENCE);
  printToFile("LP_PRESSURE: ", ParamDB->LP_PRESSURE);
  printToFile("LP_COEFF_TYPE: ", ParamDB->LP_COEFF_TYPE);
  printToFile("LP_FULL_GRADIENT_COEFF: ", ParamDB->LP_FULL_GRADIENT_COEFF);
  printToFile("LP_STREAMLINE_COEFF: ", ParamDB->LP_STREAMLINE_COEFF);
  printToFile("LP_DIVERGENCE_COEFF: ", ParamDB->LP_DIVERGENCE_COEFF);
  printToFile("LP_PRESSURE_COEFF: ", ParamDB->LP_PRESSURE_COEFF);
  printToFile("LP_FULL_GRADIENT_EXPONENT: ", ParamDB->LP_FULL_GRADIENT_EXPONENT);
  printToFile("LP_STREAMLINE_EXPONENT: ", ParamDB->LP_STREAMLINE_EXPONENT);
  printToFile("LP_DIVERGENCE_EXPONENT: ", ParamDB->LP_DIVERGENCE_EXPONENT);
  printToFile("LP_PRESSURE_EXPONENT: ", ParamDB->LP_PRESSURE_EXPONENT);
  printToFile("LP_ORDER_DIFFERENCE: ", ParamDB->LP_ORDER_DIFFERENCE);
  printToFile("LP_FULL_GRADIENT_ORDER_DIFFERENCE: ", ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE);
  printToFile("LP_STREAMLINE_ORDER_DIFFERENCE: ", ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE);
  printToFile("LP_DIVERGENCE_ORDER_DIFFERENCE: ", ParamDB->LP_DIVERGENCE_ORDER_DIFFERENCE);
  printToFile("LP_PRESSURE_ORDER_DIFFERENCE: ", ParamDB->LP_PRESSURE_ORDER_DIFFERENCE);

  printToFile("SOLVE_ADJOINT_PROBLEM: ", ParamDB->SOLVE_ADJOINT_PROBLEM);
  printToFile("SOLD_ADJOINT: ", ParamDB->SOLD_ADJOINT);
  printToFile("N_STAGES_ADJOINT: ", ParamDB->N_STAGES_ADJOINT);
  printToFile("SC_NONLIN_ITE_ADJOINT: ", ParamDB->SC_NONLIN_ITE_ADJOINT);
  printToFile("OPTIMIZATION_ITE_TYPE_ADJOINT: ", ParamDB->OPTIMIZATION_ITE_TYPE_ADJOINT);
  printToFile("BFGS_VECTORS_ADJOINT: ", ParamDB->BFGS_VECTORS_ADJOINT);
  printToFile("RELATIVE_DECREASE_ADJOINT: ", ParamDB->RELATIVE_DECREASE_ADJOINT);
  printToFile("PENALTY_ADJOINT: ", ParamDB->PENALTY_ADJOINT);
  printToFile("PENALTY_VALUE_AT_ZERO_ADJOINT: ", ParamDB->PENALTY_VALUE_AT_ZERO_ADJOINT);
  printToFile("PENALTY_SMALLEST_PARAM_FAC_ADJOINT: ", ParamDB->PENALTY_SMALLEST_PARAM_FAC_ADJOINT);
  printToFile("PENALTY_LARGEST_PARAM_FAC_ADJOINT: ", ParamDB->PENALTY_LARGEST_PARAM_FAC_ADJOINT);
  printToFile("WEIGHT_RESIDUAL_L1_ADJOINT: ", ParamDB->WEIGHT_RESIDUAL_L1_ADJOINT);
  printToFile("WEIGHT_RESIDUAL_L2_ADJOINT: ", ParamDB->WEIGHT_RESIDUAL_L2_ADJOINT);
  printToFile("WEIGHT_GRADIENT_L1_ADJOINT: ", ParamDB->WEIGHT_GRADIENT_L1_ADJOINT);
  printToFile("WEIGHT_GRADIENT_L2_ADJOINT: ", ParamDB->WEIGHT_GRADIENT_L2_ADJOINT);
  printToFile("WEIGHT_STREAM_DER_L1_ADJOINT: ", ParamDB->WEIGHT_STREAM_DER_L1_ADJOINT);
  printToFile("WEIGHT_STREAM_DER_ORTHO_L1_ADJOINT: ", ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_ADJOINT);
  printToFile("WEIGHT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT: ", ParamDB->WEIGHT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT);
  printToFile("REG_POINT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT: ", ParamDB->REG_POINT_STREAM_DER_ORTHO_L1_SQRT_ADJOINT);
  printToFile("WEIGHT_RESIDUAL_LP_ADJONT: ", ParamDB->WEIGHT_RESIDUAL_LP_ADJOINT);
  printToFile("WEIGHT_RESIDUAL_EXP_LP_ADJONT: ", ParamDB->WEIGHT_RESIDUAL_EXP_LP_ADJOINT);
  printToFile("WEIGHT_RESIDUAL_CW_ADJOINT: ", ParamDB->WEIGHT_RESIDUAL_CW_ADJOINT);
  printToFile("RESIDUAL_LP_ADJONT: ", ParamDB->RESIDUAL_LP_ADJOINT);
  printToFile("MIN_VAL_ADJOINT: ", ParamDB->MIN_VAL_ADJOINT);
  printToFile("MAX_VAL_ADJOINT: ", ParamDB->MAX_VAL_ADJOINT);
  printToFile("MIN_MAX_EXPONENT_ONE_ADJOINT: ", ParamDB->MIN_MAX_EXPONENT_ONE_ADJOINT);
  printToFile("MIN_MAX_EXPONENT_TWO_ADJOINT: ", ParamDB->MIN_MAX_EXPONENT_TWO_ADJOINT);
  printToFile("MIN_MAX_FACTOR_ONE_ADJOINT: ", ParamDB->MIN_MAX_FACTOR_ONE_ADJOINT);
  printToFile("MIN_MAX_FACTOR_TWO_ADJOINT: ", ParamDB->MIN_MAX_FACTOR_TWO_ADJOINT);
  printToFile("MIN_MAX_ADJOINT: ", ParamDB->MIN_MAX_ADJOINT);
  
  printToFile("SAVE_DATA_FILENAME: ", ParamDB->SAVE_DATA_FILENAME);
  printToFile("READ_DATA_FILENAME: ", ParamDB->READ_DATA_FILENAME);
  printToFile("POD_FILENAME: ", ParamDB->POD_FILENAME);
  printToFile("SNAP_FILENAME: ", ParamDB->SNAP_FILENAME);

  printToFile("WRITE_GRAPE: ", ParamDB->WRITE_GRAPE);
  printToFile("WRITE_GMV: ", ParamDB->WRITE_GMV);
  printToFile("WRITE_AMIRA: ", ParamDB->WRITE_AMIRA);
  printToFile("WRITE_GNU: ", ParamDB->WRITE_GNU);
  printToFile("WRITE_AMIRA: ", ParamDB->WRITE_AMIRA);
  printToFile("WRITE_CASE: ", ParamDB->WRITE_CASE);
  printToFile("WRITE_SNAPSHOTS: ", ParamDB->WRITE_SNAPSHOTS);
  printToFile("WRITE_MATLAB_MATRIX: ", ParamDB->WRITE_MATLAB_MATRIX); 
  printToFile("WRITE_MATLAB: ", ParamDB->WRITE_MATLAB);
  printToFile("SAVE_DATA: ", ParamDB->SAVE_DATA);
  printToFile("READ_DATA: ", ParamDB->READ_DATA);
  printToFile("ESTIMATE_ERRORS: ", ParamDB->ESTIMATE_ERRORS);
  printToFile("SOLVE_ADJOINT_PROBLEM: ", ParamDB->SOLVE_ADJOINT_PROBLEM);
  printToFile("COMPUTE_VORTICITY_DIVERGENCE: ", ParamDB->COMPUTE_VORTICITY_DIVERGENCE);
  printToFile("READ_GRAPE_FILE: ", ParamDB->READ_GRAPE_FILE );
  printToFile("BUILD_PODFILE:", ParamDB->BUILD_PODFILE);
  printToFile("POD_FLUCT_FIELD:", ParamDB->POD_FLUCT_FIELD);
  printToFile("POD_FLUCT_FIELD_P:", ParamDB->POD_FLUCT_FIELD_P);
  printToFile("PROJECTION_METHOD:", ParamDB->PROJECTION_METHOD);

  printToFile("P0: ", ParamDB->P0);
  printToFile("P1: ", ParamDB->P1);
  printToFile("P2: ", ParamDB->P2);
  printToFile("P3: ", ParamDB->P3);
  printToFile("P4: ", ParamDB->P4);
  printToFile("P5: ", ParamDB->P5);
  printToFile("P6: ", ParamDB->P6);
  printToFile("P7: ", ParamDB->P7);
  printToFile("P8: ", ParamDB->P8);
  printToFile("P9: ", ParamDB->P9);
  printToFile("P10: ", ParamDB->P10);
  printToFile("P11: ", ParamDB->P11);
  printToFile("P12: ", ParamDB->P12);
  printToFile("P13: ", ParamDB->P13);
  printToFile("P14: ", ParamDB->P14);
  printToFile("P15: ", ParamDB->P15);

  printToFile("PBE_P0: ", ParamDB->PBE_P0);
  printToFile("PBE_P1: ", ParamDB->PBE_P1);
  printToFile("PBE_P2: ", ParamDB->PBE_P2);
  printToFile("PBE_P3: ", ParamDB->PBE_P3);
  printToFile("PBE_P4: ", ParamDB->PBE_P4);
  printToFile("PBE_P5: ", ParamDB->PBE_P5);
  printToFile("PBE_P6: ", ParamDB->PBE_P6);
  printToFile("PBE_P7: ", ParamDB->PBE_P7);
  printToFile("PBE_P8: ", ParamDB->PBE_P8);
  printToFile("PBE_P9: ", ParamDB->PBE_P9);

  printToFile("DG_P0: ", ParamDB->DG_P0);
  printToFile("DG_P1: ", ParamDB->DG_P1);
  printToFile("DG_P2: ", ParamDB->DG_P2);
  printToFile("DG_P3: ", ParamDB->DG_P3);
  printToFile("DG_P4: ", ParamDB->DG_P4);
  printToFile("DG_P5: ", ParamDB->DG_P5);
  printToFile("DG_P6: ", ParamDB->DG_P6);
  printToFile("DG_P7: ", ParamDB->DG_P7);
  printToFile("DG_P8: ", ParamDB->DG_P8);
  printToFile("DG_P9: ", ParamDB->DG_P9);

  
  printToFile("VORTICITY THICKNESS FOR MIXING LAYER (P8): ", ParamDB->P8);


  printToFile("SC_NONLIN_ITE_ADJOINT: ", ParamDB->SC_NONLIN_ITE_ADJOINT);

  printToFile("*********** PARAMETERS FOR SADDLE POINT SOLVER ***********");
  printToFile("SC_NONLIN_ITE_TYPE_SADDLE: ", ParamDB->SC_NONLIN_ITE_TYPE_SADDLE);

  printToFile("CHAR_L0: ", ParamDB->CHAR_L0);
  printToFile("D_VISCOSITY: ", ParamDB->D_VISCOSITY);
  printToFile("SURF_TENSION: ", ParamDB->SURF_TENSION);
  printToFile("IMPACT_ANGLE: ", ParamDB->IMPACT_ANGLE);
  printToFile("Area: ", ParamDB->Area);

  // ******** parameters for VMS *********//
  printToFile("VMS_LARGE_VELOCITY_SPACE: ",  ParamDB->VMS_LARGE_VELOCITY_SPACE);
  printToFile("VMS_COARSE_MG_SMAGO: ",  ParamDB->VMS_COARSE_MG_SMAGO);
  // constants in AdaptProjectionSpace 
  printToFile("VMS_ADAPT_LOWER: ", ParamDB->VMS_ADAPT_LOWER); 
  printToFile("VMS_ADAPT_MIDDLE: ", ParamDB->VMS_ADAPT_MIDDLE); 
  printToFile("VMS_ADAPT_UPPER: ", ParamDB->VMS_ADAPT_UPPER); 
  printToFile("VMS_ADAPT_STEPS: ", ParamDB->VMS_ADAPT_STEPS); 
  printToFile("VMS_ADAPT_COMP: ", ParamDB->VMS_ADAPT_COMP); 

  printToFile("SUPERCONVERGENCE_ORDER: ", ParamDB->SUPERCONVERGENCE_ORDER);
  printToFile("WENO_TYPE: ", ParamDB->WENO_TYPE);
  
  printToFile("WRITE_SNAPSHOTS: ", ParamDB->WRITE_SNAPSHOTS);
  printToFile("DO_ROM: ", ParamDB->DO_ROM);
  printToFile("DO_ROM_P: ", ParamDB->DO_ROM_P);
  printToFile("RANK_OF_BASIS: ", ParamDB->RANK_OF_BASIS);
  printToFile("RANK_OF_BASIS_P: ", ParamDB->RANK_OF_BASIS_P);
  printToFile("POD_INNER_PRODUCT: ", ParamDB->POD_INNER_PRODUCT);
  printToFile("POD_INNER_PRODUCT_P: ", ParamDB->POD_INNER_PRODUCT_P);
  printToFile("BUILD_PODFILE: ", ParamDB->BUILD_PODFILE);
  printToFile("POD_FILENAME: ", ParamDB->POD_FILENAME);
  printToFile("SNAP_FILENAME: ", ParamDB->SNAP_FILENAME);
  printToFile("POD_FLUCT_FIELD: ", ParamDB->POD_FLUCT_FIELD);
  printToFile("POD_FLUCT_FIELD_P: ", ParamDB->POD_FLUCT_FIELD_P);
  printToFile("P_ROM_METHOD: ", ParamDB->P_ROM_METHOD);
  printToFile("PROJECTION_METHOD: ", ParamDB->PROJECTION_METHOD);

  /** write parameter for non-conforming elements */
  printToFile("NC_TYPE: ", ParamDB->NC_TYPE);

  printToFile("CHANNEL_STATISTICS2_WITH_MODEL: ", ParamDB->CHANNEL_STATISTICS2_WITH_MODEL);
  printToFile("BULK_REACTION_DISC: ", ParamDB->BULK_REACTION_DISC);
  printToFile("BULK_PB_DISC: ", ParamDB->BULK_PB_DISC);
  printToFile("BULK_PB_DISC_STAB: ", ParamDB->BULK_PB_DISC_STAB);
  printToFile("BULK_COUPLING: ", ParamDB->BULK_COUPLING);
  printToFile("BULK_GROWTH_RATE: ", ParamDB->BULK_GROWTH_RATE);
  printToFile("BULK_REACTION_MASS_LUMPING: ", ParamDB->BULK_REACTION_MASS_LUMPING);
  printToFile("BULK_REACTION_C_CUT: ", ParamDB->BULK_REACTION_C_CUT);
  printToFile("BULK_METHODS_OF_MOMENTS: ", ParamDB->BULK_METHODS_OF_MOMENTS);
  printToFile("BULK_MOM_DISC: ", ParamDB->BULK_MOM_DISC);
  printToFile("BULK_SOLD_PARAMETER_TYPE: ", ParamDB->BULK_SOLD_PARAMETER_TYPE);
  printToFile("N_CELL_LAYERS_PSD: ", ParamDB->N_CELL_LAYERS_PSD);
  printToFile("N_CELL_LAYERS_PSD_2: ", ParamDB->N_CELL_LAYERS_PSD);
  printToFile("OUTPUT_NODE_LAYER_PSD: ", ParamDB->OUTPUT_NODE_LAYER_PSD);
  printToFile("BULK_l_infty: ", ParamDB->BULK_l_infty);
  printToFile("BULK_u_infty: ", ParamDB->BULK_u_infty);
  printToFile("BULK_c_infty: ", ParamDB->BULK_c_infty);
  printToFile("BULK_c_C_infty_sat: ", ParamDB->BULK_c_C_infty_sat);
  printToFile("BULK_C_g: ", ParamDB->BULK_C_g);
  printToFile("BULK_C_nuc: ", ParamDB->BULK_C_nuc);
  printToFile("BULK_C_sat: ", ParamDB->BULK_C_sat);
  printToFile("BULK_C_2: ", ParamDB->BULK_C_2);
  printToFile("BULK_D_A: ", ParamDB->BULK_D_A);
  printToFile("BULK_D_P_0: ", ParamDB->BULK_D_P_0);
  printToFile("BULK_D_P_MAX: ", ParamDB->BULK_D_P_MAX);
  printToFile("BULK_k_g: ", ParamDB->BULK_k_g);
  printToFile("BULK_k_r: ", ParamDB->BULK_k_r);
  printToFile("BULK_k_nuc: ", ParamDB->BULK_k_nuc);
  printToFile("SSMUM_MP_X: ", ParamDB->SSMUM_MP_X);
  printToFile("SSMUM_MP_Y: ", ParamDB->SSMUM_MP_Y);
  printToFile("SSMUM_INNER_RADIUS: ", ParamDB->SSMUM_INNER_RADIUS);
  printToFile("SSMUM_OUTER_RADIUS: ", ParamDB->SSMUM_OUTER_RADIUS);
  printToFile("SSMUM_ROT_PER_SECOND: ", ParamDB->SSMUM_ROT_PER_SECOND);
  printToFile("SSMUM_MAX_CELLS_LAYERS: ", ParamDB->SSMUM_MAX_CELLS_LAYERS);
  printToFile("SSMUM_INTERPOLATION: ", ParamDB->SSMUM_INTERPOLATION);
  
  printToFile("INPUT_QUAD_RULE: ", ParamDB->INPUT_QUAD_RULE);
  
  printToFile("HEAT_TANGENTIAL_STRESS_FACTOR: ", ParamDB->HEAT_TANGENTIAL_STRESS_FACTOR);
  printToFile("HEAT_SOLID_SURFACE_FACTOR: ", ParamDB->HEAT_SOLID_SURFACE_FACTOR);
  printToFile("EQ_CONTACT_ANGLE: ", ParamDB->EQ_CONTACT_ANGLE);
  printToFile("AD_CONTACT_ANGLE: ", ParamDB->AD_CONTACT_ANGLE);
  printToFile("RE_CONTACT_ANGLE: ", ParamDB->RE_CONTACT_ANGLE);
  printToFile("DY_CONTACT_ANGLE: ", ParamDB->DY_CONTACT_ANGLE);  
  printToFile("CONTACT_ANGLE_TYPE: ", ParamDB->CONTACT_ANGLE_TYPE); 
  
}

void TDatabase::WriteTimeDB()
{
  using namespace Output; // printToFile
  
  printToFile("CURRENTTIME: ", TimeDB->CURRENTTIME);
  printToFile("CURRENTTIMESTEPLENGTH: ", TimeDB->CURRENTTIMESTEPLENGTH);
  printToFile("TIMESTEPLENGTH: ", TimeDB->TIMESTEPLENGTH);
  printToFile("MIN_TIMESTEPLENGTH: ", TimeDB->MIN_TIMESTEPLENGTH);
  printToFile("MAX_TIMESTEPLENGTH: ", TimeDB->MAX_TIMESTEPLENGTH);
  printToFile("TIMESTEPLENGTH_TOL: ", TimeDB->TIMESTEPLENGTH_TOL);
  printToFile("TIMESTEPLENGTH_CONTROL: ", TimeDB->TIMESTEPLENGTH_CONTROL);
  printToFile("TIMESTEPLENGTH_CONTROLLER: ",  TimeDB->TIMESTEPLENGTH_CONTROLLER);    
  printToFile("TIMESTEPLENGTH_PARA_KK_I: ",  TimeDB->TIMESTEPLENGTH_PARA_KK_I);
  printToFile("TIMESTEPLENGTH_PARA_KK_P: ",  TimeDB->TIMESTEPLENGTH_PARA_KK_P);
  printToFile("TIMESTEPLENGTH_PARA_KK_E: ",  TimeDB->TIMESTEPLENGTH_PARA_KK_E);
  printToFile("TIMESTEPLENGTH_PARA_KK_R: ",  TimeDB->TIMESTEPLENGTH_PARA_KK_R);
  printToFile("TIMESTEPLENGTH_PARA_KK_D: ", TimeDB->TIMESTEPLENGTH_PARA_KK_D);
  printToFile("TIMESTEPLENGTH_PARA_FAC: ",  TimeDB->TIMESTEPLENGTH_PARA_FAC);
  printToFile("TIMESTEPLENGTH_PARA_FAC_MAX: ",  TimeDB->TIMESTEPLENGTH_PARA_FAC_MAX);
  printToFile("TIMESTEPLENGTH_PARA_FAC_MIN: ",  TimeDB->TIMESTEPLENGTH_PARA_FAC_MIN);
  printToFile("TIMESTEPLENGTH_PARA_TOL: ",  TimeDB->TIMESTEPLENGTH_PARA_TOL);
  printToFile("TIMESTEPLENGTH_PARA_ATOL: ",  TimeDB->TIMESTEPLENGTH_PARA_ATOL);
  printToFile("TIMESTEPLENGTH_PARA_RTOL: ",  TimeDB->TIMESTEPLENGTH_PARA_RTOL);
  printToFile("RESET_CURRENTTIME: ", TimeDB->RESET_CURRENTTIME);
  printToFile("STEADY_STATE_TOL: ", TimeDB->STEADY_STATE_TOL);
  printToFile("SCALE_DIVERGENCE_CONSTRAINT: ", TimeDB->SCALE_DIVERGENCE_CONSTRAINT);

  printToFile("CONTROL: ", TimeDB->CONTROL);
  printToFile("CONTROL_ALPHA: ", TimeDB->CONTROL_ALPHA);
  printToFile("CONTROL_BETA: ", TimeDB->CONTROL_BETA);
  printToFile("CONTROL_GAMMA: ", TimeDB->CONTROL_GAMMA);
  printToFile("CONTROL_SAFTY: ", TimeDB->CONTROL_SAFTY);
  printToFile("CONTROL_MAXSCALE: ", TimeDB->CONTROL_MAXSCALE);
  printToFile("CONTROL_MINSCALE: ", TimeDB->CONTROL_MINSCALE);
  
  printToFile("THETA1: ", TimeDB->THETA1);
  printToFile("THETA2: ", TimeDB->THETA2);
  printToFile("THETA3: ", TimeDB->THETA3);
  printToFile("THETA4: ", TimeDB->THETA4);

  printToFile("TIME_DISC: ", TimeDB->TIME_DISC);
  printToFile("TIME_DISC2: ", TimeDB->TIME_DISC2);

  printToFile("T0: ", TimeDB->T0);
  printToFile("T1: ", TimeDB->T1);
  printToFile("T2: ", TimeDB->T2);
  printToFile("T3: ", TimeDB->T3);
  printToFile("T4: ", TimeDB->T4);
  printToFile("T5: ", TimeDB->T5);
  printToFile("T6: ", TimeDB->T6);
  printToFile("T7: ", TimeDB->T7);
  printToFile("T8: ", TimeDB->T8);
  printToFile("T9: ", TimeDB->T9);

  printToFile("STEPS_PER_IMAGE: ", TimeDB->STEPS_PER_IMAGE);
  printToFile("STEPS_PER_SNAP: ", TimeDB->STEPS_PER_SNAP);

  printToFile("RB_TYPE: ", TimeDB->RB_TYPE);
  printToFile("RB_TYPE2: ", TimeDB->RB_TYPE2);
    
  printToFile("EXTRAPOLATE_VELOCITY: ", TimeDB->EXTRAPOLATE_VELOCITY);
  printToFile("EXTRAPOLATE_PRESSURE: ", TimeDB->EXTRAPOLATE_PRESSURE);
  printToFile("EXTRAPOLATE_STEPS: ", TimeDB->EXTRAPOLATE_STEPS);
  printToFile("EXTRAPOLATE_WEIGHT: ", TimeDB->EXTRAPOLATE_WEIGHT);

  printToFile(">>> End printing old database (ParamDB) <<<");
  printToFile("");
    
//  printToFile(" n_neumann_boundary: ",ParamDB->n_neumann_boundary );
//  printToFile(" neumann_boundary_id: ", ParamDB->neumann_boundary_id);
//  printToFile("neumann_boundary_value: ", ParamDB->neumann_boundary_value);
}

TDatabase::~TDatabase()
{
  // allocate databases
  delete ParamDB;
  delete TimeDB;

  // initialize shape descriptors
  delete ShapeDB[S_Line];
  delete ShapeDB[Triangle];
  delete ShapeDB[Quadrangle];
  delete ShapeDB[Parallelogram];
  delete ShapeDB[Rectangle];
  #ifdef __3D__
    delete ShapeDB[Tetrahedron];
    delete ShapeDB[Hexahedron];
    delete ShapeDB[Brick];
  #endif
  delete [] ShapeDB;
  
  delete RefDescDB[S_Line];
  delete RefDescDB[Triangle];
  delete RefDescDB[Quadrangle];
  delete RefDescDB[Parallelogram];
  delete RefDescDB[Rectangle];
  #ifdef __3D__
    delete RefDescDB[Tetrahedron];
    delete RefDescDB[Hexahedron];
    delete RefDescDB[Brick];
  #endif
  delete RefDescDB[N_SHAPES + LineReg];
  delete RefDescDB[N_SHAPES + TriReg];
  delete RefDescDB[N_SHAPES + TriBis0];
  delete RefDescDB[N_SHAPES + TriBis1];
  delete RefDescDB[N_SHAPES + TriBis2];
  delete RefDescDB[N_SHAPES + TriBis01];
  delete RefDescDB[N_SHAPES + TriBis02];
  delete RefDescDB[N_SHAPES + TriBis10];
  delete RefDescDB[N_SHAPES + TriBis12];
  delete RefDescDB[N_SHAPES + TriBis20];
  delete RefDescDB[N_SHAPES + TriBis21];
  delete RefDescDB[N_SHAPES + QuadReg];
  delete RefDescDB[N_SHAPES + ParallReg];
  delete RefDescDB[N_SHAPES + RectReg];
  delete RefDescDB[N_SHAPES + QuadBis0];
  delete RefDescDB[N_SHAPES + QuadBis1];
  delete RefDescDB[N_SHAPES+Quad1Conf0];
  delete RefDescDB[N_SHAPES+Quad1Conf1];
  delete RefDescDB[N_SHAPES+Quad1Conf2];
  delete RefDescDB[N_SHAPES+Quad1Conf3];
  delete RefDescDB[N_SHAPES+Quad2Conf0];
  delete RefDescDB[N_SHAPES+Quad2Conf1];
  delete RefDescDB[N_SHAPES+Quad2Conf2];
  delete RefDescDB[N_SHAPES+Quad2Conf3];
  delete RefDescDB[N_SHAPES+QuadToTri0];
  delete RefDescDB[N_SHAPES+QuadToTri1];

  #ifdef __3D__
    delete RefDescDB[N_SHAPES + TetraReg];
    delete RefDescDB[N_SHAPES + TetraReg0];
    delete RefDescDB[N_SHAPES + TetraReg1];
    delete RefDescDB[N_SHAPES + TetraReg2];
    delete RefDescDB[N_SHAPES + TetraBis0];
    delete RefDescDB[N_SHAPES + TetraBis1];
    delete RefDescDB[N_SHAPES + TetraBis2];
    delete RefDescDB[N_SHAPES + TetraBis3];
    delete RefDescDB[N_SHAPES + TetraBis4];
    delete RefDescDB[N_SHAPES + TetraBis5];
    delete RefDescDB[N_SHAPES + TetraBis01];
    delete RefDescDB[N_SHAPES + TetraBis02];
    delete RefDescDB[N_SHAPES + TetraBis03];
    delete RefDescDB[N_SHAPES + TetraBis04];
    delete RefDescDB[N_SHAPES + TetraBis05];
    delete RefDescDB[N_SHAPES + TetraBis10];
    delete RefDescDB[N_SHAPES + TetraBis12];
    delete RefDescDB[N_SHAPES + TetraBis13];
    delete RefDescDB[N_SHAPES + TetraBis14];
    delete RefDescDB[N_SHAPES + TetraBis15];
    delete RefDescDB[N_SHAPES + TetraBis20];
    delete RefDescDB[N_SHAPES + TetraBis21];
    delete RefDescDB[N_SHAPES + TetraBis23];
    delete RefDescDB[N_SHAPES + TetraBis24];
    delete RefDescDB[N_SHAPES + TetraBis25];
    delete RefDescDB[N_SHAPES + TetraBis30];
    delete RefDescDB[N_SHAPES + TetraBis32];
    delete RefDescDB[N_SHAPES + TetraBis34];
    delete RefDescDB[N_SHAPES + TetraBis35];
    delete RefDescDB[N_SHAPES + TetraBis40];
    delete RefDescDB[N_SHAPES + TetraBis41];
    delete RefDescDB[N_SHAPES + TetraBis43];
    delete RefDescDB[N_SHAPES + TetraBis45];
    delete RefDescDB[N_SHAPES + TetraBis51];
    delete RefDescDB[N_SHAPES + TetraBis52];
    delete RefDescDB[N_SHAPES + TetraBis53];
    delete RefDescDB[N_SHAPES + TetraBis54];
    delete RefDescDB[N_SHAPES + TetraQuad0];
    delete RefDescDB[N_SHAPES + TetraQuad1];
    delete RefDescDB[N_SHAPES + TetraQuad2];
    delete RefDescDB[N_SHAPES + TetraQuad3];
    delete RefDescDB[N_SHAPES + HexaReg];
    delete RefDescDB[N_SHAPES + BrickReg];
  #endif
  delete [] RefDescDB;
  
  #ifdef __3D__
    //initialize mapper
    delete MapperDB[MapTriReg0];
    delete MapperDB[MapTriReg1];
    delete MapperDB[MapTriReg2];

    delete MapperDB[MapTriBis00];
    delete MapperDB[MapTriBis01];
    delete MapperDB[MapTriBis02];
    delete MapperDB[MapTriBis10];
    delete MapperDB[MapTriBis11];
    delete MapperDB[MapTriBis12];
    delete MapperDB[MapTriBis20];
    delete MapperDB[MapTriBis21];
    delete MapperDB[MapTriBis22];
    delete MapperDB[MapTriBis010];
    delete MapperDB[MapTriBis011];
    delete MapperDB[MapTriBis012];
    delete MapperDB[MapTriBis020];
    delete MapperDB[MapTriBis021];
    delete MapperDB[MapTriBis022];
    delete MapperDB[MapTriBis100];
    delete MapperDB[MapTriBis101];
    delete MapperDB[MapTriBis102];
    delete MapperDB[MapTriBis120];
    delete MapperDB[MapTriBis121];
    delete MapperDB[MapTriBis122];
    delete MapperDB[MapTriBis200];
    delete MapperDB[MapTriBis201];
    delete MapperDB[MapTriBis202];
    delete MapperDB[MapTriBis210];
    delete MapperDB[MapTriBis211];
    delete MapperDB[MapTriBis212];

    delete MapperDB[MapQuadReg0];
    delete MapperDB[MapQuadReg1];
    delete MapperDB[MapQuadReg2];
    delete MapperDB[MapQuadReg3];
  #endif
  delete [] MapperDB;

  // initialize iterators
  delete IteratorDB[It_EQ];
  delete IteratorDB[It_LE];
  delete IteratorDB[It_Finest];
  delete IteratorDB[It_EQLevel];
  delete IteratorDB[It_LELevel];
  delete IteratorDB[It_Between];
  delete IteratorDB[It_OCAF];

  #ifdef __MORTAR__
    delete IteratorDB[It_Mortar1];
    delete IteratorDB[It_Mortar2];
  #endif
  delete [] IteratorDB;
}

TParaDB::~TParaDB()
{
  // call delete on all char* which were created in 
  // TDatabase::SetDefaultParameters()
  delete [] MAPFILE;
  delete [] OUTFILE;
  delete [] SAVE_DATA_FILENAME;
  delete [] READ_DATA_FILENAME;
  delete [] SMESHFILE;
  delete [] POD_FILENAME;
  delete [] SNAP_FILENAME;
  delete [] FS_INNAME;
  delete [] FS_OUTNAME;
  delete [] MATLAB_MATRIX;
  delete [] PODFILE;
  Output::print<4>("deleted parameter database");
}


/* ************************************************************************** */
/// plays the role of the old CheckParameterConsistencyNSE
/// but with the new database.The remaining global parameters
/// can/must be removed here, progressively.
void check_parameters_consistency_NSE(ParameterDatabase& db)
{
#ifdef _MPI
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
#else
  int my_rank = 0;
#endif
  if(!db.contains("space_discretization_type"))
  {
    ErrThrow("check_parameters_consistency_NSE: you need to set the "
             "space_discretization_type");
  }

  // Newton method
  if ((TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE)&&(TDatabase::ParamDB->NSTYPE<=2))
  {
    TDatabase::ParamDB->NSTYPE+=2;
    if(my_rank==0)
      Output::info("NSE Parameter Consistency","NSTYPE changed to ", TDatabase::ParamDB->NSTYPE,
                  " because of SC_NONLIN_ITE_TYPE_SADDLE  = ",
                  TDatabase::ParamDB->SC_NONLIN_ITE_TYPE_SADDLE);
  }


  if (db["space_discretization_type"].is("sdfem") && (TDatabase::ParamDB->NSTYPE==1))
  {
      //TDatabase::ParamDB->NSTYPE = 2;
      //Output::info("NSE Parameter Consistency","NSTYPE changed from 1 to 2 because of SDFEM discretization ");
    if(my_rank==0)
        Output::info("NSE Parameter Consistency","NSTYPE 1: only reduced SDFEM, only for 2D, fixed point, not skew !!!");
  }

  if ((db["space_discretization_type"].is("sdfem")) && (TDatabase::ParamDB->NSTYPE==3))
  {
    TDatabase::ParamDB->NSTYPE = 4;
    if(my_rank==0)
      Output::warn<1>("NSE Parameter Consistency","NSTYPE changed from 3 to 4 because of SDFEM discretization ");
  }

  if ((TDatabase::ParamDB->LAPLACETYPE == 1) && (TDatabase::ParamDB->NSTYPE ==1))
  {
    TDatabase::ParamDB->NSTYPE = 3 ;
    if(my_rank==0)
      Output::warn<1>("NSE Parameter Consistency","NSTYPE changed from 1 to 3 because of LAPLACETYPE ");
  }

  if ((TDatabase::ParamDB->LAPLACETYPE == 1) && (TDatabase::ParamDB->NSTYPE ==2))
  {
    TDatabase::ParamDB->NSTYPE = 4 ;
    if(my_rank==0)
      Output::warn<1>("NSE Parameter Consistency","NSTYPE changed from 2 to 4 because of LAPLACETYPE ");
  }

  // rotational form
  Parameter nonlin_form(db["nse_nonlinear_form"]);
  if(nonlin_form.is("rotational"))
  {
    if (TDatabase::ParamDB->NSTYPE<=2)
    {
      TDatabase::ParamDB->NSTYPE+=2;
      if(my_rank==0)
      {
        Output::warn<1>("NSE Parameter Consistency",
                        "NSTYPE changed to ", TDatabase::ParamDB->NSTYPE);
        Output::warn<1>("NSE Parameter Consistency",
                        " because of nse_nonlinear_form = ", nonlin_form);
      }
    }
  }

  if (TDatabase::ParamDB->TURBULENT_VISCOSITY_TYPE==5)
  {
    if (!(db["space_discretization_type"].is("vms_projection")||
        db["space_discretization_type"].is("vms_projection_expl")))
    {
      if(my_rank==0)
      {
        Output::warn<1>("NSE Parameter Consistency","TURBULENT_VISCOSITY_TYPE = 5 only defined for projection-based VMS methods");
        Output::warn<1>("NSE Parameter Consistency","Set different TURBULENT_VISCOSITY_TYPE !!!");
        ErrThrow("BOOM!");
      }
    }
  }

  // LOCAL_PROJECTION
  if (db["space_discretization_type"].is("local_projection"))
  {
    if (TDatabase::ParamDB->LP_FULL_GRADIENT)
    {
      if (TDatabase::ParamDB->LP_STREAMLINE)
      {
        TDatabase::ParamDB->LP_STREAMLINE = 0;
        if(my_rank==0)
        Output::warn<1>("NSE Parameter Consistency","LP_STREAMLINE changed to ", TDatabase::ParamDB->LP_STREAMLINE,
                      " due to LP_FULL_GRADIENT = ", TDatabase::ParamDB->LP_FULL_GRADIENT);
      }

      if (TDatabase::ParamDB->LP_DIVERGENCE)
      {
        TDatabase::ParamDB->LP_DIVERGENCE = 0;
        if(my_rank==0)
        Output::warn<1>("NSE Parameter Consistency","LP_DIVERGENCE changed to ", TDatabase::ParamDB->LP_DIVERGENCE,
                      " due to LP_FULL_GRADIENT = ", TDatabase::ParamDB->LP_FULL_GRADIENT);
      }
    } // end LP_FULL_GRADIENT

    if (TDatabase::ParamDB->LP_DIVERGENCE)
    {
      if (TDatabase::ParamDB->NSTYPE<=2)
      {
        TDatabase::ParamDB->NSTYPE+=2;
        if(my_rank==0)
        Output::warn<1>("NSE Parameter Consistency","NSTYPE changed to ", TDatabase::ParamDB->NSTYPE,
                      " LP_DIVERGENCE = ", TDatabase::ParamDB->LP_DIVERGENCE);
      }
    }

    if(TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE == -123)
      TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

    if(TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE == -123)
      TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

    if(TDatabase::ParamDB->LP_DIVERGENCE_ORDER_DIFFERENCE == -123)
      TDatabase::ParamDB->LP_DIVERGENCE_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

    if(TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE == -123)
      TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  } // end "discretization_type" = "LOCAL_PROJECTION"
  else
  {
    // switch off all local projection terms
    TDatabase::ParamDB->LP_FULL_GRADIENT = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;

    TDatabase::ParamDB->LP_STREAMLINE = 0;
    TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
    TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;

    TDatabase::ParamDB->LP_DIVERGENCE = 0;
    TDatabase::ParamDB->LP_DIVERGENCE_COEFF = 0;
    TDatabase::ParamDB->LP_DIVERGENCE_EXPONENT = 1;

    TDatabase::ParamDB->LP_PRESSURE = 0;
    TDatabase::ParamDB->LP_PRESSURE_COEFF = 0;
    TDatabase::ParamDB->LP_PRESSURE_EXPONENT = 1;

    TDatabase::ParamDB->LP_ORDER_DIFFERENCE = 1;
    TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = 1;
    TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = 1;
    TDatabase::ParamDB->LP_DIVERGENCE_ORDER_DIFFERENCE = 1;
    TDatabase::ParamDB->LP_PRESSURE_ORDER_DIFFERENCE = 1;
  }

  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == STOKES)
    TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR = 1;
  if (TDatabase::ParamDB->FLOW_PROBLEM_TYPE == OSEEN)
  {
    TDatabase::ParamDB->INTERNAL_PROBLEM_LINEAR = 1;
    switch (TDatabase::ParamDB->NSTYPE)
    {
      case 1:
        if(my_rank==0)
            Output::warn<1>("NSE Parameter Consistency","Galerkin discretization for Oseen because of NSTYPE ",
                      TDatabase::ParamDB->NSTYPE);
        db["space_discretization_type"].set("galerkin"); // DISCTYPE = 1
        break;
      case 14:
      if(my_rank==0)
          Output::warn<1>("NSE Parameter Consistency","SUPG/PSPG/grad-div discretization for Oseen because of NSTYPE ",
                      TDatabase::ParamDB->NSTYPE);
        db["space_discretization_type"].set("supg"); // DISCTYPE = 2
        break;
      default:
      if(my_rank==0)
          Output::warn<1>("NSE Parameter Consistency","No method for Oseen implemented for NSTYPE ",
                      TDatabase::ParamDB->NSTYPE);
        exit(4711);
    }
  }
}

/* ************************************************************************** */
/// plays the role of the old SetParametersCD
/// but with the new database.The remaining global parameters
/// can/must be removed here, progressively.
void check_parameters_consistency_CD(ParameterDatabase& db, int &nonlinear_method)
{
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== KLR02_3)
    TDatabase::ParamDB->SOLD_S = 0;
  if (TDatabase::ParamDB->SOLD_PARAMETER_TYPE== LP96)
  {
    OutPut("SOLD_PARAMETER_TYPE == LP96 should be used with higher quadrature rule,"<<endl);
    OutPut("since right hand side is in general not linear !!!"<<endl);
  }

  if (!db["space_discretization_type"].is("local_projection") &&
      !db["space_discretization_type"].is("local_projection_2_level"))
  {
    // switch off all local projection terms
    TDatabase::ParamDB->LP_FULL_GRADIENT = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0;
    TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;

    TDatabase::ParamDB->LP_STREAMLINE = 0;
    TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
    TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
  }
  else
  {
    if (db["space_discretization_type"].is("local_projection"))
    {
      // check spaces and change if necessary
      switch(TDatabase::ParamDB->ANSATZ_ORDER)
      {
        case 1:
          TDatabase::ParamDB->ANSATZ_ORDER = 100;
          OutPut("ANSATZ_ORDER changed to " << TDatabase::ParamDB->ANSATZ_ORDER << endl);
          break;
        case 2:
          TDatabase::ParamDB->ANSATZ_ORDER = 201;
          OutPut("ANSATZ_ORDER changed to " << TDatabase::ParamDB->ANSATZ_ORDER << endl);
          break;
        case 3:
          TDatabase::ParamDB->ANSATZ_ORDER = 302;
          OutPut("ANSATZ_ORDER changed to " << TDatabase::ParamDB->ANSATZ_ORDER << endl);
          break;
        case 4:
          TDatabase::ParamDB->ANSATZ_ORDER = 403;
          OutPut("ANSATZ_ORDER changed to " << TDatabase::ParamDB->ANSATZ_ORDER << endl);
          break;
        case 5:
          TDatabase::ParamDB->ANSATZ_ORDER = 504;
          OutPut("ANSATZ_ORDER changed to " << TDatabase::ParamDB->ANSATZ_ORDER << endl);
          break;
        default:
          break;
      }
    }
  }

  if(TDatabase::ParamDB->LP_FULL_GRADIENT)
  {
    if(TDatabase::ParamDB->LP_STREAMLINE)
    {
      TDatabase::ParamDB->LP_STREAMLINE = 0;
      TDatabase::ParamDB->LP_STREAMLINE_COEFF = 0;
      TDatabase::ParamDB->LP_STREAMLINE_EXPONENT = 1;
      OutPut("local projection stabilisation in streamline direction ");
      OutPut("is switched off due to stabilisation of full gradient." << endl);
    }
  }

  if(TDatabase::ParamDB->LP_STREAMLINE)
  {
    if(TDatabase::ParamDB->LP_FULL_GRADIENT)
    {
      TDatabase::ParamDB->LP_FULL_GRADIENT = 0;
      TDatabase::ParamDB->LP_FULL_GRADIENT_COEFF = 0;
      TDatabase::ParamDB->LP_FULL_GRADIENT_EXPONENT = 1;
      OutPut("local projection stabilisation for gradient ");
      OutPut("is switched off due to stabilisation of streamline direction." << endl);
    }
  }

  if (TDatabase::ParamDB->LP_STREAMLINE && TDatabase::ParamDB->LP_CROSSWIND)
    nonlinear_method = 1;

  if(TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_FULL_GRADIENT_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  if(TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE == -123)
    TDatabase::ParamDB->LP_STREAMLINE_ORDER_DIFFERENCE = TDatabase::ParamDB->LP_ORDER_DIFFERENCE;

  // has to be changed for all discretizations
  // otherwise access to not allocated array in error computation
  if ((TDatabase::ParamDB->SDFEM_TYPE == 100)&&(!TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM))
  {
    TDatabase::ParamDB->SDFEM_TYPE = 2;
    OutPut("Changed Database::ParamDB->SDFEM_TYPE to " << TDatabase::ParamDB->SDFEM_TYPE
      << " since no adjoint problem is solved !!! "<<endl);
  }
  if ((TDatabase::ParamDB->SDFEM_TYPE != 100)&&(TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&
      (db["space_discretization_type"].is("sdfem")))
  {
    TDatabase::ParamDB->SDFEM_TYPE = 100;
    OutPut("Changed Database::ParamDB->SDFEM_TYPE to " << TDatabase::ParamDB->SDFEM_TYPE
      << " since adjoint problem is solved !!! "<<endl);
  }

  // SUPG, GLS
  // there is only one minor difference in method, treat them otherwise the same
  if (db["space_discretization_type"].is("gls"))
  {
    db["space_discretization_type"].set("supg");
    TDatabase::ParamDB->INTERNAL_DISC_FLAG = 1;
  }
  if ((db["space_discretization_type"].is("sdfem"))
      &&(TDatabase::ParamDB->SOLD_TYPE==0))
  {
    // this excludes some not wished side effects
    TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
  }
  if ((TDatabase::ParamDB->SOLVE_ADJOINT_PROBLEM)&&(TDatabase::ParamDB->SOLD_ADJOINT))
    TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 100;

  if  ((db["space_discretization_type"].is("sdfem"))
      &&(TDatabase::ParamDB->SOLD_PARAMETER_TYPE == FEM_TVD))
  {
    TDatabase::ParamDB->SDFEM_TYPE = 0;
    TDatabase::ParamDB->DELTA0 =  TDatabase::ParamDB->DELTA1 = 0;
    OutPut("FEM-TVD: switched stabilization off!" << endl);
  }

  if ((db["space_discretization_type"].is("sdfem"))
      &&(TDatabase::ParamDB->SOLD_TYPE))
    nonlinear_method = 1;
  if ((db["space_discretization_type"].is("sdfem"))
      &&((TDatabase::ParamDB->SDFEM_TYPE==2)
    || (TDatabase::ParamDB->SDFEM_TYPE==100)))
    SetPolynomialDegree();

  if ((db["space_discretization_type"].is("sdfem"))
      &&(TDatabase::ParamDB->SDFEM_TYPE==2))
  {
    if (TDatabase::ParamDB->CELL_MEASURE != 4)
    {
      TDatabase::ParamDB->CELL_MEASURE = 4;
      OutPut("CELL_MEASURE changed to " << TDatabase::ParamDB->CELL_MEASURE << endl);
    }
    if (TDatabase::ParamDB->DELTA0 != 1.0)
    {
      TDatabase::ParamDB->DELTA0 = 1.0;
      OutPut("DELTA0 changed to " << TDatabase::ParamDB->DELTA0 << endl);
    }
  }

  if (db["space_discretization_type"].is("cip"))
  {
    TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS = 1;
    if (TDatabase::ParamDB->CIP_TYPE < 10)
    {
      db["space_discretization_type"].set("galerkin");
      TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
    }
    else
    {
      // add SOLD term to CIP discretization
      nonlinear_method = 1;
      TDatabase::ParamDB->CIP_TYPE -=10;
      db["space_discretization_type"].set("supg");
      if (TDatabase::ParamDB->SOLD_TYPE == 0)
        TDatabase::ParamDB->SOLD_TYPE = 2;
      // switch off SUPG part
      TDatabase::ParamDB->SDFEM_TYPE = 0;
      TDatabase::ParamDB->DELTA0 = 0.0;
      TDatabase::ParamDB->DELTA1 = 0.0;
      //  TDatabase::ParamDB->SOLD_PARAMETER_TYPE = BE05_1;
      //  TDatabase::ParamDB->SOLD_TYPE = 0;
    }
  }

  if (db["space_discretization_type"].is("dg"))
  {
    db["space_discretization_type"].set("galerkin");
    TDatabase::ParamDB->INTERNAL_FACE_INTEGRALS = 2;
    if ( TDatabase::ParamDB->ANSATZ_ORDER < 10)
      TDatabase::ParamDB->ANSATZ_ORDER = -TDatabase::ParamDB->ANSATZ_ORDER-10;
    else
      // P elements on quads
      TDatabase::ParamDB->ANSATZ_ORDER = -10*TDatabase::ParamDB->ANSATZ_ORDER;
    if (TDatabase::ParamDB->ESTIMATE_ERRORS)
    {
      TDatabase::ParamDB->ESTIMATE_ERRORS = 0;
      OutPut("Error estimation does not work for DG !!!"<< endl);
    }
  }
  if  (!db["space_discretization_type"].is("sdfem"))
  {
    TDatabase::ParamDB->SOLD_TYPE = 0;
    TDatabase::ParamDB->SOLD_PARAMETER_TYPE =0;
  }
  if (db["space_discretization_type"].is("local_projection_2_level"))
  {
    nonlinear_method = 1;
    TDatabase::ParamDB->SOLD_PARAMETER_TYPE = 0;
  }
}



