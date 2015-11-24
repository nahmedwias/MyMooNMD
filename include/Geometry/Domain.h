/** ************************************************************************ 
*
* @class TDomain  
* @date  09.07.97
* @brief  contains the boundary description, the virtual cell tree and macro grid
* @author Volker Behns
* @History: collection methods (Gunar Matthies 14.10.97), 
            methods for Refine/Derefine algorithm  (Gunar Matthies 17.10.97)
            mesh class, edge generation, parallel methods and documentation (Sashikumaar Ganesan 08.08.2014)
   
************************************************************************  */

#ifndef __DOMAIN__
#define __DOMAIN__

class TDomain;

#include <BoundPart.h>
#include <Collection.h>
#include <Database.h>
#include <Iterator.h>

#ifdef __MORTAR__
struct TMortarFaceStruct
       {
         TBaseCell *Cell;
         int LocFaceNumber[2];
       };

typedef struct TMortarFaceStruct TMortarFace;
#endif

/** contains the boundary description, the virtual cell
    tree and macro grid */
class TDomain
{
  protected:
    
    /** @brief number of boundary parts */
    int N_BoundParts;    
    
    /** @brief boundary parts of domain */
    TBoundPart **BdParts;

    /** @brief number of all boundary components */
    int N_BoundComps;
    
    /** @brief start id of boundary component on each boundary part */
    int *StartBdCompID;

    /** @brief boundary part id's of all Interfaces */
    int *Interfaces;

    /** @brief number of holes */
    int N_Holes;
    
    /** @brief point in each hole */
    double *PointInHole;

    /** @brief number of regions */
    int N_Regions;
    
    /** @brief point in each region */
    double *PointInRegion;

    /** @brief array of all root cells of cell tree */
    TBaseCell **CellTree;
    
    /** @brief number of all root cells of cell tree */
    int N_RootCells;

    /** @brief number of virtuell cells on initial level */
    int N_InitVCells;

    /** @brief x coordinate of the start point (2D) */
    double StartX;
    /** @brief y coordinate of the start point (2D) */
    double StartY;
    /** @brief x length of bounding box */
    double BoundX;
    /** @brief y length of bounding box */
    double BoundY;

#ifdef __3D__
      /** @brief third coordinate of start point (3D) */
      double StartZ;
      /** @brief return number of cell in Y direction */
      int N_InitVCellsY;
      /** @brief z length of the bounding box */
      double BoundZ;
#endif

    /** @brief current refinment level */
    int RefLevel;

#ifdef __MORTAR__
      /** @brief number of mortar faces */
      int N_MortarFaces;
      /** @brief structur for mortar faces */
      TMortarFace *MortarFaces;

      /** @brief begin of each mortar face on coll */
      int *BeginMFace;
#endif
    
    friend class TTetGenMeshLoader;

#ifdef  _MPI
      /** @brief array contains the global cell number of local cells (including Halo cells) */
      int *GlobalCellIndex;

      /** @brief Number of own cells (excluding Halo cells) */
      int N_OwnCells;
#endif

  public:
    /**
     * @brief Default constructor.
     * Does only set RefLevel (refinement level) to 0.
     */
    TDomain();

    /**
     * @brief Constructor. Reads in some data.
     * @param ParamFile Path to a ParMooN parameter input textfile.
     *
     * Invokes a read-in function for the parameter file and the domain
     * description file afterwards.
     *
     * TODO This is messy in will be tidied up in the near future.
     */
    TDomain(char *ParamFile);
    
    // Methods
    /** @brief Read in initial mesh from ".GEO" file.
     *
     * @param[in] dat Input stream which contains the initial mesh information
     * in MooNMD-native ".GEO"-format (or .xGEO).
     *
     * @param[in] isXGeoFile Set true if the input is in extended GEO (.xGEO) format.
     *
     * @return Integer O on success, other numbers otherwise.
     */
    int ReadGeo(std::istream& dat, bool isXGeoFile);

#ifdef __3D__
    /** @brief Read in initial mesh from ".GEO" file as sandwich geometry.
     *
     * @param[in] dat Input stream which contains the initial mesh information
     * in MooNMD-native ".GEO"-format.
     *
     * @note The implementation of this method contains some undocumented surprises
     * ("else if(TDatabase::ParamDB->INTERNAL_PROBLEM_IDENTITY == 180)"). Handle
     * with care.
     *
     * @return Integer O on success, other numbers otherwise.
     */
    int ReadSandwichGeo(std::istream& dat);

    /** @brief make boundary parameter consistent */
    void MakeBdParamsConsistent(TCollection *coll);
    
    int CloseGrid(int level);
#endif

#ifdef __2D__
    /** @brief set boundary parameters and cell shape according to
        possibly moving vertices */
    void CorrectParametersAndShapes();
    
    /** @brief mesh genration using Triangle for give IN */
    void TriMeshGen(struct triangulateio *In);
#endif

    /** @brief read parameter file */
    int ReadParam(char *ParamFile);

    /** @brief Reads in boundary parameterization in ".PRM" file format
     *
     * @param[in] dat Input stream containing the boundary information in
     * ".PRM" style (the default MooNMD file format for domain description).
     *
     * @param[in,out] Flag Used in 3D only, set to 1 if a boundary component
     * of type TBdWall has to be constructed. Otherwise set to 0.
     *
     * @return Integer O on success, other numbers otherwise.
     */
    int ReadBdParam(std::istream& dat, int &Flag);

    /** @brief read mapping and mortar information */
    int ReadMapFile(char *MapFile, TDatabase *database);

    /** @brief get boundary part of BdCompID */
    int GetBdPartID(int BdCompID);
    /** @brief get local number of boundary component */
    int GetLocalBdCompID(int BdCompID);
    /** @brief get local number of last boundary component on part BdPartID*/
    int GetLastLocalComp(int BdPartID)
    { return StartBdCompID[BdPartID+1] - StartBdCompID[BdPartID] - 1; }
    /** @brief set start BdCompID on boundary part i */
    void SetStartBdCompID(int BdCompID, int i)
    { StartBdCompID[i] = BdCompID; }

    /** @brief get i-th boundary part */
    TBoundPart *GetBdPart(int i)
    { return BdParts[i]; }

    /** @brief get tree of cells */
    void GetTreeInfo(TBaseCell **&celltree, int &N_rootcells)
    { 
      celltree = CellTree;
      N_rootcells = N_RootCells;
    }

    /** @brief set tree of cells */
    void SetTreeInfo(TBaseCell **celltree, int N_rootcells)
    {
      CellTree = celltree;
      N_RootCells = N_rootcells;

      #ifdef  _MPI
      N_OwnCells = 0;
      #endif 
    }

    #ifdef __MORTAR__
      /** @brief set subgrid ID's on all MacroCells and generate mortar structurs */
      int SetSubGridIDs(IntFunct2D *TestFunc);
      /** @brief generate mortar structurs */
      int GenMortarStructs();
      /** @brief return number of mortar face structs */
      int GetN_MortarFace()
      { return N_MortarFaces; }
      /** @brief return mortar face struct */
      TMortarFace *GetMortarFace(int i)
      { return &(MortarFaces[i]); }

      /** @brief get a collection of all mortar cells */
      TCollection *GetMortarColl(Iterators it, int level);

      /** @brief initialize all mortar joints with FE-information */
      int InitMortarJoints(Iterators it, int level, TCollection *coll);
    #endif

    /** @brief generate initial grid using external mesh generator */
    int GenInitGrid();
    #ifdef __2D__
      /** @brief make initial 2D grid */
      int MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int N_Vertices,
                   int NVE);
      /** @brief make initial 2D grid, extended version which sets ReferenceID
       *         in cells */
      int MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int *ELEMSREF,
                   int N_Vertices, int NVE);
    #else
      /** @brief make initial 3D grid */
      int MakeGrid(double *DCORVG, int *KVERT, int *KNPR, int *ELEMSREF,
                   int N_Vertices, int NVE, int *BoundFaces, int *FaceParam,
                   int NBF, int NVpF,
                   int *Interfaceparam, int N_Interfaces);
      /** @brief make initial sandwich grid */
      int MakeSandwichGrid(double *DCORVG, int *KVERT, int *KNPR,
                           int N_Vertices, int NVE,
                           double DriftX, double DriftY, double DriftZ,
                           int N_Layers, double *Lambda);
    #endif

    /**
      * @brief Chooses in what way to construct the domain's geometry
      * and calls the corresponding method.
      *
      * The strings (or rather: char arrays) handed over to this function
      * determine how the domain description and the initial mesh are constructed.
      * The method itself does non of the initializing work but invokes those
      * functions which do.
      *
      * @param[in] PRM The description of the domain boundaries. Possibilities are
      *
      * "Default_UnitSquare" - Default domain boundary description (2D only).
      * "Default_UnitCube" - Default domain boundary description (3D only).
      *
      * If not this and not NULL, PRM is interpreted as a filepath to a .PRM file
      * and gets handed over to TDomain::ReadBdParam() to be read in.
      *
      *	@param[in] GEO The description of the initial mesh. Possibilities are:
      *
      *	"InitGrid" - Grid Generator. (2D only) Call TDomain::GenInitGrid() (creates initial grid with TRIANGLE)
      * "TwoTriangles" - Default mesh. Call TDomain::TwoTriangles. (2D only)
      * "TwoTrianglesRef" - Default mesh. Call TDomain::TwoTrianglesRef. (2D only)
      * "UnitSquare" - Default mesh. ... (2D only)
      * "UnitSquareRef" - Default mesh. ... (2D only)
      * "SquareInSquare" - Default mesh. ... (2D only)
      * "SquareInSquareRef" - Default mesh. ... (2D only)
      * "PeriodicSquares" - Default mesh. ... (2D only)
      * "PeriodicSquaresLarge" - Default mesh. ... (2D only)
      * "PeriodicTrianglesLarge" - Default mesh. ... (2D only)
      * "PeriodicRectangle_2_4" - Default mesh. ... (2D only)
      * "TestGrid3D" - Default mesh. ... (3D only)
      *
      * If none of these, the string is considered to be the path to a .GEO
      * file and thus is handed over to TDomain::ReadGeo().
      *
      * @note: The combination of the default domain and mesh descriptions
      * is neither checked nor tested. So use them carefully and be prepared for the worst!
      *
      */
    void Init(char *PRM, char *GEO);

    /** @brief write domain boundary  into a postscript file */
    int Draw(char *name, Iterators iterator, int arg);
    /** @brief write mesh into a postscript file */
    int PS(const char *name, Iterators iterator, int arg);
    /** @brief write collection into a postscript file */
    int PS(const char *name, TCollection *Coll);
    /** @brief write files for MD-Out format */
    int MD_raw(const char *name, Iterators iterator, int arg);

    /** @brief refine the grid according the cell refinement descriptors */
    int Refine();
    /** @brief refine all cells regular */
    int RegRefineAll();
    /** @brief refine all cells in subgrid ID regular */
    int RegRefineSub(int ID);
    /** @brief refine only in one direction */
    int RefineallxDirection();
    /** @brief generate a 1-regular grid */
    int Gen1RegGrid();
    /** @brief refine the finest grid according the given indicator function */
    int RefineByIndicator(DoubleFunct2D *Indicator);
    /** @brief refine the finest grid if necessary in order to get a 
        grid with conforming closures */
    int MakeConfClosure();

    /** @brief refine the finest grid according a given error estimate */
    int RefineByErrorEstimator(TCollection *Collection,double *eta_K,
                               double eta_max,double tolerance,
                               bool ConfClosure);

    /** @brief refine/derefine algorithm for a 1-regular grid, geolevel of all
        cells on the finest grid is between MinLevel and MaxLevel */
    void Refine1Reg(int MinLevel, int MaxLevel);

    /** @brief derefinemnt */
    void DeRefine();

    /** @brief convert all finest quadrangles into two triangles */
    int ConvertQuadToTri(int type);
    
    int get_ref_level() const
    { return RefLevel; }

    /** @brief produce a collection with all cells returned by iterator it */
    TCollection *GetCollection(Iterators it, int level) const;
    
    /** @brief produce a collection with all cells of a given Collection which
     *         have the given reference as ReferenceID 
     * 
     * This will give you a subcollection.
     */
    TCollection *GetCollection(TCollection *coll, int reference);
    
    /**
     * @brief get collection of cells on a certain level and with certain 
     *        reference id
     */
    TCollection *GetCollection(Iterators it, int level, int ID) const;

#ifdef  _MPI 
    /** @brief produce a own collection with all cells returned by iterator it
     *
     * CB Produces TCollection of those cells whose flag "SubDomainNumber"
     * equals ID, I think this is meant to be the rank of this process (TODO CB Change variable name ID to rank)
     * */
    TCollection *GetOwnCollection(Iterators it, int level, int ID) const;
#endif

    /** @brief produce a collection with all cells in the finest grid, sort 
        they according to their geometry level and return in Indices 
        the indices where each level starts */
    void GetSortedCollection(TCollection* &Coll, int* &Indices);

    /** @brief get bounding box parameters */
    void GetBoundBox(double &startx, double &starty,
                     double &boundx, double &boundy)
    {
      startx = StartX;
      starty = StartY;
      boundx = BoundX;
      boundy = BoundY;
    }
    
#ifdef __3D__
    /** @brief get bounding box parameters */
    void GetBoundBox(double &startx, double &starty, double &startz,
                     double &boundx, double &boundy, double &boundz)
    {
      startx = StartX;
      starty = StartY;
      startz = StartZ;
      boundx = BoundX;
      boundy = BoundY;
      boundz = BoundZ;
    }
    
    void SetBoundBox(double startx, double starty, double startz,
                     double boundx, double boundy, double boundz)
    {
      StartX = startx;
      StartY = starty;
      StartZ = startz;
      BoundX = boundx;
      BoundY = boundy;
      BoundZ = boundz;
    }
#endif
    
    // test
    #ifndef __3D__
      void TestGrid1();
      void TestGrid2();
      void TestGrid3();
      void TestMortar();
      void TestShishkin();
      void TriangleShishkin();
      void UnitSquare();
      void UnitSquareRef();
      void TwoTriangles();
      void TwoTrianglesRef();
      void SquareInSquare();
      void SquareInSquareRef();
      void SetBoundBox(double boundx, double boundy);
      void SetBoundBoxstart(double startx , double starty);
      void RefOnMortarEdge();
      void RefCardioide(double A);
      void PeriodicSquares();
      void PeriodicSquaresLarge();
      void PeriodicRectangle_2_4();
      void PeriodicTrianglesLarge();
      void QuadShishkin(double tau1, double tau2);
      void Rectangular(int dimx, int dimy);
      void TestTriaConf();
      void TestTriaConf2();
      void CheckCells();

      /**
       * @brief Initialize the domain boundary as if data/UnitSquare.PRM
       * would have been read in as .PRM file.
       *
       * This is useful for testing purposes, when the program must be
       * independent of the path from working directory to "data" directory.
       */
      void initializeDefaultUnitSquareBdry();

    #else
      void TestGrid3D();
      /**
       * @brief Initialize the domain boundary as if data/Wuerfel.PRM
       * would have been read in as .PRM file.
       *
       * This method is useful for tests, which are supposed to be independent
       * of whether an extern .PRM-file is available. To let the domain call this
       * method when initializing, pass "Default_UnitCube" as first argument
       * to the init method.
       */
      int initializeDefaultCubeBdry();

      void SetBoundBox(double boundx, double boundy, double boundz);
    #endif

    #ifdef __3D__
      int Grape(const char *name, TCollection *coll);

    #endif

    void TetrameshGen();


    #ifdef  _MPI
      void ReplaceTreeInfo(int n_cells, TBaseCell **cells, int *GLOB_cellIndex, int n_OwnCells)
       {
        if(CellTree) delete CellTree;
        N_RootCells = n_cells;
        CellTree = cells;
        GlobalCellIndex = GLOB_cellIndex;
        N_OwnCells = n_OwnCells;
       }

      void SetN_OwnCells(int n_OwnCells)
       { N_OwnCells = n_OwnCells; }

      int GetN_OwnCells()
       { return N_OwnCells; }

      int GetN_HaloCells()
       { return (N_RootCells - N_OwnCells); }
     #endif
    
#ifdef __3D__
  public: 

//      int Tetgen(const char*);

     /** @brief generate edge info in 3D mesh **/
     // TODO CB Add documentation. This method does mark certain Vertices as boundary vertices
     // (TVertex::SetAsBoundVertex(...) ) and creates edge objects, pointers to which are then
     // stored as members of certain mesh cells (TBaseCell::SetEdge(...)).
     // Does work only on currently finest level.
     int GenerateEdgeInfo();

#endif

  /** @brief adaptive refine  */
  int AdaptRefineAll();   

  /**
   * @brief Checks from the file name whether a .GEO -file will be read in or
   * a .xGEO (extended GEO) file.
   *
   * This code was moved here from the ReadGeo method.
   *
   * @param[in] GEO the filename of the .(x)GEO file.
   *
   * @note This has not been tested with an actual .xGEO-file yet.
   */
  static bool checkIfxGEO(char* GEO);
     
     
};

#endif
