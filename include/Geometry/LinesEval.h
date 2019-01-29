/******************************************************************************* 
*
* @class LineEval
* @date  19.12.18
* @brief Some tools for function evaluation over line
* @author Baptiste Moreau
* @History:
*
*******************************************************************************/

#ifndef __EVALTOOLS__
#define __EVALTOOLS__

#include <array>
#include <Collection.h>
#include "templateNames.h"
#include <ParameterDatabase.h>

class LineEval;

/**
 * @brief define all lines used for evaluation
 *
 * The single purpose of this class is to evaluate fe_functions on points which
 * are on lines. These values are then written to a file, see 
 * LinesEval::write_fe_values.
 */
template <int d>
class LinesEval
{
  protected:
    /** @brief A Database object which holds parameters which are of a certain
     *  interest to this LinesEval object. */
    ParameterDatabase db;

    /** @brief contain all lines used for evaluations */
    std::vector<LineEval> lines_for_postprocess;

    /** @brief in case of line defined from extern file, read the direction,
     * base_point and positions to define the lines
     */
    void read_position(const std::string    filename,
                       std::vector<int>&    direction,
                       std::vector<double>& position,
                       std::vector<double>& coord,
                       std::vector<int>&    n_coord);

    /** @brief in case of line defined from extern file, check if the direction,
     * base_point and positions are consistent
     */
    void check_position(const std::string         filename,
                        const std::vector<int>    direction,
                        const std::vector<double> position,
                        std::vector<int>&         n_coord);


  public:
    using FEFunction = typename Template_names<d>::FEFunction;

    /**
     * @brief default constructor
     */
    LinesEval();

    /**
     * @brief constructor, the lines are defined by database entries
     */
    LinesEval(const TDomain& domain,
              const ParameterDatabase& param_db);

    /**
     * Creates a database filled with default parameters. This database will
     * contain all necessary parameters for the LinesEval tools.
     */
    static ParameterDatabase default_lineseval_parameters();
    
    /**
     * The database passed to the constuctor should be named this or have a 
     * nested database with this name. Otherwise an exception is thrown.
     */
    constexpr static char required_database_name[] = "LineEval Database";

    /**
     * @brief write the values of the given fe functions at the points on the 
     * lines
     */
    void write_fe_values(std::vector<const FEFunction*> fe_functions,
                         double      time_step,
                         std::string file_prefix) const;
    
    /**
     * @brief return the number of defines lines
     */
    int GetLength() const;

    /**
     * @brief return the line i
     */
    LineEval GetLine(int i) const;
};


/**
 * @brief define positions on line used for evaluation
 */
class LineEval
{
  protected:
    
    /** @brief the needed collection to find the cells */
    TCollection* coll;
    
    /** @brief direction of the line */
    int direction;

    /** @brief point throuh which the line is passing */
    std::array<double, 3> base_point;

    /** @brief number of refinement points:
     *       - between two points of the extern file defining the line
     *       - or in the cells if the line is defined from a nested database */
    int n_refine;

    /** @brief flag to check if the line is defined from a nested database or
     * not, i.e. from an extern file */
    bool fromDB;

    /** @brief structure for additional informations about the points:
     * the index of the cell to which the point belongs and
     * the two intersection positions between the line and the cell,
     * in case of line defined by an extern file, only lmin_cell is used to 
     * store the position, lmax_cell is set to 0. */
    struct point_on_line
    {
      int    cell_index;
      double lmin_cell;
      double lmax_cell;

      /**
       * @brief default constructor
       */
      point_on_line();

      /**
       * @brief constructor
       */
      point_on_line(int cell_index, double lmin_cell, double lmax_cell)
      {
        this->cell_index = cell_index;
        this->lmin_cell  = lmin_cell;
        this->lmax_cell  = lmax_cell;
        
      }
    };

    /** @brief contain position informations used for evaluations along
     * the line */
    std::vector<point_on_line> line_for_postprocess;
    
    /**
     * @brief return the mean value of the function f along a line defined from
     * an extern file
     */
    double mean_value(const TFEFunction3D& f) const;


  public:

    /**
     * @brief constructor for lines defined by extern file
     * 
     * @param[in] domain    domain providing the collection
     * @param[in] direction direction of the line
     * @param[in] P[3]      point throuh which the line is passing
     * @param[in] coord     positions on the line where evaluation is to be done
     * @param[in] refine    number of refinement between two positions
     */
    LineEval(const TDomain&      domain,
             const int           direction,
             const double        P[3],
             std::vector<double> coord,
             const int           refine);

    /**
     * @brief constructor for lines defined by nested database
     *
     * @param[in] domain    domain providing the collection
     * @param[in] direction direction of the line
     * @param[in] P[3]      point throuh which the line is passing
     * @param[in] refine    number of refinement between two positions
     */
    LineEval(const TDomain& domain,
             const int      direction,
             const double   P[3],
             const int      refine);

    /** @brief return the base point which defines the line */
    std::array<double, 3> GetBasePoint() const;

    /**
     * @brief return the direction which defines the line
     */
    int GetDirection() const;
    
    /**
     * @brief return the number of points which belongs to the line (after
     * refinement)
     */
    int GetNbPoints() const;

    /**
     * @brief return the index of the cell to which the point i belongs (after
     * refinement, i.e. for i in [0 GetNbPoints()])
     */
    int GetCellIdx(int i) const;
    
    /**
     * @brief return the pointer to the cell with index cell_idx
     */
    const TBaseCell* GetCell(int cell_idx) const;

    /**
     * @brief return the position on the line of the point i (after refinement,
     * i.e. for i in [0 GetNbPoints()])
     */
    double GetPosition(int i) const;

    /**
     * @brief check if the line is defines from an extern file or from a nested
     * database
     */
    bool IsFromDB() const;

    /**
     * @brief return the mean value of the function f along the line
     */
    double space_average_value(const TFEFunction3D& f) const;

    /**
     * This operator introduces an order on the points. It will
     * compare the lmin_cell value.
     */
    friend bool operator < (const point_on_line& P, const point_on_line& Q);

    /**
     * This operator compare if the lmin_cell values are equal.
     * The lmin_cell are regarded as equal with a tolerance of 1e-16.
     */
    friend bool operator == (const point_on_line& P, const point_on_line& Q);
};

#endif
