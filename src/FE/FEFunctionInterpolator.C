/**
 * Implements class FEFunctionInterpolator declared in FEFunctioninterpolator.h
 *
 * @author Clemens Bartsch
 * @date 2016/02/05
 */

#include <FEFunctionInterpolator.h>
#include <MooNMD_Io.h>

#include <Database.h>
#include <FEDatabase2D.h>
#include <FEDatabase3D.h>

#include <FESpace1D.h>
#include <FESpace2D.h>
#include <FESpace3D.h>

#include <FEFunction2D.h>
#include <FEFunction3D.h>

#include <InterfaceJoint.h>

#include <memory>

// Experimental Macro to avoid double code.
#ifdef __2D__
#define TFESpaceXD TFESpace2D
#elif __3D__
#define TFESpaceXD TFESpace3D
#endif

FEInterpolationCheatSheet::FEInterpolationCheatSheet(
        const TFESpaceXD* old_fe_space, const TFESpaceXD* new_fe_space)
{
  size_t n_cells_new_space = new_fe_space->GetN_Cells();
  cheat_sheet= std::vector<std::vector<ContainingCells>>(n_cells_new_space);

  for(size_t c_new = 0; c_new < n_cells_new_space;++c_new)
  {
    TBaseCell* cell = new_fe_space->GetCollection()->GetCell(c_new);
#ifdef __2D__
    FE2D FEId = new_fe_space->GetFE2D(c_new, cell);
    TFE2D* Element = TFEDatabase2D::GetFE2D(FEId);
    TNodalFunctional2D* nf = Element->GetNodalFunctional2D();
    double* xi; //positions of the dof in the ref element
    double* eta;
#elif __3D__
    ErrThrow("Not yet implemented in 3D!");
#endif
    int n_points;
    nf->GetPointsForAll(n_points, xi, eta); //now we found out the number of nf points

    std::vector<ContainingCells> cheats_for_cell(n_points);
    cheat_sheet.at(c_new) = cheats_for_cell;

    BF2DRefElements RefElement = Element->GetBaseFunct2D()->GetRefElement();
    int N_Edges = 0;
    switch(RefElement)
    {
      case BFUnitSquare:
        N_Edges = 4;
        break;

      case BFUnitTriangle:
        N_Edges = 3;
        break;
    }

    RefTrans2D RefTrans;
    RefTrans2D* RefTransArray;
    RefTransArray = TFEDatabase2D::GetRefTrans2D_IDFromFE2D();
    RefTrans = RefTransArray[FEId];

    bool IsIsoparametric(false);
    if (TDatabase::ParamDB->USE_ISOPARAMETRIC)
    {
      for(int j=0; j<N_Edges ;j++)
      {
        TJoint* joint = cell->GetJoint(j);
        JointType jointtype = joint->GetType();
        if(jointtype == BoundaryEdge)
        {
          BoundTypes bdtype = ((TBoundEdge *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == InterfaceJoint)
        {
          BoundTypes bdtype = ((TInterfaceJoint *)(joint))->GetBoundComp()->GetType();
          if(bdtype != Line)
            IsIsoparametric = true;
        }
        if(jointtype == IsoInterfaceJoint ||
          jointtype == IsoBoundEdge)
          IsIsoparametric = true;
      }
    }

    if(IsIsoparametric)
    {
      switch(RefElement)
      {
        case BFUnitSquare:
          RefTrans = QuadIsoparametric;
          break;

        case BFUnitTriangle:
          RefTrans = TriaIsoparametric;
          break;
      }
    }

    TFEDatabase2D::SetCellForRefTrans(cell, RefTrans);
    //finally: get the position of the actual dof
    std::vector<double> X(MaxN_PointsForNodal2D,0.0);
    std::vector<double> Y(MaxN_PointsForNodal2D,0.0);
    double AbsDetjk[MaxN_PointsForNodal2D];

    TFEDatabase2D::GetOrigFromRef(RefTrans, n_points, xi, eta, &X.at(0), &Y.at(0), AbsDetjk);

    for(int p = 0 ; p < n_points ; ++p)
    {//now find out in which cell of the old fe function (X[p],Y[p]) lies
      int n_old_cells = old_fe_space->GetCollection()->GetN_Cells();
      for(int c_old = 0; c_old < n_old_cells; ++c_old)
      {
        if(old_fe_space->GetCollection()->GetCell(c_old)->PointInCell(X.at(p),Y.at(p)))
        {//it's in c_old - push c_old back to the list of cells which contaiin point p!
          cheat_sheet.at(c_new).at(p).push_back(c_old);
        }
      }
    }
  }
}

FEFunctionInterpolator::FEFunctionInterpolator(
    const TFESpace* into_space)
{
  // try to downcast the input pointer to figure out the dimension of the space
  // (a nice alternative would be a virtual get_dimension() method in TFESpace
  //
  // TODO: syntax for shared pointers would be
  //  if(std::dynamic_pointer_cast<TFESpace1D>(fe_space))
  //  {
  //    dimension_ = 3;
  //  }
  //
  //
  if(dynamic_cast<const TFESpace1D*>(into_space))
  {
    dimension_ = 1;
  }
  else if(dynamic_cast<const TFESpace2D*>(into_space))
  {
    dimension_ = 2;
  }
#ifdef __3D__
  else if(dynamic_cast<const TFESpace3D*>(into_space))
  {
    dimension_ = 3;
  }
#endif
  else
  {
    ErrThrow("The fe space to interpolate to must be of one of the"
        " derived types TFESpace1D, TFESpace2D or TFESpace3D!");
  }

  //store the shared_ptr to the space to interpolate to
  into_space_ = into_space;

}

/** ************************************************************************ */

void FEFunctionInterpolator::check() const
{
  //call check for casting of the stored space
  check_for_shearing();
}

/** ************************************************************************ */

TFEFunction2D FEFunctionInterpolator::interpolate(
    const TFEFunction2D& original_funct,
    std::vector<double>& values_memory,
    bool use_cheat_sheet) const
{
  //check if there is enough memory given
  if((int) values_memory.size() != into_space_->GetN_DegreesOfFreedom())
    ErrThrow("Wrong number of double memory given to FEFunctionInterpolator::interpolate.");

  double* values = &values_memory.front();
  int length = values_memory.size();
  const TFESpace2D* into_space_cast = dynamic_cast<const TFESpace2D*>(into_space_);

  TFEFunction2D interpolation(
      into_space_cast, (char*) "interpolated", (char*) "interpolated", values, length);

  if(use_cheat_sheet)
  {
    if(!cheat_sheet_)
    {//write a cheat sheet
      cheat_sheet_ = std::make_shared<FEInterpolationCheatSheet>(
          original_funct.GetFESpace2D() ,into_space_cast);
    }
    try
    {
      interpolation.Interpolate(&original_funct, cheat_sheet_.get());
    }
    catch(std::runtime_error &e)
    {
      // TODO Whenever an exception is caught here, that is a sign, that
      // the space this thing should interpolate from did change. One could think
      // about printing a warning and resetting the cheat_sheet instead of
      // ending the program here
      Output::print("Caught an exception calling interpolation.Interpolate(&original_funct, cheat_sheet_.get())");
      Output::print(e.what());
      exit(-1);
    }
  }
  else
  {
    // make use of method "Interpolate" for the moment
    // - TODO that code should be moved here entirely in time
    interpolation.Interpolate(&original_funct);
  }



  return interpolation;
}

// ////////////////////// //    private method(s)   // ////////////////////// //

void FEFunctionInterpolator::check_for_shearing() const
{
  // switch over the dimension and check if the dynamic cast
  // to the correct derived FESpace can be done
  switch (dimension_)
  {
    case 1:
      if(!dynamic_cast<const TFESpace1D*>(into_space_))
      {
        ErrThrow("In dimension 1 the cast to TFESpace1D must be possible!");
      }
      else
      {
        Output::print<1>("Correct cast to TFESpace1D.");
      }
      break;
    case 2:
      if(!dynamic_cast<const TFESpace2D*>(into_space_))
      {
        ErrThrow("In dimension 2 the cast to TFESpace2D must be possible!");
      }
      else
      {
        Output::print<1>("Correct cast to TFESpace2D.");
      }
      break;
#ifdef __3D__
    case 3:
      if(!dynamic_cast<const TFESpace3D*>(into_space_))
      {
        ErrThrow("In dimension 3 the cast to TFESpace3D must be possible!");
      }
      else
      {
        Output::print<1>("Correct cast to TFESpace3D.");
      }
      break;
#endif
    default:
      ErrThrow("What is that? Neither dimension 1, 2 nor 3?"
          " Something's terribly wrong here.")
  }
}
#undef TFESpaceXD
