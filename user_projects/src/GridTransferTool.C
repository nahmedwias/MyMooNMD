#include <FEFunction2D.h>
#include <FEFunctionInterpolator.h>
#include <FESpace2D.h>
#include <GridTransferTool.h>



// Get the barycenter of a ParMooN grid cell.
// This is of course tatlly misplaced here.
std::valarray<double> center_point(const TBaseCell& cell)
{
#ifdef __3D__
  ErrThrow("Does center point calculation work in 3D?");
#endif
  std::valarray<double> p(0.0,3);

  unsigned int n_verts = cell.GetN_Vertices();
  for(unsigned int v = 0; v < n_verts; v++)
  {
    cell.GetVertex(v)->GetX();
    p[0] += cell.GetVertex(v)->GetX();
    p[1] += cell.GetVertex(v)->GetY();
  }
  p[0] /= n_verts;
  p[1] /= n_verts;

  return p;
}

GridTransferTool::GridTransferTool(
        const TFESpace2D* to_space,
        GridTransferType type) :
        to_space_(to_space), type_(type)
{
  if (type == GridTransferType::MultiGrid)
  {
    //TODO
    ErrThrow("Multigrid type grid transfer is to be implemented");
  }
}

void GridTransferTool::transfer(
    const TFEFunction2D& input_fct, TFEFunction2D& output_fct,
    std::vector<double>& output_fct_values) const
{
  if(type_ == GridTransferType::MidPointEvaluation)
  {
    //check input
    if(output_fct.GetFESpace2D() != this->to_space_)
      ErrThrow("Space of output fe function and stored fe space do not agree.");

    // find out whether there was a call with this space before
   std::shared_ptr<FEFunctionInterpolator> intpol;
   for(auto it : interpolators_)
   {
     if(it->get_from_space() == input_fct.GetFESpace2D())
     {//found!
       Output::print<5>("Interpolator found, will be reused.");
       intpol = it;
     }
   }
   if(!intpol)
   {//have to set it up new
     Output::print<5>("New interpolator set up.");
     intpol = std::make_shared<FEFunctionInterpolator>(output_fct.GetFESpace2D());
     interpolators_.push_back(intpol);
   }
   //do the actual interpolation
   intpol->interpolate(input_fct,output_fct_values,true);

  }
  else if (type_ == GridTransferType::MultiGrid)
  {
    //TODO
    ErrThrow("Multigrid type grid transfer is to be implemented");
  }
}




