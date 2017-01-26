
void ExampleFile()
{
  OutPut("Example: 32_Variable_viscosity_Beltrami.h" << endl);
}

double get_nu()
{
  return 1;
}

constexpr bool rhs_depends_on_time = true;
constexpr bool coefficients_depend_on_time = false;

// exact solution
void Exact(double x, double y, double *values)
{
//  double t=TDatabase::TimeDB->CURRENTTIME;
  int switchviscosity = TDatabase::ParamDB->P1;
  double vmin = TDatabase::ParamDB->P2;
  double vmax = TDatabase::ParamDB->P3;

  switch(switchviscosity)
  {
    case 1: // corresponds to eps1=viscosity1
      values[0] = vmin+(vmax-vmin)*x*x*(1-x)*y*y*(1-y)*721/16;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      break;
    case 2: // corresponds to eps2=viscosity2
      values[0] = vmin+(vmax-vmin)*exp(-10e13*(pow((x-0.5),10)+pow((y-0.5),10)));
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      break;
    case 3: // corresponds to eps3=viscosity3
      values[0] = vmin+(vmax-vmin)*(1-exp(-10e13*(pow((x-0.5),10)+pow((y-0.5),10))));
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      break;
    default:
      ErrThrow("Choose parameter P1 between 1 and 3.");
  }
}

// kind of boundary condition (for FE space needed)
void BoundCondition(int BdComp, double Param, BoundCond &cond)
{
  cond = DIRICHLET;
}

// value of boundary condition
void BoundValue(int BdComp, double Param, double &value)
{
//  double t=TDatabase::TimeDB->CURRENTTIME;
  int switchviscosity = TDatabase::ParamDB->P1;
  double vmin = TDatabase::ParamDB->P2;
  double vmax = TDatabase::ParamDB->P3;

  switch(BdComp)
  {
    case 0:
      switch(switchviscosity)
      {
        case 1: // viscosity 1 at y=0
          value = vmin;
          break;
        case 2: // viscosity 2 at y=0
          value = vmin+(vmax-vmin)*exp(-10e13*(pow((Param-0.5),10)+pow((0-0.5),10)));
          break;
        case 3: // viscosity 3 at y=0
          value = vmin+(vmax-vmin)*(1-exp(-10e13*(pow((Param-0.5),10)+pow((0-0.5),10))));
          break;
        default:
          ErrThrow("Choose parameter P1 between 1 and 3.");
      }
    break;
    case 1:
      switch(switchviscosity)
      {
        case 1: // viscosity 1 at x=1
          value = vmin;
          break;
        case 2: // viscosity 2 at x=1
          value = vmin+(vmax-vmin)*exp(-10e13*(pow((1-0.5),10)+pow((Param-0.5),10)));
          break;
        case 3: // viscosity 3 at x=1
          value = vmin+(vmax-vmin)*(1-exp(-10e13*(pow((1-0.5),10)+pow((Param-0.5),10))));
          break;
        default:
          ErrThrow("Choose parameter P1 between 1 and 3.");
      }
    break;
    case 2:
      switch(switchviscosity)
      {
        case 1: // viscosity 1 at y=1
          value = vmin;
          break;
        case 2: // viscosity 2 at y=1
          value = vmin+(vmax-vmin)*exp(-10e13*(pow(((1-Param)-0.5),10)+pow((1-0.5),10)));
          break;
        case 3: // viscosity 3 at y=1
          value = vmin+(vmax-vmin)*(1-exp(-10e13*(pow(((1-Param)-0.5),10)+pow((1-0.5),10))));
          break;
        default:
          ErrThrow("Choose parameter P1 between 1 and 3.");
      }
    break;
    case 3:
      switch(switchviscosity)
      {
        case 1: // viscosity 1 at x=0
          value = vmin;
          break;
        case 2: // viscosity 2 at x=0
          value = vmin+(vmax-vmin)*exp(-10e13*(pow((0-0.5),10)+pow(((1-Param)-0.5),10)));
          break;
        case 3: // viscosity 3 at x=0
          value = vmin+(vmax-vmin)*(1-exp(-10e13*(pow((0-0.5),10)+pow(((1-Param)-0.5),10))));
          break;
        default:
          ErrThrow("Choose parameter P1 between 1 and 3.");
      }
    break;
  } // endswitch
}

// initial conditon
void InitialCondition(double x,  double y, double *values)
{
//  double t=TDatabase::TimeDB->CURRENTTIME;
  int switchviscosity = TDatabase::ParamDB->P1;
  double vmin = TDatabase::ParamDB->P2;
  double vmax = TDatabase::ParamDB->P3;

  switch(switchviscosity)
  {
    case 1: // corresponds to eps1=viscosity1
      values[0] = vmin+(vmax-vmin)*x*x*(1-x)*y*y*(1-y)*721/16;
      break;
    case 2: // corresponds to eps2=viscosity2
      values[0] = vmin+(vmax-vmin)*exp(-10e13*(pow((x-0.5),10)+pow((y-0.5),10)));
      break;
    case 3: // corresponds to eps3=viscosity3
      values[0] = vmin+(vmax-vmin)*(1-exp(-10e13*(pow((x-0.5),10)+pow((y-0.5),10))));
      break;
    default:
      ErrThrow("Choose parameter P1 between 1 and 3.");
  }
}

void BilinearCoeffs(int n_points, double *X, double *Y,
        double **parameters, double **coeffs)
{
//  double eps=1/TDatabase::ParamDB->RE_NR;
  int i;
  double *coeff, x, y;
//  double t=TDatabase::TimeDB->CURRENTTIME;
  int switchviscosity = TDatabase::ParamDB->P1;
  double vmin = TDatabase::ParamDB->P2;
  double vmax = TDatabase::ParamDB->P3;

  for(i=0;i<n_points;i++)
  {
    coeff = coeffs[i];
    x = X[i];
    y = Y[i];

    coeff[0] = 0;
    coeff[1] = 1;
    coeff[2] = 1;
    coeff[3] = 0;

    // rhs = 1.d/dx + 1.d/dy
    switch(switchviscosity)
    {
      case 1:   // rhs for viscosity 1
        coeff[4] = (721/16)*(vmax-vmin)*x*y*(-2*(-1+y)*y+x*x*(-2+3*y)+x*(2-6*y+3*y*y));
        break;
      case 2:   // rhs for viscosity 2
        coeff[4] = -(10.e14)*exp(-10e13*(pow((x-0.5),10)+pow((y-0.5),10)))*(vmax-vmin)*(pow((x-0.5),9)+pow((y-0.5),9));
        break;
      case 3:   // rhs for viscosity 3
        coeff[4] = +(10.e14)*exp(-10e13*(pow((x-0.5),10)+pow((y-0.5),10)))*(vmax-vmin)*(pow((x-0.5),9)+pow((y-0.5),9));
        break;
      default:
        ErrThrow("Choose parameter P1 between 1 and 3.");
    }
  }
}

// exact solution
void Initial(double x, double y, double *values)
{
//  double t=TDatabase::TimeDB->CURRENTTIME;
  int switchviscosity = TDatabase::ParamDB->P1;
  double vmin = TDatabase::ParamDB->P2;
  double vmax = TDatabase::ParamDB->P3;

  switch(switchviscosity)
  {
    case 1: // corresponds to eps1=viscosity1
      values[0] = vmin+(vmax-vmin)*x*x*(1-x)*y*y*(1-y)*721/16;
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      break;
    case 2: // corresponds to eps2=viscosity2
      values[0] = vmin+(vmax-vmin)*exp(-10e13*(pow((x-0.5),10)+pow((y-0.5),10)));
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      break;
    case 3: // corresponds to eps3=viscosity3
      values[0] = vmin+(vmax-vmin)*(1-exp(-10e13*(pow((x-0.5),10)+pow((y-0.5),10))));
      values[1] = 0;
      values[2] = 0;
      values[3] = 0;
      break;
    default:
      ErrThrow("Choose parameter P1 between 1 and 3.");
  }
}

