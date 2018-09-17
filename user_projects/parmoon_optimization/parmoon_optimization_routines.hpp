#include "ParameterDatabase.h"
#include "MooNMD_Io.h"
#include <nlopt.hpp>
#include <algorithm>

namespace parmoon_opt
{
  /// @brief return a database with all parameters related to optimization 
  ParameterDatabase default_optimization_database();
  
  /// @brief wrapper for nlopt.
  /// The idea is to create a class (Opt below) which implements a method
  /// 'compute_functional_and_derivative' as used below. Then, the library
  /// 'nlopt' can use this wrapper as objective functional to be minimized.
  template <class Opt>
  double evaluate_functional(unsigned n, const double *x, double *grad,
                             void *data)
  {
    Output::print(); // newline
    Output::info<5>("parmoon_opt::evaluate_functional,", "Entering wrapper "
        "from parmoon_opt namespace.");
    auto this_object = reinterpret_cast<Opt*>(data);
    return this_object->compute_functional_and_derivative(n, x, grad);
  }
  
  /// @brief print a summary of the optimization after it has finished.
  /// This includes reason for convergence/failure as well as the obtained 
  /// minimum value for the functional, the number of iterations, and the name
  /// of the used algorithm. If 'print_final_solution' is true the entire
  /// solution vector 'x' is printed.
  void print_summary(const nlopt::opt& opt, std::vector<double> x, 
                     bool print_final_solution);
  /// @brief map the result of an optimization to a string for easier printing
  std::string termination_reason(nlopt::result result);
  /// @brief set values from ParameterDatabase (especially related to stopping
  /// criteria) in the nlopt::opt object.
  void set_values_nlopt(nlopt::opt& opt, const ParameterDatabase& opt_db);
  
  template <class Opt>
  void optimize(Opt& bco, const ParameterDatabase& opt_db_in)
  {
    ParameterDatabase opt_db = default_optimization_database();
    opt_db.merge(opt_db_in, false);
    // choose algorithm and dimensionality
    unsigned n_control = bco.get_n_control();
    auto algorithm = (nlopt::algorithm)opt_db["algorithm_nlopt"].get<size_t>();
    nlopt::opt opt(algorithm, n_control);
    set_values_nlopt(opt, opt_db);
    // set the routine to evaluate the functional and derivative
    opt.set_min_objective(evaluate_functional<Opt>, &bco);
    // pick initial guess  (should this rather be in the class Opt)
    std::vector<double> x(n_control, opt_db["initial_control"].get<double>());
    //x[5] = 0.; x[6] = 0.; x[7] = 0.; x[8] = 0.; x[9] = 0.;
    
    double minf; // the minimum functional value upon return
    // call the optimization routine. Note that there is also a variant without
    // explicitly passing minf, but then there will be a copy of 'x'.
    try{ opt.optimize(x, minf); }
    catch(...){Output::print("EXCEPTION");};
    print_summary(opt, x, opt_db["print_final_solution"]);
  }
  
  //implementation for the methods defined above:
  ParameterDatabase default_optimization_database()
  {
    ParameterDatabase opt_db("optimization parameters");
    opt_db.add("algorithm_nlopt", 40u, 
               "The algorithm 'nlopt' should use. Please check the enum "
               "'algorithm' in the corresponding header file 'nlopt.hpp'.",
               0u, unsigned(nlopt::NUM_ALGORITHMS-1));
    opt_db.add("max_n_evaluations", 10u, "Maximum number of evaluations which "
               "serves as a stopping criterion. These evaluations include the "
               "functional as well as its derivative (ie, primal and adjoint "
               " solve).", 0u, 1000000u);
    opt_db.add("set_lower_and_upper_bounds", false,
               "Do or don't set the upper and lower bounds as a constant for "
               "each component of the control vector. The actual bounds are "
               "set using the parameters 'lower_bound' and 'upper_bound'");
    opt_db.add("lower_bound", 0.0, 
               "The lower bound which shall be fulfilled by each component of "
               "the control vector. To activate this number, set the parameter "
               "'set_lower_and_upper_bounds' to true.");
    opt_db.add("upper_bound", 0.0, 
               "The upper bound which shall be fulfilled by each component of "
               "the control vector. To activate this number, set the parameter "
               "'set_lower_and_upper_bounds' to true.");
    opt_db.add("ftol_rel", 0.0, 
               "Set the relative difference of the functional evaluation as a "
               "stopping criterion, i.e., something like "
               "|f_old-f_new|/|f_new|.",
               0.0, 1.e10);
    opt_db.add("ftol_abs", 0.0, 
               "Set the absolute difference of the functional evaluation as a "
               "stopping criterion, i.e., something like |f_old-f_new|.",
               0.0, 1.e10);
    opt_db.add("xtol_rel", 0.0, 
               "Set the relative difference of the control vector as a stopping"
               " criterion, i.e., something like |x_old-x_new|/|x_new|.",
               0.0, 1.e10);
    opt_db.add("xtol_abs", 0.0, 
               "Set the absolute difference of the control vector as a "
               "stopping criterion, i.e., something like |x_old-x_new|.",
               0.0, 1.e10);
    opt_db.add("stop_val", -1.e10, 
               "Stop the optimization loop if the functional value is below "
               "this value. This is useful to compare different optimization "
               "algorithms.", -1.e10, 1.e10);
    opt_db.add("n_vectors_stored", 10u, 
               "The number of vectors stored during the optimization. Some "
               "Algorithms (LBFGS, Truncated Newton) build up approximations "
               "to the Hessian. The more vectors are stored, the better the "
               "approximation, but also more memory is needed.", 0u, 1000000u);
    opt_db.add("initial_control", 0.0,
               "Choose a constant value for all components of the control "
               "vector at the start of the optimization loop.");
    opt_db.add("print_final_solution", false, 
               "Print out all components of the solution vector after the "
               "optimization algorithm finished.");
    
    return opt_db;
  }
  
  void set_values_nlopt(nlopt::opt& opt, const ParameterDatabase& opt_db)
  {
    if(opt_db["set_lower_and_upper_bounds"])
    {
      opt.set_upper_bounds(opt_db["upper_bound"].get<double>());
      opt.set_lower_bounds(opt_db["lower_bound"].get<double>());
      
      //opt.set_upper_bounds({10.,10.,10.,10.,10., 1.,1.,1.,1.,1.});
      //opt.set_lower_bounds({0.,0.,0.,0.,0.,-1.,-1.,-1.,-1.,-1.});
    }
    opt.set_ftol_rel(opt_db["ftol_rel"]);
    opt.set_ftol_abs(opt_db["ftol_abs"]);
    opt.set_xtol_rel(opt_db["xtol_rel"]);
    opt.set_xtol_abs(opt_db["xtol_abs"].get<double>());
    opt.set_maxeval(opt_db["max_n_evaluations"]);
    opt.set_stopval(opt_db["stop_val"]);
    opt.set_vector_storage(opt_db["n_vectors_stored"]);
  }
  
  void print_summary(const nlopt::opt& opt, std::vector<double> x, 
                     bool print_final_solution)
  {
    auto result = opt.last_optimize_result();
    if(result < 0)
    {
      Output::warn("optimization", "nlopt failed! ", result );
    }
    Output::print("optimization done because ", termination_reason(result));
    Output::print("found minimum ", opt.last_optimum_value(), "  after ",
                  opt.get_numevals(), " evaluations with the algorithm ",
                  opt.get_algorithm_name());
    double norm_solution = std::accumulate(x.begin(), x.end(), 0.0, 
                                           [](double a, double b) 
                                           {return a + b * b;});
    Output::print("norm of solution vector: ", std::sqrt(norm_solution));
    if(print_final_solution)
    {
      for(auto i = 0ul, n_control = x.size(); i < n_control; ++i)
      {
        Output::print("solution[", i, "] = ", x[i]);
      }
    }
  }
  
  std::string termination_reason(nlopt::result result)
  {
     switch(result)
     {
       case -1: return "generic failure"; break;
       case -2: return "invalid arguments"; break;
       case -3: return "out of memory"; break;
       case -4: return "roundoff limited"; break;
       case -5: return "forced stop"; break;
       case 1: return "success"; break;
       case 2: return "stopping value reached"; break;
       case 3: return "ftol reached"; break;
       case 4: return "xtol reached"; break;
       case 5: return "maximum number of evaluations reached"; break;
       case 6: return "maximum computing time reached"; break;
       default:
         return "unknown result " + std::to_string(result); break; 
    }
  }
  
  void print_nlopt_version()
  {
    Output::print("nlopt version: ", nlopt::version_major(), ".", 
                  nlopt::version_minor(), ".", nlopt::version_bugfix());
  }
}
