// =======================================================================
//
// Purpose:     main program for solving a stationary Brinkman equation in ParMooN
//
// Author:      Laura Blank
//
// History:     Implementation started on  11.10.2016

// =======================================================================
#include <Domain.h>
#include <Database.h>
#include <FEDatabase3D.h>
#include <LinAlg.h>
#include <Brinkman3D.h>
#include <MainUtilities.h>
#include <LocalAssembling3D.h>
#include <Example_Brinkman3D.h>
#include <ParameterDatabase.h>
#include <MooNMD_Io.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <PETScSolver.h>

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
    // Declaration of database, this is needed in every program
    TDatabase Database;
    TFEDatabase3D FEDatabase;
    
    ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
    // read in the data from the input file (argv[1]) (here: brinkman3d.dat)
    std::ifstream fs(argv[1]);
    parmoon_db.read(fs);
    fs.close();
    
    // Construct domain, thereby read in controls from the input file (argv[1]).
    TDomain domain(argv[1], parmoon_db);
    
    // Produce an outfile "... .out", this is where all output is written to (addionally to console)
    Output::set_outfile(parmoon_db["outfile"]);
    
    // Write all parameters to the outfile (not to console) for later reference
    parmoon_db.write(Output::get_outfile());
    Database.WriteParamDB(argv[0]);
    
    // Refine grid
    size_t n_ref =  domain.get_n_initial_refinement_steps();
    for(size_t i=0; i<= n_ref; i++)
    {
        if(n_ref>0) {
            domain.RegRefineAll();
            Output::print("Level:", i+1);
        }
    //}
        
        // Write grid into a Postscript file
        if(parmoon_db["output_write_ps"])
            domain.PS("Domain.ps", It_Finest, 0);
        
        // Create an example object according to the database
        Example_Brinkman3D example(parmoon_db);
        
        // Intial refinement and grabbing of grids for multigrid.
#ifdef _MPI
        int maxSubDomainPerDof = 0;
#endif
        
//        std::list<TCollection* > gridCollections = domain.refine_and_get_hierarchy_of_collections(parmoon_db
//#ifdef _MPI
//                                                         ,maxSubDomainPerDof
//#endif
//                                                         );
        
        //=========================================================================
        // Create an object of the Brinkman class
#ifdef _MPI
        Brinkman3D brinkman3d(domain, parmoon_db, example, maxSubDomainPerDof);
#else
        Brinkman3D brinkman3d(domain, parmoon_db, example);
#endif
        
        //#ifdef _MPI
        //        Brinkman3D brinkman3d(gridCollections, parmoon_db, example, maxSubDomainPerDof);
        //#else
        //        Brinkman3D brinkman3d(gridCollections, parmoon_db, example);
        //#endif
        //        Output::print<>("Database-Info:");
        parmoon_db.info();
        
        brinkman3d.assemble();
        brinkman3d.solve();
        
        // brinkman3d.solve_with_Petsc(parmoon_db);
      
        
        //        Output::print<>(TDatabase::ParamDB->VELOCITY_SPACE);
        //        Output::print<>(TDatabase::ParamDB->PRESSURE_SPACE);
        
        brinkman3d.output();
        
        //=========================================================================
    }
    
    Output::close_file();
    return 0;
    // }
} // end main
