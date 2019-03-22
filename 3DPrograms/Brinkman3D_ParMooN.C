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
#include <Brinkman3D.h>
#include <Example_Brinkman3D.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <Chrono.h> //for a stopwatch which measures time spent in program parts

// =======================================================================
// main program
// =======================================================================
int main(int argc, char* argv[])
{
    
#ifdef _MPI
    //Construct and initialise the default MPI communicator.
    MPI_Init(&argc, &argv);
#endif
    {
#ifdef _MPI
        MPI_Comm comm = MPI_COMM_WORLD;
        int my_rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        if(my_rank==0)
        {
            Output::print("<<<<< Running ParMooN: Brinkman3D Main Program >>>>>");
            Output::info("Brinkman3D", "MPI, using ", size, " processes");
        }
#else
        int my_rank = 0;
        Output::print("<<<<< Running ParMooN: Brinkman3D Main Program >>>>>");
        Output::info("Brinkman3D", "SEQUENTIAL (or OMP...)");
#endif
        
        //start a stopwatch which measures time spent in program parts
        Chrono timer;
        
        
        
        // Declaration of database, this is needed in every program
        TDatabase Database(argv[1]);
        TFEDatabase3D FEDatabase;
        
        ParameterDatabase parmoon_db = ParameterDatabase::parmoon_default_database();
        // read in the data from the input file (argv[1]) (here: brinkman3d.dat)
        std::ifstream fs(argv[1]);
        parmoon_db.read(fs);
        fs.close();
        
#ifdef _MPI
        TDatabase::ParamDB->Comm = comm;
#endif
        
        TDomain domain(parmoon_db);
        
        
        // Produce an outfile "... .out", this is where all output is written to (addionally to console)
        //open OUTFILE, this is where all output is written to (addionally to console)
        if(my_rank==0)
        {
            Output::set_outfile(parmoon_db["outfile"]);
        }
        Output::setVerbosity(parmoon_db["verbosity"]);
        
        if(my_rank==0) //Only one process should do that.
        {
            // Write all parameters to the outfile (not to console) for later reference
            parmoon_db.write(Output::get_outfile());
            Database.WriteParamDB(argv[0]);
        }
        
        
        // // Do the parameter check of the Database.
        // check_parameters_consistency_NSE(parmoon_db);
        
        
        // Refine grid
   ////     size_t n_ref =  domain.get_n_initial_refinement_steps();
   ////     for(size_t i=0; i<= n_ref; i++)
   ////     {
   ////         if(n_ref>0) {
   ////             domain.RegRefineAll();
   ////             Output::print("Level:", i+1);
   ////         }
   ////         //}
   ////
   ////        // Write grid into a Postscript file
   ////         if(parmoon_db["output_write_ps"])
   ////             domain.PS("Domain.ps", It_Finest, 0);
            
            
            
            // Create an example object according to the database
            Example_Brinkman3D example(parmoon_db);
            
            // Intial refinement and grabbing of grids for multigrid.
            std::list<TCollection* > gridCollections
             = domain.refine_and_get_hierarchy_of_collections(parmoon_db);
         Output::print(my_rank,": I CONSTRUCTED GRIDCOLLECTIONS");
        
            //        std::list<TCollection* > gridCollections = domain.refine_and_get_hierarchy_of_collections(parmoon_db
             //                                                         );
            
            //=========================================================================
            // Create an object of the Brinkman class
//            Brinkman3D brinkman3d(domain, parmoon_db, example);
           
            
                    Brinkman3D brinkman3d(gridCollections, parmoon_db, example);
            
            Output::print(my_rank,": I CONSTRuCTED BRINKMAN3D OBJECT");
        
            timer.restart_and_print("constructing Brinkman3D object");
/////            Output::print<>("Database-Info:");
/////            parmoon_db.info();
            
            
            
            brinkman3d.assemble();
        
         Output::print(my_rank,": I ASSEMBLED");
        
            timer.restart_and_print("assembling linear terms");

        
            //Timer which measures time spent in solve() method solely.
            //Chrono timer_sol;
            //timer_sol.stop();
            //timer_sol.start();

            brinkman3d.solve();
          Output::print(my_rank,": I SOLVED");
          //  timer_sol.stop();
            
            // brinkman3d.solve_with_Petsc(parmoon_db);
            
            
            //        Output::print<>(TDatabase::ParamDB->VELOCITY_SPACE);
            //        Output::print<>(TDatabase::ParamDB->PRESSURE_SPACE);


   //// }
 //       timer_sol.print_total_time("solver only");
        
        brinkman3d.output();
        timer.restart_and_print("output");
        timer.print_total_time("Brinkman3D_ParMooN program");
        
        if(my_rank==0)
            Output::print("<<<<< ParMooN Finished: Brinkman3D Main Program >>>>>");
        
        if(my_rank == 0)
            Output::close_file();
}
            //=========================================================================
#ifdef _MPI
    MPI_Finalize();
#endif
    return 0;
} // end main
