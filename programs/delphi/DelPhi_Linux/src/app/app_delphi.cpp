/**
 * @file app_delphi.cpp
 * @brief Main function to generate the executable delphicpp
 *
 * @author Chuan Li, chuanli@clemson.edu
 *
 * delphicpp is the object-oriented C++ version of the program DelPhi, which takes as input a Brookhaven
 * database coordinate file format of a molecule or equivalent data for geometrical objects and/or charge
 * distributions and calculates the electrostatic potential in and around the system, using a finite
 * difference solution to the Poisson-Boltzmann equation. This can be done for any given concentration of
 * any two different salts. The system can be composed of different parts having different dielectrics.
 *
 * =================================================================================================
 *
 * [QUICK REFERENCE TO ACCESS SVN REPOSITORY ON SERVER consus.clemson.edu] \n
 * - TO SHOW LOG HISTORY FOR THIS PROJECT:
 *   	svn log svn+ssh://username\@consus.clemson.edu/home/common/delphi_distr/delphicpp
 * - TO CHECKOUT THE LATEST REVERSION: 
 *   	svn co svn+ssh://username\@consus.clemson.edu/home/common/delphi_distr/delphicpp
 * - TO REMOVE/ADD FILES FROM A CHECKED-OUT PROJECT:
 *   	svn rm  <file/directory name>
 *   	svn add <file/directory name>
 * - TO REVIEW LOCAL CHANGES:
 *   	svn status
 * - TO COMMIT CHANGES:
 *   	export SVN_EDITOR=/usr/bin/vim
 *   	svn commit
 * - TO UPDATE YOUR LOCAL FILES:
 *   	svn update
 *
 *   (See "subversion user guide" created by S. Sarkar and C. Li for other svn options)
 *
 * =================================================================================================
 * 
 * [SUBMITTED REVERSIONS] \n
 * - r01   chuan   2013July01   The 1st reversion in SVN.
 * - r02   chuan   2013Nov25    1. Updated reading size and charge files in IO class.
 *                              2. Updated excel file. \n
 * - r03   chuan   2013Dec04    1. Completed reading/writing PDB file in IO class.
 *                              2. Completed assigning size and charge to atoms. \n
 * - r04   chuan   2013Dec05    fixed bugs in IO class.
 * - r05   chuan   2013Dec10    major changes in prototypes, delphi, main and misc.
 * - r06   chuan   2013Dec13    overhauled the subsystems of prototypes, delphi, misc, io.
 *                              1. Now the part of code before surface construction is completed.
 *                              2. Now the code requires GCC 4.3+ with option std=c++0x to compile due to the additional feature of
 *                              initializing values when creating large arrays. (see comments in makefile)
 * - r07   chuan   2013Dec24    Another set of major changes includes:
 *                              1. Converted map<string,void *> to map<string,boost::any>
 *                              2. Converted all regular and dynamically allocated multi-dimensional arrays to 1D stl::vector<T> in the map. \n
 *                              The reasons to do so are: 1) contiguous memory allocation. 2) faster 1-layer loop instead of
 *                              multi-layer loop. 3) Easier allocate/deallocate arrays
 *                              3. Started using smart pointers instead of raw pointers in  order to avoid possible memory leak.
 *                              4. Removed datadistr(:) and ndistr from the map
 *                              5. Now the class of DataContainer provides 4 template functions to access the map:
 *                            	   1) read-only by value:         getKey_constRef
 *                            	   2) read-only by pointer:       getKey_constPtr
 *                            	   3) writable access by value:   getKey_Val
 *                            	   4) writable access by pointer: getKey_Ptr
 *                              6. Reorganized the code structure.
 * - r08   chuan   2014Jan07    Re-config the code to have datacontainer and datamarshal loosely paired.
 * - n/a   chuan   2014Jan10    1. Added new piece of code to read ibnum, ibgrd, iepsmp and idebmap from forSolver.dat without Lin's
 *                              surface construction class.
 *                              2. Fixed various bugs.
 * - r09   chuan   2014Jan18    1st reversion to use Eclipse and Doxygen.
 * - r10   chuan   2014Jan18    cleaned two additional files
 * - r11   chuan   2014Feb17    completed solver class with other updates.
 * - r12   chuan   2014Feb26    Added new class of writing site and phimap.
 * - r13   chuan   2014Mar26    Added new directory of examples with 3 super-large proteins.
 * - r14   LinWang 2014Apr03    Energy Class added
 * - r15   LinLi   2014Apr03    space class and examples added
 * - r16   LinLi   2014Apr03    app_delphi.cpp for space added
 * - r17   lwang3  2014Apr07    Energy class modified: 1.Renamed variables. 2.Warning message added for some testing functions.
 * - r18   lwang3  2014Apr08    Example files cleanup.
 * - r19   chuan   2014Apr08    1. new function reset is added in datacontainer class for testing purpose.
 *                              2. modified delphi95.r16 is added to svn.
 * - r20   lwang3  2014Apr08    Fix bugs in delphi_datamarshal_getFunction.cpp.
 *                              Now it can control energy output by parameter file instead of setting flags.
 * - r21   chuan   2014Apr08    updated app_delphi.cpp to include all individual classes. F95 code have been updated accordingly.
 * - r22   LinLi   2014Apr10    Space module has been tested using showmap, now the aftersurf.dat are identical in fortran and C++
 * - r23   Chuan   2014Apr14    Modified app_delphi.cpp to show map after realizing every object.
 * - r24   LinLi   2014Apr15    changed size of ibgrd from ibmx to ibnum
 * - r25   Chuan   2014Apr21    1. Added new features for MCCE calling DelPhi
 *                              2. Added calculated energies into data container.
 *                              3. Various bug fixes.
 *                              This version is the 1st version for integrated tests.
 * - r26   chuan   2014Apr21    Fixed the bug of deblen.
 * - r27   chuan   2014Apr23    Another update to clean debug info in the code.
 * - r28   chuan   2014Apr23    Cleaned debug info in F95 code.
 * - r29   chuan   2014May02    Fixed most format and form bugs in statements reported in DelphiCPP_testscript_2A26_r28_20140424.xlsx
 * - r30   chuan   2014May10    1. Recoded reading functions using tokens.
 *                              2. Fixed most bugs in functions reported in DelphiCPP_testscript_2A26_r28_20140424.xlsx
 * - r31   chuan   2014May12    1. A few more bug fixes.
 *                              2. 2 new proteins, 1AB1.pdb and 1SF4.pdb are added.
 *                              3. Debug/ and Release/ are added.
 * - r34   chuan   2014May16    1. Removed "\r" in statements.
 *                              2. Introduced namespace "delphi" to fix ambiguous reference to "delphi_real" when compiling the program on Mac.
 * - r35   LinWang 2014May18    1. Fixed the bug of printing multiple identical messages on screen when atom radius less than ZERO.
 *                              2. Added energy_exceptions.h file and start using CWarning class for showing warning messages.
 * - r36   chuan   2014May19    1. variables named prg* were changed to vct*.
 *                              2. bug fixed for reading charge and size from pqr file
 *                              3. flag MCCE and DEBUG_MCCE are added for calling delphicpp in mcce.
 *                              4. Doxygen comments are added to files in \app and \delphi.
 * - r37   LinWang 2014May20    Added more warning messages to energy_exceptions.h
 * - r38   chuan   2014May21    A few bug fixes in CSite and CDelphiFastSOR::isFocusBndy.
 * - r39   LinWang 2014May24    The test environment setup files are added to SVN.
 * - r42   chuan   2014Jun02    1. defined the destructor of IAbstractModule to be virtual.
 *                              (Atten: pSpace.release() needs to be replaced by pSpace.reset() after fixing the problem of memory leak.)
 * - r43   chuan   2014Jun02    Added new debug flag DEBUG_OBJECT and message in constructors and destructors of all classes for
 *                              debugging objects realization.
 * - r44   chuan   2014Jun03    Merged the projects of delphicpp and mcce_delphicpp into one.
 *                              (Atten: the problem of memory leak makes mcce_delphi too expensive to run on PC)
 * - r45   chuan   2014Jun03    The following statements
 *                               - CLCSRF
 *                               - GRDCON
 *                               - LOGGRP
 *                               - MEMBRANEDATA
 *                               - PHICON
 *                               - RADPOL
 *                               - NORMC \n
 *                              and functions
 *                               - out(srf,file="filename") or out(srf,file="filename",format=BEM)
 *                               - out(frc,file="filename",format=r) or out(frc,file="filename",format=un)
 *                               - out(unpdb,file="filename") or out(unpdb,file="filename",format="whatever")
 *                               - out(unfrc,file="filename") or out(unfrc,file="filename",format="whatever")
 *                               - buffz(xxxyyyzzzxxxyyyzzz)
 *                              have been removed as discussed after the 1st round of tests.
 * - r50   chuan   2014Jun17    1. Removed shortcuts to read PQR and PQR4 files.
 *                              2. Documented \interface
 * - r51   chuan   2014Jun18    Fixed out of range and memory leaks reported by valgrind in solver and site classes.
 * - r52   chuan   2014Jun29    Took off showMap in app_delphi.
 * - r53   chuan   2014Jul02    Renamed data type real --> delphi_real and integer --> delphi_integer throughout the
 *                              code.
 * - r54   chuan   2014Jul25    A version with the flag PARALLEL_OMP to control sequential and omp-parallelized solver class.
 * - r55   LinWang 2014Aug04    Updated energy module implemented OpenMP and GCC Auto-Vectorization Optimization (-ffast-math).
 * - r56   chuan   2014Aug08    Documented \io and \misc in Doxygen.
 * - r57   LinWang 2014Aug15    1. Update CTimer function and added omp_wall_time() to adjust Timer in Delphi OpenMP program. 
 *                              2. Remove OpenMP setup and timer in energy module. 
 * - r58   LinLi   2014Aug18    All of the openMP codes are merged together.
 * - r59   LinLi   2014Aug27    ibmx bug fixed.
 * - r60   LinWang 2014Aug29    This version is for initial DelPhi C++ release beta 1. Changes are :
 *                              1. Modified /delphi/delphi_data_flash.cpp for new version.
 *                              2. Added Developer flags to print elapsed time for each class on the screen.
 *                              3. Removed -fopenmp option in makefile of the classes which do not implement OpenMP such as 
 *                                 I/O class to obtain better performance. 
 *                              4. Added a comment to makefile for compiling delphicpp for MacOSX to link OpenMP library to binary.
 */

#include "../delphi/delphi_data.h"
#include "../space/space.h"
#include "../solver/solver_fastSOR.h"
#include "../energy/energy.h"
#include "../site/site.h"

int main(int argc, char *argv[])
{
   /*
    * bool values are inserted/extracted by their textual representation: either true or false, instead of integral values.
    * This flag can be unset with the noboolalpha manipulator.
    */ 
   cout << boolalpha; 

   cerr << boolalpha;

#ifdef DEVELOPER
   cout << fixed << setprecision(7); //cout.precision(15)
#else
   cout << fixed << setprecision(3); //cout.precision(7)
#endif       

   try
   {
      CTimer * pTester =  new CTimer; // record execution time

      pTester->start();

      shared_ptr<CTimer> pTimer( new CTimer); // record execution time

      //********************************************************************************//
      //                                                                                //
      //                  read the parameter, size, charge, PDB files                   //
      //                                                                                //
      //********************************************************************************//

      //---------- a shared_ptr to an object of IDataContainer
      shared_ptr<IDataContainer> pDataContainer( new CDelphiData(argc,argv,pTimer) );

#ifdef DEBUG_DELPHI_SPACE
      cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SPACE] data container is written into file test_delphicpp_atbeginning.dat---\n\n"
           << "\033[0m";

      pDataContainer->showMap("test_delphicpp_atbeginning.dat");
#endif

#ifdef DEVELOPER
      cout << "\n\n---------- delphicpp finishes IO in "; pTester->showElapse(); cout << "\n\n";
#endif

      //********************************************************************************//
      //                                                                                //
      //    realize an object of CDelphiSpace class to construct molecular surfaces     //
      //                                                                                //
      //********************************************************************************//

      unique_ptr<IAbstractModule> pSpace( new CDelphiSpace(pDataContainer,pTimer) );

      pSpace->run();

      pSpace.reset();

#ifdef DEBUG_DELPHI_SPACE
      cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SPACE] data container is written into file test_delphicpp_aftersurf.dat---\n\n"
           << "\033[0m";

      pDataContainer->showMap("test_delphicpp_aftersurf.dat");
#endif

#ifdef DEVELOPER
      cout << "\n\n---------- delphicpp finishes SPACE class in "; pTester->showElapse(); cout << "\n\n";

      cout << endl;
#endif
    
      cout << " number of atom coordinates read  :" << right << setw(10) << pDataContainer->getKey_constRef<delphi_integer>("natom") << endl;

      if (pDataContainer->getKey_constRef<bool>("isolv"))
      {
         cout << " total number of assigned charges :" << right << setw(10) << pDataContainer->getKey_constRef<delphi_integer>("nqass") << endl;
         cout << " net assigned charge              :" << right << setw(10) << pDataContainer->getKey_constRef<delphi_real>("qnet") << endl;
         cout << " assigned positive charge         :" << right << setw(10) << pDataContainer->getKey_constRef<delphi_real>("qplus") << endl;
         cout << " centred at (gu)                  :" << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqplus").nX << " "
                                                       << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqplus").nY << " "
                                                       << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqplus").nZ << endl;
         cout << " assigned negative charge         :" << right << setw(10) << pDataContainer->getKey_constRef<delphi_real>("qmin") << endl;
         cout << " centred at (gu)                  :" << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqmin").nX << " "
                                                       << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqmin").nY << " "
                                                       << right << setw(10) << pDataContainer->getKey_constRef< SGrid<delphi_real> >("cqmin").nZ << endl;
         cout << "\nnumber of dielectric boundary points" << right << setw(10) << pDataContainer->getKey_constRef<delphi_integer>("ibnum") << endl;

         if (pDataContainer->getKey_constRef<bool>("iexun") && 0 == pDataContainer->getKey_constRef<delphi_integer>("ibnum"))
            throw CNoBndyAndDielec(pTimer); 

         //********************************************************************************//
         //                                                                                //
         //   realize an object of CDelphiFastSOR class to calculate potentials on grids   //
         //                                                                                //
         //********************************************************************************//

#ifdef DEBUG_DELPHI_SOLVER
         cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SOLVER] data container is read from file test_fromsurf.dat---\n\n" << "\033[0m";

         pDataContainer->reset("test_fromsurf.dat");

         cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SOLVER] data container is written into file test_delphicpp_beforeitr.dat---\n\n"
              << "\033[0m";

         pDataContainer->showMap("test_delphicpp_beforeitr.dat");

#endif // DEBUG_DELPHI_SOLVER

         unique_ptr<IAbstractModule> pSolver( new CDelphiFastSOR(pDataContainer,pTimer) );
         
         pSolver->run();
      
         pSolver.reset();

#ifdef DEBUG_DELPHI_SOLVER
         cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SOLVER] data container is written into file test_delphicpp_afteritr.dat---\n\n"
                    << "\033[0m";

         pDataContainer->showMap("test_delphicpp_afteritr.dat");
#endif // DEBUG_DELPHI_SOLVER

#ifdef DEVELOPER
         cout << "\n\n---------- delphicpp finishes SOLVER class in "; pTester->showElapse(); cout << "\n\n";
#endif

         //********************************************************************************//
         //                                                                                //
         //          realize an object of CDelphiEnergy class to calculate energies        //
         //                                                                                //
         //********************************************************************************//

#ifdef DEBUG_DELPHI_ENERGY
         cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_ENERGY] data container is read from file test_fromsurf.dat "
              << "and test_fromsolver.dat---\n\n" << "\033[0m";

         pDataContainer->reset("test_fromsurf.dat");

         pDataContainer->reset("test_fromsolver.dat");

         cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SOLVER] data container is written into file test_delphicpp_beforeenergy.dat---\n\n"
              << "\033[0m";

         pDataContainer->showMap("test_delphicpp_beforeenergy.dat");

#endif // DEBUG_DELPHI_ENERGY

         unique_ptr<IAbstractModule> pEnergy( new CDelphiEnergy(pDataContainer,pTimer) );

         pEnergy->run();

         pEnergy.reset();

#ifdef DEBUG_DELPHI_ENERGY
         cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_ENERGY] data container is written into file test_delphicpp_aftereng.dat---\n\n"
              << "\033[0m";

         pDataContainer->showMap("test_delphicpp_aftereng.dat");
#endif // DEBUG_DELPHI_ENERGY

#ifdef DEVELOPER
         cout << "\n\n---------- delphicpp finishes ENERGY class in "; pTester->showElapse(); cout << "\n\n";
#endif

         //********************************************************************************//
         //                                                                                //
         //               realize an object of CSite class to write site info              //
         //                                                                                //
         //********************************************************************************//

#ifdef DEBUG_DELPHI_SITE
         cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SITE] data container is read from file test_fromsurf.dat, "
              << "test_fromsolver.dat and test_fromenergy.dat---\n\n" << "\033[0m";

         pDataContainer->reset("test_fromsurf.dat");
         pDataContainer->reset("test_fromsolver.dat");
         pDataContainer->reset("test_fromenergy.dat");

         cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SOLVER] data container is written into file test_delphicpp_beforesite.dat---\n\n"
              << "\033[0m";

         pDataContainer->showMap("test_delphicpp_beforesite.dat");
#endif // DEBUG_DELPHI_SITE

         unique_ptr<CSite> pSite( new CSite(pDataContainer,pTimer) );

         if (pDataContainer->getKey_constRef<bool>("isite"))
         {
            int iisitsf = 0;
            if (pDataContainer->getKey_Ref<bool>("isitsf")) iisitsf = 1;
            pSite->writeSite(iisitsf);
         }

         if (pDataContainer->getKey_constRef<bool>("phiwrt")) pSite->writePhi();

         pSite.reset();

#ifdef DEBUG_DELPHI_SITE
         cerr << "\n\033[1;35m" << "---[DEBUG_DELPHI_SITE] data container is written into file test_delphicpp_atend.dat---\n\n"
              << "\033[0m";

         pDataContainer->showMap("test_delphicpp_atend.dat");
#endif // DEBUG_DELPHI_SITE

#ifdef DEVELOPER
         cout << "\n\n---------- delphicpp finishes SITE class in "; pTester->showElapse(); cout << "\n\n";
#endif

      }     
      
      pDataContainer.reset();
      //cout << "[+++debug+++] number of shared_ptr\'s to IDataContainer  = " << pDataContainer.use_count() << endl;

      pTimer->exit(); pTimer.reset();

      delete pTester;

   } // ---------- end of try block
   catch (CException&)
   {
      cerr << "\n\n ......... PROGRAM ABORTS WITH AN EXCEPTION AND " << CWarning::iWarningNum << " WARNING(S) ........\n\n";

      return 1;
   }

   cout << "\n\n .........  PROGRAM EXITS SUCCESSFULLY : WITH TOTAL " << CWarning::iWarningNum << " WARNING(S) ........\n\n";
 
   cout.unsetf(ios_base::floatfield); // return to cout default notation 

   return 0;

}




