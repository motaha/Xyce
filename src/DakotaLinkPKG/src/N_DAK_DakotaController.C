//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2014 Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Filename       : $N_DAK_DakotaController.C$
//
// Purpose        : This class defines a controller that Xyce can use to
//                  start a Dakota based analysis.  The object here is to
//                  allow Xyce to read a netlist with some Dakota commands
//                  embeded within it and then pass off control to Dakota
//                  to organize and run the Xyce simulations it needs.
//
// Special Notes  :
//
// Creator        : Richard L. Schiek
//
// Creation Date  : 09/07/2005
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.35 $
//
// Revision Date  : $Date: 2014/02/24 23:49:13 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#define DISABLE_DAKOTA_CONFIG_H 1

#include <N_DAK_DakotaController.h>

#define DISABLE_DAKOTA_CONFIG_H   1 // Don't need dakota's dakota_config.h 
#include <DakotaInterface.hpp>
#include <DakotaModel.hpp>
#include <DakotaResponse.hpp>
#include <DakotaVariables.hpp>
#undef  DISABLE_DAKOTA_CONFIG_H

#include <N_CIR_Xyce.h>
#include <N_ERH_ErrorMgr.h>
#include <N_DAK_DakotaInterface.h>


#include <iostream>
#include <string>

// global variable needed by Dakota for output.
int write_precision = 10;

//
// Constructor
//
N_DAK_DakotaController::N_DAK_DakotaController(int iargsIn, char *cargsIn[]):
  outputFileName( "outputDakota.txt" ),
  stdOutputFilename( "outputDakota.txt" ),
  stdErrorFilename( "errorDakota.txt" ),
  readRestartFilename(""),
  writeRestartFilename("restartOutDakota.txt"),
  restartEvals( 0 ),
  numFunctions( 0 )
{
#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "In N_DAK_DakotaController::N_DAK_DakotaController()" << std::endl;
#endif
  // save the iargsIn and cargsIn passed in to the constructor.
  // however, we want to remove the -dakota <filename> part so that 
  // we can pass this to another N_CIR_Xyce object and not have that 
  // object think we need to run Dakota.  In the copying process,
  // save the dakota input file name.
  
  iargsReduced = iargsIn - 2;  // we'll take all the args except "-dakota" and <filename>
  if( iargsReduced < 1 )
  {
    // remaingin arguments aren't enough to get at least a netlist.
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, "Not enough arguments passed in for a dakota run (i.e. runxyce <netlist> -dakota <dakota input file>)\n");
  }
  
  cargsReduced = new char * [iargsReduced];
  if( cargsReduced == NULL )
  {
    // remaingin arguments aren't enough to get at least a netlist.
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, "Couldn't allocate memory for arguments in N_DAK_DakotaController::N_DAK_DakotaController().");
  }
  
  bool dakotaFileFound = false;
  for (int i=0,k=0; i<iargsIn;++i)
  {
    // this seems odd, but it was how N_IO_CmdParse did this as well.
    if (cargsIn[i] == NULL) 
    {
      cargsReduced[k] = NULL;
      continue;
    }

    std::string tmpString(cargsIn[i]);
    if ( tmpString == "-dakota" )
    {
      // skip the copy of this 
      i++;
      // save the next item
      std::string filename(cargsIn[i]);
      dakotaInputFileName = filename;
      dakotaFileFound = true;
    }
    else
    {
      //everything else gets copied
      int size = tmpString.size()+2;
      cargsReduced[k] = new char[size];
      if( cargsReduced[k] == NULL )
      {
        // remaingin arguments aren't enough to get at least a netlist.
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, "Couldn't allocate memory for a copy of the arguments in N_DAK_DakotaController::N_DAK_DakotaController().");
      }
      for (int j=0;j<size;++j)
      {
        cargsReduced[k][j] = 0;
      }
      sprintf(cargsReduced[k], tmpString.c_str());
      k++;
    }
  }
  
  if(!dakotaFileFound)
  {
    // didn't find a dakota input file name so throw an error
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::USR_FATAL, "Could not find the dakota input filename on the command line. (i.e. -dakota <filename>)\n");
  }
  
  initializeDakota();
}

//
// destructor
//
N_DAK_DakotaController::~N_DAK_DakotaController()
{
#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "In N_DAK_DakotaController::~N_DAK_DakotaController()" << std::endl;
#endif
  // need to delete the cargsReduced we created
  for(int i=0; i<iargsReduced; i++)
  {
    if( cargsReduced[i] )
    {
      delete [] cargsReduced[i];
    }
  }
  if( cargsReduced )
  {
    delete cargsReduced;
  }

}


bool N_DAK_DakotaController::run()
{
  bool result = true;;
  result = result && executeDakotaStrategy();
  retrieveDakotaResults();
  return true;
}


//
// get the Dakota problem ready to run
//
bool N_DAK_DakotaController::initializeDakota()
{
  bool result = true;
  // to initialize Dakota, we need to call    
  result = result && initializeDakotaParallelLibrary( );
  result = result && initializeDakotaProblemDatabase( );
  result = result && initializeDakotaStrategy();
  result = result && constructApplicationInterface();

  return result;
}


//
// set up Dakota components
//
bool N_DAK_DakotaController::initializeDakotaParallelLibrary()
{
  bool result = true;
#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "In N_DAK_DakotaController::initializeDakotaParallelLibrary()" << std::endl;
#endif

  parallelLib = rcp( new Dakota::ParallelLibrary( iargsReduced, cargsReduced ) );
  
  if( Teuchos::is_null( parallelLib ) )
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, "Could not allocate Dakota::ParallelLibrary object.\n");
  }

  return result;
}


bool N_DAK_DakotaController::initializeDakotaProblemDatabase()
{
  bool result = true;

#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "In N_DAK_DakotaController::initializeDakotaProblemDatabase()" << std::endl;
#endif
  problemDatabase = rcp(new Dakota::ProblemDescDB( *parallelLib ));
  if( Teuchos::is_null(problemDatabase) )
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, "Could not allocate Dakota::ProblemDescDB object.\n");
  }
  else
  {
    problemDatabase->manage_inputs( dakotaInputFileName.c_str() );
  }
    
  return result;
}

bool N_DAK_DakotaController::initializeDakotaStrategy()
{
  bool result = true;

#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "In N_DAK_DakotaController::initializeDakotaStrategy()" << std::endl;
#endif
  problemStrategy = rcp(new Dakota::Strategy( *problemDatabase ) );
  
  if( Teuchos::is_null(problemStrategy) )
  {
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, "Could not allocate Dakota::Strategy object.\n");
  }

#ifdef Xyce_Dakota_Parallel_Debug
  numCommunicators = static_cast< int >( parallelLib->analysis_intra_communicators().size() );
  Xyce::dout() << "In N_DAK_DakotaController::initializeDakotaStrategy() numCommunicators = " << numCommunicators << std::endl;
  for( int i=0; i<numCommunicators; i++ )
  {
    Xyce::dout() << "MPI_Comm is " << (parallelLib->analysis_intra_communicators())[i] << std::endl;
  }
 #endif   
  return result;
}


bool N_DAK_DakotaController::constructApplicationInterface()
{
  bool result = true;

#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "In N_DAK_DakotaController::constructApplicationInterface()" << std::endl;
#endif

  Dakota::ModelList& models = problemDatabase->model_list();
  Dakota::ModelLIter mlIter;
  for (mlIter = models.begin(); mlIter != models.end(); mlIter++) 
  {
    Dakota::Interface& interface = (*mlIter).interface();
    // to do:  need to handle case where there is more than one element.
    //string analysisDrivers( interface.analysis_drivers().at(0) );
    if( ( interface.interface_type() == "direct" ) )
    {
      RCP<N_DAK_DakotaInterface> aDakotaInterface = rcp(new N_DAK_DakotaInterface( *problemDatabase ), false);
      
      if( Teuchos::is_null(aDakotaInterface) )
      {
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL_0, "Could not allocate N_DAK_DakotaInterface object.\n");
      }
      else
      {
        aDakotaInterface->setArguments(iargsReduced,cargsReduced);
        interface.assign_rep( aDakotaInterface.get(), false );
        //(*mlIter).reset_communicators();
        dakotaInterfacesPtrs_.push_back( aDakotaInterface );
      }
    }
  }
  
#ifdef Xyce_Dakota_Parallel_Debug
  int numCommunicators = static_cast< int >( parallelLib->analysis_intra_communicators().size() );
  Xyce::dout() << "In N_DAK_DakotaController::constructApplicationInterface() numCommunicators = " << numCommunicators << std::endl;
  for( int i=0; i<numCommunicators; i++ )
  {
    Xyce::dout() << "MPI_Comm is " << (parallelLib->analysis_intra_communicators())[i] << std::endl;
  }
#endif
  
  return result;
}

bool N_DAK_DakotaController::executeDakotaStrategy()
{
#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "In N_DAK_DakotaController::executeDakotaStrategy()" << std::endl;
#endif    
  problemStrategy->run_strategy();
  return true;
}

void N_DAK_DakotaController::retrieveDakotaResults()
{
#ifdef Xyce_Dakota_Debug
  Xyce::dout() << "In N_DAK_DakotaController::retrieveDakotaResults()" << std::endl;
#endif

  const Dakota::Variables & vars = problemStrategy->variables_results();
  const Dakota::Response & resp = problemStrategy->response_results();

}

#undef DISABLE_DAKOTA_CONFIG_H 
