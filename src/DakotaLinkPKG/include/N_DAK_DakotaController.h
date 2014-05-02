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
// Filename       : $N_DAK_DakotaController.h$
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
// Revision Number: $Revision: 1.27 $
//
// Revision Date  : $Date: 2014/02/24 23:49:13 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef N_DAK_DakotaController_h
#define N_DAK_DakotaController_h 1

// standard includes
#include <string>
#include <vector>
#include <list>
#include <iosfwd>

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// Dakota headers
#define DISABLE_DAKOTA_CONFIG_H   1 // Don't need dakota's dakota_config.h 
#include <ParallelLibrary.hpp>
#include <ProblemDescDB.hpp>
#include <DakotaStrategy.hpp>
#undef  DISABLE_DAKOTA_CONFIG_H  

// ----------   Xyce Includes   ----------
#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>

// foward declarations
#include <N_CIR_fwd.h>

class N_DAK_DakotaInterface;

class N_DAK_DakotaController
{
public:

    N_DAK_DakotaController(int iargsIn, char *cargsIn[]);
    ~N_DAK_DakotaController();

    // get the Dakota problem ready to run
    bool initializeDakota();

    // The following methods define the sequence
    // of events that occur in starting up, configuring
    // running and gathering results from a Dakota
    // analysis.

    bool initializeDakotaParallelLibrary();
    bool initializeDakotaProblemDatabase();
    bool initializeDakotaStrategy();
    bool constructApplicationInterface();
    bool executeDakotaStrategy();
    void retrieveDakotaResults();
    
    // start the Dakota controlled simulations
    bool run();

private:
    int iargsReduced;             // argumet line without -dakota <filename>
    char ** cargsReduced;         //
    
    std::string dakotaInputFileName;   // dakota input file 
    std::string outputFileName;
    std::string stdOutputFilename;
    std::string stdErrorFilename;
    std::string readRestartFilename;
    std::string writeRestartFilename;
    int restartEvals;
    int numFunctions;
    
    // dakota objects we'll need
    RCP< Dakota::ParallelLibrary > parallelLib;
    RCP< Dakota::ProblemDescDB > problemDatabase;
    RCP< Dakota::Strategy > problemStrategy;
    
    // a vector of all the DakotaInterface objects we 
    // make so that we can clean them up later
    std::vector< RCP< N_DAK_DakotaInterface > > dakotaInterfacesPtrs_;
};

#endif

