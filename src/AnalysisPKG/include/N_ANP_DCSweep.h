//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2013  Sandia Corporation
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
// Filename       : $RCSfile: N_ANP_DCSweep.h,v $
//
// Purpose        : DC sweep analysis class
//
// Special Notes  : Specify any "hidden" or subtle details of the class here.
//                  Portability details, error handling information, etc.
//
// Creator        : Richard Schiek, SNL, Electrical and Microsystem Modeling
//
// Creation Date  : 01/24/08
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.9.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:30 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_DCSweep_h
#define Xyce_N_ANP_DCSweep_h

#include <Teuchos_RefCountPtr.hpp>
using Teuchos::RefCountPtr;
using Teuchos::rcp;

// ----------   Xyce Includes   ----------
#include <N_ANP_AnalysisBase.h>

// ---------- Forward Declarations ----------
class N_ANP_AnalysisManager;

//-------------------------------------------------------------------------
// Class         : N_ANP_DCSweep
// Purpose       : Transient analysis class
// Special Notes : 
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class N_ANP_DCSweep : public N_ANP_AnalysisBase 
{
  public:
    N_ANP_DCSweep( N_ANP_AnalysisManager * anaManagerPtr ) :
    N_ANP_AnalysisBase(anaManagerPtr),
    dcLoopInitialized_( false ),
    dcLoopSize_(0)
    {
      dcParamVec_ = rcp( new vector <N_ANP_SweepParam>() ); 
    };

    virtual ~N_ANP_DCSweep( ) {};

    bool setAnalysisParams(const N_UTL_OptionBlock & paramsBlock);
    bool outputFailureStats ();

    void setParamsWithOutputMgrAdapter 
    (RefCountPtr< N_ANP_OutputMgrAdapter > & outputMgrAdapterRCPtr) 
    {
      outputMgrAdapterRCPtr->setDCParamVec( dcParamVec_ );
    };

    bool run();
    bool init();
    bool loopProcess();
    bool processSuccessfulStep();
    bool processFailedStep();
    bool finish();
    bool handlePredictor();

    // Two Level specific
    bool twoLevelStep();

    void dcSweepOutput();
    void printStepHeader();
    bool printLoopInfo(int start, int finish);

  protected:

  private:
    void initializeSolution_();
    void takeStep_ ();

  public:
  protected:
  private:
    bool dcLoopInitialized_;
    int dcLoopSize_;

    list < int > dcSweepFailures_;
    RefCountPtr< vector <N_ANP_SweepParam> > dcParamVec_;
};

#endif
