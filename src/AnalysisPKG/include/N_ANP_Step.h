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
// Filename       : $RCSfile: N_ANP_Step.h,v $
//
// Purpose        : Step analysis class
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
// Revision Number: $Revision: 1.11 $
//
// Revision Date  : $Date: 2014/02/24 23:49:12 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_ANP_Step_h
#define Xyce_N_ANP_Step_h

#include <N_ANP_fwd.h>

#include <N_ANP_AnalysisBase.h>

namespace Xyce {
namespace Analysis {

//-------------------------------------------------------------------------
// Class         : Step
// Purpose       : Step analysis class
// Special Notes : 
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 01/24/08
//-------------------------------------------------------------------------
class Step : public AnalysisBase 
{
  public:
    Step( AnalysisManager * anaManagerPtr, AnalysisBase * anaType ) :
      AnalysisBase(anaManagerPtr),
      stepLoopSize_(0),
      stepLoopIter_(0)
    {
      // this is the analysis that will be done for each step
      mainAnalysisRCPtr_ = rcp( anaType, false);

      stepParamVec_ = rcp( new std::vector <SweepParam>() );
    };

    virtual ~Step( ) {};

    bool setAnalysisParams(const N_UTL_OptionBlock & paramsBlock);

    void setParamsWithOutputMgrAdapter 
      (RefCountPtr< OutputMgrAdapter > & outputMgrAdapterRCPtr) 
    {
      outputMgrAdapterRCPtr_->setStepParamVec( stepParamVec_ );
    };

    int getStepIter () { return stepLoopIter_; }

    virtual bool run(); /* override */
    virtual bool init(); /* override */
    virtual bool loopProcess(); /* override */
    virtual bool processSuccessfulStep(); /* override */
    virtual bool processFailedStep(); /* override */
    virtual bool finish(); /* override */
    virtual bool handlePredictor() /* override */ {
      return true;
    }
    
  private:
    RefCountPtr< AnalysisBase > mainAnalysisRCPtr_;
    
    RefCountPtr< std::vector <SweepParam> > stepParamVec_;
    //
    //bool stepLoopInitialized_;    // true if the step loop has been set up.
  
    int stepLoopSize_;
    int stepLoopIter_;
};

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::Step N_ANP_Step;

#endif // Xyce_N_ANP_Step_h
