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
// Filename      : $RCSfile: N_IO_MeasureBase.C,v $
// Purpose       : Base class measure functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.15.2.3 $
// Revision Date  : $Date: 2013/12/03 23:30:12 $
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

// ---------- Standard Includes ----------

// ----------   Xyce Includes   ----------
#include<N_IO_MeasureBase.h>
#include <N_ERH_ErrorMgr.h>

//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureBase::N_IO_MeasureBase
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
N_IO_MeasureBase::N_IO_MeasureBase( const N_UTL_OptionBlock & measureBlock,  N_IO_OutputMgr &output_manager )
  : name_(""),
 mode_(""),
 type_(""),
 typeSupported_(false),
 td_(0.0),
 tdGiven_(false),
 goal_(0.0),
 weight_(0.0),
 minval_(1.0e-12),
 at_(0.0),
 from_(0.0),
 fromGiven_(false),
 to_(0.0),
 toGiven_(false),
 ymin_(0.0),
 ymax_(0.0),
 rise_(0),
 riseGiven_(false),
 fall_(0),
 fallGiven_(false),
 cross_(0),
 crossGiven_(false),
 actualRise_(0),
 isRising_(false),
 maxThresh_(0.0),
 maxThreshGiven_(false),
 minThresh_(0.0),
 minThreshGiven_(false),
 actualFall_(0),
 isFalling_(false),
 actualCross_(0),
 onValue_(0.0),
 onValueGiven_(false),
 offValue_(0.0),
 offValueGiven_(false),
 fractionToExtrema_(0.0),
 fractionToExtremaGiven_(false),
 numFreq_(10),
 gridSize_(200),
 trigOutputValueTarget_(0.0),
 trigOutputValueTargetGiven_(false),
 targOutputValueTarget_(0.0),
 targOutputValueTargetGiven_(false),
 trigFracMax_(0.0),
 trigFracMaxGiven_(false),
 targFracMax_(0.0),
 targFracMaxGiven_(false),
 numDepSolVars_(0),
 outputValueTarget_(0.0),
 outputValueTargetGiven_(false),
 lastOutputValue_(0.0),
 calculationDone_(false),
 calculationResult_(0.0),
  outputManager_(output_manager)
{
  // since many of the measure types share the use of keywords (like TD=<delay time>) we'll
  // parse those out here and store them.  We'll also pull out the out_var[= val | out_var2]
  // here since once we eliminate the keyword from the measureBlock, that's all that should left.

  // these are used to mark if we're in the trigger or target sections of a rise/fall/delay
  // measure.  This measure is different from the others in that it is essentiality
  // a compound measure which requires two objectives and thus two sets of qualifiers
  // on the objective like frac_max=, weight= etc.
  bool inTrigBlock = false;
  bool inTargBlock = false;

  list<N_UTL_Param>::iterator currentParamIt = const_cast<N_UTL_OptionBlock &>(measureBlock).begin();
  list<N_UTL_Param>::iterator endParamIt = const_cast<N_UTL_OptionBlock &>(measureBlock).end();
  while( currentParamIt != endParamIt )
  {
    string tag = currentParamIt->tag();
    if( tag == "NAME" )
    {
      name_ = currentParamIt->sVal();
    }
    else if( tag == "MODE" )
    {
      mode_ = currentParamIt->sVal();
    }
    else if( tag == "TYPE" )
    {
      type_ = currentParamIt->sVal();
      // if type is TRIG or TARG then the next N_UTL_Param in the list the
      // node or expression that is the Trigger or Target.  We'll need to
      // catch and save this.  This oddity arises because all of the other
      // measures have one objective to read while the TRIG/TARG of the Rise
      // fall, delay measure has two and they are thus named.  Also, some of
      // the qualifiers like value=, weight= and frac_max= apply to the
      // last TRIG ar TARG statement so we need to remember if we just
      // passed one of those.
      if( type_ == "TRIG" )
      {
        inTrigBlock = true;
        inTargBlock = false;
      }
      else if( type_ == "TARG" )
      {
        inTrigBlock = false;
        inTargBlock = true;
      }
    }
    else if( tag == "TD" )
    {
      td_ = currentParamIt->dVal();
      tdGiven_ = true;
    }
    else if( tag == "GOAL" )
    {
      goal_ = currentParamIt->dVal();
    }
    else if( tag == "WEIGHT" )
    {
      weight_ = currentParamIt->dVal();
    }
    else if( tag == "MAX_THRESH" )
    {
      maxThresh_ = currentParamIt->dVal();
      maxThreshGiven_ = true;
    }
    else if( tag == "MIN_THRESH" )
    {
      minThresh_ = currentParamIt->dVal();
      minThreshGiven_ = true;
    }
    else if( tag == "MINVAL" )
    {
      minval_ = currentParamIt->dVal();
    }
    else if( tag == "AT" )
    {
      at_ = currentParamIt->dVal();
    }
    else if( tag == "FROM" )
    {
      from_ = currentParamIt->dVal();
      fromGiven_ = true;
    }
    else if( tag == "TO" )
    {
      to_ = currentParamIt->dVal();
      toGiven_ = true;
    }
    else if( (tag == "IGNORE") || (tag == "YMIN") )
    {
      ymin_ = currentParamIt->dVal();
    }
    else if( tag == "YMAX" )
    {
      ymax_ = currentParamIt->dVal();
    }
    else if( tag == "RISE" )
    {
      if( currentParamIt->getType() == STR )
      {
        // user requested LAST rise in simulation
        // so measure all of them and just report the last one.
        rise_ = -1;
      }
      else
      {
        rise_ = currentParamIt->iVal();
      }
      riseGiven_ = true;
    }
    else if( tag == "FALL" )
    {
      if( currentParamIt->getType() == STR )
      {
        // user requested LAST fall in simulation
        // so measure all of them and just report the last one.
        fall_ = -1;
      }
      else
      {
        fall_ = currentParamIt->iVal();
      }
      fallGiven_ = true;
    }
    else if( tag == "CROSS" )
    {
      if( currentParamIt->getType() == STR )
      {
        // user requested LAST cross in simulation
        // so measure all of them and just report the last one.
        cross_ = -1;
      }
      else
      {
        cross_ = currentParamIt->iVal();
      }
      crossGiven_ = true;
    }
    else if( tag == "FRAC_MAX" )
    {
      if( inTrigBlock )
      {
        trigFracMax_ = currentParamIt->dVal();
        trigFracMaxGiven_ = true;
      }
      else if( inTargBlock )
      {
        targFracMax_ = currentParamIt->dVal();
        targFracMaxGiven_ = true;
      }
    }
    else if( (tag == "V") || (tag == "I") )
    {

      int nodes = currentParamIt->iVal();
      N_UTL_Param aParam;
      aParam.set( tag, nodes );

      if( inTrigBlock )
      {
        //std::cout << "In TRIG and aParam is " << aParam << std::endl;
        trig_ = aParam;
      }
      else if( inTargBlock )
      {
        //std::cout << "In TARG and aParam is " << aParam << std::endl;
        targ_ = aParam;
      }

      // at this point trig and targ hold the beginning of the argument and how
      // many nodes it depends on (typically 1 or 2)  That's all we need there
      // to later use outputVars_ and depSolVarIterVector_ to get their values.

      // here we just store the needed parts of V(a) or v(a,b) or I(device).
      // only the v(a,b) case will need an extra node in the outputVars_ array.

      numDepSolVars_++;
      // use list::insert() to keep an iterator pointing to this spot
      list<N_UTL_Param>::iterator newDepSolVarItr = outputVars_.insert( outputVars_.end(), aParam);
      depSolVarIterVector_.push_back( newDepSolVarItr );
      for( int i=0; i<nodes; i++ )
      {
        currentParamIt++;
        aParam.set( currentParamIt->tag(), currentParamIt->dVal() );
        outputVars_.push_back( aParam );
      }

    }
    else if( tag == "OBJVAL" )
    {
      if( currentParamIt->getType() == INT )
      {
        outputValueTarget_ = currentParamIt->iVal();
      }
      else
      {
        outputValueTarget_ = currentParamIt->dVal();
      }
      outputValueTargetGiven_ = true;
      if( inTrigBlock )
      {
        trigOutputValueTarget_ = outputValueTarget_;
        trigOutputValueTargetGiven_ = true;
        outputValueTargetGiven_ = false;
      }
      else if( inTargBlock )
      {
        targOutputValueTarget_ = outputValueTarget_;
        targOutputValueTargetGiven_ = true;
        outputValueTargetGiven_ = false;
      }
    }
    else if( tag == "ON" )
    {
      onValue_ = currentParamIt->dVal();
      onValueGiven_ = true;
    }
    else if( tag == "OFF" )
    {
      offValue_ = currentParamIt->dVal();
      offValueGiven_ = true;
    }
    else if( tag == "NUMFREQ" )
    {
      numFreq_ = currentParamIt->iVal();
    }
    else if( tag == "GRIDSIZE" )
    {
      gridSize_ = currentParamIt->iVal();
    }
    else
    {
      // unknown tag.  Issue a warning.
      string msg = "Unknown tag in measure statement: " + tag + " Will try to ignore.";
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::USR_WARNING, msg);
    }
    currentParamIt++;
  }

}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureBase::withinTransientWindow
// Purpose       : Checks if current time is within measurement window
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
 bool N_IO_MeasureBase::withinTransientWindow( double time )
 {
   bool retVal = true;
   if( time < td_ )
   {
     retVal = false;
   }
   return retVal;
 }

 //-----------------------------------------------------------------------------
// Function      : N_IO_MeasureBase::withinFromToWindow
// Purpose       : Checks if current time is within measurement window
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
 bool N_IO_MeasureBase::withinFromToWindow( double time )
 {
   bool retVal = true;
   if( (fromGiven_ && (time < from_ )) || (toGiven_ && (time > to_)) )
   {
     retVal = false;
   }
   return retVal;
 }

//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureBase::withinRiseFallCrossWindow
// Purpose       : Checks if current value is within measurement window
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool N_IO_MeasureBase::withinRiseFallCrossWindow( double measureVal )
{
  bool retVal = true;
  if( riseGiven_ || fallGiven_ || crossGiven_ )
  {
    retVal = false;
    // first check if we need to adjust rise/fall/cross counts.
    if( (measureVal > lastOutputValue_) && !isRising_ )
    {
      // we've started a rise
      isRising_= true;
      isFalling_ = false;
      actualRise_++;
    }
    if( (measureVal < lastOutputValue_) && !isFalling_ )
    {
      // we've started a fall
      isRising_ = false;
      isFalling_ = true;
      actualFall_++;
    }
    if( ((measureVal < 0.0) && (lastOutputValue_ > 0.0)) || ((measureVal > 0.0) && (lastOutputValue_ < 0.0)) )
    {
      // we've crossed zero
      actualCross_++;
    }

    // now check if we're in the right window
    // this could be compressed to one statement, but this looks clearer
    if( riseGiven_ && ((rise_ < 0) || (rise_ == actualRise_)))
    {
      retVal=true;
    }
    else if( fallGiven_ && ((fall_ < 0) || (fall_ == actualFall_)))
    {
      retVal=true;
    }
    else if( crossGiven_ && ((cross_< 0) || (cross_ == actualCross_)))
    {
      retVal=true;
    }
  }
  return retVal;
}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureBase::withinMinMaxThreash
// Purpose       : Check if value is withing MIN_THRESHOLD & MAX_THRESHOLD
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool N_IO_MeasureBase::withinMinMaxThreash( double value)
{
  bool returnValue = true;
  if( (minThreshGiven_ && (value < minThresh_)) )
    returnValue = false;
  if( (maxThreshGiven_ && (value > maxThresh_)) )
    returnValue = false;

  return returnValue;
}



//-----------------------------------------------------------------------------
// Function      : MeasureBase::updateOutputVars
// Purpose       : Call's the OutputMgr's getPrintValue() function to update 
//                 the objects in std::list<N_UTL_Param> outputVars_;
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 11/01/2013
//-----------------------------------------------------------------------------
void N_IO_MeasureBase::updateOutputVars(std::vector<double> & outputVarVec, const double circuitTime, RCP< N_LAS_Vector > solnVecRCP)
{

  std::list<N_UTL_Param>::iterator currentParamIt = outputVars_.begin();
  int vecIndex = 0;
  while( currentParamIt != outputVars_.end() )
  {
    if( currentParamIt->getSimContext() == UNDEFINED )
    {
      // call set param context
      outputManager_.setParamContextType_( currentParamIt );
    }
    
    if( currentParamIt->getSimContext() == NODE_OR_DEVICE_NAME )
    {
      // delete any unneeded Param objects from the outputVars_ list.
      currentParamIt = outputVars_.erase(currentParamIt); 
    }
    else
    {
      outputVarVec[ vecIndex ]  = outputManager_.getPrintValue( currentParamIt, solnVecRCP.getRawPtr() );
      vecIndex++;
      currentParamIt++;
    }
    
  }
}


//-----------------------------------------------------------------------------
// Function      : N_IO_MeasureBase::getOutputValue
// Purpose       : Call's the OutputMgr's getPrintValue() function to get sol. vars.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
double N_IO_MeasureBase::getOutputValue(list<N_UTL_Param>::iterator & paramListIt, RCP< N_LAS_Vector > solnVecRCP)
{
#ifdef Xyce_DEBUG_IO
  std::cout << "In N_IO_MeasureBase::getOutputValue() about to call getPrintValue..." << std::endl;
#endif

  if( paramListIt->getSimContext() == UNDEFINED )
  {
    // call set param context
    outputManager_.setParamContextType_( paramListIt );
  }
  double retVal = outputManager_.getPrintValue( paramListIt, solnVecRCP.getRawPtr() );
  return retVal;
}
