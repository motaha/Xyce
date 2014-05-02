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
// Filename      : $RCSfile: N_IO_MeasureBase.C,v $
// Purpose       : Base class measure functions.
// Special Notes :
// Creator       : Richard Schiek, SNL, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.33.2.4 $
// Revision Date  : $Date: 2014/03/10 21:18:46 $
// Current Owner  : $Author: rlschie $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iterator>

#include <N_IO_MeasureBase.h>
#include <N_IO_Op.h>
#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace IO {
namespace Measure {

//-----------------------------------------------------------------------------
// Function      : MeasureBase::MeasureBase
// Purpose       : Constructor
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
Base::Base( const Util::OptionBlock & measureBlock,  IO::OutputMgr &output_manager )
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
 outputValueTarget_(0.0),
 outputValueTargetGiven_(false),
 lastOutputValue_(0.0),
 calculationDone_(false),
 calculationResult_(-1.0),
 independentVarColumn_(0),
 independentVar2Column_(0),
 dependentVarColumn_(0),
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

  std::list<N_UTL_Param>::iterator currentParamIt = const_cast<N_UTL_OptionBlock &>(measureBlock).begin();
  std::list<N_UTL_Param>::iterator endParamIt = const_cast<N_UTL_OptionBlock &>(measureBlock).end();
  while( currentParamIt != endParamIt )
  {
    
    std::string tag = currentParamIt->tag();
    //Xyce::dout() << " in meaasure base setup: tag = \"" << tag << "\"" << std::endl;
    if( tag == "NAME" )
    {
      name_ = currentParamIt->stringValue();
    }
    else if( tag == "MODE" )
    {
      mode_ = currentParamIt->stringValue();
    }
    else if( tag == "TYPE" )
    {
      type_ = currentParamIt->stringValue();
      // if type is TRIG or TARG then the next Util::Param in the list the
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
      td_ = currentParamIt->getImmutableValue<double>();
      tdGiven_ = true;
    }
    else if( ( tag == "GOAL" ) || ( tag == "VALUE" ) )
    {
      goal_ = currentParamIt->getImmutableValue<double>();
    }
    else if( tag == "WEIGHT" )
    {
      weight_ = currentParamIt->getImmutableValue<double>();
    }
    else if( tag == "MAX_THRESH" )
    {
      maxThresh_ = currentParamIt->getImmutableValue<double>();
      maxThreshGiven_ = true;
    }
    else if( tag == "MIN_THRESH" )
    {
      minThresh_ = currentParamIt->getImmutableValue<double>();
      minThreshGiven_ = true;
    }
    else if( tag == "MINVAL" )
    {
      minval_ = currentParamIt->getImmutableValue<double>();
    }
    else if( tag == "AT" )
    {
      at_ = currentParamIt->getImmutableValue<double>();
    }
    else if( tag == "FROM" )
    {
      from_ = currentParamIt->getImmutableValue<double>();
      fromGiven_ = true;
    }
    else if( tag == "TO" )
    {
      to_ = currentParamIt->getImmutableValue<double>();
      toGiven_ = true;
    }
    else if( (tag == "IGNORE") || (tag == "YMIN") )
    {
      ymin_ = currentParamIt->getImmutableValue<double>();
    }
    else if( tag == "YMAX" )
    {
      ymax_ = currentParamIt->getImmutableValue<double>();
    }
    else if( tag == "RISE" )
    {
      if( currentParamIt->getType() == Xyce::Util::STR )
      {
        // user requested LAST rise in simulation
        // so measure all of them and just report the last one.
        rise_ = -1;
      }
      else
      {
        rise_ = currentParamIt->getImmutableValue<int>();
      }
      riseGiven_ = true;
    }
    else if( tag == "FALL" )
    {
      if( currentParamIt->getType() == Xyce::Util::STR )
      {
        // user requested LAST fall in simulation
        // so measure all of them and just report the last one.
        fall_ = -1;
      }
      else
      {
        fall_ = currentParamIt->getImmutableValue<int>();
      }
      fallGiven_ = true;
    }
    else if( tag == "CROSS" )
    {
      if( currentParamIt->getType() == Xyce::Util::STR )
      {
        // user requested LAST cross in simulation
        // so measure all of them and just report the last one.
        cross_ = -1;
      }
      else
      {
        cross_ = currentParamIt->getImmutableValue<int>();
      }
      crossGiven_ = true;
    }
    else if( tag == "FRAC_MAX" )
    {
      if( inTrigBlock )
      {
        trigFracMax_ = currentParamIt->getImmutableValue<double>();
        trigFracMaxGiven_ = true;
      }
      else if( inTargBlock )
      {
        targFracMax_ = currentParamIt->getImmutableValue<double>();
        targFracMaxGiven_ = true;
      }
    }
    else if( Xyce::Util::hasExpressionTag(tag) )
    {
      numDepSolVars_++;
      depSolVarIterVector_.push_back(*currentParamIt);
    }
    else if( tag == "OBJVAL" )
    {
      if( currentParamIt->getType() == Xyce::Util::INT )
      {
        outputValueTarget_ = currentParamIt->getImmutableValue<int>();
        outputValueTargetGiven_ = true;
      }
      else if( currentParamIt->getType() == Xyce::Util::DBLE )
      {
        outputValueTarget_ = currentParamIt->getImmutableValue<double>();
        outputValueTargetGiven_ = true;
      }
      else if( currentParamIt->getType() == Xyce::Util::STR )
      {
        // a bare string name that we will have to resovle
        Util::Param aParam;
        aParam.set( currentParamIt->stringValue(), 0 );
        numDepSolVars_++;
        depSolVarIterVector_.push_back(aParam);
      }
      
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
      onValue_ = currentParamIt->getImmutableValue<double>();
      onValueGiven_ = true;
    }
    else if( tag == "OFF" )
    {
      offValue_ = currentParamIt->getImmutableValue<double>();
      offValueGiven_ = true;
    }
    else if( tag == "NUMFREQ" )
    {
      numFreq_ = currentParamIt->getImmutableValue<int>();
    }
    else if( tag == "GRIDSIZE" )
    {
      gridSize_ = currentParamIt->getImmutableValue<int>();
    }
    else if( tag == "FILE" )
    {
      dataFileName_ = currentParamIt->stringValue();
    }
    else if( tag == "COMP_FUNCTION" )
    {
      comparisonFunctionName_ = currentParamIt->stringValue();
    }
    else if( tag == "INDEPVARCOL" )
    {
      independentVarColumn_ = currentParamIt->getImmutableValue<int>();
    }
    else if( tag == "INDEPVAR2COL" )
    {
      independentVar2Column_ = currentParamIt->getImmutableValue<int>();
    }
    else if( tag == "DEPVARCOL" )
    {
      dependentVarColumn_ = currentParamIt->getImmutableValue<int>();
    }
    else if( tag == "DEFAULT_VAL" )
    {
      calculationResult_ = currentParamIt->getImmutableValue<double>();
    }
    else if( tag[0]=='V' || tag[0]=='I' || tag[0]=='N' ) 
    {
      // this if clause must come last because we are only checking the 
      // first letter and don't with to get confused with kewords 
      // that happen to start with V, I or N
      int nodes = currentParamIt->getImmutableValue<int>();
      Util::Param aParam;
      aParam.set( tag, nodes );

      if( inTrigBlock )
      {
        trig_ = aParam;
      }
      else if( inTargBlock )
      {
        targ_ = aParam;
      }

      // at this point trig and targ hold the beginning of the argument and how
      // many nodes it depends on (typically 1 or 2)  That's all we need there
      // to later use outputVars_ and depSolVarIterVector_ to get their values.

      // here we just store the needed parts of V(a) or v(a,b) or I(device).
      // only the v(a,b) case will need an extra node in the outputVars_ array.

      numDepSolVars_++;

      depSolVarIterVector_.push_back(aParam);
      for( int i=0; i<nodes; i++ )
      {
        currentParamIt++;
        aParam.set( currentParamIt->tag(), currentParamIt->getImmutableValue<double>() );
        depSolVarIterVector_.push_back( aParam );
      }
    }

    else
    {
      Xyce::Report::UserWarning() << "Unknown tag in measure statement: " << tag << ", ignoring";
    }
    currentParamIt++;
  }

}


void
Base::fixupMeasureParameters() 
{
  makeOps(outputManager_, depSolVarIterVector_.begin(), depSolVarIterVector_.end(), std::back_inserter(outputVars_));

  prepareOutputVariables();
}


//-----------------------------------------------------------------------------
// Function      : MeasureBase::withinTransientWindow
// Purpose       : Checks if current time is within measurement window
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
 bool Base::withinTransientWindow( double time )
 {
   bool retVal = true;
   if( tdGiven_ && (time < td_) )
   {
     retVal = false;
   }
   return retVal;
 }

 //-----------------------------------------------------------------------------
// Function      : MeasureBase::withinFromToWindow
// Purpose       : Checks if current time is within measurement window
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
 bool Base::withinFromToWindow( double time )
 {
   bool retVal = true;
   if( (fromGiven_ && (time < from_ )) || (toGiven_ && (time > to_)) )
   {
     retVal = false;
   }
   return retVal;
 }

//-----------------------------------------------------------------------------
// Function      : MeasureBase::withinRiseFallCrossWindow
// Purpose       : Checks if current value is within measurement window
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool Base::withinRiseFallCrossWindow( double measureVal, double crossVal )
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
    if( (((measureVal-crossVal) < 0.0) && ((lastOutputValue_-crossVal) > 0.0)) 
     || (((measureVal-crossVal) > 0.0) && ((lastOutputValue_-crossVal) < 0.0)) )
    {
      // we've crossed measureVal-crossVal == 0 
      actualCross_++;
      //Xyce::dout() << name_ << " actualCross = " << actualCross_ << " with lastOutputValue_ = " << lastOutputValue_ << std::endl;
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
    lastOutputValue_=measureVal;
  }
  return retVal;
}


//-----------------------------------------------------------------------------
// Function      : MeasureBase::withinMinMaxThreash
// Purpose       : Check if value is withing MIN_THRESHOLD & MAX_THRESHOLD
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
bool Base::withinMinMaxThreash( double value)
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
// Purpose       : Call's the OutputMgr's getValue() function to update 
//                 the objects in std::list<N_UTL_Param> outputVars_;
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 11/01/2013
//-----------------------------------------------------------------------------
void Base::updateOutputVars(std::vector<double> & outputVarVec, const double circuitTime, const N_LAS_Vector *solnVec, 
  const N_LAS_Vector *stateVec, const N_LAS_Vector * storeVec, const N_LAS_Vector *imaginaryVec )

{
  int vecIndex = 0;
  for (std::vector<Util::Operator *>::const_iterator it = outputVars_.begin(); it != outputVars_.end(); ++it)
  {
    outputVarVec[vecIndex] = getValue(outputManager_.getCommPtr()->comm(), *(*it), solnVec, imaginaryVec, stateVec, storeVec ).real();
    vecIndex++;
  }
}


//-----------------------------------------------------------------------------
// Function      : MeasureBase::getOutputValue
// Purpose       : Call's the OutputMgr's getPrgetImmutableValue<int>() function to get sol. vars.
// Special Notes :
// Scope         : public
// Creator       : Richard Schiek, Electrical and Microsystem Modeling
// Creation Date : 03/10/2009
//-----------------------------------------------------------------------------
double Base::getOutputValue(Xyce::Util::Operator *op, const N_LAS_Vector *solnVec, const N_LAS_Vector *stateVec, 
  const N_LAS_Vector * storeVec, const N_LAS_Vector *imaginaryVec )
{
  double retVal = getValue(outputManager_.getCommPtr()->comm(), *op, solnVec, imaginaryVec, stateVec, storeVec ).real();
  return retVal;
}

} // namespace IO
} // namespace IO
} // namespace Xyce
