// ----------------------------------------------------------------------------
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
// ----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Resistor.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 02/28/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.229.2.3 $
//
// Revision Date  : $Date: 2014/03/06 23:33:43 $
//
// Current Owner  : $Author: tvrusso $
//----------------------------------------------------------------------------
#include <Xyce_config.h>

#include <N_DEV_Resistor.h>

#include <N_DEV_DeviceOptions.h>
#include <N_DEV_ExternData.h>
#include <N_DEV_Message.h>

#include <N_LAS_Matrix.h>

namespace Xyce {
namespace Device {
namespace Resistor {


///
/// Common Jacobian Stamp for all Resistor devices.
/// Because all resistors have identical Jacobian stamps, this data is
/// declared static and is shared by all resistor instances.
/// 
std::vector<std::vector<int> > Instance::jacStamp;

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::initializeJacobianStamp
// Purpose       : 
// Special Notes : Initialization of jacobian stamp moved from constructor
//                 in revision 1.227 of N_DEV_Resistor.C by David Baur.
//                 The code itself code was written by R. Hoekstra, 9/27/2002
// Scope         : private
// Creator       : David Baur
// Creation Date : 2/11/2014
//-----------------------------------------------------------------------------
///
/// @brief Common Jacobian stamp initializer for all Resistor devices.
///
/// The Jacobian stamp is a sparse-matrix representation of the pattern
/// of non-zero elements that a device will put into the Jacobian matrix.
///
/// The Jacobian stamp is used by the Topology package to determine indices
/// into the full Jacobian matrix for elements that correspond to this 
/// device.
///
/// There is one row of the Jacobian stamp for each equation associated with
/// a device.  The number of elements in a row are the number of non-zero
/// elements in that row of the device's contribution to the Jacobian.
/// The values of the elements are numbers local to the device that
/// represent the column in which the non-zero element appears.
///
/// For the resistor, there are two external nodes (the positive and negative
/// terminals of the device).  The positive node is referred to as the 0th
/// node of the device, and the negative node the 1st node.
/// Considering positive current flow from the positive node to the negative
/// node, current out of the positive node is \f$(V_+-V_-)*G\f$, and current
/// out of the negative node is \f$-(V_+-V_-)*G\f$.  Thus, the Jacobian
/// matrix contribution for the resistor is:
/// \f[
/// \left[\begin{array}{cc}
///  G& -G\\
/// -G& G
/// \end{array}
/// \right] \f]
///
/// This is a dense Jacobian with two rows.  The first row is the row
/// for the positive node's KCL, the second row is the row for the
/// negative node KCL.  Each row has two non-zero elements.  The
/// columns correspond to the nodes of the device: the first column is
/// the positive node, the second the negative node.  The element of
/// the jacobian is the dependence of the equation associated with the
/// row on the variable associated with the column.
///
/// The Jacobian stamp therefore has two rows, and each row has two elements.
/// In this trivial device (because the matrix is small and fully dense),
/// the stamp values are 0 ("positive node") for the first nonzero in 
/// each row, and 1 ("negative node") for the second nonzero.
///
void Instance::initializeJacobianStamp()
{
  if (jacStamp.empty())
  {
    jacStamp.resize(2);
    jacStamp[0].resize(2);
    jacStamp[1].resize(2);
    jacStamp[0][0] = 0;
    jacStamp[0][1] = 1;
    jacStamp[1][0] = 0;
    jacStamp[1][1] = 1;
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Traits::loadInstanceParameters
// Purpose       : 
// Special Notes : The addPar calls here were refactored and moved here
//                 from the instance constructor.  Those addPars had been
//                 in place from 2/4/2005.
// Scope         : private
// Creator       : David Baur
// Creation Date : 1/27/2014
//-----------------------------------------------------------------------------
///
/// Loads the parameter definition into the instance parameter map.
///
/// @param p     instance parameter map
///
/// Each parameter supported by the resistor device instance is
/// defined via the addPar call.  The minimum information required is
/// the name of the parameter, the default value, and a member pointer
/// to a variable in the instance class to which the value will be
/// assigned.
///
/// Additional information such as documentation for the parameter, units
/// (used only in output of tables for documentation), whether a special
/// boolean should be set if the parameter is given, etc. may be specified
/// using the various set* calls of the Xyce::Device::Descriptor class.
///
/// It is important to note that since parameters defined by addPar are
/// defined as metadata separate from any instance, it is not possible to
/// establish interrelationships between parameter defaults here.  Parameter
/// defaults specified in addPar are always constants.  If a device requires
/// that parameter defaults depend on values of other parameters, this has to
/// be done in the instance constructor.  Examples of such parameters in this 
/// device arethe "DTEMP" and "W" parameters, which are set to special defaults
/// at device instantiation time.  Defaults for those parameters in the addPar
/// calls (and hence in the LaTeX tables in the reference guide produced from
/// this data) are misleading.
///
void Traits::loadInstanceParameters(ParametricData<Resistor::Instance> &p)
{
  p.addPar("R", 1000.0, &Resistor::Instance::R)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_OHM)
    .setDescription("Resistance");
  p.addPar("L", 0.0, &Resistor::Instance::length)
    .setUnit(U_METER)
    .setDescription("Length");
  p.addPar("W", 0.0, &Resistor::Instance::width)
    .setUnit(U_METER)
    .setDescription("Width");
  p.addPar("TEMP", 0.0, &Resistor::Instance::temp)
    .setExpressionAccess(ParameterType::TIME_DEP)
    .setUnit(U_DEGC)
    .setDescription("Device temperature");

  p.addPar("TC1", 0.0, &Resistor::Instance::tempCoeff1)
    .setGivenMember(&Resistor::Instance::tempCoeff1Given)
    .setUnit(U_DEGCM1)
    .setDescription("Linear Temperature Coefficient");
  p.addPar("TC2", 0.0, &Resistor::Instance::tempCoeff2)
    .setGivenMember(&Resistor::Instance::tempCoeff2Given)
    .setUnit(U_DEGCM2)
    .setDescription("Quadratic Temperature Coefficient");
  p.makeVector("TC", 2); // Allow TC to be entered as a vector (TC=1,2)

  p.addPar("DTEMP", 0.0, &Resistor::Instance::dtemp)
    .setGivenMember(&Resistor::Instance::dtempGiven)
    .setUnit(U_DEGC)
    .setDescription("Device Temperature -- For compatibility only. Parameter is NOT used");
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Traits::loadModelParameters
// Purpose       : 
// Special Notes : The addPar calls here were refactored and moved here
//                 from the model constructor.  Those addPars had been
//                 in place from 2/4/2005.
// Scope         : private
// Creator       : David Baur
// Creation Date : 1/27/2014
//-----------------------------------------------------------------------------
///
/// Loads the parameter definition into the model parameter map.
///
/// @param p     model parameter map
///
/// @see Xyce::Device::Resistor::Traits::loadInstanceParameters
///
/// Resistors may optionally be given a model name to enable a semiconductor
/// resistor model with a sheet resistance.  The resistance of an instance
/// is then determined by the length and width given on the instance line.
/// This loadModelParameters method defines the parameters that may be 
/// specified on model cards associated with semiconductor resistors.
///
void Traits::loadModelParameters(ParametricData<Resistor::Model> &p)
{
  // Create parameter definitions for parameter member variables
  p.addPar("TC1", 0.0, &Resistor::Model::tempCoeff1)
    .setUnit(U_DEGCM1)
    .setDescription("Linear Temperature Coefficient");
  p.addPar("TC2", 0.0, &Resistor::Model::tempCoeff2)
    .setUnit(U_DEGCM2)
    .setDescription("Quadratic Temperature Coefficient");
  p.addPar("RSH",   0.0, &Resistor::Model::sheetRes)
    .setUnit(U_OHM)
    .setDescription("Sheet Resistance");
  p.addPar("DEFW",  1.e-5, &Resistor::Model::defWidth)
    .setUnit(U_METER)
    .setDescription("Default Instance Width");
  p.addPar("NARROW",0.0, &Resistor::Model::narrow)
    .setUnit(U_METER)
    .setDescription("Narrowing due to side etching");
  p.addPar("TNOM",  0.0, &Resistor::Model::tnom)
    .setUnit(U_DEGC)
    .setDescription("Parameter Measurement Temperature");
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::Instance
// Purpose       : Instance constructor
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter
// Creation Date : 3/16/2000
//-----------------------------------------------------------------------------
///
/// Construct a resistor instance.
///
/// @param[in] configuration  Device configuration and traits.
/// @param[in] instance_block Instance information from parser.
/// @param[in] model Resistor Model to which we should add this instance.
/// @param[in] factory_block  Device options defined in netlist.
///
/// @note The parameter member variable values from the initializers
/// are immediately replaced with the values from the parameter
/// definitions by the setDefaultParams() function.

Instance::Instance(
  const Configuration & configuration,
  const InstanceBlock & instance_block,
  Model &               model,
  const FactoryBlock &  factory_block)
  : DeviceInstance(instance_block, configuration.getInstanceParameters(), factory_block),
    model_(model),
    R(0.0),
    length(0.0),
    width(0.0),
    temp(factory_block.deviceOptions_.temp.getImmutableValue<double>()),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    dtemp(0.0),
    tempCoeff1Given(false),
    tempCoeff2Given(false),
    dtempGiven(false),
    G(0.0),
    i0(0.0),
    li_Pos(-1),
    li_Neg(-1),
    li_store_dev_i(0),
    APosEquPosNodeOffset(-1),
    APosEquNegNodeOffset(-1),
    ANegEquPosNodeOffset(-1),
    ANegEquNegNodeOffset(-1)
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  ,
    f_PosEquPosNodePtr(0),
    f_PosEquNegNodePtr(0),
    f_NegEquPosNodePtr(0),
    f_NegEquNegNodePtr(0)
#endif
{
  // Initialize DeviceInstance values
  numIntVars   = 0;    // Initialize number if internal nodes in DeviceInstance
  numExtVars   = 2;    // Initialize number if external nodes in DeviceInstance
  numStateVars = 0;    // Initialize number if state variables in DeviceInstance
  setNumStoreVars(0);  // Initialize number if store variables in DeviceInstance
  numLeadCurrentStoreVars = 1; // Initialize number if lead current variables 
                               //in DeviceInstance

  initializeJacobianStamp();

  // Set params to constant default values from parameter definition
  setDefaultParams();

  // Set params according to instance line and constant defaults from metadata
  setParams(instance_block.params);

  // Set any non-constant parameter defaults:
  if (!given("TEMP"))
    temp = factory_block.deviceOptions_.temp.getImmutableValue<double>();
  if (!given("W"))
    width = model_.defWidth;

  // Get temperature values from model is not given in instance
  if (!tempCoeff1Given)
    tempCoeff1 = model_.tempCoeff1;
  if (!tempCoeff2Given)
    tempCoeff2 = model_.tempCoeff2;

  // Calculate any parameters specified as expressions
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors
  if (!given("R"))
  {
    if (model_.given("RSH") && given("L") && (model_.sheetRes != 0) &&
        (length != 0))
    {
      R = model_.sheetRes * (length - model_.narrow)
          / (width - model_.narrow);
    }
    else
    {
      R = 1000.0;
      UserWarning0(*this) << "Resistance is set to 0, setting to the default, " << R << " ohms";
    }
  }

  // Process the parameters to complete initialization
  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::processParams
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/3/02
//-----------------------------------------------------------------------------
/// Process parameters.
///
/// @return true on success
///
/// In general, the processParams method is intended as a place for complex
/// manipulation of parameters that must happen if temperature or other
/// external changes happen.  In the resistor device it does nothing other
/// than call updateTemperature.
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   6/03/02
bool Instance::processParams()
{
  // now set the temperature related stuff.
  return updateTemperature(temp);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::registerLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
///
/// Register local IDs
///
/// Register the local internal and external node IDs.
///
/// @param intLIDVecRef internal local IDs from topology package
/// @param extLIDVecRef external local IDs from topology package
/// 
/// Instantiation (calling the device constructor) of the device
/// sets up variables numIntVars and numExtVars, the numbers of internal and
/// external variables associated with the device.  This information is then
/// used by the Topology package to assign locations in the solution vector
/// (and all other vectors of the same shape) for those variables.
/// The "Local IDs" (LIDs) of these locations are provided by Topology
/// so the device can know where to load its data.
///
/// This method saves the LIDs from Topology and associates each one with
/// a particular local name for the internal or external variable.  They 
/// are then used when we load the F and Q vectors.
///
/// The resistor device has no internal variables, so this method makes no use 
/// of the intLIDVecRef array.
///
/// @author Robert Hoekstra, SNL, Parallel Computational Sciences
/// @date   6/12/02
void Instance::registerLIDs(
  const std::vector<int> & intLIDVecRef,
  const std::vector<int> & extLIDVecRef)
{
  AssertLIDs(intLIDVecRef.size() == numIntVars);
  AssertLIDs(extLIDVecRef.size() == numExtVars);

  // Copy the local ID lists.
  intLIDVec = intLIDVecRef;                           // Set the internal local IDs in DeviceInstance
  extLIDVec = extLIDVecRef;                           // Set the external local IDs in DeviceInstance

  li_Pos = extLIDVec[0];
  li_Neg = extLIDVec[1];

  if (DEBUG_DEVICE && getDeviceOptions().debugLevel > 0)
  {
    dout() << getName() << " LIDs"
      //<< Util::push << std::endl
      << std::endl
           << "li_Pos_ = " << li_Pos << std::endl
           << "li_Neg_ = " << li_Neg << std::endl
           //<< Util::pop << std::endl;
           << std::endl;
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::registerStateLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/12/02
//-----------------------------------------------------------------------------
///
/// Register the local state IDs
///
/// @note The resistor does not have any state vars, so this function
/// does nothing.
///
/// @param staLIDVecRef State variable local IDs
///
/// In general, devices may declare at construction that they require storage
/// locations in the "state vector."  Topology assigns locations in the 
/// state vector and returns "Local IDs" (LIDs) for devices to use for their
/// state vector entries.  If a device uses state vector entries, it
/// uses the registerStateLIDs method to save the local IDs for later use.
/// 
/// @author Robert Hoekstra, SNL, Parallel Computational Sciences
/// @date   06/12/02

void Instance::registerStateLIDs(const std::vector<int> & staLIDVecRef)
{
  AssertLIDs(staLIDVecRef.size() == numStateVars);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::registerStoreLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
/// Register the local store IDs
///
/// In addition to state vector, Xyce maintains a separate datastructure
/// called a "store" vector.  As with other such vectors, the device
/// declares at construction time how many store vector entries it needs,
/// and later Topology assigns locations for devices, returning LIDs.
///
/// These LIDs are stored in this method for later use.
///
/// The Resistor device uses exactly one "store vector" element, where
/// it keeps the "lead current" that may be used on .PRINT lines as
/// "I(R1)" for the current through resistor R1.
///
/// @param stoLIDVecRef Store variable local IDs
///
/// @author Richard Schiek, Electrical Systems Modeling
/// @date   12/18/2012

void Instance::registerStoreLIDs(const std::vector<int> & stoLIDVecRef)
{
  AssertLIDs(stoLIDVecRef.size() == getNumStoreVars());

  if (loadLeadCurrent)
  {
    li_store_dev_i = stoLIDVecRef[0];
  }
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::getStoreNameMap
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Richard Schiek, Electrical Systems Modeling
// Creation Date : 12/18/2012
//-----------------------------------------------------------------------------
/// Populates and returns the store name map.
///
/// If the DeviceInstance::storeNameMap is empty, populate it.
///
/// @return reference to the DeviceInstance::storeNameMap
///
/// For the purpose of lead currents, store vector elements must be given 
/// names that can be used to locate lead currents at print time.
/// When a netlist attempts to print, say, "I(R1)" the output code looks for
/// an entry in the store vector named R1:DEV_I.
/// 
/// This method does the assignment of names to store vector elements.
///
/// @author Richard Schiek, Electrical Systems Modeling
/// @date   12/18/2012

std::map<int,std::string> & Instance::getStoreNameMap()
{
  // set up the internal name map, if it hasn't been already.
  if (loadLeadCurrent && storeNameMap.empty ())
  {
    // change subcircuitname:devicetype_deviceName to
    // devicetype:subcircuitName:deviceName
    std::string modName(getName());
    spiceInternalName(modName);
    std::string tmpstr;
    tmpstr = modName+":DEV_I";
    storeNameMap[ li_store_dev_i ] = tmpstr;
  }

  return storeNameMap;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::registerJacLIDs
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 08/27/01
//-----------------------------------------------------------------------------
///
/// Register the Jacobian local IDs
///
/// @param jacLIDVec Jacobian local Ids
///
/// @see Xyce::Device::Resistor::Instance::initializeJacobianStamp
///
/// Having established local IDs for the solution variables, Topology must
/// also assign local IDs for the elements of the Jacobian matrix.
///
/// For each non-zero element that was identified in the jacobianStamp,
/// Topology will assign a Jacobian LID.  The jacLIDVec will have the 
/// same structure as the JacStamp, but the values in the jacLIDVec will
/// be offsets into the row of the sparse Jacobian matrix corresponding
/// to the non-zero identified in the stamp.
/// 
/// These offsets are stored in this method for use later when we load
/// the Jacobian.
///
/// @author Robert Hoekstra, SNL, Parallel Computational Sciences
/// @date   08/27/01

void Instance::registerJacLIDs(const std::vector< std::vector<int> > & jacLIDVec)
{
  // Let DeviceInstance do its work.
  DeviceInstance::registerJacLIDs(jacLIDVec);

  // Store offsets of the components of the Jacobian of this instance
  APosEquPosNodeOffset = jacLIDVec[0][0];
  APosEquNegNodeOffset = jacLIDVec[0][1];
  ANegEquPosNodeOffset = jacLIDVec[1][0];
  ANegEquNegNodeOffset = jacLIDVec[1][1];
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::setupPointers
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/30/08
//-----------------------------------------------------------------------------
///
/// Setup direct access pointer to solution matrix and vectors.
///
/// @see Xyce::Device::Resistor::Instance::registerJacLIDs
///
/// As an alternative to the row offsets defined in registerJacLIDs, it 
/// is also possible to obtain direct pointers of the Jacobian elements.
///
/// This method uses the offsets obtained in registerJacLIDs to retrieve
/// the pointers.
///
/// In the resistor device the pointers to the matrix are only saved
/// (and are only used for matrix access) if
/// Xyce_NONPOINTER_MATRIX_LOAD is NOT defined at compile time.  For
/// some devices the use of pointers instead of array indexing can be
/// a performance enhancement.
///
/// Use of pointers in this device is disabled by defining
/// Xyce_NONPOINTER_MATRIX_LOAD at compile time.  When disabled, array
/// indexing with the offsets from registerJacLIDs is used in
/// the matrix load methods.
///
/// @author Eric Keiter, SNL
/// @date   11/30/08

void Instance::setupPointers()
{
#ifndef Xyce_NONPOINTER_MATRIX_LOAD
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  f_PosEquPosNodePtr = &(dFdx[li_Pos][APosEquPosNodeOffset]);
  f_PosEquNegNodePtr = &(dFdx[li_Pos][APosEquNegNodeOffset]);
  f_NegEquPosNodePtr = &(dFdx[li_Neg][ANegEquPosNodeOffset]);
  f_NegEquNegNodePtr = &(dFdx[li_Neg][ANegEquNegNodeOffset]);
#endif
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::updateIntermediateVars
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 3/05/2004
//-----------------------------------------------------------------------------
///
/// Update the intermediate variables
///
/// @return true on success
///
/// The bulk of any device's computation is handled in the instance class'
/// updateIntermediateVars method.  For the resistor, this is
/// merely the computation of the current through the resistor given the
/// voltage difference between its terminals.
///
/// Intermediate variables computed here are used in the methods that
/// load data into the F, Q, dFdX and dQdX data structures.
///
/// @note This method is called by the updatePrimaryState
/// method. Since the Resistor class reimplements the "Master"
/// "loadState" function that loads the contributions from all
/// resistor devices in a single loop, THIS FUNCTION IS NOT ACTUALLY
/// USED!  
///
/// @see Xyce::Device::Resistor::Instance::loadDAEFVector Xyce::Device::Resistor::Instance::loadDAEQVector Xyce::Device::Resistor::Instance::loadDAEdFdx Xyce::Device::Resistor::Instance::loadDAEdQdx Xyce::Device::Resistor::Instance::updatePrimaryState
/// @author Eric R. Keiter, Dept. 9233.
/// @date   3/05/04
///
bool Instance::updateIntermediateVars()
{
  double * solVec = extData.nextSolVectorRawPtr;
  double v_pos = solVec[li_Pos];
  double v_neg = solVec[li_Neg];
  i0 = (v_pos-v_neg)*G;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::updatePrimaryState
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/29/01
//-----------------------------------------------------------------------------
///
/// Update the state variables.
///
/// @return true on success
///
/// This function is the function that is called from the device manager
/// each time a load of vectors and Jacobian is required.  Its first
/// job is to call updateIntermediateVars.
///
/// After calling updateIntermediateVars, the updatePrimaryState method
/// may load state vector elements as needed.
/// 
/// The resistor device has no state vector elements, so all this method does
/// in the resistor instance class is call updateIntermediateVars.
///
/// There is no longer a "secondary" state.  The "primary" in
/// this method's name is purely historical.
///
/// @note This method is called by the default implementation of the
/// loadState master function. Since the Resistor class reimplements
/// the "Master" "loadState" function that loads the contributions
/// from all resistor devices in a single loop, THIS FUNCTION IS NOT
/// ACTUALLY USED, NOR is the updateIntermediateVars method it calls!
/// The updatePrimaryState method is only called when a device class
/// does not re-implement the master class.  This can be a source of
/// confusion when attempting to modify the Resistor device, or any
/// other device That reimplements the Master classes.
///
/// @see Xyce::Device::Resistor::Master::updateState
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   01/29/01
///
bool Instance::updatePrimaryState()
{
  return updateIntermediateVars();
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::loadDAEFVector
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 1/24/03
//-----------------------------------------------------------------------------
///
/// Load the DAE F vector.
///
/// @return true on success
///
/// The Xyce DAE formulation solves the differential-algebraic
/// equations \f$F(x)+dQ(x)/dt=0\f$ These are vector equations
/// resulting from the modified nodal analysis of the circuit.
/// 
/// This method loads the F-vector contributions for a single resistor
/// instance.
///
/// In this method, the offsets defined in registerLIDs are used to
/// store the device's F contributions into the F vector.
///
/// The Q vector is used for devices that store charge or magnetic
/// energy, and since the resistor does none of that, the F vector
/// contributions are the whole of the resistor's contribution to the
/// full set of DAEs.
///
/// @note This method is called by the default implementation of the
/// loadDAEVectors master function. Since the Resistor class
/// reimplements the "Master" "loadDAEVectors" function that loads the
/// contributions from all resistor devices in a single loop, THIS
/// FUNCTION IS NOT ACTUALLY USED.  The loadDAEFVector method is only
/// called when a device class does not re-implement the master class.
/// This can be a source of confusion when attempting to modify the Resistor
/// device, or any other device that reimplements the Master classes.
///
/// @see Xyce::Device::Resistor::Master::loadDAEVectors
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   01/24/03
///
bool Instance::loadDAEFVector()
{
  double * fVec = extData.daeFVectorRawPtr;
  fVec[li_Pos] += i0;
  fVec[li_Neg] += -i0;

  if( loadLeadCurrent )
  {
    double * stoVec = extData.nextStoVectorRawPtr;
    stoVec[li_store_dev_i] = i0;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::loadDAEdFdx
// Purpose       : Loads the F-vector contributions for a single
//                 resistor  instance.
// Special Notes : 
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 03/05/04
//-----------------------------------------------------------------------------
///
/// Load the DAE the derivative of the F vector with respect to the
/// solution vector x, dFdx
///
/// Loads the contributions for a single resistor instance to the 
/// dFdx matrix (the F contribution to the Jacobian).
///
/// This method uses the Jacobian LIDs (row offsets) that were stored by
/// registerJacLIDs.
///
/// @see Xyce::Device::Resistor::Instance::registerJacLIDs
///
/// @note This method is called by the default implementation of the
/// loadDAEMatrices master function. Since the Resistor class
/// reimplements the "Master" "loadDAEMatrices" function that loads the
/// contributions from all resistor devices in a single loop, THIS
/// FUNCTION IS NOT ACTUALLY USED.  The loadDAEdFdx method is only
/// called when a device class does not re-implement the master class.
/// This can be a source of confusion when attempting to modify the Resistor
/// device, or any other device that reimplements the Master classes.
///
/// @see Xyce::Device::Resistor::Master::loadDAEMatrices
///
/// @return true on success
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   03/05/04
bool Instance::loadDAEdFdx()
{
  N_LAS_Matrix & dFdx = *(extData.dFdxMatrixPtr);
  dFdx[li_Pos][APosEquPosNodeOffset] += G;
  dFdx[li_Pos][APosEquNegNodeOffset] -= G;
  dFdx[li_Neg][ANegEquPosNodeOffset] -= G;
  dFdx[li_Neg][ANegEquNegNodeOffset] += G;
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Instance::updateTemperature
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Tom Russo, Component Information and Models
// Creation Date : 02/27/01
//-----------------------------------------------------------------------------
///
/// Update the parameters that depend on the temperature of the device
///
/// @param temp_tmp temperature
///
/// Xyce has a number of mechanisms that allow temperature to be changed
/// after a device has been instantiated.  These include .STEP loops over
/// temperature.  When temperature is changed, any device that has parameters
/// that depend on temperature must be updated.  That updating happens here.
///
/// The Resistor device supports temperature-dependent resistance through its
/// TC1 (linear dependence) and TC2 (quadratic dependence) parameters.
/// If these parameters are specified, the resistance must be updated.
///
/// @return true on success
///
/// @author Tom Russo, Component Information and Models
/// @date   02/27/01
bool Instance::updateTemperature(const double & temp_tmp)
{
  bool bsuccess = true;
  double difference, factor;

  if (temp_tmp != -999.0)
    temp = temp_tmp;
  difference = temp - model_.tnom;
  factor = 1.0 + tempCoeff1*difference + tempCoeff2*difference*difference;

  if (R*factor != 0.0)
    G = 1.0/(R * factor);
  else
    G = 0.0;

  return bsuccess;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Model::processParams
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 6/03/02
//-----------------------------------------------------------------------------
///
/// Process model parameters
///
/// @return true on success
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   6/03/02
bool Model::processParams()
{
  return true;
}

//----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Model::processInstanceParams
// Purpose       : 
// Special Notes : 
// Scope         : public
// Creator       : Dave Shirely, PSSI
// Creation Date : 03/23/06
//----------------------------------------------------------------------------
///
/// Process the instance parameters of instance owned by this model
///
/// This method simply loops over all instances associated with this
/// model and calls their processParams method.
///
/// @return true
///
/// @author Dave Shirely, PSSI
/// @date   03/23/06

bool Model::processInstanceParams()
{
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    (*it)->processParams();
  }

  return true;
}
//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Model::N_DEV_ResistorModel
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/16/00
//-----------------------------------------------------------------------------
///
/// Construct a resistor model from a "model block" that was created
/// by the netlist parser.
///
/// @param configuration
/// @param model_block
/// @param factory_block
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   5/16/00
Model::Model(
  const Configuration & configuration,
  const ModelBlock &    model_block,
  const FactoryBlock &  factory_block)
  : DeviceModel(model_block, configuration.getModelParameters(), factory_block),
    tempCoeff1(0.0),
    tempCoeff2(0.0),
    sheetRes(0.0),
    defWidth(10e-6),
    narrow(0.0),
    tnom(getDeviceOptions().tnom)
{
  // Set params to constant default values.
  setDefaultParams();

  // Set params according to .model line and constant defaults from metadata.
  setModParams(model_block.params);

  // Set any non-constant parameter defaults.
  if (!given("TNOM"))
    tnom = getDeviceOptions().tnom;

  // Calculate any parameters specified as expressions.
  updateDependentParameters();

  // calculate dependent (ie computed) params and check for errors.
  processParams();
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Model::~N_DEV_ResistorModel
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
///
/// Destroy this model.
///
/// Also destroys all instances that use this model.
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   3/16/00
Model::~Model()
{
  // Destory all owned instances
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    delete (*it);
  }
}

//-----------------------------------------------------------------------------
// Function      : N_DEV_ResistorModel::printOutInstances
// Purpose       : debugging tool.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 4/03/00
//-----------------------------------------------------------------------------
///
/// Print instances associated with this model.
///
/// Used only for debugging
///
/// @param os output stream
///
/// @return reference to output stream
///
/// @author Eric Keiter, SNL, Parallel Computational Sciences
/// @date   4/03/00
///
std::ostream &Model::printOutInstances(std::ostream &os) const
{
  os << std::endl;
  os << "Number of Resistor Instances: " << instanceContainer.size() << std::endl;
  os << "    name     model name  Parameters" << std::endl;

  int i = 0;
  for (InstanceVector::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
  {
    os << "  " << i << ": " << (*it)->getName() << "\t";
    os << getName();
    os << "\t\tR(Tnom) = " << (*it)->R;
    os << "\tG(T) = " << (*it)->G;
    os << std::endl;
    ++i;
  }

  os << std::endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : Model::forEachInstance
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 2/4/2014
//-----------------------------------------------------------------------------
/// Apply a device instance "op" to all instances associated with this
/// model
/// 
/// @param[in] op Operator to apply to all instances.
/// 
/// 
void Model::forEachInstance(DeviceInstanceOp &op) const /* override */ 
{
  for (std::vector<Instance *>::const_iterator it = instanceContainer.begin(); it != instanceContainer.end(); ++it)
    op(*it);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Master::updateState
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// Update state for all resistor instances, regardless of model.
///
/// @param solVec solution vector
/// @param staVec state vector
/// @param stoVec store vector
///
/// @return true on success
///
/// @note Because the resistor device re-implements the base-class
/// Master::updateState, the Instance::updatePrimaryState method is never
/// called, nor is the Instance::updateIntermediateVars method.  This method
/// replaces those, and does the same work but inside a loop over all
/// resistor instances.
///
/// @see Xyce::Device::Resistor::Instance::updatePrimaryState
/// @author Eric Keiter, SNL
/// @date   11/26/08

bool Master::updateState(double * solVec, double * staVec, double * stoVec)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ri = *(*it);

    double v_pos = solVec[ri.li_Pos];
    double v_neg = solVec[ri.li_Neg];

    // Load RHS vector element for the positive circuit node KCL equ.
    ri.i0 = (v_pos-v_neg)*ri.G;
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Master::loadDAEVectors
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// Load DAE vectors of all resistor instances, regardless of model
///
/// @param solVec solution vector
/// @param fVec f vector
/// @param qVec q vector
/// @param storeLeadF store lead current f vector
/// @param storeLeadQ store lead current q vector
///
/// @return true on success
///
/// @note Because the resistor device re-implements the base-class
/// Master::loadDAEVectors, the Instance::loadDAEFVector method is
/// never called.  This method replaces those, and does the same work
/// but inside a loop over all resistor instances.
///
/// @see Xyce::Device::Resistor::Instance::loadDAEFVector
///
/// @author Eric Keiter, SNL
/// @date   11/26/08

bool Master::loadDAEVectors (double * solVec, double * fVec, double *qVec,  double * storeLeadF, double * storeLeadQ)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ri = *(*it);
    fVec[ri.li_Pos] += ri.i0;
    fVec[ri.li_Neg] += -ri.i0;
    if( ri.loadLeadCurrent )
    {
      storeLeadF[ri.li_store_dev_i] = ri.i0;
    }
  }
  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Master::loadDAEMatrices
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL
// Creation Date : 11/26/08
//-----------------------------------------------------------------------------
///
/// Load DAE matrices for all resistor instances, regardless of model
///
/// @param dFdx matrix of derivatives of F vector with respect to solution
/// @param dQdx matrix of derivatives of Q vector with respect to solution
///
/// @return true on success
///
/// @note Because the resistor device re-implements the base-class
/// Master::loadDAEMatrices, the Instance::loadDAEdFdx method is
/// never called.  This method replaces those, and does the same work
/// but inside a loop over all resistor instances.
///
/// @see Xyce::Device::Resistor::Instance::loadDAEdFdx
///
/// @author Eric Keiter, SNL
/// @date   11/26/08

bool Master::loadDAEMatrices(N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx)
{
  for (InstanceVector::const_iterator it = getInstanceBegin(); it != getInstanceEnd(); ++it)
  {
    Instance & ri = *(*it);

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    *(ri.f_PosEquPosNodePtr) += ri.G;
    *(ri.f_PosEquNegNodePtr) -= ri.G;
    *(ri.f_NegEquPosNodePtr) -= ri.G;
    *(ri.f_NegEquNegNodePtr) += ri.G;
#else
    dFdx[ri.li_Pos][ri.APosEquPosNodeOffset] += ri.G;
    dFdx[ri.li_Pos][ri.APosEquNegNodeOffset] -= ri.G;
    dFdx[ri.li_Neg][ri.ANegEquPosNodeOffset] -= ri.G;
    dFdx[ri.li_Neg][ri.ANegEquNegNodeOffset] += ri.G;
#endif
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::Traits::factory
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 
//-----------------------------------------------------------------------------
///
/// Create a new instance of the Resistor device.
///
/// @param configuration
/// @param factory_block
///
Device *
Traits::factory(const Configuration &configuration, const FactoryBlock &factory_block)
{
  return new Master(configuration, factory_block, factory_block.solverState_, factory_block.deviceOptions_);
}

//-----------------------------------------------------------------------------
// Function      : Xyce::Device::Resistor::registerDevice
// Purpose       : 
// Special Notes :
// Scope         : public
// Creator       : David Baur
// Creation Date : 
//-----------------------------------------------------------------------------
///
/// Define how to use the device in a netlist.
///
/// This method is called from the Xyce::Device::registerOpenDevices
/// function, which in turn is called by the device manager.
///
/// The device is declared here to be an "R" device, which may optionally
/// take a model card of type "R".  This device will correspond to model
/// level 1 of resistor models.
void
registerDevice()
{
  Config<Traits>::addConfiguration()
    .registerDevice("r", 1)
    .registerModelType("r", 1);
}

} // namespace Resistor
} // namespace Device
} // namespace Xyce
