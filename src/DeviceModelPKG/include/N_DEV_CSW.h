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
// Filename       : $RCSfile: N_DEV_CSW.h,v $
//
// Purpose        : Switch base classes
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
// Revision Number: $Revision: 1.15.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_CSW_h
#define Xyce_N_DEV_CSW_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Device.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

#include <N_DEV_Param.h>

namespace Xyce {
namespace Device {

class CSWModel;
class MatrixLoadData;

//-----------------------------------------------------------------------------
// Class         : CSWInstance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class CSWInstance : public DeviceInstance
{
  public:
    CSWInstance(InstanceBlock &IB,
                     MatrixLoadData & mlData1,
                     SolverState &ss1,
                     ExternData  &ed1,
                     DeviceOptions & do1);

    CSWInstance(const CSWInstance &right);

    ~CSWInstance();

    void registerLIDs( const vector<int> & intLIDVecRef,
	               const vector<int> & extLIDVecRef );
    void registerStateLIDs( const vector<int> & staLIDVecRef );

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    bool updateIntermediateVarsBlock ();
    bool updatePrimaryStateBlock ();
    bool updateSecondaryStateBlock ();

  protected:
  private:

  public:
  protected:
  private:
    // user specified parameters
    double R; // resistance (ohms)
    bool ON; // whether switch is on or off initially

    // derived parameters
    double G; // conductance (1.0/ohms)

    // places to store node voltages
    double v_pos;
    double v_neg;

    // double for current state of switch
    double SW_STATE;

    // and a state variable to save the SW_STATE
    double switch_state;
    bool switch_state_initialized;

    int li_switch_state;

    // iterator reference to the switch model which owns this instance
    vector<CSWModel*>::iterator SWM_iter;

    // local indices (offsets)
    int li_Pos;
    int li_Neg;
    int li_PosCntl;
    int li_NegCntl;

    // Pointer variables corresponding to the above declared indices.
    int APosEquPosNodeOffset;
    int APosEquNegNodeOffset;
    int ANegEquPosNodeOffset;
    int ANegEquNegNodeOffset;

    static vector< vector<int> > jacStamp;

    friend class CSW;
    friend class CSWModel;
};

//-----------------------------------------------------------------------------
// Class         : CSWModel
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class CSWModel : public DeviceModel
{
  public:
    CSWModel(const CSWModel &right);
    CSWModel(const ModelBlock &MB,
                        SolverState & ss1,
                     DeviceOptions & do1);

    ~CSWModel();

  virtual std::ostream &printOutInstances(std::ostream &os) const;

  protected:
  private:

  public:
  protected:
  private:
    double RON;
    double ROFF;
    double VON;
    double VOFF;
    double VdInv;   // the inverse of (VON-VOFF) or 1e-12, if too small.
    double Lm; //log mean of resistances
    double Lr; //log ratio of resistor values

    vector<CSWInstance*> instanceContainer;

    friend class CSW;
    friend class CSWInstance;
};

//-----------------------------------------------------------------------------
// Class         : CSW
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class CSW : public Device
{
  public:
    ~CSW();

    static Device * factory(SolverState & ss1,
                     DeviceOptions & do1);
    string & deviceName();

    DeviceModel * addModel( const ModelBlock & MB );
    DeviceInstance * addInstance( InstanceBlock & IB,
                                        MatrixLoadData & mlData1,
                                        SolverState &ss1,
                                        ExternData  &ed1,
                                        DeviceOptions & do1);

    void printOutModels();

  DeviceModel * getModelPointer( const std::string & modelNameRef,
                                           bool & bsuccess );

  bool deleteInstance (const std::string & tmpName);

  protected:

  private:
    CSW(SolverState & ss1,
                     DeviceOptions & do1);
    DeviceModel * addModel_internal( const ModelBlock & MB,
                              vector<CSWModel*>::iterator & SWM_iter );

    CSW( const CSW & right );

  public:

  protected:

  private:
    vector<CSWModel*> modelContainer;

    friend class CSWModel;
    friend class CSWInstance;
};

//-----------------------------------------------------------------------------
// Function      : CSW::factory
// Purpose       : Factory function for the CSW class.
//
// Special Notes : ERK.  10/16/2005.  This used to be a singleton (ie a
//                 static pointer was returned) but had to be changed
//                 so that the library version of Xyce would work
//                 correctly.
//
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline Device * CSW::factory(SolverState & ss1,
                                        DeviceOptions & do1)
{
   Device * devptr = new CSW(ss1,do1);
   return devptr;
}

//-----------------------------------------------------------------------------
// Function      : CSW::deviceName
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline string & CSW::deviceName()
{
  return name;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::CSWInstance N_DEV_CSWInstance;
typedef Xyce::Device::CSWModel N_DEV_CSWModel;
typedef Xyce::Device::CSW N_DEV_CSW;

#endif
