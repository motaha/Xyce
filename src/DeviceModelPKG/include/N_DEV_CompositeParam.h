//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2013  Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
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
// Filename       : $RCSfile: N_DEV_CompositeParam.h,v $
//
// Purpose        : This file contains the device entity base class.
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 05/05/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.16.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------
#ifndef Xyce_N_DEV_CompositeParam_h
#define Xyce_N_DEV_CompositeParam_h

// ---------- Standard Includes ----------
#include <string>
#include <vector>
#include <list>
#include <map>

// ------------- Xyce Includes ------------
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_Xyce.h>
#include <N_DEV_Param.h>
#include <N_DEV_Pars.h>

// ---------- Forward Declarations ----------

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : CompositeParam
//
// Purpose       : Base class for vector-composite classes.  This class mostly
//                 contains parameter processing similar to that of the
//                 DeviceEntity class.
//
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/05/05
//-----------------------------------------------------------------------------
class CompositeParam : public ParameterBase
{
  public:
    typedef Descriptor Pars;
    typedef Descriptor Parameter;
    typedef ParametricData<void>::ParameterMap ParameterMap;

    CompositeParam();
    virtual ~CompositeParam();

  private:
    CompositeParam(const CompositeParam &);
    CompositeParam &operator=(const CompositeParam &);

  public:
    virtual void processParams() {}

    void setDefaultParams();
    void setParams( const std::string &, const Param & );
    bool given(const string &parameter_name);

    std::ostream &printParams(std::ostream &os);

    const std::string &getName() const {
      return name_;
    }

  protected:
    virtual const ParametricData<void> &getMyParametricData() const = 0;

    const ParameterMap *getPMap() const {
      return reinterpret_cast<const ParameterMap *>(&getMyParametricData().getMap());
    }

  private:
    std::string           name_;
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::CompositeParam N_DEV_CompositeParam;

#endif

