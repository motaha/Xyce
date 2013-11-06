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
// Filename       : $RCSfile: N_DEV_CCVS.h,v $
//
// Purpose        : Current-controlled voltage source.
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
// Revision Number: $Revision: 1.11.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_CCVS_h
#define Xyce_N_DEV_CCVS_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Source.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_Device.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : CCVSInstance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class CCVSInstance : public DeviceInstance
{
  public:
    CCVSInstance(const CCVSInstance &right);
    ~CCVSInstance();

  protected:
  private:

  private:
    friend class CCVS;
    friend class CCVSModel;

};

//-----------------------------------------------------------------------------
// Class         : CCVSModel
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class CCVSModel : public DeviceModel
{
  public:
    CCVSModel(const CCVSModel &right);
    ~CCVSModel();

    const CCVSInstance * get_the_CCVSInstance () const;
    void set_the_CCVSInstance (CCVSInstance * value);

  protected:
  private:

  private:
    CCVSInstance *the_CCVSInstance;

    friend class CCVS;
    friend class CCVSInstance;

};

//-----------------------------------------------------------------------------
// Class         : CCVS
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------

class CCVS : public Device
{
  public:
    CCVS(const CCVS &right);
    ~CCVS();


    const CCVSModel * get_the_CCVSModel () const;
    void set_the_CCVSModel (CCVSModel * value);

    static Device * factory(SolverState & ss1,
                    DeviceOptions & do1);

  protected:
  private:
    CCVS(SolverState & ss1,
                    DeviceOptions & do1);

  private:
    CCVSModel *the_CCVSModel;

    friend class CCVSInstance;
    friend class CCVSModel;
};

// Class CCVSInstance

// Class CCVSModel


//-----------------------------------------------------------------------------
// Function      : CCVSModel::get_the_CCVSInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline const CCVSInstance * CCVSModel::get_the_CCVSInstance () const
{
  return the_CCVSInstance;
}
//-----------------------------------------------------------------------------
// Function      : CCVSModel::set_the_CCVSInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline void CCVSModel::set_the_CCVSInstance (CCVSInstance * value)
{
  the_CCVSInstance = value;
}

// Class CCVS


//-----------------------------------------------------------------------------
// Function      : CCVS::get_the_CCVSModel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline const CCVSModel * CCVS::get_the_CCVSModel () const
{
  return the_CCVSModel;
}

//-----------------------------------------------------------------------------
// Function      : CCVS::set_the_CCVSModel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline void CCVS::set_the_CCVSModel (CCVSModel * value)
{
  the_CCVSModel = value;
}

//-----------------------------------------------------------------------------
// Function      : CCVS::factory
// Purpose       :
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
inline Device * CCVS::factory(SolverState & ss1,
                    DeviceOptions & do1)
{
   Device * devptr = new CCVS(ss1,do1);
   return devptr;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::CCVSInstance N_DEV_CCVSInstance;
typedef Xyce::Device::CCVSModel N_DEV_CCVSModel;
typedef Xyce::Device::CCVS N_DEV_CCVS;

#endif
