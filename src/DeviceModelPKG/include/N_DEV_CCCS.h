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
// Filename       : $RCSfile: N_DEV_CCCS.h,v $
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
// Revision Number: $Revision: 1.11.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:36 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_CCCS_h
#define Xyce_N_DEV_CCCS_h

// ----------   Xyce Includes   ----------
#include <N_DEV_Source.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_Device.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : CCCSInstance
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class CCCSInstance : public DeviceInstance
{
  public:
      CCCSInstance(const CCCSInstance &right);
      ~CCCSInstance();

  protected:
  private:

  private:

    friend class CCCS;
    friend class CCCSModel;

};


//-----------------------------------------------------------------------------
// Class         : CCCSModel
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class CCCSModel : public DeviceModel
{
  public:
    CCCSModel(const CCCSModel &right);

    ~CCCSModel();

    CCCSModel & operator=(const CCCSModel &right);

    int operator==(const CCCSModel &right) const;

    int operator!=(const CCCSModel &right) const;

    const CCCSInstance * get_the_CCCSInstance () const;
    void set_the_CCCSInstance (CCCSInstance * value);


  protected:

  private:

  private:

    CCCSInstance *the_CCCSInstance;

    friend class CCCS;
    friend class CCCSInstance;

};

//-----------------------------------------------------------------------------
// Class         : CCCS
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
class CCCS : public Device
{
  public:

      CCCS(const CCCS &right);

      ~CCCS();

      CCCS & operator=(const CCCS &right);

      int operator==(const CCCS &right) const;

      int operator!=(const CCCS &right) const;


      const CCCSModel * get_the_CCCSModel () const;
      void set_the_CCCSModel (CCCSModel * value);

      static Device * factory(SolverState & ss1,
                    DeviceOptions & do1);

  protected:

  private:
      CCCS(SolverState & ss1,
                    DeviceOptions & do1);

  private:

    CCCSModel *the_CCCSModel;

    friend class CCCSModel;
    friend class CCCSInstance;
};

// Class CCCSInstance

// Class CCCSModel

//-----------------------------------------------------------------------------
// Function      : Models::get_the_CCCSInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline const CCCSInstance * CCCSModel::get_the_CCCSInstance () const
{
  return the_CCCSInstance;
}

//-----------------------------------------------------------------------------
// Function      : Models::set_the_CCCSInstance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline void CCCSModel::set_the_CCCSInstance (CCCSInstance * value)
{
  the_CCCSInstance = value;
}

// Class CCCS

//-----------------------------------------------------------------------------
// Function      : CCCS::get_the_CCCSModel
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline const CCCSModel * CCCS::get_the_CCCSModel () const
{
  return the_CCCSModel;
}

//-----------------------------------------------------------------------------
// Function      : CCCS::set_the_CCCSModel
// Purpose       : Factory function for the CCCS class.
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
inline void CCCS::set_the_CCCSModel (CCCSModel * value)
{
  the_CCCSModel = value;
}

//-----------------------------------------------------------------------------
// Function      : CCCS::factory
// Purpose       : Factory function for the CCCS class.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
inline Device * CCCS::factory(SolverState & ss1,
                    DeviceOptions & do1)
{
   Device * devptr = new CCCS(ss1,do1);
   return devptr;
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::CCCSInstance N_DEV_CCCSInstance;
typedef Xyce::Device::CCCSModel N_DEV_CCCSModel;
typedef Xyce::Device::CCCS N_DEV_CCCS;

#endif
