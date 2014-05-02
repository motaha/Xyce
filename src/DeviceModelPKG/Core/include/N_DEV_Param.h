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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Param.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric Keiter
//
// Creation Date  : 5/22/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.32.2.1 $
//
// Revision Date  : $Date: 2014/02/26 20:16:30 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef  N_DEV_PARAM_H
#define  N_DEV_PARAM_H

#include <iosfwd>
#include <string>
#include <vector>
#include <map>

#include <N_DEV_fwd.h>
#include <N_UTL_NoCase.h>
#include <N_UTL_Param.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : N_DEV_Param
// Purpose       :
// Special Notes :
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/10/01
//-----------------------------------------------------------------------------
class Param : public Util::Param
{
public:
  Param()
    : Util::Param(),
      isGiven_(false),
      isDefault_(false)
  {}

  template <class T>
  Param(const std::string &tag, const T &value, bool is_given = false)
    : Util::Param(tag, value),
      isGiven_(is_given),
      isDefault_(false)
  {}

  Param(const Param &rhsParam)
    : Util::Param(rhsParam),
      isGiven_(rhsParam.isGiven_),
      isDefault_(rhsParam.isDefault_)
  {}

  Param &operator=(const Param &rhsParam) 
  {
    Util::Param::operator=(rhsParam);
    isGiven_ = rhsParam.isGiven_;
    isDefault_ = rhsParam.isDefault_;

    return *this;
  }

  virtual ~Param()
  {}

  void setGiven(bool is_given) 
  {
    isGiven_ = is_given;
  }

  void setDefault(bool is_default) 
  {
    isDefault_ = is_default;
  }

  bool given() const 
  {
    return isGiven_;
  }

  bool default_val() const 
  {
    return isDefault_;
  }


  virtual Packable * instance() const /* override */ 
  {
    return new Param();
  }

  virtual int packedByteCount() const /* override */ ;
  virtual void pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const /* override */ ;
  virtual void unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm ) /* override */ ;

private:
  bool isGiven_;
  bool isDefault_;
};

inline void setParamValue(Param &param, const Param &from_param) 
{
  param.setVal(static_cast<const Util::Param &>(from_param));
}

inline void setParam(Param &param, const std::string &tag, const Param &from_param) 
{
  param.set(tag, static_cast<const Util::Param &>(from_param));
}

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Param N_DEV_Param;

#endif
