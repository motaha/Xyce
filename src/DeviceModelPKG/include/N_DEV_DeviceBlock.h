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
// Filename       : $RCSfile: N_DEV_DeviceBlock.h,v $
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
// Revision Number: $Revision: 1.20.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------


#ifndef Xyce_N_DEV_DeviceBlock_h
#define Xyce_N_DEV_DeviceBlock_h

// ---------- Standard Includes ----------

#include <string>
#include <vector>
#include <iosfwd>

// ----------   Xyce Includes   ----------

#include <N_DEV_Param.h>
#include <N_UTL_Misc.h>
#include <N_UTL_Xyce.h>
#include <N_UTL_Packable.h>

// ----------   Other Includes   ----------

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Class         : ModelBlock
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------
class ModelBlock : public Packable
{
  public:
    ModelBlock(const std::string &name_ = "", const std::string &type_ = "", int level_ = 1);
    ModelBlock(const ModelBlock & right);
    ~ModelBlock();

    ModelBlock & operator=(const ModelBlock & right);

    void clear ();

//Packing Utils
    Packable * instance() const;
    int packedByteCount() const;

    void pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const;
    void unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm );

    int operator==(const ModelBlock &right) const;
    int operator!=(const ModelBlock &right) const;

    string name;
    string type;
    int level;
    vector<Param> params;

    string netlistFileName_;
    int lineNumber_;

    friend ostream& operator<<(ostream& os, const ModelBlock & mb);
};

//-----------------------------------------------------------------------------
// Class         : InstanceBlock
// Purpose       :
// Special Notes :
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/31/00
//-----------------------------------------------------------------------------
class InstanceBlock : public Packable
{

public:
  InstanceBlock(const std::string & n_str = string() );
  InstanceBlock(const InstanceBlock & right);
  ~InstanceBlock();

  const std::string &getName() const {
    return name_;
  }
  void setName(const std::string &name) {
    name_ = name;
  }

  const std::string &getModelName() const {
    return modelName_;
  }
  void setModelName(const std::string &modelName) {
    modelName_ = modelName;
  }

  InstanceBlock & operator=(InstanceBlock & right);

  void clear ();

  //Packing Utils
  Packable * instance() const;
  int packedByteCount() const;

  void pack( char * buf, int bsize, int & pos, N_PDS_Comm * comm ) const;
  void unpack( char * pB, int bsize, int & pos, N_PDS_Comm * comm );

  int operator==(InstanceBlock &right) const;
  int operator!=(InstanceBlock &right) const;

private:
  string name_;
  string modelName_;

public:
  vector<Param> params;

  int iNumNodes;
  int numIntVars;
  int numExtVars;
  int numStateVars;

  bool modelFlag;
  bool sourceFlag;
  bool bsourceFlag;
  bool offFlag;
  bool off;

  string netlistFileName_;
  int lineNumber_;

  friend ostream& operator<<(ostream& os, const InstanceBlock & ib);
};

} // namespace Device
} // namespace Xyce

typedef Xyce::Device::InstanceBlock N_DEV_InstanceBlock;
typedef Xyce::Device::ModelBlock N_DEV_ModelBlock;

#endif
