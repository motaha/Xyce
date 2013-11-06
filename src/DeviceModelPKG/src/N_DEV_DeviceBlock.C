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
// Filename       : $RCSfile: N_DEV_DeviceBlock.C,v $
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
// Revision Number: $Revision: 1.23.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


#include <iostream>

#include <N_DEV_DeviceBlock.h>
#include <N_ERH_ErrorMgr.h>

#include <N_PDS_Comm.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : ModelBlock::ModelBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
ModelBlock::ModelBlock(const std::string & name_, const std::string &type_, int level_)
  : name(name_),
    level(level_),
    type(type_),
    lineNumber_(0)
{
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::ModelBlock
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ModelBlock::ModelBlock(const ModelBlock &right)
  : name(right.name),
    type(right.type),
    level(right.level),
    netlistFileName_(right.netlistFileName_),
    lineNumber_(right.lineNumber_),
    params(right.params)
{
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::~ModelBlock
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ModelBlock::~ModelBlock()
{
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::operator=
// Purpose       : "=" operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
ModelBlock & ModelBlock::operator=(const ModelBlock &right)
{
  name   = right.name;
  type   = right.type;
  level =  right.level;

  netlistFileName_ = right.netlistFileName_;
  lineNumber_ = right.lineNumber_;

  params = right.params;

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::operator==
// Purpose       : "==" operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
int ModelBlock::operator==(const ModelBlock &right) const
{
  return name == right.name;
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::operator!=
// Purpose       : "!=" operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
int ModelBlock::operator!=(const ModelBlock &right) const
{
  return name != right.name;
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 4/06/00
//-----------------------------------------------------------------------------
ostream& operator<<(ostream & os, const ModelBlock & mb)
{
  vector<Param>::const_iterator it_pL, end_pL;

  os << "Model Block" << endl;
  os << "Model:  " << mb.name << endl;
  os << " type:  " << mb.type << endl;
  os << " Level: " << mb.level << endl;
  os << " netlistFileName_: " << mb.netlistFileName_ << endl;
  os << " lineNumber_: " << mb.lineNumber_ << endl;
  os << " Tagged Params" << endl;
  os << " -------------" << endl;

  it_pL=mb.params.begin();
  end_pL=mb.params.end();
  for ( ; it_pL != end_pL; ++it_pL)
  {
    os << it_pL->tag() << "\t" << it_pL->sVal() << endl;
  }

  os << " -------------" << endl;
  os << endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::clear
// Purpose       : empties out the block.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
void ModelBlock::clear()
{
  name = "";
  type = "";
  level = 0;

  netlistFileName_ = "";
  lineNumber_ = 0;

  params.clear();

}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::instance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/28/01
//-----------------------------------------------------------------------------
Packable * ModelBlock::instance() const
{
  return new ModelBlock();
}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::packedByteCount
// Purpose       : Counts bytes needed to pack block
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
int ModelBlock::packedByteCount() const
{

  int byteCount = 0;

  int size, length, i;

  vector<Param>::const_iterator it_tpL;

  //----- count name
  length = name.length();
  byteCount += sizeof(int);
  byteCount += length;

  //----- count type
  length = type.length();
  byteCount += sizeof(int);
  byteCount += length;

  //----- count level
  byteCount += sizeof(int);

  //----- count params
  size = params.size();
  byteCount += sizeof(int);
  it_tpL = params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    byteCount += it_tpL->packedByteCount();
  }

  //----- count netlistFileName_
  length = netlistFileName_.length();
  byteCount += sizeof(int);
  byteCount += length;

  //----- count lineNumber_
  byteCount += sizeof(int);

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::pack
// Purpose       : Packs ModelBlock into char buffer using MPI_PACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
void ModelBlock::pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const
{

  int size, length;
  int i;
  vector<Param>::const_iterator it_tpL;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  //----- pack name
  length = name.length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( name.c_str(), length, buf, bsize, pos );

  //----- pack type
  length = type.length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( type.c_str(), length, buf, bsize, pos );

  //----- pack level
  comm->pack(&level, 1, buf, bsize, pos );

  //----- pack params
  size = params.size();
  comm->pack(&size, 1, buf, bsize, pos );
  it_tpL = params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    it_tpL->pack( buf, bsize, pos, comm );
  }

  //----- pack netlistFileName_
  length = netlistFileName_.length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( netlistFileName_.c_str(), length, buf, bsize, pos );

  //----- packlineNumber_
  comm->pack(&lineNumber_, 1, buf, bsize, pos );

#ifdef Xyce_DEBUG_TOPOLOGY
  cout << "Packed " << pos << " bytes for ModelBlock: " <<
    name << endl;
#endif
#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    string msg = "Predicted pos does not match actual pos in ModelBlock::pack";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, msg );
  }
#endif

}

//-----------------------------------------------------------------------------
// Function      : ModelBlock::unpack
// Purpose       : Unpacks ModelBlock from char buffer using MPI_UNPACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
void ModelBlock::unpack(char * pB, int bsize,int & pos, N_PDS_Comm * comm)
{

  int size, length;
  int i;

  //----- unpack name
  comm->unpack( pB, bsize, pos, &length, 1 );
  name = string( (pB+pos), length);
  pos += length;

  //----- unpack type
  comm->unpack( pB, bsize, pos, &length, 1 );
  type = string( (pB+pos), length);
  pos += length;

  //----- unpack level
  comm->unpack( pB, bsize, pos, &level, 1 );

  //----- unpack params
  comm->unpack( pB, bsize, pos, &size, 1 );
  params.clear();
  Param dp;
  for( i = 0; i < size; ++i )
  {
    dp.unpack( pB, bsize, pos, comm );
    params.push_back( dp );
  }

  //----- unpack netlistFileName_
  comm->unpack( pB, bsize, pos, &length, 1 );
  netlistFileName_ = string( (pB+pos), length);
  pos += length;

  //----- unpack lineNumber_
  comm->unpack( pB, bsize, pos, &lineNumber_, 1 );

#ifdef Xyce_DEBUG_TOPOLOGY
  cout << "Unpacked " << pos << " bytes for ModelBlock: " <<
    name << endl;
#endif

}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::InstanceBlock
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/18/00
//-----------------------------------------------------------------------------
InstanceBlock::InstanceBlock (const std::string & n_str)
  : name_(n_str),
    modelName_(),
    iNumNodes(0),
    numIntVars(0),
    numExtVars(0),
    numStateVars(0),
    modelFlag(0),
    sourceFlag(0),
    bsourceFlag(0),
    offFlag(0),
    off(0),
    lineNumber_(0)
{
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::InstanceBlock
// Purpose       : copy constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
InstanceBlock::InstanceBlock (const InstanceBlock &right)
  : name_      (right.name_),
    modelName_(right.modelName_),
    iNumNodes (right.iNumNodes),
    numIntVars(right.numIntVars),
    numExtVars(right.numExtVars),
    numStateVars(right.numStateVars),
    modelFlag (right.modelFlag),
    sourceFlag(right.sourceFlag),
    bsourceFlag(right.bsourceFlag),
    offFlag   (right.offFlag),
    off       (right.off),
    netlistFileName_  (right.netlistFileName_),
    lineNumber_       (right.lineNumber_),
    params    (right.params)
{
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::~InstanceBlock
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
InstanceBlock::~InstanceBlock ()
{
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::operator=
// Purpose       : "=" operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
InstanceBlock & InstanceBlock::operator=(InstanceBlock &right)
{
  name_      = right.name_;
  modelName_ = right.modelName_;
  iNumNodes = right.iNumNodes;
  numIntVars= right.numIntVars;
  numExtVars= right.numExtVars;
  numStateVars= right.numStateVars;
  modelFlag = right.modelFlag;
  sourceFlag= right.sourceFlag;
  bsourceFlag= right.bsourceFlag;
  offFlag   = right.offFlag;
  off       = right.off;
  netlistFileName_  = right.netlistFileName_;
  lineNumber_       = right.lineNumber_;
  params    = right.params;

  return *this;
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::operator==
// Purpose       : "==" operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
int InstanceBlock::operator==(InstanceBlock &right) const
{
  return name_ == right.name_;
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::operator!=
// Purpose       : "!=" operator
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/16/00
//-----------------------------------------------------------------------------
int InstanceBlock::operator!=(InstanceBlock &right) const
{
  return name_ != right.name_;
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::clear
// Purpose       : empties out the block.
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 5/01/00
//-----------------------------------------------------------------------------
void InstanceBlock::clear ()
{
  name_ = "";
  modelName_  = "";
  iNumNodes  = 0;
  numIntVars = 0;
  numExtVars = 0;
  numStateVars = 0;
  modelFlag  = 0;
  sourceFlag = 0;
  bsourceFlag = 0;
  offFlag    = 0;
  off        = 0;
  netlistFileName_  = "";
  lineNumber_       = 0;

  params.clear();
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::operator<<
// Purpose       : "<<" operator
// Special Notes :
// Scope         : public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 4/06/00
//-----------------------------------------------------------------------------
ostream& operator<<(ostream & os, const InstanceBlock & ib)
{
  vector<Param>::const_iterator it_tpL, end_tpL;

  os << "Instance Block" << endl;
  os << "Name:    " << ib.name_ << endl;
  os << " Model:  " << ib.getModelName() << endl;
  os << " # Nodes: " << ib.iNumNodes << endl;
  os << " # Int Vars: " << ib.numIntVars << endl;
  os << " # Ext Vars: " << ib.numExtVars << endl;
  os << " # State Vars: " << ib.numStateVars << endl;
  os << " modelFlag: " << ib.modelFlag << endl;
  os << " sourceFlag: " << ib.sourceFlag << endl;
  os << " bsourceFlag: " << ib.bsourceFlag << endl;
  os << " offFlag: " << ib.offFlag << endl;
  os << " off: " << ib.off << endl;
  os << " netlistFileName_: " << ib.netlistFileName_ << endl;
  os << " lineNumber_: " << ib.lineNumber_ << endl;
  os << " Tagged Params" << endl;
  os << " -------------" << endl;

  it_tpL=ib.params.begin();
  end_tpL=ib.params.end();
  for ( ; it_tpL != end_tpL; ++it_tpL)
  {
    os << it_tpL->tag() << "\t" << it_tpL->sVal() << endl;
  }

  os << " -------------" << endl;
  os << endl;

  return os;
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::instance
// Purpose       :
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/28/01
//-----------------------------------------------------------------------------
Packable * InstanceBlock::instance() const
{
  return new InstanceBlock();
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::packedByteCount
// Purpose       : count bytes needed to pack block
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 6/13/00
//-----------------------------------------------------------------------------
int InstanceBlock::packedByteCount() const
{

  int byteCount = 0;

  int size, length, i;
  vector<Param>::const_iterator it_tpL;

  //----- count name
  length = getName().length();
  byteCount += sizeof(int);
  byteCount += length * sizeof(char);

  //----- count getModelName()
  length = getModelName().length();
  byteCount += sizeof(int);
  byteCount += length * sizeof(char);

  //----- count params
  size = params.size();
  byteCount += sizeof(int);
  it_tpL = params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    byteCount += it_tpL->packedByteCount();
  }

  //----- count iNumNodes
  byteCount += sizeof(int);

  //----- count numIntVars
  byteCount += sizeof(int);

  //----- countnumExtVars
  byteCount += sizeof(int);

  //----- count numStateVars
  byteCount += sizeof(int);

  //----- count modelFlag
  byteCount += sizeof(int);

  //----- count sourceFlag
  byteCount += sizeof(int);

  //----- count bsourceFlag
  byteCount += sizeof(int);

  //----- count offFlag
  byteCount += sizeof(int);

  //----- pack off
  byteCount += sizeof(int);

  //----- count netlistFileName_
  length = netlistFileName_.length();
  byteCount += sizeof(int);
  byteCount += length * sizeof(char);

  //----- count lineNumber_
  byteCount += sizeof(int);

  return byteCount;

}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::pack
// Purpose       : Packs InstanceBlock into char buffer using MPI_PACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
void InstanceBlock::pack(char * buf, int bsize, int & pos, N_PDS_Comm * comm) const
{

  int size, length;
  int i;
  vector<Param>::const_iterator it_tpL;
#ifdef Xyce_COUNT_PACKED_BYTES
  int predictedPos = pos+packedByteCount();
#endif

  //----- pack name
  length = getName().length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( name_.c_str(), length, buf, bsize, pos );

  //----- pack getModelName()
  length = getModelName().length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( getModelName().c_str(), length, buf, bsize, pos );

  //----- pack params
  size = params.size();
  comm->pack(&size, 1, buf, bsize, pos );
  it_tpL = params.begin();
  for (i = 0; i < size; ++i, ++it_tpL)
  {
    it_tpL->pack( buf, bsize, pos, comm );
  }

  //----- pack iNumNodes
  comm->pack(&iNumNodes, 1, buf, bsize, pos );

  //----- pack numIntVars
  comm->pack(&numIntVars, 1, buf, bsize, pos );

  //----- pack numExtVars
  comm->pack(&numExtVars, 1, buf, bsize, pos );

  //----- pack numStateVars
  comm->pack(&numStateVars, 1, buf, bsize, pos );

  //----- pack modelFlag
  i = modelFlag;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack sourceFlag
  i = sourceFlag;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack bsourceFlag
  i = bsourceFlag;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack offFlag
  i = offFlag;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack off
  i = off;
  comm->pack(&i, 1, buf, bsize, pos );

  //----- pack name
  length = netlistFileName_.length();
  comm->pack(&length, 1, buf, bsize, pos );
  comm->pack( netlistFileName_.c_str(), length, buf, bsize, pos );

  //----- pack lineNumber_
  comm->pack(&lineNumber_, 1, buf, bsize, pos );

#ifdef Xyce_DEBUG_TOPOLOGY
  cout << "Packed " << pos << " bytes for InstanceBlock: " <<
    getName() << endl;
#endif
#ifdef Xyce_COUNT_PACKED_BYTES
  if (pos != predictedPos)
  {
    string msg = "Predicted pos does not match actual pos in InstanceBlock::pack";
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, msg );
  }
#endif
}

//-----------------------------------------------------------------------------
// Function      : InstanceBlock::unpack
// Purpose       : Unpacks InstanceBlock from char buffer using MPI_UNPACK
// Special Notes :
// Scope         : public
// Creator       : Rob Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 5/24/00
//-----------------------------------------------------------------------------
void InstanceBlock::unpack(char * pB, int bsize, int & pos, N_PDS_Comm * comm)
{

  int size, length;
  int i;

  //----- unpack name
  comm->unpack( pB, bsize, pos, &length, 1 );
  name_ = string( (pB+pos), length);
  pos += length;

  //----- unpack getModelName()
  comm->unpack( pB, bsize, pos, &length, 1 );
  modelName_ = string( (pB+pos), length);
  pos += length;

  //----- unpack params
  comm->unpack( pB, bsize, pos, &size, 1 );
  params.clear();
  Param dp;
  for( i = 0; i < size; ++i )
  {
    dp.unpack( pB, bsize, pos, comm );
    params.push_back( dp );
  }

  //----- unpack iNumNodes
  comm->unpack( pB, bsize, pos, &iNumNodes, 1 );

  //----- unpack numIntVars
  comm->unpack( pB, bsize, pos, &numIntVars, 1 );

  //----- unpack numExtVars
  comm->unpack( pB, bsize, pos, &numExtVars, 1 );

  //----- unpack numStateVars
  comm->unpack( pB, bsize, pos, &numStateVars, 1 );

  //----- unpack modelFlag
  comm->unpack( pB, bsize, pos, &i, 1 );
  modelFlag = ( i != 0 );

  //----- unpack sourceFlag
  comm->unpack( pB, bsize, pos, &i, 1 );
  sourceFlag = ( i != 0 );

  //----- unpack bsourceFlag
  comm->unpack( pB, bsize, pos, &i, 1 );
  bsourceFlag = ( i != 0 );

  //----- unpack offFlag
  comm->unpack( pB, bsize, pos, &i, 1 );
  offFlag = ( i != 0 );

  //----- unpack off
  comm->unpack( pB, bsize, pos, &i, 1 );
  off = ( i != 0 );

  //----- unpack netlistFileName_
  comm->unpack( pB, bsize, pos, &length, 1 );
  netlistFileName_ = string( (pB+pos), length);
  pos += length;

  //----- unpack lineNumber_
  comm->unpack( pB, bsize, pos, &lineNumber_, 1 );

#ifdef Xyce_DEBUG_TOPOLOGY
  cout << "Unpacked " << pos << " bytes for InstanceBlock: " <<
    getName() << endl;
#endif

}

} // namespace Device
} // namespace Xyce
