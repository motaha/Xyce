//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright (C) 2002-2011  Sandia Corporation
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
// Filename       : $RCSfile: N_TOP_fwd.h,v $
//
// Purpose        : Forward declarations
//
// Special Notes  : Forward declaring everything as a class breaks if the implementation of the type changes (like during
//                  templatization)
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.2.2.1 $
//
// Revision Date  : $Date: 2014/03/03 18:29:28 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_TOP_fwd_h
#define Xyce_N_TOP_fwd_h

namespace Xyce {

class NodeID;

namespace Topo {

class CktGraph;
class CktGraphBasic;
class CktGraphCreator;
class CktGraphCreatorBasic;
class CktGraphSupport;
class CktNode;
class CktNodeCreator;
class CktNode_Ckt;
class CktNode_Dev;
class CktNode_V;
class Directory;
class Indexor;
class InsertionTool;
class Graph;
class GraphNode;
class Manager;
class Node;
class NodeBlock;
class NodeDevBlock;
class TopoLSUtil;
class Topology;
class System;

} // namespace Topo
} // namespace Xyce

typedef Xyce::Topo::CktGraph N_TOP_CktGraph;
typedef Xyce::Topo::CktGraphBasic N_TOP_CktGraphBasic;
typedef Xyce::Topo::CktGraphCreator N_TOP_CktGraphCreator;
typedef Xyce::Topo::CktGraphCreatorBasic N_TOP_CktGraphCreatorBasic;
typedef Xyce::Topo::CktGraphSupport N_TOP_CktGraphSupport;
typedef Xyce::Topo::CktNode N_TOP_CktNode;
typedef Xyce::Topo::CktNodeCreator N_TOP_CktNodeCreator;
typedef Xyce::Topo::CktNode_Ckt N_TOP_CktNode_Ckt;
typedef Xyce::Topo::CktNode_Dev N_TOP_CktNode_Dev;
typedef Xyce::Topo::CktNode_V N_TOP_CktNode_V;
typedef Xyce::Topo::Directory N_TOP_Directory;
typedef Xyce::Topo::Indexor N_TOP_Indexor;
typedef Xyce::Topo::Manager N_TOP_Manager;
typedef Xyce::Topo::Node N_TOP_Node;
typedef Xyce::Topo::NodeBlock N_TOP_NodeBlock;
typedef Xyce::Topo::NodeDevBlock N_TOP_NodeDevBlock;
typedef Xyce::Topo::TopoLSUtil N_TOP_TopoLSUtil;
typedef Xyce::Topo::Topology N_TOP_Topology;

typedef Xyce::NodeID NodeID;

#endif // Xyce_N_UTL_fwd_h
