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
// Filename       : $RCSfile: N_LAS_Builder.h,v $
//
// Purpose        : Builder for LAS objects, hides parallel map stuff
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 01/22/01
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.26 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  Xyce_LAS_BUILDER_H
#define  Xyce_LAS_BUILDER_H

// ---------- Standard Includes ----------

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>
// to eliminate RCP warnings, putting N_PDS_Manager header here.
#include <N_PDS_Manager.h>

// ---------- Forward Declarations ----------

class N_PDS_ParMap;

class N_LAS_QueryUtil;
class N_LAS_System;

class N_LAS_Vector;
class N_LAS_MultiVector;
class N_LAS_Matrix;

class Epetra_MapColoring;
class Epetra_CrsGraph;

using Teuchos::RCP;
using Teuchos::rcp;

//-----------------------------------------------------------------------------
// Class         : N_LAS_Builder
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 1/22/02
//-----------------------------------------------------------------------------
class N_LAS_Builder
{

public:

  // Default Constructor
  N_LAS_Builder()
  : pdsMgr_(0),
    lasQueryUtil_(0),
    lasSystem_(0)
  {}

  // Destructor
  virtual ~N_LAS_Builder();

  // Registration methods for necessary utilities
  virtual bool registerPDSManager(N_PDS_Manager * PDS_Manager);

  virtual bool registerPDSManager(RCP<N_PDS_Manager> PDS_Manager);

  virtual bool registerQueryUtil(N_LAS_QueryUtil * LAS_QUtil)
  { return (lasQueryUtil_ = LAS_QUtil); }

  virtual bool registerSystem(N_LAS_System * system)
  { return (lasSystem_ = system); }

  // Vector and Matrix creators which use N_LAS_QueryUtil and N_PDS_ParMap
  // attributes to transparently construct proper objects for this linear
  // system

  // Multivector factory with num vectors and initial value
  virtual N_LAS_MultiVector * createMultiVector( const int numVectors = 1, const double value = 0.0 ) const;
  // State Multivector factory with num vectors and initial value
  virtual N_LAS_MultiVector * createStateMultiVector( const int numVectors = 1, const double value = 0.0 ) const;
  // Store Multivector factory with num vectors and initial value
  virtual N_LAS_MultiVector * createStoreMultiVector( const int numVectors = 1, const double value = 0.0 ) const;
  // Vector factory with initial value
  virtual N_LAS_Vector * createVector( const double value = 0.0 ) const;
  // State-vector factory
  virtual N_LAS_Vector * createStateVector( const double value = 0.0 ) const;
  // Store-vector factory
  virtual N_LAS_Vector * createStoreVector( const double value = 0.0 ) const;

  // Matrix factory
  virtual N_LAS_Matrix * createMatrix( const double initialValue = 0.0 ) const;

  //Coloring Assoc with Variable Types in Solution Vector
  virtual Epetra_MapColoring * createSolnColoring() const;

  //Coloring needed for imposing .IC and .NODESET
  virtual Epetra_MapColoring * createInitialConditionColoring() const;

  virtual bool generateParMaps();

  virtual bool generateGraphs();

  virtual RCP<const Epetra_Map> getSolutionMap() const;
  virtual RCP<const Epetra_Map> getStateMap() const;
  virtual RCP<const Epetra_Map> getStoreMap() const;

  virtual void getSolutionMaps(
      RCP<N_PDS_ParMap>& soln_map,
      RCP<N_PDS_ParMap>& soln_ognd_map
      ) const;
  virtual void getStateMap(
      RCP<N_PDS_ParMap>& state_map
      ) const;
  virtual void getStoreMap(
      RCP<N_PDS_ParMap>& store_map
      ) const;
  virtual void getdQdxGraphs(
      RCP<Epetra_CrsGraph>& dQdx_graph,
      RCP<Epetra_CrsGraph>& dQdx_ognd_graph
      ) const;
  virtual void getdFdxGraphs(
      RCP<Epetra_CrsGraph>& dFdx_graph,
      RCP<Epetra_CrsGraph>& dFdx_ognd_graph
      ) const;

private:

  RCP<N_PDS_Manager> pdsMgr_;
  N_LAS_QueryUtil * lasQueryUtil_;
  N_LAS_System * lasSystem_;
};

#endif


