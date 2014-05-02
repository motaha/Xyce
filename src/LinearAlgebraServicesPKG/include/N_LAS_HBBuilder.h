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
// Filename       : $RCSfile: N_LAS_HBBuilder.h,v $
//
// Purpose        : Builder for HB specific linear objects
//
// Special Notes  :
//
// Creator        : Robert Hoekstra, SNL, Parallel Computational Sciences
//
// Creation Date  : 03/11/04
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.6 $
//
// Revision Date  : $Date: 2014/02/24 23:49:22 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef  Xyce_LAS_HBBUILDER_H
#define  Xyce_LAS_HBBUILDER_H

// ---------- Standard Includes ----------

#include <string>
#include <vector>
#include <Teuchos_RCP.hpp>

// ----------   Xyce Includes   ----------

#include <N_UTL_Xyce.h>

#include <N_LAS_Builder.h>

// ---------- Forward Declarations ----------

class N_LAS_Vector;
class N_LAS_BlockVector;
class N_LAS_MultiVector;
class N_LAS_Matrix;

class N_MPDE_Discretization;

class Epetra_MapColoring;
class Epetra_Map;
class Epetra_CrsGraph;

class N_PDS_ParMap;

//-----------------------------------------------------------------------------
// Class         : N_LAS_HBBuilder
// Purpose       :
// Special Notes :
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 3/11/04
//-----------------------------------------------------------------------------
class N_LAS_HBBuilder : public N_LAS_Builder
{

 public:

  // Default Constructor
  N_LAS_HBBuilder( const int Size,
                Teuchos::RCP<const N_MPDE_Discretization> Disc );

  // Destructor
  ~N_LAS_HBBuilder() {}

  // Vector and Matrix creators

  // Vector factory with initial value
  N_LAS_Vector * createVector( double initialValue = 0.0 ) const;

  // State Vector factory with initial value
  N_LAS_Vector * createStateVector( double initialValue = 0.0 ) const;

  // Store Vector factory with initial value
  N_LAS_Vector * createStoreVector( double initialValue = 0.0 ) const;

  // HB time-domain block vector creation:
  Teuchos::RCP<N_LAS_BlockVector> createTimeDomainBlockVector() const;
  Teuchos::RCP<N_LAS_BlockVector> createTimeDomainStateBlockVector() const;
  Teuchos::RCP<N_LAS_BlockVector> createTimeDomainStoreBlockVector() const;
  
  Teuchos::RCP<N_LAS_BlockVector> createExpandedRealFormBlockVector() const;
  Teuchos::RCP<N_LAS_BlockVector> createExpandedRealFormTransposeBlockVector() const;
  Teuchos::RCP<N_LAS_BlockVector> createExpandedRealFormTransposeStateBlockVector() const;
  Teuchos::RCP<N_LAS_BlockVector> createExpandedRealFormTransposeStoreBlockVector() const;

  // Matrix factory
  N_LAS_Matrix * createMatrix( double initialValue = 0.0 ) const
  { return createDAEFullMatrix(initialValue); }

  // DAE Jacobians
  N_LAS_Matrix * createDAEdQdxMatrix( double initialValue = 0.0 ) const;
  N_LAS_Matrix * createDAEdFdxMatrix( double initialValue = 0.0 ) const;
  N_LAS_Matrix * createDAEFullMatrix( double initialValue = 0.0 ) const;

  bool generateMaps( const Teuchos::RCP<N_PDS_ParMap>& BaseMap, 
                     const Teuchos::RCP<N_PDS_ParMap>& oBaseMap );

  bool generateStateMaps( const Teuchos::RCP<N_PDS_ParMap>& BaseStateMap );
  bool generateStoreMaps( const Teuchos::RCP<N_PDS_ParMap>& BaseStoreMap );

  bool generateGraphs( const Epetra_CrsGraph & BaseFullGraph );

  Teuchos::RCP<const Epetra_Map> getSolutionMap() const;
  Teuchos::RCP<const Epetra_Map> getStateMap() const;
  Teuchos::RCP<const Epetra_Map> getStoreMap() const;

  // Return the base map for each block in the expanded maps (a.k.a. time-domain maps)
  Teuchos::RCP<const N_PDS_ParMap> getBaseSolutionMap() const
  { return BaseMap_; }

  Teuchos::RCP<const N_PDS_ParMap> getBaseStateMap() const
  { return BaseStateMap_; }

  Teuchos::RCP<const N_PDS_ParMap> getBaseStoreMap() const
  { return BaseStoreMap_; }
  
  // Return GID offset for blocks for construction of Loader
  int getHBOffset()
  { return offset_; }

  int getHBStateOffset()
  { return stateOffset_; }

  int getHBStoreOffset()
  { return storeOffset_; }

private:

  const int numHarmonics_;
  int numSolVariables_, numStateVariables_;
  int numStoreVariables_;
  std::vector<std::vector<int> > blockPattern_;
  Teuchos::RCP<const N_MPDE_Discretization> Disc_;
  
  int offset_, stateOffset_;
  int storeOffset_;

  // HB maps for block vectors:
  // numBlocks = 2*(number of harmonics), numElem = number of solution variables
  Teuchos::RCP<N_PDS_ParMap> HBExpandedRealFormBVMap_; 
  Teuchos::RCP<N_PDS_ParMap> HBExpandedRealFormStateBVMap_;
  Teuchos::RCP<N_PDS_ParMap> HBExpandedRealFormStoreBVMap_;
  
  // numBlocks = number of solution variables, numElem = 2*(number of harmonics) 
  // We don't need a special map here, its the same as the non-tranpose

  Teuchos::RCP<N_PDS_ParMap> BaseMap_, oBaseMap_;

  Teuchos::RCP<N_PDS_ParMap> BaseStateMap_;
  Teuchos::RCP<N_PDS_ParMap> BaseStoreMap_;

  Teuchos::RCP<Epetra_CrsGraph> BaseFullGraph_;
  Teuchos::RCP<Epetra_CrsGraph> HBFullGraph_;

  Teuchos::RCP<N_PDS_ParMap> HBMap_, oHBMap_;

  Teuchos::RCP<N_PDS_ParMap> HBStateMap_;
  Teuchos::RCP<N_PDS_ParMap> HBStoreMap_;

};

#endif


