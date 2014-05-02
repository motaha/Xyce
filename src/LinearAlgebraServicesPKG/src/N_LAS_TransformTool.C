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
// Filename       : $RCSfile: N_LAS_TransformTool.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Robert Hoeksra, SNL, Parallel Computational Sciences
//
// Creation Date  : 4/2/03
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.59 $
//
// Revision Date  : $Date: 2014/02/24 23:49:23 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------

#include <EpetraExt_LPTrans_From_GraphTrans.h>
#include <EpetraExt_LPTrans_From_MatrixTrans.h>

#include <EpetraExt_CrsSingletonFilter_LinearProblem.h>

#include <EpetraExt_SolverMap_LinearProblem.h>

#ifdef Xyce_PARALLEL_MPI
#ifdef Xyce_USE_ISORROPIA
#include <EpetraExt_AmesosBTFGlobal_LinearProblem.h>
#include <EpetraExt_Isorropia_CrsGraph.h>
#endif
#endif

#ifdef Xyce_RCM
#include <EpetraExt_SymmRCM_CrsGraph.h>
#endif

#ifdef Xyce_AMD
#include <EpetraExt_AMD_CrsGraph.h>
#endif

#include <EpetraExt_AmesosBTF_CrsGraph.h>

#include <EpetraExt_Scale_LinearProblem.h>

#ifdef Xyce_TRILINOS_DEV
#include <EpetraExt_AmesosAMDGlobal_CrsGraph.h>
#endif

#include <EpetraExt_Reindex_LinearProblem.h>

#include <Teuchos_ScalarTraits.hpp>

// ----------   Xyce Includes   ----------

#include <N_LAS_TransformTool.h>

#include <N_UTL_OptionBlock.h>

#include <N_ERH_ErrorMgr.h>

//-----------------------------------------------------------------------------
// Function      : N_LAS_TransformTool::operator()
// Purpose       :
// Special Notes :
// Scope         : Public
// Creator       : Robert Hoekstra, SNL, Parallel Computational Sciences
// Creation Date : 04/2/03
//-----------------------------------------------------------------------------
Teuchos::RCP<N_LAS_Transform>
N_LAS_TransformTool::
operator()( const N_UTL_OptionBlock & options )
{
  int    rcm       = 0;
  int    rcmTestWidth = 0;

  int    sFilter   = 0;

  bool   scale     = false;
  int    lScale    = 0;
  int    rScale    = 0;
  double expScale  = 1.0;
  int    iterScale = 1;
  bool   btf       = false;
  int    globalbtf = 0;
  double globalbtftol = Teuchos::ScalarTraits<double>::eps();
  bool   globalbtfverb = false;
#ifdef Xyce_TRILINOS_DEV
  int    globalamd = 0;
  bool   globalamdverb = false;
#endif

  int    partition = 1;   // Zoltan via Isorropia partitioning
  std::string partition_type="GRAPH";

  // AMD should not be used in serial if KLU is being used.
#ifdef Xyce_PARALLEL_MPI
  bool   amd       = true;
#else
  bool   amd       = false;
#endif
  bool   amd_verbose = false;

  bool   reindex   = true;
  bool   solverMap = true;

  typedef std::list<N_UTL_Param>::const_iterator ParamListCIter;
  ParamListCIter iterL = options.getParams().begin();
  ParamListCIter endL  = options.getParams().end();
  for( ; iterL != endL; ++iterL )
  {
    std::string tag = iterL->uTag();

    if     ( tag == "TR_RCM" )               rcm = iterL->getImmutableValue<int>();
    else if( tag == "TR_RCM_TEST_WIDTH" )    rcmTestWidth= iterL->getImmutableValue<int>();
    else if( tag == "TR_SINGLETON_FILTER" )  sFilter = iterL->getImmutableValue<int>();
    else if( tag == "TR_SCALE" )             scale = iterL->getImmutableValue<int>();
    else if( tag == "TR_SCALE_LEFT" )        lScale = iterL->getImmutableValue<int>();
    else if( tag == "TR_SCALE_RIGHT" )       rScale = iterL->getImmutableValue<int>();
    else if( tag == "TR_SCALE_EXP" )         expScale = iterL->getImmutableValue<double>();
    else if( tag == "TR_SCALE_ITER" )        iterScale = iterL->getImmutableValue<int>();
    else if( tag == "TR_BTF" )               btf = iterL->getImmutableValue<int>();
    else if( tag == "TR_GLOBAL_BTF" )        globalbtf = iterL->getImmutableValue<int>();
    else if( tag == "TR_GLOBAL_BTF_DROPTOL" )globalbtftol = iterL->getImmutableValue<double>();
    else if( tag == "TR_GLOBAL_BTF_VERBOSE" )globalbtfverb = iterL->getImmutableValue<int>();
#ifdef Xyce_TRILINOS_DEV
    else if( tag == "TR_GLOBAL_AMD" )        globalamd = iterL->getImmutableValue<int>();
    else if( tag == "TR_GLOBAL_AMD_VERBOSE" )globalamdverb = iterL->getImmutableValue<int>();
#endif
    else if( tag == "TR_PARTITION" )         partition = iterL->getImmutableValue<int>();
#ifdef Xyce_USE_ISORROPIA
    else if( tag == "TR_PARTITION_TYPE" )    partition_type = iterL->usVal();
#endif
    else if( tag == "TR_AMD" )               amd = iterL->getImmutableValue<int>();
    else if( tag == "TR_AMD_VERBOSE" )       amd_verbose = iterL->getImmutableValue<int>();
    else if( tag == "TR_REINDEX" )           reindex = iterL->getImmutableValue<int>();
    else if( tag == "TR_SOLVERMAP" )         solverMap = iterL->getImmutableValue<int>();
  }

// Sort out which partitioning to use.
#ifndef Xyce_PARALLEL_MPI
  // If we are not running in parallel partitioning is not necessary.
  partition = 0;
#else
  // Turn off partitioning if Isorropia is not enabled.
#ifndef Xyce_USE_ISORROPIA
  if (partition == 1) {
    partition = 0;
  }
#endif
#endif // Xyce_PARALLEL_MPI

#ifndef Xyce_AMD
  amd = false;
#endif

#ifdef Xyce_VERBOSE_LINEAR
  Xyce::lout() << "Linear Transforms" << std::endl
               << "-----------------" << std::endl;
  if( sFilter )
    Xyce::lout() << "Singleton Filter" << std::endl;
#ifdef Xyce_TRILINOS_DEV
  if( globalamd )
    Xyce::lout() << "Global AMD" << std::endl;
#endif
  if( globalbtf )
    Xyce::lout() <<  "Global BTF" << std::endl;
  if( scale )
    Xyce::lout() << "Scaling" << std::endl;
#ifdef Xyce_USE_ISORROPIA
  if( partition == 1 ) {
    Xyce::lout() << "Isorropia Partitioning (" << partition_type << ") " << std::endl;
  }
#endif
  if( rcm )
    Xyce::lout() << "RCM" << std::endl;
  if( amd )
    Xyce::lout() << "AMD" << std::endl;
  if( btf )
    Xyce::lout() << "BTF" << std::endl;
  if( reindex )
    Xyce::lout() << "Reindexing" << std::endl;
  if( solverMap )
    Xyce::lout() << "Column Remapping" << std::endl;
  Xyce::lout() << "-----------------" << std::endl;
#endif

  Teuchos::RCP<N_LAS_Transform> CompTrans = Teuchos::rcp( new N_LAS_Transform() );

  bool TransFlag = false;

  // First remove any dense rows or columns using singleton filtering
  if( sFilter > 0 )
  {
    CompTrans->addTransform(
        dynamic_cast<EpetraExt::SameTypeTransform<Epetra_LinearProblem>*>
	(new EpetraExt::LinearProblem_CrsSingletonFilter(true)) );
    TransFlag = true;

#ifdef Xyce_PARALLEL_MPI
    if (partition_type == "GRAPH") {
      // Need to reindex because Isorropia does not do that at this time,
      // and symmetrization using different GID indices for row and column maps
      // can cause an issue.
      CompTrans->addTransform( new EpetraExt::LinearProblem_SolverMap() );

      CompTrans->addTransform( new EpetraExt::LinearProblem_Reindex( 0 ) );
    }
#endif

  }

  // Permute the global graph using AMD ordering (partitioning should be used in parallel)
#ifdef Xyce_TRILINOS_DEV
  if ( globalamd )
  {
    EpetraExt::AmesosAMDGlobal_CrsGraph * globalAMDTrans = new EpetraExt::AmesosAMDGlobal_CrsGraph();
    EpetraExt::LinearProblem_GraphTrans * globalAMD_LPTrans =
      new EpetraExt::LinearProblem_GraphTrans(
      *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(globalAMDTrans)) );
    CompTrans->addTransform( globalAMD_LPTrans );
    Teuchos::set_extra_data( Teuchos::rcp( globalAMDTrans ), "globalAMDTrans", Teuchos::inOutArg(CompTrans) );
    TransFlag = true;
  }
#endif

  // Permute the global graph using BTF ordering and a block partitioning (partitioning should *not* be used in parallel)
  if( globalbtf )
  {
#ifdef Xyce_PARALLEL_MPI
#ifdef Xyce_USE_ISORROPIA
    std::string balanceType = "linear";
    if (globalbtf == 2)
      balanceType = "isorropia";
    EpetraExt::AmesosBTFGlobal_LinearProblem * GlobalBTF_LPTrans =
        new EpetraExt::AmesosBTFGlobal_LinearProblem( globalbtftol, balanceType, false, globalbtfverb );
    CompTrans->addTransform( GlobalBTF_LPTrans );
    TransFlag = true;
#endif // Xyce_USE_ISORROPIA
#endif // Xyce_PARALLEL_MPI
  }

  // Partition the global graph using Isorropia
#ifdef Xyce_PARALLEL_MPI
  if( partition )
  {
    if (partition == 1) {
#ifdef Xyce_USE_ISORROPIA
      // Create the parameter list and fill it with values.
      Teuchos::ParameterList paramlist;
      Teuchos::ParameterList& sublist = paramlist.sublist("ZOLTAN");
      if (partition_type == "HYPERGRAPH") {
        paramlist.set("PARTITIONING METHOD", "HYPERGRAPH");
        sublist.set("CHECK_GRAPH", "0");
        sublist.set("LB_APPROACH", "PARTITION");
      }
      else if (partition_type == "GRAPH") {
        paramlist.set("PARTITIONING METHOD", "GRAPH");
        sublist.set("CHECK_GRAPH", "0");
        sublist.set("LB_APPROACH", "PARTITION");
        sublist.set("GRAPH_PACKAGE", "PARMETIS");
        sublist.set("PARMETIS_METHOD", "PARTKWAY");
      }
      else {
        N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Selected Partitioning Type is Invalid!\n" );
      }

#ifdef Xyce_VERBOSE_LINEAR
      sublist.set("DEBUG_LEVEL", "2" );
#endif
      EpetraExt::Isorropia_CrsGraph * ITrans = new EpetraExt::Isorropia_CrsGraph( paramlist );
      EpetraExt::LinearProblem_GraphTrans * I_LPTrans =
        new EpetraExt::LinearProblem_GraphTrans(
        *(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(ITrans)) );
      CompTrans->addTransform( I_LPTrans );
      Teuchos::set_extra_data( Teuchos::rcp( ITrans ), "ITrans", Teuchos::inOutArg(CompTrans) );
      TransFlag = true;
#endif
    }
    else {
      N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "Selected Partitioning Method is Invalid!\n" );
    }
  }
#endif

  // Perform RCM ordering on the local graph
  if( rcm )
  {
#ifdef Xyce_RCM
    bool brute = false;
    if( rcm == 2 ) brute = true;
    EpetraExt::CrsGraph_SymmRCM * RCMTrans = 0;
    if( rcmTestWidth != 0 )
      RCMTrans = new EpetraExt::CrsGraph_SymmRCM( brute, rcmTestWidth );
    else
      RCMTrans = new EpetraExt::CrsGraph_SymmRCM( brute );
    EpetraExt::LinearProblem_GraphTrans * RCM_LPTrans =
	new EpetraExt::LinearProblem_GraphTrans(
	*(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(RCMTrans)) );
    CompTrans->addTransform( RCM_LPTrans );
    Teuchos::set_extra_data( Teuchos::rcp( RCMTrans ), "RCMTrans", Teuchos::inOutArg(CompTrans) );
    TransFlag = true;
#else
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "TR_RCM not Enabled!\n" );
#endif
  }

  // Perform AMD ordering on the local graph
  if( amd )
  {
#ifdef Xyce_AMD
    EpetraExt::CrsGraph_AMD * AMDTrans = 0;
    AMDTrans = new EpetraExt::CrsGraph_AMD( amd_verbose );
    EpetraExt::LinearProblem_GraphTrans * AMD_LPTrans =
	new EpetraExt::LinearProblem_GraphTrans(
	*(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(AMDTrans)) );
    CompTrans->addTransform( AMD_LPTrans );
    Teuchos::set_extra_data( Teuchos::rcp( AMDTrans ), "AMDTrans", Teuchos::inOutArg(CompTrans) );
    TransFlag = true;
#else
    N_ERH_ErrorMgr::report( N_ERH_ErrorMgr::DEV_FATAL, "TR_AMD not Enabled!\n" );
#endif
  }

  // Perform BTF ordering on the local graph
  if( btf )
  {
    EpetraExt::AmesosBTF_CrsGraph * BTFTrans = new EpetraExt::AmesosBTF_CrsGraph();
    EpetraExt::LinearProblem_GraphTrans * BTF_LPTrans =
	new EpetraExt::LinearProblem_GraphTrans(
	*(dynamic_cast<EpetraExt::StructuralSameTypeTransform<Epetra_CrsGraph>*>(BTFTrans)) );
    CompTrans->addTransform( BTF_LPTrans );
    Teuchos::set_extra_data( Teuchos::rcp( BTFTrans ), "BTFTrans", Teuchos::inOutArg(CompTrans) );
    TransFlag = true;
  }

  // Scale the linear problem
  if ( scale )
  {
    CompTrans->addTransform(
	new EpetraExt::LinearProblem_Scale(
		static_cast<EpetraExt::LinearProblem_Scale::ScaleType>(lScale),
		static_cast<EpetraExt::LinearProblem_Scale::ScaleType>(rScale),
		expScale,
		iterScale ) );
    TransFlag = true;
  }

  if( solverMap && TransFlag )
  {
    CompTrans->addTransform( new EpetraExt::LinearProblem_SolverMap() );
    TransFlag = true;
  }

  if( reindex && TransFlag )
  {
    CompTrans->addTransform( new EpetraExt::LinearProblem_Reindex( 0 ) );
    TransFlag = true;
  }

  if( TransFlag ) return CompTrans;
  else
  {
    return Teuchos::null;
  }

}

