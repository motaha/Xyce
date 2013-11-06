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
// Filename      : $RCSfile: N_ANP_ModelEvaluator.C,v $
// Purpose       : This file contains the functions which define the 
//                 EpetraExt::ModelEvaluator interface for Xyce
// Special Notes :
// Creator       : Todd Coffey, 1414
// Creation Date : early 2009
//
// Revision Information:
// ---------------------
// Revision Number: $Revision: 1.11.2.2 $
// Revision Date  : $Date: 2013/10/03 17:23:31 $
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>


// ----------   Xyce Includes   ----------
#include <N_ANP_ModelEvaluator.h>
#include <N_LAS_Vector.h>
#include <N_LAS_BlockVector.h>
#include <N_LAS_BlockMatrix.h>
#include <N_PDS_ParMap.h>
// ----------   Trilinos Includes   ----------
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Operator.h>

using Teuchos::rcp;
using Teuchos::is_null;


RCP<N_LAS_BlockVector> convertEpetraToNLASBlockVectorView(
    const RCP<const Epetra_Vector>& vec,
    const Epetra_Map& map
    ) 
{
  // Add check that the size of vec is the same as 2*(size of map) 
  int Nv = vec->GlobalLength();
  int Nm = map.NumGlobalElements();
  TEUCHOS_TEST_FOR_EXCEPTION( Nv != 2*Nm, std::logic_error,
      "Error!  The size of the vector must be twice the size of the map!"
      );
  Epetra_Vector* nonconst_vec = const_cast<Epetra_Vector*>(&*vec);
  RCP<N_LAS_BlockVector> nlasBlockVec = rcp(new N_LAS_BlockVector( nonconst_vec, map, 2, false ));
  return nlasBlockVec;
}


RCP<N_LAS_Vector> convertEpetraToNLASVectorView(
    const RCP<const Epetra_Vector>& vec
    ) 
{
  Epetra_Vector* nonconst_vec = const_cast<Epetra_Vector*>(&*vec);
  RCP<N_LAS_Vector> nlasVec = rcp(new N_LAS_Vector( nonconst_vec, false ));
  return nlasVec;
}


N_ANP_ModelEvaluator::N_ANP_ModelEvaluator()
  :isInitialized_(false)
{
}


void N_ANP_ModelEvaluator::initialize(int iargs, char* cargs[])
{
  xycePtr_ = rcp(new N_CIR_Xyce());
  xycePtr_->initialize(iargs,cargs);
  xycePtr_->initializeTransientModel();
  this->setupInOutArgs_();
  this->setupMapsAndGraphs_(); // must be called after xycePtr_->initialize
  //voltLimFVector_ = rcp(new N_LAS_Vector(*x_map_,*x_ognd_map_));
  //voltLimQVector_ = rcp(new N_LAS_Vector(*x_map_,*x_ognd_map_));

  x_gnd_ = rcp(new N_LAS_Vector(*x_map_,*x_ognd_map_));
  z_gnd_ = rcp(new N_LAS_Vector(*x_map_,*x_ognd_map_));
  xdot_gnd_ = rcp(new N_LAS_Vector(*x_map_,*x_ognd_map_));
  zdot_gnd_ = rcp(new N_LAS_Vector(*x_map_,*x_ognd_map_));
  f_0_gnd_ = rcp(new N_LAS_Vector(*x_map_,*x_ognd_map_));
  f_0_gnd_.release(); // Xyce deletes this pointer in N_TIA_DataStore
  f_1_gnd_ = rcp(new N_LAS_Vector(*x_map_,*x_ognd_map_));
  f_1_gnd_.release(); // Xyce deletes this pointer in N_TIA_DataStore
  eVec_ = rcp(new N_LAS_Vector(*x_map_,*x_ognd_map_));
  eVec_->putScalar(1.0);

  dQdx_gnd_matrix_ = rcp(new N_LAS_Matrix( &*dQdx_ognd_graph_, &*dQdx_graph_ ));
  dQdx_gnd_matrix_.release(); // Xyce deletes this pointer in N_TIA_DataStore
  dFdx_gnd_matrix_ = rcp(new N_LAS_Matrix( &*dFdx_ognd_graph_, &*dFdx_graph_ ));
  dFdx_gnd_matrix_.release(); // Xyce deletes this pointer in N_TIA_DataStore

  isInitialized_ = true;
}


N_ANP_ModelEvaluator::~N_ANP_ModelEvaluator()
{
  //xycePtr_->finalize();       
}


std::vector<std::string> N_ANP_ModelEvaluator::getVariableNames()
{
  return xycePtr_->getVariableNames();
}


void N_ANP_ModelEvaluator::setupInOutArgs_()
{
  Np_ = 2;
  Ng_ = 3;
  EpetraExt::ModelEvaluator::InArgsSetup inArgs;
  inArgs.setSupports(IN_ARG_t,true);
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.setSupports(IN_ARG_x_dot,true);
  inArgs.setSupports(IN_ARG_alpha,true);
  inArgs.setSupports(IN_ARG_beta,true);
  inArgs.set_Np(Np_);
  inArgs_ = inArgs;

  EpetraExt::ModelEvaluator::OutArgsSetup outArgs;
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_Np_Ng(Np_,Ng_);
  outArgs_ = outArgs;
}


void N_ANP_ModelEvaluator::setupMapsAndGraphs_()
{
  xycePtr_->getMapsAndGraphs(x_map_,x_ognd_map_,s_map_,store_map_,dQdx_graph_,dQdx_ognd_graph_,dFdx_graph_,dFdx_ognd_graph_);

  // create block map for conversion between Xyce and Trilinos
  Epetra_Map& x_map = *x_map_->petraMap();
  int BlockSize = x_map.NumGlobalElements();
  int Size = 2; // two blocks
  int MaxGID = x_map.MaxAllGID();
  int offset = 1;
  while( offset <= MaxGID ) offset *= 10;
  int BaseIndex = x_map.IndexBase();
  const Epetra_Comm & Comm  = x_map.Comm();
  vector<int> blockGIDs;
  {
    blockGIDs.resize(Size*BlockSize);
    vector<int> BaseGIDs(BlockSize);
    x_map.MyGlobalElements( &BaseGIDs[0] );
    for( int i = 0; i < Size; ++i )
    {
      for( int j = 0; j < BlockSize; ++j )
      {
        blockGIDs[i*BlockSize+j] = BaseGIDs[j] + offset*i;
      }
    }
  }
  //std::cout << "N_ANP_ModelEvaluator::setupMapsAndGraphs_: blockGIDs = " << std::endl;
  //for (int n = 0 ; n<blockGIDs.size() ; ++n) {
  //  std::cout << "blockGIDs[" << n << "] = " << blockGIDs[n] << std::endl;
  //}
  blockMap_ = rcp(new Epetra_Map( BlockSize*Size,
                            BlockSize*Size,
                            &blockGIDs[0],
                            BaseIndex,
                            Comm 
                            )
                );

  // Augment the base graph to include the diagonal entries.
  dFdx_graph_with_diagonal_ = 
    rcp(new Epetra_CrsGraph( 
          Copy,
          dynamic_cast<Epetra_BlockMap&>(x_map),
          0
          )
       );
  // Copy graph non-zeros from dQdx and dFdx graphs:
  // W = [ dFdx, I ]
  //     [ dQdx, I ]
  int MaxIndices = dFdx_graph_->MaxNumIndices();
  vector<int> Indices;
  Indices.resize(MaxIndices+1);
  int NumIndices;
  int BaseRow;
  for( int j = 0; j < BlockSize; ++j )
  {
    BaseRow = x_map.GID(j);
    dFdx_graph_->ExtractGlobalRowCopy( BaseRow, MaxIndices, NumIndices, &Indices[0] );
    // Here we add the diagonal:
    Indices[NumIndices] = BaseRow;
    NumIndices++;
    dFdx_graph_with_diagonal_->InsertGlobalIndices( BaseRow, NumIndices, &Indices[0] );
  }
  dFdx_graph_with_diagonal_->TransformToLocal();
  //std::cout << "-----------------------------------------------" << std::endl;
  //std::cout << "dFdx_graph_with_diagonal = " << std::endl;
  //dFdx_graph_with_diagonal_->Print(std::cout);
  //std::cout << "-----------------------------------------------" << std::endl;
}


EpetraExt::ModelEvaluator::InArgs N_ANP_ModelEvaluator::createInArgs() const
{
  return inArgs_;
}


EpetraExt::ModelEvaluator::OutArgs N_ANP_ModelEvaluator::createOutArgs() const
{
  return outArgs_;
}


RCP<const Epetra_Map> N_ANP_ModelEvaluator::get_x_map() const
{
  return blockMap_;
}


RCP<const Epetra_Map> N_ANP_ModelEvaluator::get_f_map() const
{
  return blockMap_;
}


RCP<const Epetra_Map> N_ANP_ModelEvaluator::get_p_map(int p) const
{
  TEUCHOS_ASSERT( ((0 <= p) && (p < Np_)) );
  return Teuchos::rcpFromRef(*s_map_->petraMap());
}


RCP<const Epetra_Map> N_ANP_ModelEvaluator::get_g_map(int p) const
{
  TEUCHOS_ASSERT( ((0 <= p) && (p < Ng_)) );
  if (p == 0) { // s
    return Teuchos::rcpFromRef(*s_map_->petraMap());
  }
  else if (p == 1) { // voltLimQ
    return blockMap_;
  }
  else if (p == 2) { // voltLim F
    return blockMap_;
  }
  return Teuchos::null; // Should never get here.
}


RCP<const Epetra_Map> N_ANP_ModelEvaluator::get_small_x_map() const
{
  return Teuchos::rcpFromRef(*x_map_->petraMap());
}


void N_ANP_ModelEvaluator::evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const
{
  TEUCHOS_TEST_FOR_EXCEPTION( !isInitialized_, std::logic_error,
      "Error!  Please call initialize before evalModel"
      );
  // Get data out of inArgs
  double t = inArgs.get_t();
  RCP<const Epetra_Vector> x_in = inArgs.get_x().assert_not_null();
  RCP<const Epetra_Vector> xdot_in = inArgs.get_x_dot().assert_not_null();
  RCP<Epetra_Vector> s_out = outArgs.get_g(0);
  RCP<Epetra_Vector> store_out = outArgs.get_g(0); // ???

  RCP<N_LAS_BlockVector> xyce_x    = convertEpetraToNLASBlockVectorView(x_in,*(x_map_->petraMap()));
  RCP<N_LAS_BlockVector> xyce_xdot = convertEpetraToNLASBlockVectorView(xdot_in,*(x_map_->petraMap()));

  // xyce_x = (x,z)
  N_LAS_Vector & x = xyce_x->block(0);
  N_LAS_Vector & z = xyce_x->block(1);
  N_LAS_Vector & xdot = xyce_xdot->block(0);
  N_LAS_Vector & zdot = xyce_xdot->block(1);

  // The base N_LAS_Vector map has an extra element at the end for ground and
  // we're choosing for the moment to copy into this type of data structure
  // rather than fix our blockMap.
  x_gnd_->update(1.0,x,0.0); // x_gnd_ = x
  z_gnd_->update(1.0,z,0.0); // z_gnd_ = z
  xdot_gnd_->update(1.0,xdot,0.0); // xdot_gnd_ = xdot
  zdot_gnd_->update(1.0,zdot,0.0); // zdot_gnd_ = zdot

  N_LAS_Vector & x_ref = *x_gnd_;
  N_LAS_Vector & z_ref = *z_gnd_;
  //N_LAS_Vector & xdot_ref = *xdot_gnd_;
  N_LAS_Vector & zdot_ref = *zdot_gnd_;

  if (!Teuchos::is_null(s_out)) {
    // We're loading the state vectors
    RCP<N_LAS_Vector>      xyce_s_out    = convertEpetraToNLASVectorView(s_out);
    N_LAS_Vector & s_out = *xyce_s_out;

    RCP<N_LAS_Vector>      xyce_store_out    = convertEpetraToNLASVectorView(store_out);
    N_LAS_Vector & store_out = *xyce_store_out;

    xycePtr_->evalTransientModelState(
      t,
      & x_ref,
      & s_out,
      & store_out
      );

  } else {
    // We're loading f and W
    double alpha = inArgs.get_alpha();
    double beta = inArgs.get_beta();
    RCP<const Epetra_Vector> s_in = inArgs.get_p(0).assert_not_null();
    RCP<const Epetra_Vector> sdot_in = inArgs.get_p(1).assert_not_null();

    RCP<Epetra_Vector> f_out = outArgs.get_f().assert_not_null();
    RCP<Epetra_Operator> W_out = outArgs.get_W().assert_not_null();
    RCP<N_LAS_BlockMatrix> bMat;
    {
      // Pull out dFdx and dQdx matrices:
      RCP<Epetra_CrsMatrix> crsMat = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out,true);
      std::string label = "N_LAS_BlockMatrix";
      bMat = Teuchos::get_extra_data<RCP<N_LAS_BlockMatrix> >(crsMat,label);
      //N_LAS_Matrix& dFdx_tmp = bMat->block(0,0);
      //N_LAS_Matrix& dQdx_tmp = bMat->block(1,0);
      
			// Copy dFdx_tmp into dFdx_gnd_matrix_
      dFdx_gnd_matrix_->put(0.0);
			// Commented out
      //dFdx_gnd_matrix_->add(dFdx_tmp);
      
			// Copy dQdx_tmp into dQdx_gnd_matrix_
      dQdx_gnd_matrix_->put(0.0);
			// Commented out
      //dQdx_gnd_matrix_->add(dQdx_tmp);

      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "Before evalTransientModel:" << std::endl;
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "dQdx_gnd_matrix_ = " << std::endl;
      //dQdx_gnd_matrix_->printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
      //std::cout << "dFdx_gnd_matrix_ = " << std::endl;
      //dFdx_gnd_matrix_->printPetraObject();
      //std::cout << "-----------------------------------------------" << std::endl;
    }
    RCP<Epetra_Vector> voltLimQVector = outArgs.get_g(1);
    RCP<Epetra_Vector> voltLimFVector = outArgs.get_g(2);

    RCP<N_LAS_Vector>      xyce_s        = convertEpetraToNLASVectorView(s_in);
    RCP<N_LAS_Vector>      xyce_sdot     = convertEpetraToNLASVectorView(sdot_in);
    RCP<N_LAS_BlockVector> xyce_f        = convertEpetraToNLASBlockVectorView(f_out,*(x_map_->petraMap()));
    RCP<N_LAS_Vector>      xyce_store    = convertEpetraToNLASVectorView(s_in);
    // 08/05/09 tscoffe:  This is a hack to avoid having to deal with voltage
    // limiting for the moment.
    RCP<N_LAS_Vector>      xyce_voltLimQ;
    if (Teuchos::is_null(voltLimQVector)) {
      if (Teuchos::is_null(tempVoltLimQVector_)) {
        tempVoltLimQVector_ = rcp(new N_LAS_Vector(*x_map_,*x_ognd_map_));
      }
      xyce_voltLimQ = tempVoltLimQVector_;
    } 
    else {
      xyce_voltLimQ = convertEpetraToNLASVectorView(voltLimQVector);
    }
    RCP<N_LAS_Vector>      xyce_voltLimF;
    if (Teuchos::is_null(voltLimFVector)) {
      if (Teuchos::is_null(tempVoltLimFVector_)) {
        tempVoltLimFVector_ = rcp(new N_LAS_Vector(*x_map_,*x_ognd_map_));
      }
      xyce_voltLimF = tempVoltLimFVector_;
    }
    else {
      xyce_voltLimF = convertEpetraToNLASVectorView(voltLimFVector);
    }
    N_LAS_Vector & s = *xyce_s;
    N_LAS_Vector & store = *xyce_store;
    N_LAS_Vector & sdot = *xyce_sdot;
    N_LAS_Vector & f_0 = xyce_f->block(0);
    N_LAS_Vector & f_1 = xyce_f->block(1);
    //f_0_gnd_->update(1.0,f_0,0.0); // f_0_gnd_ = f_0
    //f_1_gnd_->update(1.0,f_1,0.0); // f_1_gnd_ = f_1
    f_0_gnd_->putScalar(0.0);
    f_1_gnd_->putScalar(0.0);
    N_LAS_Vector & f_0_ref = *f_0_gnd_;
    N_LAS_Vector & f_1_ref = *f_1_gnd_;
    N_LAS_Vector & voltLimQ_ref = *xyce_voltLimQ;
    N_LAS_Vector & voltLimF_ref = *xyce_voltLimF;

    // Eval model through AnalysisInterface
    // F(t,x,xdot,s,sdot) = [ zdot + f, q - z ]^T
    xycePtr_->evalTransientModel(
      t,
      & x_ref,
      & x_ref,
      & x_ref,
      & s,
      & s,
      & s,
      & sdot,
      & store,
      & store,
      & store,
      & store,   // this should e storLeadCurrQComp !
      & f_1_ref, // f_1 = q
      & f_0_ref, // f_0 = f
      & voltLimF_ref,
      & voltLimQ_ref,
      &* dQdx_gnd_matrix_,
      &* dFdx_gnd_matrix_
      );

    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "After evalTransientModel:" << std::endl;
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "dQdx_gnd_matrix_ = " << std::endl;
    //dQdx_gnd_matrix_->printPetraObject();
    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "dFdx_gnd_matrix_ = " << std::endl;
    //dFdx_gnd_matrix_->printPetraObject();
    //std::cout << "-----------------------------------------------" << std::endl;

    // Assemble outbound residual:
    // f_0 = zdot + f
    // f_1 = q - z
    f_0_ref.update(1.0,zdot_ref,+1.0);
    f_1_ref.update(-1.0,z_ref,1.0);

    // Now we copy the data back into our client's Epetra_Vectors
    f_0.update(1.0,f_0_ref,0.0); // f_0 = f_0_ref
    f_1.update(1.0,f_1_ref,0.0); // f_1 = f_1_ref

    {
      // Assemble the outbound Jacobian:
      // W = [ beta*dFdx_gnd_matrix_, alpha*I ]
      //     [ beta*dQdx_gnd_matrix_, -beta*I ]
      N_LAS_Matrix & ul = bMat->block(0,0);
      ul.put(0.0);
      ul.add(*dFdx_gnd_matrix_);
      ul.scale(beta);
      N_LAS_Matrix & ll = bMat->block(1,0);
      ll.put(0.0);
      ll.add(*dQdx_gnd_matrix_);
      ll.scale(beta);
      N_LAS_Matrix & ur = bMat->block(0,1);
      N_LAS_Matrix & lr = bMat->block(1,1);
      ur.put(0.0);
      lr.put(0.0);
      ur.replaceDiagonal(*eVec_);
      ur.scale(alpha);
      lr.replaceDiagonal(*eVec_);
      lr.scale(-beta);
    }

    //std::cout << "-----------------------------------------------" << std::endl;
    //std::cout << "outgoing:  bMat(Statefull): = " << std::endl;
    //bMat->printPetraObject();
    //std::cout << "-----------------------------------------------" << std::endl;
    
  }

}


RCP<Epetra_Operator> N_ANP_ModelEvaluator::create_W() const
{
  // Define which blocks of the blockMatrix are non-zero (ALL):
  int Size = 2;
  vector< vector<int> > Cols;
  Cols.resize(Size);
  for (int i=0 ; i<Size ; ++i) {
    Cols[i].resize(Size);
    for (int j=0 ; j<Size ; ++j) {
      Cols[i][j] = j;
    }
  }
  // Create the N_LAS_BlockMatrix:
  RCP<N_LAS_BlockMatrix> bMat = 
    N_LAS_blockMatrix( 
        Cols,
        *blockMap_,
        *dFdx_graph_with_diagonal_
        );
  bMat.release(); // will leak memory now.
  // Get an Epetra_CrsMatrix out of this
  Epetra_CrsMatrix & crsMat = bMat->epetraObj();
  // We attach the RCP<N_LAS_BlockMatrix> to the RCP object for the Epetra_Operator, so its not deleted.
  RCP<Epetra_CrsMatrix> outMatrix = rcp(&crsMat,false); // do not delete the underlying epetraObj()
  std::string label = "N_LAS_BlockMatrix";
  Teuchos::set_extra_data( bMat, label, Teuchos::outArg(outMatrix) );

  return outMatrix;
}


bool N_ANP_ModelEvaluator::isInitialized() const 
{ 
    return isInitialized_; 
}

