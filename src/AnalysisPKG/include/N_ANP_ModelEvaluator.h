//-------------------------------------------------------------------------
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
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_ANP_ModelEvaluator.h,v $
//
// Purpose        : 
//                  
//                  
//
// Special Notes  : 
//                  
//
// Creator        : 
//
// Creation Date  : 
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.12 $
//
// Revision Date  : $Date: 2014/02/24 23:49:12 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------
#ifndef N_ANP_MODEL_EVALUATOR_H
#define N_ANP_MODEL_EVALUATOR_H

#undef HAVE_LIBPARMETIS
#include <EpetraExt_ModelEvaluator.h>
#include <N_CIR_Xyce.h>

// Forward Declarations
class Epetra_Map;
class Epetra_Vector;
class Epetra_CrsGraph;
class Epetra_Operator;
class N_LAS_Vector;
class N_LAS_BlockVector;

using Teuchos::RCP;

namespace Xyce {
namespace Analysis {

RCP<N_LAS_BlockVector> convertEpetraToNLASBlockVectorView(
    const RCP<const Epetra_Vector>& vec,
    const RCP<Epetra_Map>& map
    );
RCP<N_LAS_Vector> convertEpetraToNLASVectorView(
    const RCP<const Epetra_Vector>& vec
    );


class ModelEvaluator : public EpetraExt::ModelEvaluator {
public:
  ModelEvaluator();
  virtual ~ModelEvaluator();

  void initialize(int iargs, char* cargs[]);

  /** \name Overridden from EpetraExt::ModelEvaluator . */
  //@{
  Teuchos::RCP<const Epetra_Map> get_x_map() const;
  Teuchos::RCP<const Epetra_Map> get_f_map() const;
  Teuchos::RCP<const Epetra_Map> get_p_map(int p) const;
  Teuchos::RCP<const Epetra_Map> get_g_map(int p) const;
  Teuchos::RCP<const Epetra_Map> get_small_x_map() const;
  Teuchos::RCP<Epetra_Operator> create_W() const;
  EpetraExt::ModelEvaluator::InArgs createInArgs() const;
  EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
  void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;
  //@}

  std::vector<std::string> getVariableNames(); 

  bool isInitialized() const; 

private:
  // Private member functions:
  void setupInOutArgs_();
  void setupMapsAndGraphs_();

  // Private data:
  bool isInitialized_;
  RCP<N_CIR_Xyce> xycePtr_;
  EpetraExt::ModelEvaluator::InArgs inArgs_;
  EpetraExt::ModelEvaluator::OutArgs outArgs_;
  RCP<N_PDS_ParMap> x_map_;
  RCP<N_PDS_ParMap> x_ognd_map_;
  RCP<N_PDS_ParMap> s_map_;
  RCP<N_PDS_ParMap> store_map_;
  RCP<Epetra_CrsGraph> dQdx_graph_;
  RCP<Epetra_CrsGraph> dQdx_ognd_graph_;
  RCP<Epetra_CrsGraph> dFdx_graph_;
  RCP<Epetra_CrsGraph> dFdx_graph_with_diagonal_;
  RCP<Epetra_CrsGraph> dFdx_ognd_graph_;
  int Np_;
  int Ng_;
  mutable RCP<N_LAS_Vector> tempVoltLimFVector_;
  mutable RCP<N_LAS_Vector> tempVoltLimQVector_;
  RCP<N_LAS_Vector> eVec_; // Used for filling diagonal on identity matrices
  RCP<N_PDS_ParMap> blockMap_;

  // Used for Xyce loads where overlap_gnd map is required, so we copy the data in & out.
  RCP<N_LAS_Vector> x_gnd_;
  RCP<N_LAS_Vector> z_gnd_;
  RCP<N_LAS_Vector> xdot_gnd_;
  RCP<N_LAS_Vector> zdot_gnd_;
  RCP<N_LAS_Vector> f_0_gnd_;
  RCP<N_LAS_Vector> f_1_gnd_;
  RCP<N_LAS_Matrix> dQdx_gnd_matrix_;
  RCP<N_LAS_Matrix> dFdx_gnd_matrix_;

};

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::ModelEvaluator N_ANP_ModelEvaluator;

#endif // N_ANP_MODEL_EVALUATOR_H

