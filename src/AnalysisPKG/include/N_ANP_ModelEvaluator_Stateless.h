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
// Filename       : $RCSfile: N_ANP_ModelEvaluator_Stateless.h,v $
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
// Revision Number: $Revision: 1.7 $
//
// Revision Date  : $Date: 2014/02/24 23:49:12 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef N_ANP_MODEL_EVALUATOR_STATELESS_H
#define N_ANP_MODEL_EVALUATOR_STATELESS_H

#include <EpetraExt_ModelEvaluator.h>
#include <N_ANP_ModelEvaluator.h>

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

// This defines the class ModelEvaluator_Stateless derived from
// EpetraExt:ModelEvaluator
class ModelEvaluator_Stateless : public EpetraExt::ModelEvaluator {

  public:
    // (Default) Constructor of the class  
    ModelEvaluator_Stateless();

    // Destructor of the class 
    virtual ~ModelEvaluator_Stateless();

    void set_XyceModelEvaluator(const RCP<Analysis::ModelEvaluator>& xyceME); 

    /** \name Overridden from EpetraExt::ModelEvaluator . */
    //@{
    Teuchos::RCP<const Epetra_Map> get_x_map() const;
    Teuchos::RCP<const Epetra_Map> get_f_map() const;
    Teuchos::RCP<const Epetra_Map> get_p_map(int p) const;
    Teuchos::RCP<const Epetra_Map> get_g_map(int p) const;
    Teuchos::RCP<Epetra_Operator> create_W() const;
    EpetraExt::ModelEvaluator::InArgs createInArgs() const;
    EpetraExt::ModelEvaluator::OutArgs createOutArgs() const;
    void evalModel( const InArgs& inArgs, const OutArgs& outArgs ) const;
    //@}

  private:
    // Private member functions:
    void setupInOutArgs_();

    // Private data:
    bool isInitialized_;

    //RCP<N_CIR_Xyce> xycePtr_;
    RCP<Analysis::ModelEvaluator> xyceME_;

    EpetraExt::ModelEvaluator::InArgs inArgs_;
    EpetraExt::ModelEvaluator::OutArgs outArgs_;
    int Np_; // 0
    int Ng_; // 0
    RCP<Epetra_Vector> tempStateVector_;
    RCP<Epetra_Vector> tempStateDotVector_;
    RCP<Epetra_Vector> tempVoltLimFVector_;
    RCP<Epetra_Vector> tempVoltLimQVector_;
    RCP<Epetra_Vector> tempFVector_;
    RCP<Epetra_Operator> tempWOperator_;
};

// Nonmember constructor
RCP<ModelEvaluator_Stateless> N_ANP_modelEvaluator_Stateless();

// Nonmember constructor
RCP<ModelEvaluator_Stateless> N_ANP_modelEvaluator_Stateless(
    const RCP<N_ANP_ModelEvaluator>& xyceME);

} // namespace Analysis
} // namespace Xyce

typedef Xyce::Analysis::ModelEvaluator_Stateless N_ANP_ModelEvaluator_Stateless;

#endif // N_ANP_MODEL_EVALUATOR_STATELESS_H
