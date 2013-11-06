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
// Filename       : $RCSfile: N_DEV_ROM.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Heidi Thornquist, SNL
//
// Creation Date  : 12/11/09
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.23.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_ROM_h
#define Xyce_N_DEV_ROM_h

// ----------   Xyce Includes   ----------
#include <N_DEV_DeviceTemplate.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>

namespace Xyce {
namespace Device {
namespace ROM {

// ---------- Forward Declarations ----------
class Model;
class Instance;

//-----------------------------------------------------------------------------
// Class         : Instance
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
//-----------------------------------------------------------------------------

class Instance : public DeviceInstance
{
  friend class ParametricData<Instance>;
  friend class Model;
  friend class Master;

  public:
  static ParametricData<Instance> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

    Instance(InstanceBlock & IB,
			    Model & Citer,
                            MatrixLoadData & mlData1,
                            SolverState &ss1,
                            ExternData  &ed1,
                            DeviceOptions & do1);


    ~Instance();

private:
  Instance(const Instance &);
  Instance &operator=(const Instance &);

public:
    // Additional Public Declarations
    void registerLIDs( const vector<int> & intLIDVecRef,
                       const vector<int> & extLIDVecRef );
    void registerStateLIDs( const vector<int> & staLIDVecRef );

    map<int,string> & getIntNameMap ();

    const vector< vector<int> > & jacobianStamp() const;
    void registerJacLIDs( const vector< vector<int> > & jacLIDVec );

    bool processParams (string param = "");
    bool updateTemperature(const double & temp_tmp);

    bool updateIntermediateVars () { return true; };
    bool updatePrimaryState ();

    bool setIC ();

    bool loadDeviceMask ();

    // load functions, residual:
    bool loadDAEQVector ();
    bool loadDAEFVector ();

    // load functions, Jacobian:
    bool loadDAEdQdx ();
    bool loadDAEdFdx ();

    void setupPointers();

    void varTypes( vector<char> & varTypeVec );

  public:
      // iterator reference to the resistor model which owns this instance.
  // Getters and setters
  Model &getModel() {
    return model_;
  }

private:

  Model &       model_;         //< Owning model

      // Data Members for Class Attributes
      bool isCSparse;
      bool isGSparse;

      // User-specified parameters:
      bool maskROMVars;
      int usePortDesc;
      int numROMVars;
      string baseFileName;   // base file name for reduced-order model files
      vector<double> Chat;  // Reduced-order model V'*C*V
      vector<int> Chat_colIdx, Chat_rowPtr;  // Chat structures if stored in CSR format
      vector<double> Ghat;  // Reduced-order model V'*G*V
      vector<int> Ghat_colIdx, Ghat_rowPtr;  // Ghat structures if stored in CSR format
      vector<int> CG_colIdx, CG_rowPtr;  // Union of Chat and Ghat maps stored in CSR format
      vector<double> Bhat;  // Reduced-order model V'*B
      vector<double> Lhat;  // Reduced-order model L'*V
      vector<double> Qhat;  // Workspace Qhat = Chat * xhat
      vector<double> Fhat;  // Workspace Fhat = [Iq - Lhat'* xhat; Ghat*xhat - Bhat*up]
      vector<double> i_ip;  // Storage for Iq

      // Two-level stamps (BNB)
      vector<double> Jstamp; // stamp for Jacobian
      vector<double> Fstamp; // stamp for F
      vector<double> G2;     // intermediate variable
      vector<double> C2;     // intermediate varaible
      vector<double> A2;     // intermediate varaible
      vector<double> A2last;
      vector<double> G2p;    // intermediate varaible
      vector<double> Gp2;    // intermediate varaible
      vector<double> A2sol;  // intermediate variable
      double dt, dt_last, alph, alph_last, coef, coefLast;
      double currentOrder, usedOrder;
      int lastTimeStepNumber;
      vector<int> ipiv_A2;  // for LAPACK math

      //local id's (offsets)
      vector<int> li_ROM;  // Interior variables
      vector<int> li_state; // Internal state

      // Offsets for Jacobian
      vector<int> AEqu_up_NodeOffset;
      vector<int> AEqu_ip_NodeOffset;
      vector< vector<int> > AEqu_NodeOffset;
      vector<int> ROMEqu_Lt_NodeOffset;
      vector<int> ROMEqu_B_NodeOffset;
      vector<int> ROMEqu_GpC_NodeOffset;
      // Offsets for sparse C and C in Jacobian
      vector<int> ROMEqu_C_NodeOffset;
      vector<int> ROMEqu_G_NodeOffset;

      // Pointers for Jacobian
      vector<double *> fEqu_up_NodePtr;
      vector<double *> fEqu_ip_NodePtr;
      vector<double *> fEqu_un_NodePtr; // BNB

      vector<double *> qROMEqu_Chat_VarsPtrs;
      vector<double *> fROMEqu_Ghat_VarsPtrs;
      vector<double *> fROMEqu_Lhat_VarsPtrs;
      vector<double *> fROMEqu_Bhat_VarsPtrs;

      vector< vector<int> > jacStamp;
};

//-----------------------------------------------------------------------------
// Class         : Model
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
//-----------------------------------------------------------------------------
class Model  : public DeviceModel
{
  typedef std::vector<Instance *> InstanceVector;

  friend class ParametricData<Model>;
  friend class Instance;
  friend class Master;

  public:
  static ParametricData<Model> &getParametricData();

  virtual const ParametricData<void> &getMyParametricData() const {
    return getParametricData();
  }

      Model    (const ModelBlock & MB,
                                     SolverState & ss1,
                                     DeviceOptions & do1);
      ~Model   ();

private:
  Model();
  Model(const Model &);
  Model &operator=(const Model &);

public:
  virtual std::ostream &printOutInstances(std::ostream &os) const;
      bool processParams (string param = "");
      bool processInstanceParams (string param = "");

public:
  InstanceVector &getInstanceVector() {
    return instanceContainer;
  }

  const InstanceVector &getInstanceVector() const {
    return instanceContainer;
  }

private:
    vector<Instance*> instanceContainer;

private:

};

//-----------------------------------------------------------------------------
// Class         : Master
// Purpose       :
// Special Notes :
// Creator       : Heidi Thornquist, SNL, Parallel Computational Sciences
// Creation Date : 12/11/09
//-----------------------------------------------------------------------------
class Master : public Xyce::Device::DeviceTemplate<Model, Instance>
{
  public:
    Master (
      const std::string &dn,
      const std::string &cn,
      const std::string &dmName,
           LinearDevice linearDev,
           SolverState & ss1,
           DeviceOptions & do1)
      : Xyce::Device::DeviceTemplate<Model, Instance>(
           dn, cn, dmName, linearDev, ss1, do1)
    {

    }

    virtual bool updateState (double * solVec, double * staVec, double * stoVec);

    void printMatrix (std::string vname, double * Matrix, int Nrows, int Ncols); // BNB

    // load functions:
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ);
    virtual bool loadDAEMatrices (N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx);

    friend class Instance;
    friend class Model;
};

} // namespace ROM
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::ROM::Instance N_DEV_ROMInstance;
typedef Xyce::Device::ROM::Model N_DEV_ROMModel;
typedef Xyce::Device::ROM::Master N_DEV_ROMMaster;

#endif // Xyce__h

