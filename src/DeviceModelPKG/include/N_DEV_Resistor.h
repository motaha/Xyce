//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2013  Sandia Corporation
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
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
// Filename       : $RCSfile: N_DEV_Resistor.h,v $
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
// Revision Number: $Revision: 1.92.2.4 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-----------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Resistor_h
#define Xyce_N_DEV_Resistor_h

#include <N_DEV_fwd.h>
#include <N_DEV_DeviceInstance.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceTemplate.h>

namespace Xyce {
namespace Device {
namespace Resistor {

class Model;
class Instance;

/**
 * Resistor device instance.
 *
 * An instance is created for each occurance of the device in the netlist.
 *
 * It contains "unique" resistor information - ie stuff that will be true of only one resistor in the circuit, such as
 * the nodes to which it is connected.  A resistor is connected to only two circuit nodes.
 *
 * This class does not directly contain information about its node indices. It contains indices into the 3 parts (A, dx,
 * and b) of the matrix problem A*dx = b, and also x.  A is the Jacobian matrix, dx is the update to the solution
 * std::vector x, and b is the right hand side function std::vector.  These indices are global, and determined by
 * topology during the initialization stage of execution.
 *
 */
class Instance : public DeviceInstance
{
    friend class ParametricData<Instance>;            ///< Allow ParametricData to changes member values
    friend class Model;                               ///< Don't force a lot of pointless getters
    friend class Master;                              ///< Don't force a lot of pointless getters

  public:
    static ParametricData<Instance> &getParametricData(); ///< Singleton holding configuration and parameters

    Instance(
      InstanceBlock &   instance_block,               ///< Instance information from parser
      Model &           model,                        ///< Resistor model to add this instance to
      MatrixLoadData &  matrix_load_data,             ///< Solution matrix load data
      SolverState &     solver_state,                 ///< Solution sover state
      ExternData &      extern_data,                  ///< Solution external data
      DeviceOptions &   device_options);              ///< Device options defined in netlist

    /** 
     * Destroys this instance 
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Mon Aug 12 08:36:23 2013
     */
    ~Instance()
    {}

  private:
    Instance(const Instance &);
    Instance &operator=(const Instance &);

  public:
    /** 
     * Gets the parametric data for the resistor instance class.
     *
     * @return const reference to the resistor instance ParametricData singleton 
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Mon Aug 12 08:36:29 2013
     */
    virtual const ParametricData<void> &getMyParametricData() const {
      return getParametricData();
    }

    /** 
     * Gets the resistor model that owns this instance. 
     *
     * @return reference to the owning Resistor::Model
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Mon Aug 12 08:36:37 2013
     */
    Model &getModel() {
      return model_;
    }

    virtual void registerLIDs(const std::vector<int> & intLIDVecRef, const std::vector<int> & extLIDVecRef) /* override */;
    virtual void registerStateLIDs(const std::vector<int> & staLIDVecRef) /* override */;
    virtual void registerStoreLIDs(const std::vector<int> & stoLIDVecRef) /* override */;
    virtual void registerJacLIDs(const std::vector< std::vector<int> > & jacLIDVec) /* override */;
    virtual std::map<int, std::string> &getStoreNameMap() /* override */;

    virtual bool processParams(string param = "") /* override */;
    virtual bool updateTemperature(const double & temp_tmp) /* override */;
    virtual bool updateIntermediateVars() /* override */;
    virtual bool updatePrimaryState() /* override */;

    /** 
     * Return Jacobian stamp that informs topology of the layout of the resistor jacobian. 
     *
     * The Jacobian stamp describes the shape of the Jacobian to the Topology subsystem.  The Topology subsystem, in
     * turn, returns the offsets into the matrix and solution vectors where this instance data is located.
     *
     * @return const reference to a std::vector of std::vector of integers describing Jacobian stamp shape
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Mon Aug 12 08:36:45 2013
     */
    virtual const std::vector< std::vector<int> > &jacobianStamp() const  /* override */ {
      return jacStamp;
    }

    virtual bool loadDAEFVector() /* override */;
    virtual bool loadDAEdFdx() /* override */;

    /** 
     * Load residual vector
     *
     * @return true on success
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Mon Aug 12 08:36:51 2013
     */
    virtual bool loadDAEQVector() /* override */ {
      return true;
    }

    /** 
     * Load Jacobian
     *
     * @return true on success
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Mon Aug 12 08:36:58 2013
     */
    virtual bool loadDAEdQdx() /* override */ {
      return true;
    }

    virtual void setupPointers() /* override */;

  private:
    static std::vector< std::vector<int> >  jacStamp; ///< All Resistors have a common Jacobian Stamp

    Model &     model_;                 ///< Owning model

    // User-specified parameters:
    double      R;                      ///< Resistance (ohms)

    // These are for the semiconductor resistor
    double      length;                 ///< Resistor length.
    double      width;                  ///< Resistor width.
    double      temp;                   ///< Temperature of this instance

    // Temperature dependence parameters, these can override values specified in the model
    double      tempCoeff1;             ///< First order temperature coeff.
    double      tempCoeff2;             ///< Second order temperature coeff.
    double      dtemp;                  ///< Externally specified device temperature.
                                        ///<   NOT used, only here for compatibility in parsing
                                        ///<   netlist from simulators that support it

    // Flags used to tell if the user has specified one of these values on the command line.
    bool        tempCoeff1Given;        ///< First order temperation value was given in netlist
    bool        tempCoeff2Given;        ///< Second order temperature coeff was given in netlist
    bool        dtempGiven;             ///< Externally specified device temperature was given in netlist

    // Derived parameters:
    double      G;                      ///< Conductance(1.0/ohms)
    double      i0;                     ///< Current(ohms)

    int         li_Pos;                 ///< Index for Positive Node
    int         li_Neg;                 ///< Index for Negative Node
    int         li_store_dev_i;         ///< Index for Lead Current

    // Offset variables corresponding to the above declared indices.
    int         APosEquPosNodeOffset;   ///< Column index into force matrix of Pos/Pos conductance
    int         APosEquNegNodeOffset;   ///< Column index into force matrix of Pos/Neg conductance
    int         ANegEquPosNodeOffset;   ///< Column index into force matrix of Neg/Pos conductance
    int         ANegEquNegNodeOffset;   ///< Column index into force matrix of Neg/Neg conductance

#ifndef Xyce_NONPOINTER_MATRIX_LOAD
    // Pointers for Jacobian
    double *    f_PosEquPosNodePtr;
    double *    f_PosEquNegNodePtr;
    double *    f_NegEquPosNodePtr;
    double *    f_NegEquNegNodePtr;
#endif
};


/**
 * Resistor model
 *
 */
class Model : public DeviceModel
{
    friend class ParametricData<Model>;               ///< Allow ParametricData to changes member values
    friend class Instance;                            ///< Don't force a lot of pointless getters
    friend class Master;                              ///< Don't force a lot of pointless getters

  public:
    typedef std::vector<Instance *> InstanceVector;

    static ParametricData<Model> &getParametricData();

    Model(const ModelBlock &model_block, SolverState &solver_state, DeviceOptions &device_options);
    ~Model();

  private:
    Model();
    Model(const Model &);
    Model &operator=(const Model &);

  public:
    /** 
     * Gets the parametric data for the resistor model class.
     *
     * @return const reference to the resistor model ParametricData singleton 
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Mon Aug 12 08:36:29 2013
     */
    virtual const ParametricData<void> &getMyParametricData() const {
      return getParametricData();
    }

    /** 
     * Get the instance vector for all resistors owned by this model.
     *
     * @return reference to InstanceVector containing all resistors owned by this model
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Mon Aug 12 09:10:00 2013
     */
    InstanceVector &getInstanceVector() {
      return instanceContainer;
    }

    /** 
     * Get the instance vector for all resistors owned by this model.
     *
     * @return const reference to InstanceVector containing all resistors owned by this model
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Mon Aug 12 09:10:00 2013
     */
    const InstanceVector &getInstanceVector() const {
      return instanceContainer;
    }

    virtual std::ostream &printOutInstances(std::ostream &os) const;

    virtual bool processParams(string param = "") /* override */;
    virtual bool processInstanceParams(string param = "") /* override */;

  private:
    InstanceVector      instanceContainer;            ///< List of owned resistor intances

  private:
    // Semiconductor resistor parameters
    double      tempCoeff1;     ///< First order temperature coefficient
    double      tempCoeff2;     ///< Second order temperature coefficient
    double      sheetRes;       ///< Sheet resistance
    double      defWidth;       ///< Default width
    double      narrow;         ///< Narrowing due to side etching
    double      tnom;           ///< Parameter measurement temperature
};


/**
 * Resistor master
 *
 * 
 *
 */
class Master : public DeviceTemplate<Model, Instance>
{
    friend class Instance;                            ///< Don't force a lot of pointless getters
    friend class Model;                               ///< Don't force a lot of pointless getters

  public:
     /**
      * Construct a Resistor Device.
      *
      * @param device_name
      * @param class_name
      * @param default_model_name
      * @param model_type
      * @param linear_device
      * @param solver_state
      * @param device_options
      */
    Master(
      const std::string & device_name,
      const std::string & class_name,
      const std::string & default_model_name,
      LinearDevice        linear_device,
      SolverState &       solver_state,
      DeviceOptions &     device_options)
      : DeviceTemplate<Model, Instance>(device_name, class_name, default_model_name, linear_device, solver_state, device_options)
    {}

    virtual bool updateState(double * solVec, double * staVec, double * stoVec) /* override */;
    virtual bool loadDAEVectors (double * solVec, double * fVec, double * qVec, double * storeLeadF, double * storeLeadQ) /* override */;
    virtual bool loadDAEMatrices(N_LAS_Matrix & dFdx, N_LAS_Matrix & dQdx) /* override */;
};

} // namespace Resistor
} // namespace Device
} // namespace Xyce

typedef Xyce::Device::Resistor::Instance N_DEV_ResistorInstance;
typedef Xyce::Device::Resistor::Model N_DEV_ResistorModel;
typedef Xyce::Device::Resistor::Master N_DEV_ResistorMaster;

#endif
