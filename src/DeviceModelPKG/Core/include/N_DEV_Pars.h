//-----------------------------------------------------------------------------
// Copyright Notice
//
//   Copyright 2002 Sandia Corporation. Under the terms
//   of Contract DE-AC04-94AL85000 with Sandia Corporation, the U.S.
//   Government retains certain rights in this software.
//
//    Xyce(TM) Parallel Electrical Simulator
//    Copyright(C) 2002-2011  Sandia Corporation
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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_Pars.h,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Netlist device and model parameter management
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.8.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:37 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_DEV_Pars_h
#define Xyce_N_DEV_Pars_h

#include <algorithm>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <N_DEV_Configuration.h>
#include <N_DEV_Const.h>
#include <N_DEV_Param.h>
#include <N_DEV_Units.h>
#include <N_DEV_fwd.h>
#include <N_ERH_ErrorMgr.h>
#include <N_UTL_NoCase.h>

/**
 * ParameterType::ExprAccess is enumeration of parameter usage types and masks
 */
namespace ParameterType {
  enum ExprAccess
  {
    NO_DEP    = 0x0,              ///< Parameter can only be set to a constant from netlist
    TIME_DEP  = 0x1,              ///< Parameter may be specified as time dependent expression from netlist
    SOLN_DEP  = 0x2,              ///< Parameter may be specified as a solution dependent expression from netlist
    NO_INPUT  = 0x4,              ///< Parameter cannot be initialized from netlist // THIS IS SCHEDULED FOR TERMINATION
    LOG_T_DEP = 0x8,              ///< Parameter uses temperature interpolation based on log of value
    MIN_RES   = 0x10,             ///< Parameter is subject to being set to minimum lead resistance
    MIN_CAP   = 0x20              ///< Parameter is subject to being set to minimum junction capacitance
  };
};

using namespace ParameterType;


namespace Xyce {
namespace Device {

/**
 * @addtogroup xyce_device_parameters_detail
 *
 * @brief
 *
 * Device models, device instances and composite parameters are collectively called an
 * <b>entity</b>, which are all derived from the ParameterBase tag class, and manage the mapping
 * from a parameter's string name to a member variable in each created entity object.  The mapping
 * is refered to as parameter binding and the object's bound member variable is refered to as the
 * parameter member variable.  The object may have an additional member variable for each binding
 * that indicates if the value was specified in the netlist device and is known as the given member
 * variable.
 *
 * Each class of entity has a ParametricData<T> class, where T is the entity class, and a singleton
 * object of that class.  This singleton manages the parameter binding.  A parameter descriptor,
 * Descriptor, has the parameter member variable pointer for entity object of type T that sets and
 * retrieves the current value of an entity's member variable using the parameters string name.
 *
 * The Descriptor maintains other information about the parameter.  It has a serial number for each
 * parameter and a boolean member variable pointer that determines if the parameter has been given
 * in the netlist.  The given value map, GivenValueMap, is a mapping from entity object pointer and
 * serial number to a boolean given member varaible indicating if the parameter was provided by the
 * netlist.
 *
 * Any entity may also wish to restore a parameter member variable value to the value originally set
 * during initialization.  The Descriptor maintains an index into an original value map to store and
 * retrieve the parameters original value.  The original value map, OriginalValueMap, is a mapping
 * from entity object pointer and original value index to the parameter's origin value.  If the
 * original value index is -1, then no original value information is maintained.
 *
 * A parameter may also be vectorized, allowing it to maintain a vector values.  The Descriptor
 * maintains that vector length.
 *
 * A parameter has specific usages as described by the ParameterType::ExprAccess enumeration.
 *
 * A parameter may also serves as an aggregation of parametric values which is stored as a composite
 * parameter, CompositeParam, entity.
 *
 * And a parameter has units, a documentation category and a documentation description.
 *
 * @todo Replace the s_originalValueMap and s_givenValueMap with a non-static implementation.
 *
 * Entity - DeviceModel, DeviceInstance or CompositeParam
 *
 * Given - flag and possibly and entity member variable that indicates if a parameter was given in
 * the netlist.
 *
 * OriginalValue - parameters value after initialization
 *
 */

typedef std::map<std::string, CompositeParam *> CompositeMap;

/**
 * Tag for all classes that manage parameters AKA entities
 */
struct ParameterBase
{};

template<class C>
class ParametricData;

template<class T>
class Entry;

template <class T>
inline std::ostream &printEntry(std::ostream &os, const Entry<T> &entry);

template <class T>
inline std::ostream &printEntry(std::ostream &os, const Entry<std::vector<T> > &entry);

template <>
inline std::ostream &printEntry(std::ostream &os, const Entry<std::string> &entry);

template <>
inline std::ostream &printEntry(std::ostream &os, const Entry<bool> &entry);

template <>
inline std::ostream &printEntry(std::ostream &os, const Entry<CompositeMap> &entry);

void typeMismatch(const std::type_info &from_type, const std::type_info &to_type);

/**
 * Class Entry<void> defines the parameter binding value entry interface.
 *
 * This defines the interface to check the data type of an Entry<T> object and to print the value.
 * Type specific Entry classes inherit from this class.  The entry_cast<T>() function is used to
 * cast an object of this type to the derived Entry<T> class safely.
 *
 */
template<>
class Entry<void>
{
  protected:
    /**
     * Constructs the Entry base class
     *
     * @note The construct is protected so it may only be constructed by Entry<T> classes.
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 08:32:36 2013
     */
    Entry()
    {}

  public:
    /**
     * Destroys Entry
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 08:34:16 2013
     */
    virtual ~Entry()
    {}

  private:
    Entry(const Entry &);
    Entry &operator=(const Entry &);

  public:
    /**
     * Returns the type_info of the data type being stored in the entry.
     *
     * @return const reference to the type_info of the data type being stored in the entry.
     */
    virtual const std::type_info &type() const = 0;

    /**
     * Prints the value of the entry to the output stream
     *
     * @param os output stream to write to
     *
     * @return reference to the output stream
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:28:02 2013
     */
    std::ostream &print(std::ostream &os) const {
      doPrint(os);

      return os;
    }

  private:
    /**
     * Prints the value of the entry to the output stream
     *
     * @param os output stream to write to
     *
     * @return reference to the output stream
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:29:29 2013
     */
    virtual std::ostream &doPrint(std::ostream &os) const = 0;
};

/**
 * Class Entry<T> defines the parameter member variable access for parameter member variable of type T
 *
 * The pointer to the parameter member variable and the default value to set are contained in this
 * class.
 */
template<class T>
class Entry : public Entry<void>
{
  public:
    /**
     * Constructs an Entry.
     *
     * Initializes the parameter member variable pointer of type T.  Note the member variable
     * pointer are actually offsets into an object where the data for the member exists.  So, by
     * storing this pointer the value of any object of type can be retrieved.
     *
     * The default value ot type T is also contained here.
     *
     * @param member_value_ptr parameter member variable pointer
     * @param default_value Default value of parameter
     */
    Entry(T ParameterBase::*member_value_ptr, const T &default_value = T())
      : memberValuePtr_(member_value_ptr),
        defaultValue_(default_value)
    {}

    /**
     * Destroys the entry
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:30:29 2013
     */
    virtual ~Entry()
    {}

  private:
    Entry(const Entry &);
    Entry &operator=(const Entry &);

  public:
    /**
     * Returns type_info of this entry.
     *
     * @return const reference to type_info of this entry
     */
    virtual const std::type_info &type() const {
      return typeid(T);
    }

  private:
    /**
     * Prints the value of the entry to the output stream
     *
     * @param os output stream to write to
     *
     * @return reference to the output stream
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:29:29 2013
     */
    virtual std::ostream &doPrint(std::ostream &os) const {
      return printEntry(os, *this);
    }

  public:
    /**
     * Return the member pointer to the data member variable that holds the value associated with
     * this parameter
     *
     * @return member pointer to the data member variable
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:31:21 2013
     */
    T ParameterBase::*getMemberPtr() const {
      return memberValuePtr_;
    }

    /**
     * Return the value of the entity's parameter
     *
     * @param entity device class or device instance
     *
     * @return const reference to the value of the entity's parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:32:58 2013
     */
    const T &getValue(const ParameterBase &entity) const {
      return entity.*memberValuePtr_;
    }

    /**
     * Return the value of the entity's parameter
     *
     * @param entity device class or device instance
     *
     * @return reference to the value of the entity's parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:32:58 2013
     */
    T &getValue(ParameterBase &entity) const {
      return entity.*memberValuePtr_;
    }

    /**
     * Sets the value of the entity's parameter
     *
     * @param entity device class or device instance
     * @param value value to set the entity's parameter to
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:32:58 2013
     */
    void setValue(ParameterBase &entity, const T &value) const {
      entity.*memberValuePtr_ = value;
    }

    /**
     * Return the default value of the parameter
     *
     * All parameters provide a default value when created.  The parameter is set to this value and
     * the given flag is cleared on construction.  If the value is provided by the net list, the
     * parameter's value is set accordingly and given flag is set to true.
     *
     * @return const reference to the parameter's default value
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:35:18 2013
     */
    const T &getDefaultValue() const {
      return defaultValue_;
    }

    /**
     * Sets the parameter's default value
     *
     * All parameters provide a default value when created.  The parameter is set to this value and
     * the given flag is cleared on construction.  If the value is provided by the net list, the
     * parameter's value is set accordingly and given flag is set to true.
     *
     * @param value default value of the parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:37:50 2013
     */
    void setDefaultValue(const T &value) {
      defaultValue_ = value;
    }

  private:
    T                           defaultValue_;        ///< Default value of parameter
    T ParameterBase::*          memberValuePtr_;      ///< Member pointer containing value
};

/**
 * Casts the entry to type T
 *
 * If the entry if not of the specified type, typeMismatch() is calls to emit a fatal error.
 *
 * @param T type to cast to
 * @param entry entry to cast
 *
 * @return const reference to Entry<T>
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Thu Jul 25 13:21:09 2013
 */
template<class T>
const Entry<T> &entry_cast(const Entry<void> &entry) {
  if (entry.type() != typeid(T))
    typeMismatch(entry.type(), typeid(T));

  return static_cast<const Entry<T> &>(entry);
}

/**
 * Casts the entry to type T
 *
 * If the entry if not of the specified type, typeMismatch() is calls to emit a fatal error.
 *
 * @param T type to cast to
 * @param entry entry to cast
 *
 * @return reference to Entry<T>
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Thu Jul 25 13:21:09 2013
 */
template<class T>
Entry<T> &entry_cast(Entry<void> &entry) {
  if (entry.type() != typeid(T))
    typeMismatch(entry.type(), typeid(T));

  return static_cast<Entry<T> &>(entry);
}

/**
 * Class Descriptor describes the parameters stored in the ParametricData parameter map.
 *
 * The descriptor contains the Entry for the parameter member variable.  Each parameter is assigned
 * a serialNumber_ when it is created and is used by manage the given value map.  The parameter may
 * have a given member variable as well.  The parameter can be marked to store its original
 * initialized value and the originalValueIndex_ is used to store this value in the s_originalValueMap.
 *
 * If the parameter is vectorized, the index to this parameter is in vec_.
 *
 * The expressionAccess_ is used to indicate the parameter's usage.
 *
 * The Descriptor also contains the units, catagory and description for documentation generation.
 * serial number, its original value index if the initiazes value needs to be restored, its usage
 * (ExprAccess), units,
 */
class Descriptor
{
  public:
    /**
     * Constructs Descriptor
     *
     * @param entry the parameter member variable pointer, type and default value
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:38:43 2013
     */
    Descriptor(Entry<void> *const entry)
      : serialNumber_(0),
        originalValueIndex_(-1),
        vec_(0),
        expressionAccess_(ParameterType::NO_DEP),
        entry_(entry),
        unit_(U_NONE),
        category_(CAT_NONE),
        description_(""),
        compositeParametricData_(0),
        given_(static_cast<bool ParameterBase::*>(0))
    {}

    /**
     * Destroy Descriptor
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 15:42:55 2013
     */
    virtual ~Descriptor() {
      delete entry_;
    }

  private:
    Descriptor(const Descriptor &descriptor);
    Descriptor &operator=(const Descriptor &descriptor);

  public:
    /**
     * Tests entry data type
     *
     * @return true if entry is of type T
     *
     * @date   Wed Aug  7 10:31:04 2013
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     */
    template <class T>
    bool isType() const {
      if (entry_)
        return entry_->type() == typeid(T);

      return typeid(T) == typeid(int);
    }

    /**
     * Sets the original value index used to store and retrieve values from the OriginalValueMap
     *
     * @param index index of the original value in the OriginalValueMap
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:37:03 2013
     */
    void setOriginalValueIndex(int index) {
      originalValueIndex_ = index;
    }

    /**
     * Returns the index of the origina value in the OriginalValueMap
     *
     * @return original value index
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:19:59 2013
     */
    int getOriginalValueIndex() const {
      return originalValueIndex_;
    }

    /**
     * Sets the serial number used to store and retrieve given boolean from the GivenValueMap
     *
     * @param serial_number serial number of the parameter in the GivenValueMap
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:38:18 2013
     */
    void setSerialNumber(int serial_number) {
      serialNumber_ = serial_number;
    }

    /**
     * Gets the serial number used to store and retireve given boolean fromt he GivenValueMap
     *
     * @return parameter serial number
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:40:30 2013
     */
    int getSerialNumber() const {
      return serialNumber_;
    }

    /**
     * Sets the expression access which describe the usage of the parameter
     *
     * @param expression_access usage of the paraemeter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:41:56 2013
     */
    void setExpressionAccess(ExprAccess expression_access) {
      expressionAccess_ = expression_access;
    }

    /**
     * Gets the expression access which describes the usage of the paramter
     *
     * @return usage of the parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:41:26 2013
     */
    ExprAccess getExpressionAccess() const {
      return expressionAccess_;
    }

    /**
     * Sets the units of the parameter, only used to document the parameter
     *
     * @param unit units of the parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:43:32 2013
     */
    void setUnit(ParameterUnit unit) {
      unit_ = unit;
    }

    /**
     * Gets the units of the parameter, only used to document the parameter
     *
     * @return units of the parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:44:16 2013
     */
    ParameterUnit getUnit() const {
      return unit_;
    }

    /**
     * Sets the category of the parameter, only used to document the parameter
     *
     * @param category category of the parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:45:09 2013
     */
    void setCategory(ParameterCategory category) {
      category_ = category;
    }

    /**
     * Gets the category of the parameter, only used to document the parameter
     *
     * @return categort of the parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:45:40 2013
     */
    ParameterCategory getCategory() const {
      return category_;
    }

    /**
     * Sets the description of the parameter, only used to document the parameter
     *
     * @param description description of the parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:46:21 2013
     */
    void setDescription(const std::string &description) {
      description_ = description;
    }

    /**
     * Gets the description of the parameter, only used to document the parameter
     *
     * @return description of the parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:46:46 2013
     */
    const std::string &getDescription() const {
      return description_;
    }

    /**
     * Sets the composite parametric
     *
     * A composite parameter is a named aggregation of parameters
     *
     * @param composite_parametric_data composite parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:47:19 2013
     */
    void setCompositeParametricData(const ParametricData<void> *composite_parametric_data) {
      compositeParametricData_ = composite_parametric_data;
    }

    /**
     * Return the composite parameter
     *
     * A composite parameter is a named aggregation of parameters
     *
     * @return const pointer to the composite parameter
     *
     * @note I think the template parameter is really only ever CompositeParam, though sometimes its
     * currently void.
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:09:01 2013
     */
    template <class U>
    const ParametricData<U> *getCompositeParametricData() const {
      return static_cast<const ParametricData<U> *>(compositeParametricData_);
    }

    /**
     * sets the vector length of an vectorized parameter
     *
     * @param index
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:50:05 2013
     */
    void setVec(int index) {
      vec_ = index;
    }

    /**
     * Gets the vector length of a vectorized parameter
     *
     * @return length of the vectorized parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 14:51:10 2013
     */
    int getVec() const {
      return vec_;
    }

    /**
     * Gets the entry object of the parameter
     *
     * @return const reference to the Entry<void>
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 15:29:16 2013
     */
    const Entry<void> &getEntry() const {
      return *entry_;
    }

    /**
     * Gets the entry object of the parameter
     *
     * @return reference to the Entry<void>
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 15:29:16 2013
     */
    Entry<void> &getEntry() {
      return *entry_;
    }

public:
    /**
     * Returns the value of the parameter for the entity
     *
     * @param entity device class or device instance
     *
     * @return reference to the value of the parameter for the entity
     *
     * @date   Wed Aug  7 10:36:10 2013
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     */
    template <class T>
    const T &value(const ParameterBase &entity) const {
      const Entry<T> &entry = entry_cast<T>(*entry_);

      return entry.getValue(entity);
    }

    /**
     * Returns the value of the parameter for the entity
     *
     * @param entity device class or device instance
     *
     * @return reference to the value of the parameter for the entity
     *
     * @date   Wed Aug  7 10:36:10 2013
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     */
    template <class T>
    T &value(ParameterBase &entity) const {
      const Entry<T> &entry = entry_cast<T>(*entry_);

      return entry.getValue(entity);
    }

    /**
     * Returns the parameter member variable pointer of the enrtry
     *
     * @return parameter member pointer of type T
     *
     * @date   Wed Aug  7 10:54:58 2013
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     */
    template <class T>
    T ParameterBase::*getMemberPtr() const {
      const Entry<T> &entry = entry_cast<T>(*entry_);

      return entry.getMemberPtr();
    }

    /**
     * Tests if parameter has a given data member
     *
     * Parameters may provide a boolean member variable that is set true if the netlist provides
     * the value.
     *
     * @return true if the parameter has a boolean member variable to set if it is provided by the netlist
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:03:55 2013
     */
    bool hasGivenMember() const {
      return given_ != static_cast<bool ParameterBase::*>(0);
    }

    /**
     * Sets the boolean member variable to set if the netlist provides the value.
     *
     * Parameters may provide a boolean member variable that is set true if the netlist provides
     * the value.
     *
     * @param given boolean member variable to be set
     *
     * @date   Wed Aug  7 11:02:38 2013
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     */
    template<class U>
    void setGivenMember(bool U::*given) {
      given_ = static_cast<bool ParameterBase::*>(given);
    }

    /**
     * Tests if the parameter has been given by the netlist.
     *
     * @return true if the parameter was given by the netline
     *
     * @date   Wed Aug  7 11:01:42 2013
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     */
    bool getGiven(ParameterBase &entity) const {
      return entity.*given_;
    }

    /**
     * Sets the given state of the parameter to value
     *
     * The parameter's boolean member variable that is set true if exists and the netlist provides
     * the value.
     *
     * @param entity device class or device instance
     * @param value true if provided by netlist
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Wed Aug  7 11:09:00 2013
     */
    void setGiven(ParameterBase &entity, bool value) const {
      if (hasGivenMember()) {
        entity.*given_ = value;
      }
    }

  private:
    int                                 serialNumber_;          ///< Unique identifier of descriptor
    int                                 originalValueIndex_;    ///< Index where original value is store usually -1, meaning none saved.
    int                                 vec_;                   ///< If > 0 specifies a vector of params.(eg: if = 3 then IC becomes IC1, IC2, IC3)
    ExprAccess                          expressionAccess_;      ///< Flags for parameter attributes, such as whether can be input by user, may depend on time, etc.

    Entry<void> * const                 entry_;                 ///< Pointer to entry which contains the value

    ParameterUnit                       unit_;                  ///< Unit designator for documentation
    ParameterCategory                   category_;              ///< Category designator for documentation
    std::string                         description_;           ///< Description of parameter for documentation
    const ParametricData<void> *        compositeParametricData_;

    bool ParameterBase::*               given_;                 ///< Pointer to given bool, usually NULL.
};


/**
 * Prints the entry default value to the output stream
 *
 * @param os output stream
 * @param entry entry to print
 *
 * @return reference to the output stream
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Fri Aug  9 15:32:26 2013
 */
template <class T>
inline std::ostream &printEntry(std::ostream &os, const Entry<T> &entry) {
  os << entry.getDefaultValue();

  return os;
}

/**
 * Prints the entry default values of a vectorized parameter
 *
 * @param os output stream
 * @param entry entry to print
 *
 * @return reference to the output stream
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Fri Aug  9 15:33:24 2013
 */
template <class T>
inline std::ostream &printEntry(std::ostream &os, const Entry<std::vector<T> > &entry) {
  for (typename std::vector<T>::const_iterator it = entry.getDefaultValue().begin(); it != entry.getDefaultValue().end(); ++it)
    os << (*it) << std::endl;

  return os;
}

/**
 * Prints the entry default string value, within single quotes
 *
 * @param os output stream
 * @param entry entry to print
 *
 * @return reference to the output stream
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Fri Aug  9 15:34:18 2013
 */
template <>
inline std::ostream &printEntry(std::ostream &os, const Entry<std::string> &entry) {
  os << "'" << entry.getDefaultValue() << "'";

  return os;
}

/**
 * Prints the entry default boolean value, printed as true or false
 *
 * @param os output stream
 * @param entry entry to print
 *
 * @return reference to the output stream
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Fri Aug  9 15:34:18 2013
 */
template <>
inline std::ostream &printEntry(std::ostream &os, const Entry<bool> &entry) {
  os << (entry.getDefaultValue() ? "true" : "false");

  return os;
}

/**
 * Prints the entry composite value as newline terminated list of colon separated name, value pairs
 *
 * @param os output stream
 * @param entry entry to print
 *
 * @return reference to the output stream
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Fri Aug  9 15:34:18 2013
 */
template <>
inline std::ostream &printEntry(std::ostream &os, const Entry<CompositeMap> &entry) {
  for (CompositeMap::const_iterator it = entry.getDefaultValue().begin(); it != entry.getDefaultValue().end(); ++it)
    os << (*it).first << ": " << (*it).second << std::endl;

  return os;
}

/**
 * Gets the default value of the parameter
 *
 * @param descriptor descriptor of the parameter
 *
 * @return const reference to the default value of the entry of the descriptor
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Fri Aug  9 15:37:33 2013
 */
template <class T>
inline const T &getDefaultValue(const Descriptor &descriptor) {
  const Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  return entry.getDefaultValue();
}

/**
 * Sets the default value of the parameter
 *
 * @param descriptor descriptor of the parameter
 * @param t default value
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Fri Aug  9 15:38:49 2013
 */
template <class T>
inline void setDefaultValue(Descriptor &descriptor, const T &t) {
  Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  entry.setDefaultValue(t);
}

/**
 * Returns the value of the parameter for the entity
 *
 * @param entity device class or device instance
 * @param descriptor descriptor of the parameter
 *
 * @return const reference to the value of the entry of the descriptor
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Fri Aug  9 15:39:40 2013
 */
template <class T>
inline const T &value(const ParameterBase &entity, const Descriptor &descriptor) {
  const Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  return entry.getValue(entity);
}

/**
 * Returns the value of the parameter for the entity
 *
 * @param entity device class or device instance
 * @param descriptor descriptor of the parameter
 *
 * @return reference to the value of the entry of the descriptor
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Fri Aug  9 15:39:40 2013
 */
template <class T>
inline T &value(ParameterBase &entity, const Descriptor &descriptor) {
  const Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  return entry.getValue(entity);
}


/**
 * Gets the value of the parameter for the entity
 *
 * @param entity device class or device instance
 * @param descriptor descriptor of the parameter
 *
 * @return const reference to the value of the entry of the descriptor
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Fri Aug  9 15:39:40 2013
 */
template <class T, class U>
inline const T &getValue(const ParameterBase &entity, const Descriptor &descriptor) {
  const Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  return entry.getValue(entity);
}

/**
 * Sets the value of the parameter
 *
 * @param entity device class or device instance
 * @param descriptor descriptor of the parameter
 * @param value value
 *
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 * @date   Fri Aug  9 15:42:13 2013
 */
template <class T, class U>
inline void setValue(ParameterBase &entity, const Descriptor &descriptor, const T &value) {
  const Entry<T> &entry = entry_cast<T>(descriptor.getEntry());

  entry.setValue(entity, value);
}

/**
 * Class ParametricData<void> manages the configuration information and the parameter binding map
 *
 * Parametric data associated with a device instance, device model or composite parameter
 *
 * The Parametric data class manages the mapping of parameter string names to descriptors and the
 * general configuration information associated with a device model.
 *
 * To restore original values during perturbation, the originalValueCount_ and serialNumber_ members
 * maintain counts of original values to be stored and of parameters declared.
 *
 * @todo configTable_ needs to be moved elsewhere someday as it is only used for device instance
 * entities  It was added here temporarily to make it available without having to construct one ala
 * the parameter binding map.
 *
 * @date   Tue Aug  6 13:10:21 2013
 * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
 */
template<>
class ParametricData<void>
{
  public:
    typedef std::map<std::string, Descriptor *, LessNoCase> ParameterMap;

    /**
     * Constructs a ParametricData object
     *
     * @date   Tue Aug  6 13:10:21 2013
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     */
    ParametricData()
      : map_(),
        serialNumberCount_(0),
        originalValueCount_(0),
        configTable_()
    {}

    /**
     * Destroys a ParametricData object
     *
     * @date   Tue Aug  6 13:13:47 2013
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     */
    virtual ~ParametricData() {
      for (ParameterMap::iterator it = map_.begin(); it != map_.end(); ++it)
        delete (*it).second;
    }

  private:
    ParametricData(const ParametricData &parametric_data);
    ParametricData &operator=(const ParametricData &parametric_data);

  public:
    /** 
     * Gets the parameter binding map map
     *
     * @return reference to the parameter binding map
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 15:50:06 2013
     */
    ParameterMap &getMap() {
      return map_;
    }

    /** 
     * Gets the parameter binding map map
     *
     * @return const reference to the parameter binding map
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 15:50:06 2013
     */
    const ParameterMap &getMap() const {
      return map_;
    }

    /** 
     * Gets the configuation table.
     *
     * @return reference to the configuation table
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:11:37 2013
     */
    Configuration &getConfigTable() {
      return configTable_;
    }

    /** 
     * Gets the configuation table.
     *
     * @return const reference to the configuation table
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:11:37 2013
     */
    const Configuration &getConfigTable() const {
      return configTable_;
    }

    /** 
     * Sets the number of nodes in the configuation table 
     *
     * @param node_count number of nodes
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:13:14 2013
     */
    void setNumNodes(int node_count) {
      configTable_.numNodes = node_count;
    }

    /** 
     * Gets the number of nodes in the configuation table 
     *
     * @return number of nodes
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:13:14 2013
     */
    int getNumNodes() const {
      return configTable_.numNodes;
    }

    /** 
     * Sets the number of optional nodes in the configuation table 
     *
     * @param node_count number of optional nodes
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:13:14 2013
     */
    void setNumOptionalNodes(int node_count) {
      configTable_.numOptionalNodes = node_count;
    }

    /** 
     * Gets the number of optional nodes in the configuation table 
     *
     * @return number of optional nodes
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:13:14 2013
     */
    int getNumOptionalNodes() const {
      return configTable_.numOptionalNodes;
    }

    /** 
     * Sets the number of fill nodes in the configuation table 
     *
     * @param node_count number of fill nodes
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:13:14 2013
     */
    void setNumFillNodes(int node_count) {
      configTable_.numFillNodes = node_count;
    }

    /** 
     * Gets the number of fill nodes in the configuation table 
     *
     * @return number of fill nodes
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:13:14 2013
     */
    int getNumFillNodes() const {
      return configTable_.numFillNodes;
    }

    /** 
     * Sets the model required flag 
     *
     * @param model_required 1 if the model is required
     *
     * @todo This otta be a boolean
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:17:23 2013
     */
    void setModelRequired(int model_required) {
      configTable_.modelRequired = model_required;
    }

    /** 
     * Gets the model required falg 
     *
     * @return 1 if the model is required
     *
     * @todo This otta be a boolean
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:19:42 2013
     */
    int getModelRequired() const {
      return configTable_.modelRequired;
    }

    /** 
     * Sets the primary parameter
     *
     * This name is applied to the netlist as is specified is the first primary name is missing.
     *
     * @param parameter_name name of the primary parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:20:18 2013
     */
    void setPrimaryParameter(const std::string &parameter_name) {
      configTable_.primaryParameter = parameter_name;
    }

    /** 
     * Gets the primary parameter
     *
     * This name is applied to the netlist as is specified is the first primary name is missing.
     *
     * @return name of the primary parameter
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:20:18 2013
     */
    const std::string &getPrimaryParameter() const {
      return configTable_.primaryParameter;
    }

    /** 
     * Adds the model type to the list of model types of this device instance 
     *
     * @param model_type_name name of the model type
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:22:13 2013
     */
    void addModelType(const std::string &model_type_name) {
      std::vector<std::string>::iterator it = std::find_if(configTable_.modelTypes.begin(), configTable_.modelTypes.end(), EqualNoCaseOp(model_type_name));

      // Add if not found
      if (it == configTable_.modelTypes.end())
        configTable_.modelTypes.push_back(model_type_name);
    }

    /** 
     * Returns the list of model types
     *
     * @return const reference to the list of model type
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:23:04 2013
     */
    const std::vector<std::string> &getModelTypes() const {
      return configTable_.modelTypes;
    }

  protected:
    /** 
     * Adds the parameter to the parameter binding map 
     *
     * @param name parameter name
     * @param descriptor descriptor created for the parameter
     * @param parameter_data_class typeinfo to get the class name for diagnostics
     *
     * @author David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
     * @date   Fri Aug  9 16:23:45 2013
     */
    void addDescriptor(const std::string &name, Descriptor *descriptor, const std::type_info &parameter_data_class);

  protected:
    ParameterMap  map_;                 ///< Mapping from parameter name to descriptor
    int           serialNumberCount_;   ///< Count of parameters
    int           originalValueCount_;  ///< Count of original values
    Configuration configTable_;         ///< Device configuration information
};

void checkExprAccess(const std::string &name, ParameterType::ExprAccess &expr_access, const std::type_info &parameter_data_class);

/**
 * Manages parameter binding for class C.
 *
 * @note There is a bewildering set of addPar() functions to take in to account nearly any
 * conceivable combination of parameter.  This should either be refactored to allow a cascade of
 * setters or the parameter management change entirely.
 */
template<class C>
class ParametricData : public ParametricData<void>
{
  public:
    ParametricData()
      : ParametricData<void>()
    {}

    virtual ~ParametricData()
    {}

  private:
    ParametricData(const ParametricData &parametric_data);
    ParametricData &operator=(const ParametricData &parametric_data);

  public:
    template<class T, class U>
    void addPar(const char *parName, T def, bool orig, ParameterType::ExprAccess depend, T U::*varPtr, bool U::*givenPtr, ParameterUnit unit, ParameterCategory category, const char *description)
    {
      Descriptor *descriptor = new Descriptor(new Entry<T>(static_cast<T ParameterBase::*>(varPtr), def));

      checkExprAccess(parName, depend, typeid(C));

      descriptor->setExpressionAccess(depend);
      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);
      descriptor->setGivenMember(givenPtr);
      descriptor->setUnit(unit);
      descriptor->setCategory(category);
      descriptor->setDescription(description);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class T, class U>
    void addPar(const char *parName, T def, bool orig, ParameterType::ExprAccess depend, T U::*varPtr, void *noGivenPtr, ParameterUnit unit, ParameterCategory category, const char *description)
    {
      Descriptor *descriptor = new Descriptor(new Entry<T>(static_cast<T ParameterBase::*>(varPtr), def));

      checkExprAccess(parName, depend, typeid(C));

      descriptor->setExpressionAccess(depend);
      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);
      descriptor->setUnit(unit);
      descriptor->setCategory(category);
      descriptor->setDescription(description);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class T, class U>
    void addPar(const char *parName, T def, bool orig, ParameterType::ExprAccess depend, T U::*varPtr, ParameterUnit unit, ParameterCategory category, const char *description)
    {
      Descriptor *descriptor = new Descriptor(new Entry<T>(static_cast<T ParameterBase::*>(varPtr), def));

      checkExprAccess(parName, depend, typeid(C));

      descriptor->setExpressionAccess(depend);
      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);
      descriptor->setUnit(unit);
      descriptor->setCategory(category);
      descriptor->setDescription(description);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class U>
    void addPar(const char *parName, const char *def, bool orig, ParameterType::ExprAccess depend, std::string U::*varPtr, bool U::*givenPtr, ParameterUnit unit, ParameterCategory category, const char *description)
    {
      Descriptor *descriptor = new Descriptor(new Entry<std::string>(static_cast<std::string ParameterBase::*>(varPtr), def));

      checkExprAccess(parName, depend, typeid(C));

      descriptor->setExpressionAccess(depend);
      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);
      descriptor->setGivenMember(givenPtr);
      descriptor->setUnit(unit);
      descriptor->setCategory(category);
      descriptor->setDescription(description);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class U>
    void addPar(const char *parName, const char *def, bool orig, ParameterType::ExprAccess depend, std::string U::*varPtr, void *noGivenPtr, ParameterUnit unit, ParameterCategory category, const char *description)
    {
      Descriptor *descriptor = new Descriptor(new Entry<std::string>(static_cast<std::string ParameterBase::*>(varPtr), def));

      checkExprAccess(parName, depend, typeid(C));

      descriptor->setExpressionAccess(depend);
      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);
      descriptor->setUnit(unit);
      descriptor->setCategory(category);
      descriptor->setDescription(description);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class U>
    void addPar(const char *parName, const char *def, bool orig, ParameterType::ExprAccess depend, std::string U::*varPtr, ParameterUnit unit, ParameterCategory category, const char *description)
    {
      Descriptor *descriptor = new Descriptor(new Entry<std::string>(static_cast<std::string ParameterBase::*>(varPtr), def));

      checkExprAccess(parName, depend, typeid(C));

      descriptor->setExpressionAccess(depend);
      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);
      descriptor->setUnit(unit);
      descriptor->setCategory(category);
      descriptor->setDescription(description);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class T, class U>
    void addPar(const char *parName, const T &def, bool orig, ParameterType::ExprAccess depend, T U::*varPtr, bool U::*givenPtr)
    {
      Descriptor *descriptor = new Descriptor(new Entry<T>(static_cast<T ParameterBase::*>(varPtr), def));

      checkExprAccess(parName, depend, typeid(C));

      descriptor->setExpressionAccess(depend);
      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);
      descriptor->setGivenMember(givenPtr);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class T, class U>
    void addPar(const char *parName, const T &def, bool orig, ParameterType::ExprAccess depend, T U::*varPtr, void *noGivenPtr)
    {
      Descriptor *descriptor = new Descriptor(new Entry<T>(static_cast<T ParameterBase::*>(varPtr), def));

      checkExprAccess(parName, depend, typeid(C));

      descriptor->setExpressionAccess(depend);
      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class T, class U>
    void addPar(const char *parName, const T &def, bool orig, ParameterType::ExprAccess depend, T U::*varPtr)
    {
      Descriptor *descriptor = new Descriptor(new Entry<T>(static_cast<T ParameterBase::*>(varPtr), def));

      checkExprAccess(parName, depend, typeid(C));

      descriptor->setExpressionAccess(depend);
      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class T, class U>
    void addPar(const char *parName, T def, bool orig, T U::*varPtr, bool U::*givenPtr)
    {
      Descriptor *descriptor = new Descriptor(new Entry<T>(static_cast<T ParameterBase::*>(varPtr), def));

      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);
      descriptor->setGivenMember(givenPtr);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class T, class U>
    void addPar(const char *parName, T def, bool orig, T U::*varPtr)
    {
      Descriptor *descriptor = new Descriptor(new Entry<T>(static_cast<T ParameterBase::*>(varPtr), def));

      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class U>
    void addPar(const char *parName, const char *def, bool orig, std::string U::*varPtr, bool U::*givenPtr)
    {
      Descriptor *descriptor = new Descriptor(new Entry<std::string>(static_cast<std::string ParameterBase::*>(varPtr), def));

      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);
      descriptor->setGivenMember(givenPtr);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class U>
    void addPar(const char *parName, const char *def, bool orig, std::string U::*varPtr)
    {
      Descriptor *descriptor = new Descriptor(new Entry<std::string>(static_cast<std::string ParameterBase::*>(varPtr), def));

      if (orig)
        descriptor->setOriginalValueIndex(originalValueCount_++);
      descriptor->setSerialNumber(serialNumberCount_++);

      addDescriptor(parName, descriptor, typeid(C));
    }

    template<class U, class V>
    void addComposite(const char *comp_name, const ParametricData<U> &composite_pars, CompositeMap V::*composite_map)
    {
      Descriptor *descriptor = new Descriptor(new Entry<CompositeMap>(static_cast<CompositeMap ParameterBase::*>(composite_map)));

      descriptor->setSerialNumber(serialNumberCount_++);
      descriptor->setExpressionAccess(ParameterType::NO_DEP);
      descriptor->setUnit(U_INVALID);
      descriptor->setCategory(CAT_INVALID);
      descriptor->setCompositeParametricData(&composite_pars);

      addDescriptor(comp_name, descriptor, typeid(C));
    }

    void makeVector (const std::string &cname, int len) {
      for (int i = 1; i <= len; ++i) {
        std::ostringstream oss;
        oss << cname << i;
        std::string param = oss.str();

        ParameterMap::iterator it = map_.find(param);
        if (it == map_.end()) {
          std::string msg(" C::makeVector: parameter: ");
          msg += param;
          msg += " does not exist";
          N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
        }

        Descriptor &descriptor = *(*it).second;
        descriptor.setVec(i);
      }
    }
};

double getOriginalValue(ParameterBase *entity, int original_value_index);
void setOriginalValue(ParameterBase *entity, int original_value_index, double value);

bool wasValueGiven(ParameterBase *entity, int serial_number);
void setValueGiven(ParameterBase *entity, int serial_number, bool value);

} // namespace Device
} // namespace Xyce

#endif // Xyce_N_DEV_Pars_h
