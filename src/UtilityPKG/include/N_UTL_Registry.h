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
// Filename       : $RCSfile: N_UTL_Registry.h,v $
//
// Purpose        : Manage functions in a keyed registry
//
// Special Notes  : 
//
// Creator        : David G. Baur  Raytheon  Sandia National Laboratories 1355 <dgbaur@sandia.gov>
//
// Creation Date  : 2013/04/18 18:01:27
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.5.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:52 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#ifndef Xyce_N_UTL_UserPlugin_h
#define Xyce_N_UTL_UserPlugin_h

#include <memory>
#include <map>
// #include <unordered_map>
#include <vector>
#include <string>
#include <ostream>
#include <sstream>
#include <typeinfo>
#include <stdexcept>

/**
 * @file
 *
 * Historically, the term User Subroutine has referred to Fortran subroutines which are
 * generally called from within a procedure, allowing the user to create or select a
 * mathematical calculation to be applied to set of arguments.
 *
 * The purpose of this package is add User Plugins and User Functions as well as provide an
 * analogous implementation of the traditional User Subroutine.
 *
 * A User Plugin is a C++ class in which the application developer designs and implements
 * a base class that end users can create or select derivative classes to before the
 * desired operations.
 *
 * A User Subroutine is a function with a specific calling signature allowing any
 * registered subroutine with that signature to be called by its registered key.
 *
 * A User Function is a functional calculation which accepts one or more independent const
 * variables and returns a single dependent variable.  These variables may be of scalar,
 * vector or object quantities.  These functions can be created or selected by the end
 * user.
 *
 */

namespace Xyce {
namespace Plugin {

template<class T>
inline const T &demangle(const T &t) {
  return t;
}

/**
 * Class <b>Registry</b> serves as a singleton for holding templatized
 * createInstance and UserSubroutine function pointers and pointer to class factory objects.
 * The registry is simply a mapping of key pairs to a void pointer.  The key pair
 * consists of the base class key and the derived class keys.  And, since the only legal
 * way to get things into the registry is via the and UserSubroutine classes, there is no
 * signature checking performed at this level
 *
 * There should never be a need to instantiate this singleton as the UserSubroutine
 * registration process should perform that function when necessary.
 */

template <typename KeyType, class Cmp = std::less<KeyType> >
class Registry
{
public:
  typedef std::pair<const std::type_info *, KeyType> Key_;

  struct Cmp_ : std::binary_function<Key_, Key_, bool>
  {
    bool operator()(const Key_ &lhs, const Key_ &rhs ) const {
      if (*lhs.first == *rhs.first)
        return cmp(lhs.second, rhs.second);
      else
        return (*lhs.first).before(*rhs.first);
    }

    Cmp cmp;
  };

  typedef std::map<std::pair<const std::type_info *, KeyType>, void *, Cmp_> Map;

  /**
   * @brief Typedef <b>KeyPair</b> is the derived class key.
   *
   */
  typedef KeyType Key;
  typedef std::vector<Key> KeyVector;

  /**
   * @brief Typedef <b>KeyPair</b> is the derived class key.
   *
   */
  typedef typename Map::key_type KeyPair;

  /**
   * Creates a new <b>Registry</b> instance.
   *
   */
  Registry()
  {}

  /**
   * @brief Creates a new <b>Registry</b> instance and registers it, and more
   * importantly the derived class factory with the specified key pair.
   *
   * @param key_pair           a <b>KeyPair</b> gives the key to the registry
   *                            entry.
   */
  explicit Registry(const KeyPair &key_pair) {
    registerIt(key_pair, this);
  }

  /**
   * Destructor <b>~Registry</b> is virtual to fake polymorphism so that the registry class can
   * utilize additional compiler/runtime checks that the registered class is indeed a class factory.
   *
   */
  virtual ~Registry()
  {}

   /**
   * @brief Member function <b>registerDL</b> opens a dynamic library and optionally executes a "C"
   * registration function.
   *
   * If function is specified, and not zero length, the function must exist.  If the function name is
   * specified with zero length, the function "dl_register" is executed if it exists.  If the
   * function name is not found, if platform specific fortran suffix is appended and the function is
   * searched again.
   *
   * If no function name is specified, no registration function is executed.
   *
   * NOTE: Loading C++ sharable objects results in static object construction upon load.
   *
   * @param so_path             a <b>char</b> const pointer to the path of the
   *                            shareable object file.  If the path does not contains a
   *                            '/' character, the file is searched for through the
   *                            LD_LIBRARY_PATH envirnment variable.
   *
   * @param function_key       a <b>char</b> const pointer to the name of a registration function
   *                            which will be called immediately after loading the sharable object.
   *
   */
  static void registerDL(const char *so_path, const char *function_key = 0) {
#ifdef XYCE_DLOPEN_ENABLED
    slibout.m(Slib::LOG_PLUGIN) << "Loading dynamic library " << so_path << stk::diag::dendl;
    void *dl = dlopen(so_path, RTLD_NOW);
    if (!dl){
      throw std::runtime_error(dlerror());
    }

    if (function_key) {
      std::string s = std::strlen(function_key) ? function_key : "dl_register";

      dl_register_t f = (dl_register_t) dlsym(dl, s.c_str());
      if (!f) {
        s = s + XYCE_FORTRAN_SUFFIX;

        f = (dl_register_t) dlsym(dl, s.c_str());
      }

      if (f) {
        slibout.m(Slib::LOG_PLUGIN) << "Executing dynamic library " << so_path << " function " << s << "()" << stk::diag::dendl;
        (*f)();
      }
      else {
        if (std::strlen(function_key)) {
          std::ostringstream str;
          str << "Registration function " << function_key << " not found in " << so_path;
          throw std::runtime_error(str.str().c_str());
        }
      }
    }

#else
    throw std::runtime_error("Dynamic linkage is not supported on this platform");
#endif
  }

  template <typename T>
  static T getsym(const char *sym) {
#ifdef XYCE_DLOPEN_ENABLED
    void *s = NULL;
    void *dl = dlopen(NULL, RTLD_LAZY);
    if (dl) {
      s = dlsym(dl, sym);
      dlclose(dl);
    }

    return s;
#else
    return NULL;
#endif
  }

  /**
   * @brief Member function <b>registerIt</b> registers a key pair with a void
   * pointer.
   *
   * If the key pair already exists within the registry, a std::invalid_argument exception
   * is thrown with the offending key pair called out.
   *
   * @param key_pair           a <b>KeyPair</b> const reference to the key pair
   *                            to be registered.
   *
   * @param func_ptr            a <b>void</b> pointer to the function to be
   *                            registered.
   *
   * @throws                    a <b>std::invalid_argument</b> exception is thrown
   *                            if there is an instance creation function already
   *                            registered for the derived class.
   *
   */
  void registerIt(const KeyPair &key_pair, void *func_ptr) {
    // slibout.m(Slib::LOG_PLUGIN) << "Registering " << key_pair.second
    //                             << " of type " << demangle(key_pair.first->key())
    //                             << " at " << func_ptr << stk::diag::dendl;

    typename Map::const_iterator registry_entry = map.find(key_pair);
    if (registry_entry != map.end() && (*registry_entry).second != func_ptr) {
      std::ostringstream strout;
      strout << "Function with signature " << demangle((*registry_entry).first.first->name());
      strout << " and derived key '";
      strout << (*registry_entry).first.second;
      strout << "' already registered to create function at address " << (*registry_entry).second;
      throw std::invalid_argument(strout.str());
    }
    map[key_pair] = func_ptr;
  }


  /**
   * @brief Member function <b>getPluginPtr</b> find the function with the key pair
   * specified.
   *
   * If the key pair does not exist within the registry, a std::invalid_argument
   * exception is thrown with the offending key pair called out.
   *
   * @param key_pair           a <b>KeyPair</b> const reference to the key pair
   *                            to be retrieved.
   *
   * @throws                    a <b>std::invalid_argument</b> exception is thrown
   *                            if the function is not found.
   *
   */
  void *getPluginPtr(const KeyPair &key_pair) const {
    void *creator_function = getFuncPtr(key_pair);
    if (creator_function)
      return creator_function;
    else {
      std::ostringstream strout;

      strout << "User plugin creator function with base class '" << demangle(key_pair.first->name())
             << "' and derived class name '" << key_pair.second
             << "' not found in registry";

      throw std::invalid_argument(strout.str());
    }
  }


  /**
   * @brief Member function <b>getFunctionPtr</b> find the function with the key pair
   * specified.
   *
   * If the key pair does not exist within the registry, a std::invalid_argument
   * exception is thrown with the offending key pair called out.
   *
   * @param key_pair           a <b>KeyPair</b> const reference to the key pair
   *                            to be retrieved.
   *
   * @throws                    a <b>std::invalid_argument</b> exception is thrown
   *                            if the function is not found.
   *
   */
  void *getFunctionPtr(const KeyPair &key_pair) const {
    void *creator_function = getFuncPtr(key_pair);
    if (creator_function)
      return creator_function;
    else {
      std::ostringstream strout;

      strout << "User subroutine " << key_pair.second << "\n"
             << " with signature " << demangle(key_pair.first->name()) << "\n"
             << " not found in registry";

      throw std::invalid_argument(strout.str());
    }
  }


  /**
   * @brief Member function <b>getFuncPtr</b> returns the function pointer with the
   * specfied <it>key_pair</i>.
   *
   * @param key_pair           a <b>KeyPair</b> const reference to the registered
   *                            key pair.
   *
   * @returns                   a <b>void</b> function pointer with the specfied
   *                            <it>key_pair</i>.
   */
  Registry *getFactoryPtr(const KeyPair &key_pair) const {
    Registry *creator_function = (Registry *) getFuncPtr(key_pair);
    if (creator_function)
      return creator_function;
    else {
      std::ostringstream strout;
      strout << "Registry does not contain function with signature " << demangle(key_pair.first->name())
             << " and derived name '" << key_pair.second << "'";
      throw std::invalid_argument(strout.str());
    }
  }

  /**
   * @brief Member function <b>getFuncPtr</b> returns the function pointer with the
   * specfied <em>key_pair</em>.
   *
   * @param key_pair           a <b>KeyPair</b> const reference to the registered
   *                            key pair.
   *
   * @returns                   a <b>void</b> function pointer with the specfied
   *                            <em>key_pair</em>.
   */
  void *getFuncPtr(const KeyPair &key_pair) const {
    typename Map::const_iterator registry_entry = map.find(key_pair);
    return registry_entry == map.end() ? NULL : (*registry_entry).second;
  }

  /**
   * @brief Member function <b>getDerivedkeys</b> returns keys assocaited with the
   * function pointers of the specified type.
   *
   * @param type                a <b>std::type_info</b> const reference to typeid to
   *                            retrieve the derived names.
   *
   * @returns                   a <b>std::vector<str::string></b> value of the derived
   *                            names.
   */
  KeyVector getDerivedKey(const std::type_info &type) const {
    KeyVector derived_keys;

    for (typename Map::const_iterator it = map.begin(); it != map.end(); ++it)
      if (*(*it).first.first == type)
        derived_keys.push_back((*it).first.second);

    return derived_keys;
  }


  /**
   * Member template function <b>create</b> creates an instance of the desired
   * object by providing the factory responsible for generating that object type.  The
   * create factory is retrieved using the base class name specified by the factory base
   * class and the specified derived name.
   *
   * The derived factory is responsible for implementing the appropriate operator()
   * functions to constuct the derived object.
   *
   * @param derived_key        a <b>Key</b> const reference to the derived
   *                            object's name
   *
   * @return                    a <b>T</b> reference to the creation factory object.
   */
  template<class T>
  static T &create(const Registry<Key> &registry, const Key &derived_key) {
    return static_cast<T &>(*registry.getFactoryPtr(KeyPair(&typeid(T), derived_key)));
  }

  /**
   * @brief Member function <b>dump</b> dumps the registry.
   *
   * @param os                  a <b>std::ostream</b> reference to dump the registry
   *                            to.
   *
   * @return                    a <b>std::ostream</b> reference to <em>os</em>.
   */
  std::ostream &verbose_print(std::ostream &os) const {
    for (typename Map::const_iterator it = map.begin(); it != map.end(); ++it)
      os << (*it).first.second << " of type " << demangle((*it).first.first->name()) << " at " << (*it).second << std::endl;
    return os;
  }

private:
  Map           map;
};


/**
 * Template class <b>UserPlugin</b> is a template used for the association of base
 * and derived classed to be registered and created via the UserPlugin mechanism.  The
 * template traits enforces the signature matching of the base class constructor, derived class
 * constructor, derived class creator static function and the usage of the creator static
 * function.
 *
 * The registration requires a unique base class name for each base class type.  And,
 * since typeid is not reliable for that implementation, each base class is required to
 * implement a traits class with a Base typedef which specifies the base class, a
 * Signature typedef which specifies the signature of the create function and a static
 * function named getUserPluginCreatorName() which returns a const std::string reference to
 * the base classes name.  This name must be unique across users of the registry.  There
 * is no way to enforce this programatically, so it is recommended that a application
 * prefix be attached to the base class name.
 *
 */
template <class Registry, class Creator, typename S = Creator *(*)()>
class UserPlugin
{
public:
  typedef S Signature;                                  ///< Creator signature
  typedef typename Registry::Key Key;
  typedef typename Registry::KeyPair KeyPair;

private:
  UserPlugin();                                         ///< Not implemented
  UserPlugin(const UserPlugin&);                        ///< Not implemented
  UserPlugin &operator=(const UserPlugin&);             ///< Not implemented

public:
  /**
   * @brief Member function <b>registerCreator</b> registers the base class name and
   * specified derived class name with the specified creator function.
   *
   * The base class name is determined by the BaseTraits template argument's
   * getUserPluginCreatorName() static member function.  The signature is defined by the
   * BaseTraits template argument's Signature typedef.
   *
   * @param derived_key        a <b>Key</b> const reference to the derived
   *                            class name.
   *
   * @param function            a <b>signature</b> function pointer to the creator
   *                            function.
   *
   */
  static void registerCreator(Registry &registry, const Key &derived_key, Signature function) {
    registry.registerIt(KeyPair(&typeid(Signature), derived_key), (void *) function);
  }

  /**
   * @brief Member function <b>create</b> returns the createInstance() function
   * associated with the specified base class and derived_key.
   *
   * @param derived_key        a <b>Key</b> const reference to the derived
   *                            classes name.
   *
   * @throws                    a <b>std::invalid_argument</b> exception is thrown
   *                            if there is no instance creation function registered for
   *                            the specified name of the base class.
   *
   * @return                    a <b>Signature</b> function pointer to the instance
   *                            create function.
   */
  static Signature create(const Registry &registry, const Key &derived_key) {
    Signature creator_function = (Signature) registry.getPluginPtr(KeyPair(&typeid(Signature), derived_key));

    return (*creator_function);
  }

  /**
   * @brief Member function <b>exists</b> returns true if class of the type
   * specified by derived_key exists in BaseClass.
   *
   * @param derived_key        a <b>Key</b> const reference to the derived
   *                            classes name.
   *
   * @return                    a <b>bool</b> of true if class of the type specified by derived_key exists in BaseClass.
   */
  static bool exists(const Registry &registry, const Key &derived_key) {
    return registry.getFuncPtr(KeyPair(&typeid(Signature), derived_key)) != NULL;
  }

  /**
   * @brief Class template <b>Register</b> registers the <i>createInstance()</i>
   * function with the <i>derived_key</i> on object creation.
   *
   * @param DerivedClass        a <b>class</b> which specifies the derived class
   *                            which holds the <i>createInstance()</i> function.
   *
   */
  template <class DerivedClass>
  class Register
  {
  public:
    typedef DerivedClass XDerivedClass;

    /**
     * @brief Creates a new <b>Register</b> instance.  Upon creation, the
     * <i>DerivedClass::createInstance()</i> instance creation function is registered
     * with the <i>derived_key</i>.
     *
     * @param derived_key      a <b>Key</b> const reference to the derived
     *                          class' name.
     *
     */
    explicit Register(Registry &registry, const Key &derived_key)
      : m_function(DerivedClass::createInstance)
    {
      UserPlugin<Key, Creator, Signature>::registerCreator(registry, derived_key, m_function);
    }

    /**
     * @brief Creates a new <b>Register</b> instance.  Upon creation, the
     * <i>DerivedClass::createInstance()</i> instance creation function is registered
     * with the <i>derived_key</i>.
     *
     * @param derived_key      a <b>Key</b> const reference to the derived
     *                          class' name.
     *
     */
    Register(Registry &registry, const Key &derived_key, Signature create_instance)
      : m_function(create_instance)
    {
      UserPlugin<Registry, Creator, Signature>::registerCreator(registry, derived_key, m_function);
    }

  private:
    Signature           m_function;                     ///< Place to hold function pointer
  };
};


/**
 * Class template <b>UserSubroutine</b> is a template used for the association of
 * user function to be registered and acquired via the UserSubroutine mechanism.  The
 * template traits enforces the signature matching of the user function and the its
 * usage.
 *
 * The registration requires a unique function name and is required to implement a
 * traits class with a Signature typedef which specifies the signature of the user
 * function and a static function named getUserSubroutineName() which returns a const
 * std::string reference to the user function's name.  This name must be unique across users
 * of the registry.  There is no way to enforce this programatically, so it is recommended
 * that a application prefix be attached to the function name.
 *
 */
template <class Registry, class S>
class UserSubroutine
{
private:
  UserSubroutine();                                     ///< Not implemented
  UserSubroutine(const UserSubroutine&);                ///< Not implemented
  UserSubroutine &operator=(const UserSubroutine&);     ///< Not implemented

public:
  typedef S Signature;                                  ///< Subroutine call signature
  typedef typename Registry::Key Key;
  typedef typename Registry::KeyPair KeyPair;

  /**
   * @brief Member function <b>registerFunction</b> registers the user function's
   * name with the specified user function function pointer.
   *
   * The base class name is determined by the BaseClass template argument's
   * getCreatorName() static member function.  The signature is defined
   * UserSubroutineTraits template argument's Signature typedef.
   *
   * @param function_key       a <b>Key</b> const reference to the user
   *                            function's name.
   *
   * @param function            a <b>Signature</b> function pointer to the user
   *                            function.
   *
   */
  inline static void registerFunction(Registry &registry, const Key &function_key, Signature *function) {
    registry.registerIt(KeyPair(&typeid(Signature), function_key), (void *) function);
  }

  /**
   * @brief Member function <b>execute</b> returns the user function function
   * associated with the specified signature and derived_key.
   *
   * @param function_key       a <b>Key</b> const reference to the user
   *                            function's name.
   *
   * @throws                    a <b>std::invalid_argument</b> exception is thrown
   *                            if there is no user function registered for the
   *                            specified name.
   *
   * @return                    a <b>Signature</b> user function.
   */
  static Signature *execute(Registry &registry, const Key &function_key) {
    Signature *user_function = (Signature *) registry.getFunctionPtr(KeyPair(&typeid(Signature), function_key));

    return (*user_function);
  }

  /**
   * @brief Member function <b>execute</b> returns the user function function
   * pointer associated with the specified signature and derived_key.
   *
   * @param function_key       a <b>Key</b> const reference to the user
   *                            function's name.
   *
   * @throws                    a <b>std::invalid_argument</b> exception is thrown
   *                            if there is no user function registered for the
   *                            specified name.
   *
   * @return                    a <b>Signature</b> user function pointer.
   */
  static Signature *getFunction(const Registry &registry, const Key &function_key) {
    Signature *user_function = (Signature *) registry.getFunctionPtr(KeyPair(&typeid(Signature), function_key));

    return user_function;
  }

  /**
   * @brief Member function <b>exists</b> returns true if user function specified by
   * derived_key exists.
   *
   * @param function_key       a <b>Key</b> const reference to the user
   *                            function's name.
   *
   * @return                    a <b>bool</b> of true if user function specified
   *                            signature and <i>function_key</i> exists in BaseClass.
   */
  static bool exists(const Registry &registry, const Key &derived_key) {
    return registry.getFuncPtr(KeyPair(&typeid(Signature), derived_key)) != NULL;
  }

  /**
   * @brief Class template <b>Register</b> registers the user function function
   * pointer with the <i>function_key</i> on object creation.
   *
   */
  class Register
  {
  public:
    /**
     * @brief Creates a new <b>Register</b> instance.  Upon creation, the
     * <i>func_ptr()</i> function is registered with the <i>function_key</i>.
     *
     * @param function_key     a <b>Key</b> const reference to the user
     *                          function's name.
     *
     */
    Register(Registry &registry, const Key &function_key, Signature *function)
      : m_function(function)
    {
      UserSubroutine<Registry, Signature>::registerFunction(registry, function_key, *m_function);
    }

  private:
    Signature *                 m_function;                     ///< Holder of the function pointer
  };
};

// template <typename Key>
// void *Registry<Key>::getsym<void *>(const char *sym);

// template <typename Key, typename T>
// inline T Registry<Key>::getsym(const char *sym) {
//   return static_cast<T>(getsym<void *>(sym));
// }

} // namespace Plugin
} // namespace Xyce

typedef std::type_info *type_info_func();

#endif // Xyce_N_UTL_UserPlugin_h
