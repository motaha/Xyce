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

//-------------------------------------------------------------------------
// Filename       : $RCSfile: N_DEV_CompositeParam.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Dave Shirley, PSSI
//
// Creation Date  : 05/05/05
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.16.2.2 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>
#include <string>

// ----------   Xyce Includes   ----------

#include <N_DEV_CompositeParam.h>
#include <N_DEV_Const.h>
#include <N_DEV_Param.h>

namespace Xyce {
namespace Device {

//-----------------------------------------------------------------------------
// Function      : CompositeParam::CompositeParam
// Purpose       : constructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/06/05
//-----------------------------------------------------------------------------
CompositeParam::CompositeParam()
  : name_    ("")
{}

//-----------------------------------------------------------------------------
// Function      : CompositeParam::~CompositeParam
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/06/05
//-----------------------------------------------------------------------------
CompositeParam::~CompositeParam()
{
}

//-----------------------------------------------------------------------------
// Function      : CompositeParam::setDefaultParams
// Purpose       : Set parameters according to default values
// Special Notes : Normally called before setParams
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/06/05
//-----------------------------------------------------------------------------
void CompositeParam::setDefaultParams( )
{
// First, allocate and zero out original and given vals
  for (ParameterMap::const_iterator par_i =(*getPMap()).begin() ; par_i !=(*getPMap()).end() ; ++par_i)
  {
    Pars &p = static_cast<Pars &>(*(*par_i).second);

    if (p.hasGivenMember())
      p.setGiven(*this, false);

    if (p.isType<double>())
      setValue<double, CompositeParam>(*this, p, getDefaultValue<double>(p));
    else if (p.isType<bool>())
      setValue<bool, CompositeParam>(*this, p, getDefaultValue<bool>(p));
    else if (p.isType<int>())
      setValue<int, CompositeParam>(*this, p, getDefaultValue<int>(p));
    else if (p.isType<long>())
      setValue<long, CompositeParam>(*this, p, getDefaultValue<long>(p));
    else if (p.isType<std::string>())
      setValue<std::string, CompositeParam>(*this, p, getDefaultValue<std::string>(p));
    else if (p.isType<std::vector<double> >())
      value<std::vector<double> >(*this, p).clear();
    else if (p.isType<std::vector<std::string> >())
      value<std::vector<std::string> >(*this, p).clear();
  }
}

//-----------------------------------------------------------------------------
// Function      : CompositeParam::outputParams
// Purpose       : Write out parameter set(for debugging)
// Special Notes :
// Scope         : protected
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/06/05
//-----------------------------------------------------------------------------
std::ostream &CompositeParam::printParams(std::ostream &os)
{
  for (ParameterMap::const_iterator it_parameter = getPMap()->find("NAME"); it_parameter != getPMap()->end() ; )
  {
    os << endl << "   " <<(*it_parameter).first << " ";
    const Pars &p = static_cast<const Pars &>(*(*it_parameter).second);
    if (p.isType<double>())
      os << "(double) : " << p.value<double>(*this);
    else if (p.isType<bool>())
      os << "(bool) : " << p.value<bool>(*this);
    else if (p.isType<int>())
      os << "(int) : " << p.value<int>(*this);
    else if (p.isType<long>())
      os << "(long) : " << p.value<long>(*this);
    else if (p.isType<std::string>())
      os << "(string) : " << p.value<std::string>(*this);
    else if (p.isType<std::vector<std::string> >())
    {
      const std::vector<std::string> &string_vector = p.value<std::vector<std::string> >(*this);

      for (std::vector<std::string>::const_iterator it = string_vector.begin(); it != string_vector.end(); ++it)
        os << "  " << *it;
    }
    else if (p.isType<std::vector<double> >())
    {
      const std::vector<double> &string_vector = p.value<std::vector<double> >(*this);

      for (std::vector<double>::const_iterator it = string_vector.begin(); it != string_vector.end(); ++it)
        os << "  " << *it;
    }

    if ((*it_parameter).first == "NAME")
      it_parameter =(*getPMap()).begin();
    else
      it_parameter++;

    if (it_parameter !=(*getPMap()).end() &&(*it_parameter).first == "NAME")
      it_parameter++;
  }
}

//-----------------------------------------------------------------------------
// Function      : CompositeParam::setParams
// Purpose       : Set parameter according to input Param.  Used to
//                 set instance or model parameter sets to netlist values
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/06/05
//-----------------------------------------------------------------------------

void CompositeParam::setParams( const string & pName, const Param & ndParam )
{
  ParameterMap::const_iterator p_i = getPMap()->find(pName);
  if (p_i != getPMap()->end())
  {
    const Pars &p = static_cast<const Pars &>(*(*p_i).second);
    if (p.hasGivenMember())
    {
      if (ndParam.given())
      p.setGiven(*this, true);
      else if (p.getGiven(*this))
        return;
    }
    setValueGiven(this, p.getSerialNumber(), ndParam.given());
    if (ndParam.given() || ndParam.default_val())
    {
      if (p.isType<double>())
      {
        p.value<double>(*this) = ndParam.dVal();
        if (pName == "TNOM" || pName == "TEMP")
          p.value<double>(*this) += CONSTCtoK;
        if (p.getOriginalValueIndex() >= 0)
          setOriginalValue(this, p.getOriginalValueIndex(), p.value<double>(*this));
      }
      else if (p.isType<std::string>())
      {
        p.value<std::string>(*this) = ndParam.sVal();
      }
      else if (p.isType<int>())
      {
        p.value<int>(*this) = ndParam.iVal();
        if (p.getOriginalValueIndex() >= 0)
          setOriginalValue(this, p.getOriginalValueIndex(), static_cast<double>(p.value<int>(*this)));
      }
      else if (p.isType<long>())
      {
        p.value<long>(*this) = ndParam.lVal();
        if (p.getOriginalValueIndex() >= 0)
          setOriginalValue(this, p.getOriginalValueIndex(), static_cast<double>(p.value<long>(*this)));
      }
      else if (p.isType<bool>())
      {
        p.value<bool>(*this) = (ndParam.dVal() != 0.0);
        if (p.getOriginalValueIndex() >= 0)
        {
          if (p.value<bool>(*this))
            setOriginalValue(this, p.getOriginalValueIndex(), 1.0);
          else
            setOriginalValue(this, p.getOriginalValueIndex(), 0.0);
        }
      }
      else if (p.isType<std::vector<double> >())
      {
       (p.value<std::vector<double> >(*this)).push_back(ndParam.dVal());
      }
      else if (p.isType<std::vector<std::string> >())
      {
        p.value<std::vector<std::string> >(*this).push_back(ndParam.sVal());
      }
      else
      {
        string msg("CompositeParam::setParams: unknown type");
        N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
      }
    }
  }
  else
  {
    string msg("CompositeParam::setParams: undefined parameter: ");
    msg += pName;
    msg += "\nThis parameter is in metadata, but not recognized in constructor";
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }
}

//-----------------------------------------------------------------------------
// Function      : CompositeParam::given
// Purpose       : Return whether param was given
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 05/06/05
//-----------------------------------------------------------------------------

bool CompositeParam::given(const string &parameter_name )
{
  ParameterMap::const_iterator it =(*getPMap()).find(parameter_name);
  if (it == (*getPMap()).end())
  {
    string msg("CompositeParam::Given: unrecognized param: ");
    msg += parameter_name;
    N_ERH_ErrorMgr::report(N_ERH_ErrorMgr::DEV_FATAL, msg);
  }

  return wasValueGiven(this, (*it).second->getSerialNumber());
}

} // namespace Device
} // namespace Xyce
