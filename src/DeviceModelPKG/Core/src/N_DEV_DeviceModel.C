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
// Filename       : $RCSfile: N_DEV_DeviceModel.C,v $
//
// Purpose        :
//
// Special Notes  :
//
// Creator        : Eric R. Keiter, SNL, Parallel Computational Sciences
//
// Creation Date  : 04/03/00
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.67 $
//
// Revision Date  : $Date: 2014/02/24 23:49:15 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>

#include <iostream>
#include <set>

#include <N_DEV_Const.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Param.h>
#include <N_DEV_Message.h>
#include <N_DEV_DeviceMgr.h>
#include <N_DEV_Configuration.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Device {

const char *modelEntityType = "model";

//-----------------------------------------------------------------------------
// Function      : DeviceModel::DeviceModel
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceModel::DeviceModel(
  const ModelBlock &            model_block,
  ParametricData<void> &        parametric_data,
  const FactoryBlock &          factory_block)
  : DeviceEntity(modelEntityType, model_block.name, parametric_data, factory_block.solverState_, factory_block.deviceOptions_, model_block.netlistFileName_, model_block.lineNumber_),
    type_(model_block.type),
    level_(model_block.level),
    temperatureModel(""),
    doseModel(""),
    iModel(TEMP),
    iMethod(QUAD),
    base_temp(CONSTREFTEMP)
{}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::setModParams
// Purpose       : Set up parameter fits from model line
// Special Notes :
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/26/05
//-----------------------------------------------------------------------------
void DeviceModel::setModParams(const std::vector<Param> &params)
{
  std::vector<int> m_start;
  std::set<std::string> pname;
  std::map<std::string, int> ptype;

  int param_index = 0;
  m_start.push_back(param_index);
  for (std::vector<Param>::const_iterator mp = params.begin(); mp != params.end(); ++mp)
  {
    ++param_index;
    if ((*mp).tag() == "INDEPENDENT;PARAM")
    {
      m_start.push_back(param_index);
      pname.clear();
    }
    else
    {
      if (pname.find((*mp).tag()) != pname.end())
      {
        UserError0(*this) << "Duplicate specification of parameter " << (*mp).tag();
        return;
      }
      pname.insert((*mp).tag());
      if (m_start.size() == 1)
      {
        ptype[(*mp).tag()] = (*mp).getType();
      }
    }
  }
  if (m_start.size() == 1)
  {
    setParams(params);
  }
  else
  {
    m_start.push_back(param_index + 1);

    // An interpolation method is present, first figure out what it
    // is and make sure that all models agree on the method.
    std::string tmod("");
    std::string dmod("");
    std::string modName("");

    for (int i = 0; i < m_start.size()-1; ++i)
    {
      for (int j = m_start[i]; j < m_start[i+1]-1; ++j)
      {
        if (params[j].tag() == "TEMPMODEL" && params[j].stringValue() != "NONE")
        {
          if (i == 0)
          {
            tmod = params[j].stringValue();
          }
          else
          {
            if (tmod != params[j].stringValue())
            {
              UserError0(*this) << "Inconsistent or missing TEMPMODEL parameter, " << params[j].stringValue() << " specified here, " << tmod << " specified perviously";
              return;
            }
          }
        }
        if (params[j].tag() == "DOSEMODEL" && params[j].stringValue() != "NONE")
        {
          if (i == 0)
          {
            dmod = params[j].stringValue();
          }
          else
          {
            if (dmod != params[j].stringValue())
            {
              UserError0(*this) << "Inconsistent or missing DOSEMODEL parameter, " << params[j].stringValue() << " specified here, " << dmod << " specified perviously";
              return;
            }
          }
        }
      }
    }

    if (tmod == "" && dmod == "")
    {
      UserError0(*this) << "Duplicate model specification implies parameter interpolation, TEMPMODEL or DOSEMODEL parameters must specified";
      return;
    }
    else if (tmod != "" && dmod != "")
    {
      UserError0(*this) << "Only one of TEMPMODEL or DOSEMODEL parameters may be specified";
      return;
    }
    else if (tmod != "")
    {
      modName = tmod;
      iModel = TEMP;
    }
    else if (dmod != "")
    {
      modName = dmod;
      iModel = DOSE;
    }

    if (modName == "QUADRATIC")
    {
      iMethod = QUAD;
      if (m_start.size() != 4)
      {
        UserError0(*this) << "Three model specifications required for QUADRATIC fit";
        return;
      }
      fit.resize(3);
    }
    else if (modName == "PWL")
    {
      iMethod = PWL;
      fit.resize(m_start.size() - 1);
    }
    else
    {
      UserError0(*this) << "Only QUADRATIC or PWL interpolation method is supported";
      return;
    }

    // First find the params that vary between the specified models:
    std::map<std::string,double> basePars;
    std::vector<double> t;
    std::string par;

    for (int i = 0 ;i < m_start.size()-1; ++i)
    {
      for (int j = m_start[i]; j < m_start[i+1]-1; ++j)
      {
        par = params[j].tag();
        ParameterMap::const_iterator parIt = getParameterMap().find(par);
        if (parIt == getParameterMap().end())
          DevelFatal0(*this).in("DeviceModel::setModParams") << "Parameter " << par << " not found";

        const Descriptor &nPar = *(*parIt).second;
        if (params[j].given() && nPar.isType<double>())
        {
          if (iModel == TEMP)
          {
            if (par == "TNOM")
            {
              t.push_back(params[j].getImmutableValue<double>());
            }
          }
          else if (iModel == DOSE)
          {
            if (par == "DOSE")
            {
              t.push_back(params[j].getImmutableValue<double>());
            }
          }
          if (i == 0)
          {
            if (params[j].getType() != Util::EXPR)
            {
              basePars[par] = params[j].getImmutableValue<double>();
            }
          }
          else
          {
            if (ptype[par] == Util::EXPR)
            {
              UserWarning0(*this) << "Non-constant expression for parameter " << par << ", it will not interpolated";
            }
            else
            {
              if (basePars.find(par) == basePars.end())
              {
                UserError0(*this) << "Unknown parameter " << params[j].tag() <<  " in temperature compensation .MODEL statement";
                return;
              }
              if (basePars[par] != params[j].getImmutableValue<double>() && params[j].given())
              {
                if (fitMap.find(par) == fitMap.end())
                {
                  fitMap[par] = fit[0].size();
                  fitParams.push_back(nPar.getMemberPtr<double>());
                  fit[0].push_back(basePars[par]);
                }
              }
            }
          }
        }
      }
    }

    if (t.size() != m_start.size()-1)
    {
      UserError0(*this) << (iModel == TEMP ? "TNOM" : "DOSE") << " not specified in all .MODEL statements";
      return;
    }

    for (int i = 1; i < t.size(); ++i)
    {
      for (int j = 0; j < i; ++j)
      {
        if (t[i] == t[j])
        {
          UserError0(*this) << "Identical " << (iModel == TEMP ? "TNOM" : "DOSE") << " values in .MODEL statements";
          return;
        }
      }
    }

    // Now, collect the values to use for the fits:
    int nFit = fitMap.size();
    //int nSet = t.size();
    std::vector<std::vector<double> > temp(nFit);
    std::vector<std::vector<double> > val(nFit);
    oldParams.resize(nFit);

    for (int i = 1; i < fit.size(); ++i)
    {
      fit[i].resize(nFit);
    }
    min_par.resize(nFit);
    max_par.resize(nFit);
    parType.resize(nFit);
    std::map<std::string, int>::iterator fm = fitMap.begin();
    std::map<std::string, int>::iterator fm_end = fitMap.end();

    for ( ; fm != fm_end ; ++fm)
    {
      par = (*fm).first;
      ParameterMap::const_iterator parIt = getParameterMap().find(par);
      if (parIt == getParameterMap().end())
        throw std::runtime_error(std::string("Parameter ") + par + " not found");
      {
        const Descriptor &nPar = *(*parIt).second;

        if (nPar.getExpressionAccess() & ParameterType::LOG_T_DEP)
        {
          parType[fitMap[par]] = LOG_FIT;
        }
        else
        {
          parType[fitMap[par]] = LINEAR_FIT;
        }
      }
    }

    base_temp = t[0];
    if (iModel == TEMP)
    {
      base_temp += CONSTCtoK;
    }

    for (int i = 0; i < nFit; ++i)
    {
      temp[i].push_back(t[0]);
      val[i].push_back(fit[0][i]);
      if (parType[i] == LOG_FIT)
      {
        if (val[i][0] <= 0)
        {
          UserError0 message(*this);
          message << "Negative parameter value for log interpolation based parameter ";
          fm = fitMap.begin();
          for ( ; fm != fm_end ; ++fm)
          {
            if ((*fm).second == i)
            {
              par = (*fm).first;
              break;
            }
          }
          message << par;
          return;
        }
        val[i][0] = log(val[i][0]);
      }
      min_par[i] = val[i][0];
      max_par[i] = val[i][0];
    }

    for (int i = 1; i < m_start.size()-1; ++i)
    {
      for (int j = m_start[i]; j < m_start[i+1]-1; ++j)
      {
        if (params[j].getType() == Util::DBLE && params[j].given())
        {
          par = params[j].tag();
          std::map<std::string, int>::iterator fm1 = fitMap.find(par);
          if (fm1 != fitMap.end())
          {
            int k = fm1->second;
            temp[k].push_back(t[i]);
            double p = params[j].getImmutableValue<double>();
            if (parType[k] == LOG_FIT)
            {
              if (p <= 0)
              {
                UserError0 message(*this);
                message << "Negative parameter value for log interpolation based parameter ";
                fm = fitMap.begin();
                for ( ; fm != fm_end ; ++fm)
                {
                  if ((*fm).second == k)
                  {
                    par = (*fm).first;
                    break;
                  }
                }
                message << par;
                return;
              }
              p = log(p);
            }
            val[k].push_back(p);
            if (p > max_par[k])
              max_par[k] = p;
            if (p < min_par[k])
              min_par[k] = p;
          }
        }
      }
    }

    // Finally, do the actual fits:
    if (fitMap.size() > 0)
    {
      if (iMethod == QUAD)
      {
        std::map<std::string, int>::iterator f;
        for (f=fitMap.begin() ; f!=fitMap.end() ; ++f)
        {
          int i = (*f).second;
          if (temp[i].size() == 2)
          {
            fit[0][i] = val[i][0];
            fit[1][i] = (val[i][1] - val[i][0])/(temp[i][1] - temp[i][0]);
            fit[2][i] = 0;
          }
          else if (temp[i].size() == 3)
          {
            fit[0][i] = val[i][0];
            double x1,x2,y1,y2;
            x1 = temp[i][1] - temp[i][0];
            y1 = val[i][1];
            x2 = temp[i][2] - temp[i][0];
            y2 = val[i][2];
            fit[2][i] = (x2*y1-x1*y2-fit[0][i]*(x2-x1))/(x2*x1*x1-x1*x2*x2);
            fit[1][i] = (y1-fit[2][i]*x1*x1-fit[0][i])/x1;
          }
          else
          {
            DevelFatal0(*this).in("setModParams") << "Internal error in DeviceModel, illegal number of fit points for parameter " << (*f).first;
          }
          if ((*f).first == "TNOM")
          {
            fit[0][i] += CONSTCtoK;
          }
        }
      }
      else if (iMethod == PWL)
      {
        int nT = fit.size();
        std::map<double,int> tOrd;
        for (int i = 0; i < nT; ++i)
        {
          tOrd[t[i]] = 0;
        }
        int i = 0;
        std::map<double,int>::iterator tOrd_i = tOrd.begin();
        std::map<double,int>::iterator tOrd_end = tOrd.end();
        for ( ; tOrd_i != tOrd_end; ++tOrd_i)
        {
          if (iModel == TEMP)
          {
            base.push_back((*tOrd_i).first+CONSTCtoK);
          }
          else
          {
            base.push_back((*tOrd_i).first);
          }
          (*tOrd_i).second = i++;
        }
        std::map<std::string, int>::iterator f;
        std::map<std::string, int>::iterator f_end;
        std::vector<bool> p(nT,false);
        f=fitMap.begin();
        f_end=fitMap.end();
        for ( ; f!=f_end; ++f)
        {
          i = (*f).second;
          for (int j = 0 ;j < nT; ++j)
          {
            p[j] = false;
          }
          for (int j = 0; j < temp[i].size() ; ++j)
          {
            int ind = tOrd[temp[i][j]];
            p[ind] = true;
            fit[ind][i] = val[i][j];
          }
          for (int j = 0; j < nT ; ++j)
          {
            if (!p[j])
            {
              int k_lo = j;
              int k_hi = j;
              while (k_lo >= 0 && !p[k_lo])
              {
                k_lo--;
              }
              while (k_hi <nT && !p[k_hi])
              {
                ++k_hi;
              }
              if (k_lo == -1)
              {
                if (k_hi < nT)
                {
                  fit[j][i] = fit[k_hi][i];
                }
                else
                {
                  DevelFatal0(*this).in("DeviceModel::setModParams") <<"DeviceModel::setModParams: Internal error forming PWL interpolation";
                }
              }
              else
              {
                if (k_hi < nT)
                {
                  double frac = (base[j]-base[k_lo])/(base[k_hi]-base[k_lo]);
                  fit[j][i] = fit[k_hi][i]*frac+fit[k_lo][i]*(1-frac);
                }
                else
                {
                  fit[j][i] = fit[k_lo][i];
                }
              }
            }
            if ((*f).first == "TNOM")
            {
              fit[j][i] += CONSTCtoK;
            }
          }
        }
      }
    }

    // params.resize(m_start[1]-1);
    // setParams(params);
    std::vector<Param> remaining_params(&params[0], &params[m_start[1] - 1]);
    setParams(remaining_params);
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::~DeviceModel
// Purpose       : destructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceModel::~DeviceModel()
{
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::saveParams
// Purpose       : save existing param values before fitting params to temperature
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/29/05
//-----------------------------------------------------------------------------
void DeviceModel::saveParams()
{
  int nFit = fitMap.size();
  int i;

  if (nFit == 0)
  {
    return;
  }

  for (i=0 ; i<nFit ; ++i)
  {
    oldParams[i] = this->*(fitParams[i]);
  }
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::interpolateTNOM
// Purpose       : interpolate param values to a specified temperature
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/29/05
//-----------------------------------------------------------------------------
bool DeviceModel::interpolateTNOM(double t)
{
  if (iModel != TEMP)
  {
    return false;
  }

  return interpolate(t);
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::interpolateDOSE
// Purpose       : interpolate param values to a specified temperature
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/11/05
//-----------------------------------------------------------------------------
bool DeviceModel::interpolateDOSE(double d)
{
  if (iModel != DOSE)
  {
    return false;
  }

  return interpolate(d);
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::interpolated
// Purpose       : returns true if an interpolation is in effect
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 1/31/06
//-----------------------------------------------------------------------------
bool DeviceModel::interpolated()
{
  return (fitMap.size() > 0);
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::interpolate
// Purpose       : interpolate param values to a specified temperature
// Special Notes :
// Scope         : private
// Creator       : Dave Shirley, PSSI
// Creation Date : 11/11/05
//-----------------------------------------------------------------------------
bool DeviceModel::interpolate(double t)
{
  int nFit = fitMap.size();
  int i, j, k_hi, k_lo;
  double del;
  double frac, p;

  if (nFit == 0)
  {
    return false;
  }

  if (iMethod == QUAD)
  {
    del = t - base_temp;
    //    for (i=0 ; i<nFit ; ++i)
    std::map<std::string,int>::iterator fp;
    std::map<std::string,int>::iterator fm_begin=fitMap.begin();
    std::map<std::string,int>::iterator fm_end=fitMap.end();
    for (fp=fm_begin; fp != fm_end; fp++)
    {
      i=fp->second;
      p = (fit[2][i]*del + fit[1][i])*del + fit[0][i];

      if (p > max_par[i] && i>0)
      {
        this->*(fitParams[i]) = max_par[i];
      }
      else if (p < min_par[i] && i>0)
      {
        this->*(fitParams[i]) = min_par[i];
      }
      else
      {
        this->*(fitParams[i]) = p;
      }
    }
  }
  else if (iMethod == PWL)
  {
    del = t;
    k_hi = 0;
    for (j=0 ; j<fit.size() ; ++j)
    {
      if (base[j] >= del)
      {
        break;
      }
      k_hi = j+1;
    }
    if (k_hi == 0)
    {
      frac = 0;
    }
    else if (k_hi == fit.size())
    {
      k_hi = fit.size()-1;
      frac = 1;
    }
    else
    {
      k_lo = k_hi-1;
      frac = (del-base[k_lo])/(base[k_hi]-base[k_lo]);
    }
    if (frac == 1)
    {
      for (i=0 ; i<nFit ; ++i)
      {
        this->*(fitParams[i]) = fit[k_hi][i];
      }
    }
    else
    {
      for (i=0 ; i<nFit ; ++i)
        this->*(fitParams[i]) = fit[k_hi][i]*frac + fit[k_lo][i]*(1-frac);
    }
  }
  for (i=0 ; i<nFit ; ++i)
  {
    if (parType[i] == LOG_FIT)
    {
      this->*(fitParams[i]) = exp(this->*(fitParams[i]));
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
// Function      : DeviceModel::restoreParams
// Purpose       : restore previously saved param values
// Special Notes :
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/29/05
//-----------------------------------------------------------------------------
void DeviceModel::restoreParams()
{
  int nFit = fitMap.size();
  int i;

  if (nFit == 0)
  {
    return;
  }

  for (i=0 ; i<nFit ; ++i)
  {
    this->*(fitParams[i]) = oldParams[i];
  }
}

} // namespace Device
} // namespace Xyce
