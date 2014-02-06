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
// Revision Number: $Revision: 1.53.2.3 $
//
// Revision Date  : $Date: 2013/10/03 17:23:38 $
//
// Current Owner  : $Author: tvrusso $
//-------------------------------------------------------------------------

#include <Xyce_config.h>


// ---------- Standard Includes ----------
#include <iostream>
#include <set>

// ----------   Xyce Includes   ----------
#include <N_DEV_Const.h>
#include <N_DEV_DeviceModel.h>
#include <N_DEV_DeviceBlock.h>
#include <N_DEV_DeviceOptions.h>
#include <N_DEV_SolverState.h>
#include <N_DEV_Param.h>

#include <N_ERH_ErrorMgr.h>

namespace Xyce {
namespace Device {


//-----------------------------------------------------------------------------
// Function      : DeviceModel::DeviceModel
// Purpose       : model block constructor
// Special Notes :
// Scope         : public
// Creator       : Eric Keiter, SNL, Parallel Computational Sciences
// Creation Date : 3/30/00
//-----------------------------------------------------------------------------
DeviceModel::DeviceModel(
  const ModelBlock &    model_block,
  SolverState &         solver_state,
  DeviceOptions &       device_options)
  : DeviceEntity(solver_state, device_options, model_block.name, model_block.netlistFileName_, model_block.lineNumber_),
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
// ERK:  5/27/08.  The argument params cannot be a reference because the data
//       passed to it (generally MB.params) is const.  I tried making
//       it both const, and a reference, but the resize of params at the
//       very end of this function violates const.
//
//       As a result, this function implicitly depends on the copying
//       the params object into a temporary object, and doing that copy
//       correctly.
//
// Scope         : public
// Creator       : Dave Shirley, PSSI
// Creation Date : 8/26/05
//-----------------------------------------------------------------------------
void DeviceModel::setModParams(vector<Param> params)
{
  vector<Param>::const_iterator mp;
  vector<Param>::const_iterator mp_begin = params.begin();
  vector<Param>::const_iterator mp_end = params.end();
  int i, j, k, k_lo, k_hi;
  vector<int> m_start;
  set<string> pname;
  map<string,int> ptype;

  i = 0;
  m_start.push_back(i);
  for (mp = mp_begin ; mp != mp_end ; ++mp)
  {
    ++i;
    if ((*mp).tag() == "INDEPENDENT;PARAM")
    {
      m_start.push_back(i);
      pname.clear();
    }
    else
    {
      if (pname.find((*mp).tag()) != pname.end())
      {
        string msg = "Duplicate specification of parameter: ";
        msg += (*mp).tag();
        msg += " in model: ";
        msg += getName();
        std::ostringstream oss;
        oss << "Error in " << netlistLocation() << "\n" << msg;
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, oss.str());
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
    return;
  }
  m_start.push_back(i+1);

			 // An interpolation method is present, first figure out what it
			 // is and make sure that all models agree on the method.
  string tmod("");
  string dmod("");
  string modName("");

  for (i=0 ; i<m_start.size()-1 ; ++i)
  {
    for (j=m_start[i] ; j<m_start[i+1]-1 ; ++j)
    {
      if (params[j].tag() == "TEMPMODEL" && params[j].sVal() != "NONE")
      {
        if (i == 0)
        {
          tmod = params[j].sVal();
        }
        else
        {
          if (tmod != params[j].sVal())
          {
            string msg = "Inconsistent or missing TEMPMODEL for model: "+getName();
            std::ostringstream oss;
            oss << "Error in " << netlistLocation() << "\n" << msg;
            N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, oss.str());
          }
        }
      }
      if (params[j].tag() == "DOSEMODEL" && params[j].sVal() != "NONE")
      {
        if (i == 0)
        {
          dmod = params[j].sVal();
        }
        else
        {
          if (dmod != params[j].sVal())
          {
            string msg = "Inconsistent or missing DOSEMODEL for model: "+getName();
            std::ostringstream oss;
            oss << "Error in " << netlistLocation() << "\n" << msg;
            N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, oss.str());
          }
        }
      }
    }
  }

  if (tmod == "" && dmod == "")
  {
    string msg = "Neither TEMPMODEL or DOSEMODEL specified for model: "+getName();
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, oss.str());
  }
  else if (tmod != "" && dmod != "")
  {
    string msg = "Both TEMPMODEL and DOSEMODEL specified for model: "+getName();
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, oss.str());
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
    if (m_start.size() > 4)
    {
      string msg =
        "Mode than three model specifications given for model: "+getName();
      std::ostringstream oss;
      oss << "Error in " << netlistLocation() << "\n" << msg;
      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, oss.str());
    }
    fit.resize(3);
  }
  else if (modName == "PWL")
  {
    iMethod = PWL;
    fit.resize(m_start.size()-1);
  }
  else
  {
    string msg =
      "Only QUADRATIC or PWL interpolation method is supported for model: "+getName();
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, oss.str());
  }

   // First find the params that vary between the specified models:
  map<string,double> basePars;
  vector<double> t;
  string par;

  for (i=0 ; i<m_start.size()-1 ; ++i)
  {
    for (j=m_start[i] ; j<m_start[i+1]-1 ; ++j)
    {
      par = params[j].tag();
      ParameterMap::const_iterator parIt = (*getPMap()).find(par);
      if (parIt == (*getPMap()).end())
	throw std::runtime_error(string("Parameter ") + par + " not found");

      {
	const DeviceEntity::Pars &nPar = static_cast<const DeviceEntity::Pars &>(*(*parIt).second);
	if (params[j].given() && nPar.isType<double>())
	{
	  if (iModel == TEMP)
	  {
	    if (par == "TNOM")
	    {
	      t.push_back(params[j].dVal());
	    }
	  }
	  else if (iModel == DOSE)
	  {
	    if (par == "DOSE")
	    {
	      t.push_back(params[j].dVal());
	    }
	  }
	  if (i == 0)
	  {
	    if (params[j].getType() != EXPR)
	    {
	      basePars[par] = params[j].dVal();
	    }
	  }
	  else
	  {
	    if (ptype[par] == EXPR)
	    {
	      string msg = "Non-constant expression for parameter: ";
	      msg += par + " not interpolated in model: "+getName();
	      N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_WARNING, msg);
	    }
	    else
	    {
	      if (basePars.find(par) == basePars.end())
	      {
		string msg = "Unknown param: " + params[j].tag() + " in model: ";
		msg += getName() + " temperature compensation .model statement";
                std::ostringstream oss;
                oss << "Error in " << netlistLocation() << "\n" << msg;
                N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL, oss.str());
	      }
	      if (basePars[par] != params[j].dVal() && params[j].given())
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
  }

  if (t.size() != m_start.size()-1)
  {
    string msg;
    if (iModel == TEMP)
    {
      msg = "TNOM";
    }
    else if (iModel == DOSE)
    {
      msg = "DOSE";
    }
    msg += " not specified in all .model statements for model: "+getName();
    std::ostringstream oss;
    oss << "Error in " << netlistLocation() << "\n" << msg;
    N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, oss.str());
  }

  for (i=1 ; i<t.size() ; ++i)
  {
    for (j=0 ; j<i ; ++j)
    {
      if (t[i] == t[j])
      {
        string msg = "identical ";
        if (iModel == TEMP)
        {
          msg += "TNOM";
        }
        else if (iModel == DOSE)
        {
          msg += "DOSE";
        }
        msg += " values in .model statements for model: "+getName();
        std::ostringstream oss;
        oss << "Error in " << netlistLocation() << "\n" << msg;
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, oss.str());
      }
    }
  }

   // Now, collect the values to use for the fits:
  int nFit = fitMap.size();
			   //int nSet = t.size();
  vector<vector<double> > temp(nFit);
  vector<vector<double> > val(nFit);
  oldParams.resize(nFit);

  for (i=1 ; i<fit.size() ; ++i)
  {
    fit[i].resize(nFit);
  }
  min_par.resize(nFit);
  max_par.resize(nFit);
  parType.resize(nFit);
  map<string, int>::iterator fm = fitMap.begin();
  map<string, int>::iterator fm_end = fitMap.end();

  for ( ; fm != fm_end ; ++fm)
  {
    par = (*fm).first;
    DeviceEntity::ParameterMap::const_iterator parIt = (*getPMap()).find(par);
    if (parIt == (*getPMap()).end())
      throw std::runtime_error(string("Parameter ") + par + " not found");
    {
      const DeviceEntity::Pars &nPar = static_cast<const DeviceEntity::Pars &>(*(*parIt).second);

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

  for (i=0 ; i<nFit ; ++i)
  {
    temp[i].push_back(t[0]);
    val[i].push_back(fit[0][i]);
    if (parType[i] == LOG_FIT)
    {
      if (val[i][0] <= 0)
      {
        string msg = "Negative parameter value for log interpolation based parameter: ";
        fm = fitMap.begin();
        for ( ; fm != fm_end ; ++fm)
        {
          if ((*fm).second == i)
          {
            par = (*fm).first;
            break;
          }
        }
        msg += par;
        std::ostringstream oss;
        oss << "Error in " << netlistLocation() << "\n" << msg;
        N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, oss.str());
      }
      val[i][0] = log(val[i][0]);
    }
    min_par[i] = val[i][0];
    max_par[i] = val[i][0];
  }

  double p;
  for (i=1 ; i<m_start.size()-1 ; ++i)
  {
    for (j=m_start[i] ; j<m_start[i+1]-1 ; ++j)
    {
      if (params[j].getType() == DBLE && params[j].given())
      {
        par = params[j].tag();
        map<string, int>::iterator fm1 = fitMap.find(par);
        if (fm1 != fitMap.end())
        {
          k = fm1->second;
          temp[k].push_back(t[i]);
          p = params[j].dVal();
          if (parType[k] == LOG_FIT)
          {
            if (p <= 0)
            {
              string msg = "Negative parameter value for log interpolation based parameter: ";
              fm = fitMap.begin();
              for ( ; fm != fm_end ; ++fm)
              {
                if ((*fm).second == k)
                {
                  par = (*fm).first;
                  break;
                }
              }
              msg += par;
              std::ostringstream oss;
              oss << "Error in " << netlistLocation() << "\n" << msg;
              N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::USR_FATAL_0, oss.str());
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
      map<string, int>::iterator f;
      for (f=fitMap.begin() ; f!=fitMap.end() ; ++f)
      {
        i = (*f).second;
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
          string msg =
            "Internal error in DeviceModel, illegal number "
	    "of fit points for parameter: ";
          msg += (*f).first;
          N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0, msg);
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
      map<double,int> tOrd;
      for (i=0 ; i<nT ; ++i)
      {
        tOrd[t[i]] = 0;
      }
      i = 0;
      map<double,int>::iterator tOrd_i = tOrd.begin();
      map<double,int>::iterator tOrd_end = tOrd.end();
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
      map<string, int>::iterator f;
      map<string, int>::iterator f_end;
      vector<bool> p(nT,false);
      f=fitMap.begin();
      f_end=fitMap.end();
      for ( ; f!=f_end; ++f)
      {
        i = (*f).second;
        for (j=0 ; j<nT ; ++j)
        {
          p[j] = false;
        }
        for (j=0 ; j<temp[i].size() ; ++j)
        {
          int ind = tOrd[temp[i][j]];
          p[ind] = true;
          fit[ind][i] = val[i][j];
        }
        for (j=0 ; j<nT ; ++j)
        {
          if (!p[j])
          {
            k_lo = j;
            k_hi = j;
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
                string msg =
                  "DeviceModel::setModParams: Internal error "
                  "forming PWL interpolation";
                N_ERH_ErrorMgr::report ( N_ERH_ErrorMgr::DEV_FATAL_0, msg);
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

  params.resize(m_start[1]-1);
  setParams(params);
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
    map<string,int>::iterator fp;
    map<string,int>::iterator fm_begin=fitMap.begin();
    map<string,int>::iterator fm_end=fitMap.end();
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
