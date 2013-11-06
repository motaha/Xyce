#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <set>

using namespace std;
 
#include "N_CIR_Xyce.h"

int main() 
{
  string netlist("c7i.cir");
  N_CIR_Xyce *xycePtr;
  vector< string > dacNames;
  vector< string >::iterator d_i, d_end;
  map< string, vector< pair<double,double> >* > updates;
  map< string, vector< pair<double,double> > > dummy;
  double time, step_time, actual_step_time;
  int i;
  bool accept;

  cout << "Starting stand alone test of dcop in Xyce for mixed signal simulation" << endl;

  xycePtr = new N_CIR_Xyce();
  if (xycePtr) 
  {
    cout << "Xyce instance successfully created" << endl;
  }
     
  int iargs = 2;
  char* cargs[2]; 
  cargs[0] = (char*) malloc(5);
  strcpy(cargs[0], "Xyce");
  cargs[1] = (char*) malloc(netlist.size()+1);
  strcpy(cargs[1], netlist.c_str());
  if ( ! xycePtr->initialize(iargs, cargs) )
  { 
    cout << "Failed to initialize Xyce for netlist " << netlist << endl;
    return 2;  
  } 
  free(cargs[0]); 
  free(cargs[1]);

  // Get the names of all DACs in the analog circuit.
  dacNames.clear();
  if ( ! xycePtr->getDACDeviceNames(dacNames) )
  {
    cout << "Failed to get the DAC device names for netlist " << netlist << endl;
    return 3;
  } 

  vector< pair<double,double> > inp_update;

  inp_update.push_back(pair<double,double>(0.0, 3.3));

  updates[string("XINPUT:Y%DAC%DA_INVD0")] = &inp_update;

  if (updates.size() != dacNames.size()) 
  {
    cout << "updates size = " << updates.size() << " and number of DACs is " << dacNames.size() << endl;
  }
  d_i = dacNames.begin();
  d_end = dacNames.end();
  for ( ; d_i != d_end ; ++d_i ) 
  {
    cout << "Checking DAC named: '" << *d_i << "'" << endl;
    if (updates.find(*d_i) == updates.end())
      cout << "DAC " << *d_i << " is not initialized" << endl;
    else
      cout << "DAC " << *d_i << " is initialized" << endl;
  }

  int i1, i2;
  double transitionTime = 0;
  double transitionVal = 3.3;

  xycePtr->updateTimeVoltagePairs (updates);
  cout << "updateTimeVoltagePairs called" << endl;
  step_time = 1e-6;
  actual_step_time = 0;
  bool stat;
  dummy.clear();
  time = xycePtr->getTime();
  stat = xycePtr->provisionalStep(step_time, actual_step_time, dummy);
  xycePtr->acceptProvisionalStep();
  time = 0;
  i = 0;
  while (time < 3e-9) 
  {
    dummy.clear();
    cout << "Preparing to call provisionalStep at time = " << time << endl;
    stat = xycePtr->provisionalStep(step_time, actual_step_time, dummy);
    cout << "From provisional step, stat = " << stat << ", actual step time = " << actual_step_time << endl;
    if (stat != 1)
    {
      cout << "status from provisionalStep = fail!" << endl;
      break;
    }
#if 1
    i1 = static_cast<int>(time/1e-9);
    i2 = static_cast<int>((time+actual_step_time)/1e-9);
    if (time > transitionTime && i1 != i2) 
    {
      cout << "Rejected" << endl;
      transitionTime += 1e-9;
      transitionVal = 3.3 - transitionVal; 
      inp_update.clear();
      inp_update.push_back(pair<double,double>(transitionTime, transitionVal));
      xycePtr->updateTimeVoltagePairs (updates);
      cout << "updateTimeVoltagePairs called for transition to " << transitionVal <<
              " at time = " << transitionTime << endl;
      xycePtr->rejectProvisionalStep();
    }
    else 
    {
#endif
      xycePtr->acceptProvisionalStep();
      cout << "Accepted" << endl;
      time += actual_step_time;
#if 1
    }
#endif
  }
  cout << "Simulation Complete" << endl;
  return 0;
}

