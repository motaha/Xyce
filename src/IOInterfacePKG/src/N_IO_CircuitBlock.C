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

//-----------------------------------------------------------------------------
// Filename       : $RCSfile: N_IO_CircuitBlock.C,v $
//
// Purpose        : Define the circuit level containers for holding netlist
//                  circuit data and the associated circuit level methods.
//
// Special Notes  :
//
// Creator        : Lon Waters, SNL
//
// Creation Date  : 09/06/2001
//
// Revision Information:
// ---------------------
//
// Revision Number: $Revision: 1.279.2.2 $
//
// Revision Date  : $Date: 2014/03/10 16:15:21 $
//
// Current Owner  : $Author: dgbaur $
//-----------------------------------------------------------------------------

#include <Xyce_config.h>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <functional>

#include <N_IO_CircuitBlock.h>
#include <N_IO_DeviceBlock.h>
#include <N_IO_OutputMgr.h>
#include <N_IO_ParameterBlock.h>
#include <N_DEV_DeviceInterface.h>
#include <N_IO_DistributionTool.h>
#include <N_IO_CmdParse.h>

namespace Xyce {
namespace IO {

//-----------------------------------------------------------------------------
// Class         : CircuitBlockData
// Purpose       : Chesire cat pattern for CircuitBlock, this class
//                 contains the private data and methods of CircuitBlock.
//
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : ??
//-----------------------------------------------------------------------------
class CircuitBlockData
{
  public:

  //Original constructor
    CircuitBlockData(
        CmdParse & cp,
        CircuitContext & cc,
        CircuitMetadata & md,
        std::map<std::string,int> & mn,
        std::map<std::string,FileSSFPair> & ssfm,
        bool & gPI,
        const std::vector< std::pair< std::string, std::string> > & externalNetlistParams
        );

    // New constructor, to accept boolean remove parameters.
    CircuitBlockData(
        CmdParse & cp,
        CircuitContext & cc,
        CircuitMetadata & md,
        std::map<std::string,int> & mn,
        std::map<std::string,FileSSFPair> & ssfm,
        bool & gPI,
        const std::vector< std::pair< std::string, std::string> > & externalNetlistParams,
        bool removeCvar,
        bool removeDvar,
        bool removeIvar,
        bool removeLvar,
        bool removeMvar,
        bool removeQvar,
        bool removeRvar,
        bool removeVvar,
        bool replgndvar
        );

    // Destructor.
    ~CircuitBlockData();

    // Pointer to CircuitBlock that "owns" a given instance
    // of this class.
    CircuitBlock* circuitBlockPtr_;

    // CircuitBlock private data
    std::string title_;    // For top level circuit, given by first line of netlist.
    std::string name_;     // For subcircuits.

    std::ifstream netlistIn_;

    SpiceSeparatedFieldTool* ssfPtr_;

    std::streampos fileStartPosition_;
    std::streampos fileEndPosition_;
    int lineStartPosition_;
    int lineEndPosition_;

    CmdParse & commandLine_;

    std::map<std::string,int> & modelNames_;

    CircuitBlock* mainCircuitPtr_;

    DistributionTool* distToolPtr_;

    std::map<std::string,FileSSFPair> & ssfMap_;

    bool & globalParamsInserted_;

    std::vector< std::pair< std::string, std::string> > externalNetlistParams_;

    DeviceBlock device_;

    // These are added to check allow removal of "redundant" devices (where
    // all device nodes are the same).
    bool remove_redundant_C_;     //capacitors
    bool remove_redundant_D_;     //diodes
    bool remove_redundant_I_;     //independent current sources
    bool remove_redundant_L_;     //inductors
    bool remove_redundant_M_;     //MOSFETS
    bool remove_redundant_Q_;     //BJTs
    bool remove_redundant_R_;     //resistors
    bool remove_redundant_V_;     //independent voltage sources

    // This is added to turn on an option where we replace any occurrence of
    // "GND", "GND!", or "GROUND" with "0" (so that all four terms are
    // synonyms)
    bool replace_ground_;

    std::vector<SpiceSeparatedFieldTool::StringToken> parsedLine_;

    ParameterBlock tmpModel;

    std::vector<std::string> xmlBufSave_;

    Xyce::Topo::InsertionTool* insertionToolPtr_;
  Device::DeviceInterface* devIntPtr_;

    std::vector<std::string> includeFiles_;


    // CircuitBlock private methods.

    // Read and parse the XML based circuit metadata.
    std::vector<std::string> metadataBufs_;

    //This function preprocesses the netlist file to provide the user the
    //option of removing "redundant" devices (devices where all the nodes are
    //the same.  The info gathered here affects the phase 1 and phase 2 parse.
    bool parsePreprocess(const std::string & netlistFileName);

    //This function will reproduce a copy of the netlist file under the name
    //netlistfilename_copy.cir (to be used, in the end, to produce netlist files
    //which contain large resistors connecting "dangling" nodes to ground.
    void produceUnflattenedNetlist(const std::string & netlistFileName);

    // Handle a netlist line, determine the line type and take the
    // appropriate action.
    bool handleLinePass1( bool & result, std::map<std::string,int> & devMap, std::map<std::string,int> & par,
                          std::map<std::string,int> & fun, CircuitContext::ModelMap & modMap,
                          std::map<std::string,int> & sub, const std::string &libSelect, std::string &libInside );

    bool getLinePass2(std::vector<SpiceSeparatedFieldTool::StringToken>& line,
                      const std::string &libSelect, std::string &libInside);

    bool removeTwoTerminalDevice(const char linetype,
                                 const ExtendedString & node1,
                                 const ExtendedString & node2);

    bool removeThreeTerminalDevice(const char linetype,
                                   const ExtendedString & node1,
                                   const ExtendedString & node2,
                                   const ExtendedString & node3);

    // Handle a line for Mutual Inductance Pass
    bool getLinePassMI();

    // Handle a netlist .include or .lib line, return the include file, and lib strings.
    void handleIncludeLine(
        std::vector<SpiceSeparatedFieldTool::StringToken> const& parsedLine,
        const ExtendedString &, std::string& includefile, std::string &libSelect, std::string &libInside);

    // Handle a netlist .endl line.
    void handleEndlLine(
        std::vector<SpiceSeparatedFieldTool::StringToken> const& parsedLine,
        const std::string &libSelect, 
        std::string &libInside);

    // Post process a DC sweep if one was specified.
    bool handleDCSweep();

    // Post process a STEP sweep if one was specified.
    bool handleSTEPSweep();

    // Post process the mutual inductors in the current circuit.
    bool handleMutualInductances();

    // Post process a mutual inductor in the current circuit.
    bool handleMutualInductance( DeviceBlock & device );

    // Parse the given include file adding the contents
    // to the current CircuitBlock.
    bool parseIncludeFile(std::string const& includeFile, std::string const& libSelect,
         std::map<std::string,int> & devMap, std::map<std::string,int> & par, std::map<std::string,int> & fun,
                          CircuitContext::ModelMap & modMap, std::map<std::string,int> & sub);

    // Parse the given include file for 2nd pass
    bool parseIncludeFile2(std::string const& includeFiles,
            const std::string &libSelect);

    // Retrieve separate IC= data from line or external file and temporarily
    // store in CircuitBlock
    void handleInitCond(
     std::vector<SpiceSeparatedFieldTool::StringToken> const& parsedLine );

    // Fully parse and instantiate a single device.
    bool instantiateDevice( DeviceBlock & device,
        std::string & prefix, const std::map<std::string,std::string>& nodeMap,
        const std::string &libSelect, const std::string &libInside);

    // Expand a subcircuit instance by adding the devices and
    // device models that compose the subcircuit to the main
    // (top level) circuit. Prepend device names and nodes with
    // subcircuitPrefix.
    bool expandSubcircuitInstance(DeviceBlock & subcircuitInstance,
            const std::string &libSelect, const std::string &libInside);

    //
    bool getFileLocation( std::string const& path,
                          std::string const& file,
                          std::string& location);
    bool getFileLocation( std::string const& path,
                          std::string const& file,
                          char separator,
                          std::string& location);
  private:

    // Copy constructor.
    CircuitBlockData( CircuitBlockData const& rhsCBD );
    CircuitBlockData& operator=( const CircuitBlockData& rhsCBD );

};

// ----------   Static Initializations   ----------

//--------------------------------------------------------------------------
// Function      : CircuitBlock::CircuitBlock
// Purpose       : Constructor
// Special Notes : This constructor is the one used for subcircuiting.
//                 Changed on 10/10/2007 by KRS to accept boolean remove
//                 variables.
// Creator       : Lon Waters
// Creation Date : 09/02/2001
//--------------------------------------------------------------------------
CircuitBlock::CircuitBlock(
    std::string const& fileName,
    std::vector<SpiceSeparatedFieldTool::StringToken> const& parsedInputLine,
    CmdParse & cp,
    CircuitMetadata & md,
    std::map<std::string,int> & mn,
    std::map<std::string,FileSSFPair> & ssfm,
    CircuitContext & cc,
    OutputMgr * outputMgrPtr,
    int & uc,
    bool & gPI,
    CircuitBlock* mainCircPtr,
    DistributionTool* dtPtr,
    Xyce::Topo::InsertionTool* itPtr,
    Device::DeviceInterface* diPtr,
    std::map<std::string, RCP<Device::InstanceBlock> > & dNames,
    std::set<std::string> & nNames,
    const std::vector< std::pair< std::string, std::string> > & externalNetlistParams,
    bool removeCvar,
    bool removeDvar,
    bool removeIvar,
    bool removeLvar,
    bool removeMvar,
    bool removeQvar,
    bool removeRvar,
    bool removeVvar,
    bool replgndvar
    )
  : netlistFileName(fileName),
    parentCircuitPtr(NULL),
    metadata(md),
    circuitContext(cc),
    outputMgrPtr_(outputMgrPtr),
    useCount(uc),
    nodeNames_(nNames),
    deviceNames_(dNames),
    data_(new CircuitBlockData(cp,cc,md,mn,ssfm,gPI, externalNetlistParams,
                    removeCvar,removeDvar,
            removeIvar,removeLvar,removeMvar,
            removeQvar,removeRvar,removeVvar,
            replgndvar)),
    commandLine_(cp),
    externalNetlistParams_(externalNetlistParams),
    numXMLBufReceived(0),
    cmdChecked(false),
    netlistSave(true),
    devProcessedNumber(0)
{
  useCount = 1;
  data_->circuitBlockPtr_ = this;
  data_->parsedLine_ = parsedInputLine;
  data_->mainCircuitPtr_ = mainCircPtr;
  data_->distToolPtr_ = dtPtr;
  data_->insertionToolPtr_ = itPtr;
  data_->devIntPtr_ = diPtr;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::CircuitBlock
// Purpose       : Constructor
// Special Notes : Use of this constructor coincides with the top-level or
//                 main circuit, the pointer to the main circuit is
//                 set here.  This constructor is only called from
//                 IO_NetlistImportTool.  It is never called from inside of
//                 CircuitBlock.
//
// Creator       : Lon Waters
// Creation Date : 09/02/2001
//--------------------------------------------------------------------------
CircuitBlock::CircuitBlock(
    const std::string & netlistFileNameIn,
    CmdParse & cp,
    CircuitMetadata & md,
    std::map<std::string,int> & mn,
    std::map<std::string,FileSSFPair> & ssfm,
    CircuitContext & cc,
    OutputMgr * outputMgrPtr,
    int & uc,
    bool & gPI,
    std::map<std::string,RCP<Device::InstanceBlock> > & dNames,
    std::set<std::string> & nNames,
    const std::vector< std::pair< std::string, std::string> > & externalNetlistParams
    )
  : netlistFileName(netlistFileNameIn),
    parentCircuitPtr(NULL),
    metadata(md),
    circuitContext(cc),
    outputMgrPtr_(outputMgrPtr),
    useCount(uc),
    nodeNames_(nNames),
    deviceNames_(dNames),
    data_(new CircuitBlockData(cp,cc,md,mn,ssfm,gPI,externalNetlistParams)),
    commandLine_(cp),
    externalNetlistParams_(externalNetlistParams),
    numXMLBufReceived(0),
    cmdChecked(false),
    netlistSave(true),
    devProcessedNumber(0)
{
  useCount = 1;
  data_->circuitBlockPtr_ = this;
  data_->mainCircuitPtr_ = this;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::~CircuitBlock
// Purpose       : Destructor
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 09/02/2001
//--------------------------------------------------------------------------
CircuitBlock::~CircuitBlock()
{
  if (useCount == 1)
  {
    // The original instance of the class is being destroyed.
    std::map< std::string, CircuitBlock * >::iterator itcbt = circuitBlockTable_.begin();
    for ( ; itcbt != circuitBlockTable_.end(); ++itcbt )
    {
      delete itcbt->second;
    }
    circuitBlockTable_.clear();

    // Delete the blocks pointed to by the device, node and model tables.
    std::map<std::string, Topo::NodeDevBlock*>::iterator deviceTableIter;
    std::map<std::string, Topo::NodeDevBlock*>::iterator dtStart = deviceTable.begin();
    std::map<std::string, Topo::NodeDevBlock*>::iterator dtEnd = deviceTable.end();

    for (deviceTableIter = dtStart; deviceTableIter != dtEnd; ++deviceTableIter)
    {
      if( deviceTableIter->second != 0 ) delete deviceTableIter->second;
    }

    std::map<std::string, Topo::NodeBlock*>::iterator nodeTableIter;
    std::map<std::string, Topo::NodeBlock*>::iterator ntStart = nodeTable.begin();
    std::map<std::string, Topo::NodeBlock*>::iterator ntEnd = nodeTable.end();

    for (nodeTableIter = ntStart; nodeTableIter != ntEnd; ++nodeTableIter)
    {
      if( nodeTableIter->second != 0 ) delete nodeTableIter->second;
    }

    std::map<std::string, Device::ModelBlock*>::iterator modelTableIter;
    std::map<std::string, Device::ModelBlock*>::iterator mtStart = modelTable.begin();
    std::map<std::string, Device::ModelBlock*>::iterator mtEnd = modelTable.end();

    for (modelTableIter = mtStart; modelTableIter != mtEnd; ++modelTableIter)
    {
      if( modelTableIter->second != 0 ) delete modelTableIter->second;
    }

    if (parentCircuitPtr == NULL)
    {
      //Destroy the SSFs in data_
      std::map<std::string,FileSSFPair>::iterator iterSSF = data_->ssfMap_.begin();
      std::map<std::string,FileSSFPair>::iterator endSSF = data_->ssfMap_.end();
      for( ; iterSSF != endSSF; ++iterSSF )
      {
        if( iterSSF->second.first != 0 )
        {
          iterSSF->second.first->close();
          delete iterSSF->second.first;
        }
        if( iterSSF->second.second != 0 )
          delete iterSSF->second.second;
      }
    }
  }
  else
  {
    // A copy of the original is being destroyed.
    useCount--;
  }

  // Delete the private data (hidden by Chesire cat pattern).
  delete data_;
  data_ = NULL;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::CircuitBlockData
// Purpose       : Constructor
// Special Notes : Modified on 10/10/2007 by KRS to default boolean remove
//                 parameters to false.
// Creator       : Lon Waters
// Creation Date : 09/02/2001
//--------------------------------------------------------------------------
CircuitBlockData::CircuitBlockData(
      CmdParse & cp,
      CircuitContext & cc,
      CircuitMetadata & md,
      std::map<std::string,int> & mn,
      std::map<std::string,FileSSFPair> & ssfm,
      bool & gPI,
      const std::vector< std::pair< std::string, std::string> > & externalNetlistParams )
    : title_(""),
      name_(""),
      ssfPtr_(0),
      fileStartPosition_(0),
      fileEndPosition_(0),
      lineStartPosition_(1),
      lineEndPosition_(1),
      commandLine_(cp),
      modelNames_(mn),
      mainCircuitPtr_(NULL),
      distToolPtr_(NULL),
      ssfMap_(ssfm),
      globalParamsInserted_(gPI),
      externalNetlistParams_(externalNetlistParams),
      device_(cc,md),
      remove_redundant_C_(false),
      remove_redundant_D_(false),
      remove_redundant_I_(false),
      remove_redundant_L_(false),
      remove_redundant_M_(false),
      remove_redundant_Q_(false),
      remove_redundant_R_(false),
      remove_redundant_V_(false),
      replace_ground_(false)

{
}


//--------------------------------------------------------------------------
// Function      : CircuitBlockData::CircuitBlockData
// Purpose       : Constructor which explicitly takes in boolean remove
//                 arguments
// Special Notes :
// Creator       : Keith Santarelli
// Creation Date : 10/10/2007
//--------------------------------------------------------------------------
CircuitBlockData::CircuitBlockData(
      CmdParse & cp,
      CircuitContext & cc,
      CircuitMetadata & md,
      std::map<std::string,int> & mn,
      std::map<std::string,FileSSFPair> & ssfm,
      bool & gPI,
      const std::vector< std::pair< std::string, std::string> > & externalNetlistParams,
      bool removeCvar,
      bool removeDvar,
      bool removeIvar,
      bool removeLvar,
      bool removeMvar,
      bool removeQvar,
      bool removeRvar,
      bool removeVvar,
      bool replgndvar)
    : title_(""),
      name_(""),
      ssfPtr_(0),
      fileStartPosition_(0),
      fileEndPosition_(0),
      lineStartPosition_(1),
      lineEndPosition_(1),
      commandLine_(cp),
      modelNames_(mn),
      mainCircuitPtr_(NULL),
      distToolPtr_(NULL),
      ssfMap_(ssfm),
      globalParamsInserted_(gPI),
      externalNetlistParams_( externalNetlistParams ),
      device_(cc,md),
      remove_redundant_C_(removeCvar),
      remove_redundant_D_(removeDvar),
      remove_redundant_I_(removeIvar),
      remove_redundant_L_(removeLvar),
      remove_redundant_M_(removeMvar),
      remove_redundant_Q_(removeQvar),
      remove_redundant_R_(removeRvar),
      remove_redundant_V_(removeVvar),
      replace_ground_(replgndvar)
{
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::~CircuitBlockData
// Purpose       : Destructor
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 09/02/2001
//--------------------------------------------------------------------------
CircuitBlockData::~CircuitBlockData( )
{
  netlistIn_.close();
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::getName
// Purpose       : Get the subcircuit name.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 02/05/2002
//--------------------------------------------------------------------------
std::string const& CircuitBlock::getName() const
{
  return data_->name_;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::registerDistributionTool
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 07/23/2003
//----------------------------------------------------------------------------
void CircuitBlock::registerDistributionTool(DistributionTool* dtPtr)
{
  if( dtPtr == NULL )
  {
    Report::DevelFatal0() << "Distribution Tool failed to register with CircuitBlock";
  }
  else
  {
    data_->distToolPtr_ = dtPtr;
  }
}


//----------------------------------------------------------------------------
// Function       : CircuitBlock::registerInsertionTool
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/28/2003
//----------------------------------------------------------------------------
void CircuitBlock::registerInsertionTool(
    Xyce::Topo::InsertionTool* insertionToolPtr)
{
  if( insertionToolPtr == NULL )
  {
    Report::DevelFatal0() << "Insertion Tool failed to register with CircuitBlock";
  }
  else
  {
    data_->insertionToolPtr_ = insertionToolPtr;
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::registerDeviceInterface
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/29/2003
//----------------------------------------------------------------------------
void CircuitBlock::registerDeviceInterface(
  Device::DeviceInterface* devIntPtr)
{
  if( devIntPtr == NULL )
  {
    Report::DevelFatal0() << "Device Interface failed to register with CircuitBlock";
  }
  else
  {
    data_->devIntPtr_ = devIntPtr;
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitBlockData::getFileLocation
// Purpose        :
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 09/26/2003
//----------------------------------------------------------------------------
bool CircuitBlockData::getFileLocation( std::string const& path,
                                             std::string const& file,
                                             std::string& location )
{
#ifdef WIN32
  return getFileLocation(path,file,';',location);
#else
  return getFileLocation(path,file,':',location);
#endif
}

bool CircuitBlockData::getFileLocation( std::string const& path,
                                             std::string const& file,
                                             char separator,
                                             std::string& location)
{
  location = "";

  bool failure = false;
  bool success = true;

  std::list<char const*> dirList;
  // List of directories in the path.

  // Construct a list of directory names.
  char const* start = path.c_str();
  char* end = NULL;

  while ( ( end = const_cast<char*>(strchr(start, separator)) ) != NULL )
  {
    *end = '\0';
    dirList.push_back(start);
    start = end+1;
  }
  dirList.push_back(start);

  // see if we're dealing with a full path already.

  // look for current directory file as well.

  FILE* in = fopen(file.c_str(), "r");
  if ( in != NULL )
  {
    fclose(in);
    location = "";
    return success;
  }

  // Extract the full file name.
  std::list<char const*>::const_iterator iter = dirList.begin();
  std::list<char const*>::const_iterator iter_end = dirList.end();
  for ( iter = dirList.begin(); iter != iter_end;
        ++iter )
  {
    std::string dirname(*iter);
    std::string fullname(dirname);
    fullname += "/";
    fullname += file;

    // Try to open the file for reading. If the file can be opened
    // for reading, the file exists.

    // It is possible that the file exists and is not readable.
    // That case can be taken care of only through other system
    // specific calls. For the sake portability, this case can
    // be ignored for the time being.

    FILE* in = fopen(fullname.c_str(), "r");
    if ( in != NULL )
    {
      location = dirname;
      location += "/";
      fclose(in);
      return success;
    }
  }

  return failure;
}


//--------------------------------------------------------------------------
// Function      : CircuitBlock::readline
// Purpose       : Line-terminator-agnostic istream::getline() workalike for
//               : reading a single line from input stream
//
// Special Notes : file ptr is advanced
//               : line is truncated at 256
//               : line terminator is extracted but not stored
//               : null char is appended
//--------------------------------------------------------------------------
void CircuitBlock::readline( std::istream & in, char * line )
{
  int pos = 0;

  // read up to 256 chars into line
  while( in.good() && pos < 256 )
  {
    // read a char and append to line
    in.get( line[pos] );

    if( line[pos] == '\r' )
    {
      if( in.peek() == '\n' )
      {
        // extract the CR LF pair and discard
        in.get();
      }

      // flag end of line
      line[pos] = '\n';
    }

    if( in.eof() || line[pos] == '\n' )
    {
      // replace the eol/eof with the null char
      line[pos] = '\0';

      // force while exit
      pos = 256;
    }

    // append to string
    ++pos;
  }
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::parseNetlistFilePass1
// Purpose       : Top level entry for parsing netlist This is where the
//                 library context is initialized.
//
// Special Notes :
//
// Creator       : Dave Shirley, PSSI
//
// Creation Date : 10/08/2009
//--------------------------------------------------------------------------
bool CircuitBlock::parseNetlistFilePass1( )
{
  std::string libSelect, libInside;
  return parseNetlistFilePass1(libSelect, libInside);
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::parseNetlistFilePass1
// Purpose       : Parse the netlist file. Netlists parsing is a two-phase
//                 operation since the models are needed when determining
//                 how to handle the device lines. In the first phase,
//                 the netlist is read, the type of each line is determined,
//                 and an object is instantiated corresponding to the type
//                 of the line. If the line is not a device line it can
//                 be fully treated in this phase. Device lines are completed
//                 in the second phase.
//
// Special Notes :
//
// Creator       : Lon Waters
//
// Creation Date : 09/02/2001
//--------------------------------------------------------------------------
bool CircuitBlock::parseNetlistFilePass1( const std::string &libSelect, std::string libInside )
{
  CircuitContext *myContext = circuitContext.getCurrentContextPtr();
  bool result = true;

  // If this is the parent circuit, open the netlist file register an
  // instance of SpiceSeparatedFieldTool for the netlist input stream and
  // get the title line. Note: this instance of SpiceSeparatedFieldTool
  // will be shared by subcircuits (if any) that recursively call this
  // method.

  if ( parentCircuitPtr == NULL )
  {
    // Open the netlist file.  Using binary to avoid issues with compiler/plat
    // *fstream differences in implementation
    data_->netlistIn_.open( netlistFileName.c_str(), std::ios::in | std::ios::binary );

    if ( !data_->netlistIn_.is_open() )
    {
      Report::UserError0() << "Could not find netlist file " << netlistFileName;
      return false;
    }

    data_->ssfPtr_ = new
      SpiceSeparatedFieldTool(data_->netlistIn_, netlistFileName, externalNetlistParams_);
    data_->ssfMap_[netlistFileName] = FileSSFPair((std::ifstream*)0, data_->ssfPtr_);

    char ch_input[256];
    readline( data_->netlistIn_, ch_input );
    data_->title_ = ch_input;
    outputMgrPtr_->setTitle(data_->title_);

    // Increment the line number in the SSF object for the netlist to
    // account for the title line.
    data_->ssfPtr_ ->changeCursorLineNumber(1);

    //Added 10/8/07, KRS.  Adding a preprocessing phase to determine which, if
    //any, devices will be selected to remove "redundancy" (i.e., if "resistor"
    //is selected as an option in a .PREPROCESS statement, Xyce will ignore
    //all resistors for which both nodes are the same, i.e.
    // "R1 3 3 1" will not be added to the netlist.  This first phase just
    //extracts the parameters (diodes, capacitors, resistors, etc.) for which
    //we want Xyce to eliminate redundancy.
    //
    // 12/10/07 KRS:  Also using this preprocess phase to detect if we want to
    // add resistors of specified value between ground and "dangling" nodes
    data_->parsePreprocess(netlistFileName);

    //Now reset the location/line number (parsepreprocess brings us to end of
    //file).
    data_->ssfPtr_->setLocation(data_->fileStartPosition_);
    data_->ssfPtr_->setLineNumber( data_->lineStartPosition_ );
    data_->netlistIn_.clear();
    data_->netlistIn_.seekg(0, std::ios::beg);
    readline( data_->netlistIn_, ch_input );
    data_->ssfPtr_->changeCursorLineNumber( 1 );

  }

  if (DEBUG_IO) {
    if ( parentCircuitPtr == NULL )
    {
      Xyce::dout() << "Pass 1 parsing of netlist file: " << netlistFileName << std::endl;
    }
  }

  std::map<std::string,int> par, fun, sub;
  for (;;) {
    bool line_parsed = true;

    if (data_->handleLinePass1( line_parsed, myContext->devMap, par, fun, myContext->modMap, sub, libSelect, libInside ) )
      result = result && line_parsed;
    else
      break;
  }
  if (!result)
    return result;

  // if K lines found, collect coupled inductance data for parsing later
  if( ( parentCircuitPtr == NULL ) && ( !(rawMIs.empty()) ) )
  {
    std::multimap< CircuitContext *, DeviceBlock >::iterator mm =
     rawMIs.begin();

    if (DEBUG_IO)
      Xyce::dout() << "Total K lines found:  " << rawMIs.size() << std::endl;

    for( ; mm != rawMIs.end(); ++mm )
    {
      circuitContext.setContext( ( *mm ).first );

      //The mutual inductance might be an expression involving parameters,
      // which are normally handled in pass 2.  We therefore have to do
      // a special resolve before extractData for K lines?
      std::vector<Device::Param> junkSubcircuitParams;


      circuitContext.resolve(junkSubcircuitParams);

      // Parses the K line
      //( ( *mm ).second ).extractData( &circuitContext, &metadata );
      ( ( *mm ).second ).extractData( );

      // Add mutual inductance to circuit context
      circuitContext.addMutualInductance( ( *mm ).second );

      if (DEBUG_IO)
        Xyce::dout() << "In Pass 1:  adding: "
                     << ( ( *mm ).second ).getName() << " with model "
                     << ( ( *mm ).second ).getModelName() << " to "
                     << circuitContext.getCurrentContextName() << std::endl;

      circuitContext.restorePreviousContext();
    }

    // temporary blocks no longer needed; free memory
    rawMIs.clear();
  }

  myContext->modMap.clear();
  myContext->devMap.clear();

  if ( parentCircuitPtr == NULL )
  {
    // Finish resolving Mutual Inductances if necessary
    // This means finding all the assoc. inductors and including their inductance in MI.

    // Note that the parsing of mutual inductances includes converting coupled sets of
    // inductors into devices like Y%MIL% (linear coupling) and Y%MIN% (non-linear
    // coupling) that contain the coupling plus the inductors.  This creates an extra
    // pass between pass1 and pass2 which occurs on proc zero in parallel.

    if( circuitContext.totalMutualInductanceCount() )
    {
      parseMutualInductances();
    }

    // Get the total device count for the circuit and set its value in the
    // distribution tool.

    int count = circuitContext.getTotalDeviceCount();
    data_->distToolPtr_->deviceCount( count );

    // Broadcast the circuit context to all procs.
    data_->distToolPtr_->circuitContext(&circuitContext);

    // Resolve current context parameters.
    std::vector<Device::Param> params;
    result = circuitContext.resolve(params);
    if (!result)
return result;

    // resolve any functions and parameters in expression on the print line
    // at this point (or for that matter any functions/parameters in
    // the optionsTable data (so ".OP" ".OPTIONS" ".DCOP"  ".OUTPUT"
    // ".PRINT"".TRAN"  ".STEP"  ".RESULT" ".OBJECTIVE" ".IC" ".DCVOLT"
    //  ".NODESET" ".SAVE" ".LOAD" ".MPDE" ".HB"  ".AC" ".MEASURE" ".MEAS")
    // This could happen later, as in the actual classes that handle the above
    // functions, however at this stage we have all the contextual information
    // to resolve this without duplicating code elsewhere.
    resolveExpressionsInOptionBlocks();

    // Broadcast the circuit options to all procs.
    data_->distToolPtr_->circuitOptions(optionsTable);
  }

  if (DEBUG_IO) {
    if ( parentCircuitPtr == NULL )
    {
      Xyce::dout() << "Done with pass 1 netlist file parsing" << std::endl;
    }
  }

  return true; // Only get here on success.
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::parseNetlistFilePass2
// Purpose        : Second pass over the netlist. This phase primarily
//                  handles devices.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 02/21/2003
//----------------------------------------------------------------------------
bool CircuitBlock::parseNetlistFilePass2()
{
  bool result = false;

  if (DEBUG_IO) {
    if (parentCircuitPtr == NULL)
    {
      Xyce::dout() << "Pass 2 parsing of netlist file: " <<
        netlistFileName << std::endl;
    }
  }

  data_->distToolPtr_->setFileName(netlistFileName);
  // Set the start location of the circuit or subcircuit in its
  // associated file.
  data_->ssfPtr_->setLocation(data_->fileStartPosition_);
  data_->ssfPtr_->setLineNumber( data_->lineStartPosition_ );

  // If this is the main circuit, skip over the title line and continue.
  if ( parentCircuitPtr == NULL )
  {
    std::vector<SpiceSeparatedFieldTool::StringToken> line;
    char ch_input[256];
    data_->netlistIn_.clear();
    data_->netlistIn_.seekg(0, std::ios::beg);
    readline( data_->netlistIn_, ch_input );
    data_->ssfPtr_->changeCursorLineNumber( 1 );
  }
  else
  {
    Report::DevelFatal().in("CircuitBlock::parseNetlistFilePass2") << "Called from child context";
  }

  std::string libSelect, libInside;
  // Instantiate all devices in the current context.
  result = instantiateDevices(libSelect, libInside);
  if (!result)
    return result;

  //KRS, 7/9/08:  This next chunk of code dealing with the MIs used to be in
  //getLinePass2, but it was only encountered if a ".END" statement was
  //explicitly in the netlist file (and, hence, caused horrible failures for
  //netlists with mutual inductors that didn't have an explicit ".END"
  //statement in them.  Moved the code here so that it gets run regardless of
  //whether a ".END" statement is in the netlist file.

  // send MIs if present in top circuit
  if( data_->circuitBlockPtr_->circuitContext.haveMutualInductances() )
  {
    // Distribute each Y%MI?%name device line, end of circuit
    int n = data_->circuitBlockPtr_->circuitContext.getNumMILines();
    for( int i = 0; i < n; ++i )
    {
      // parse locally if distool does not distribute; normally the
      // CircuitBlock::instantiateDevices() performs this step
      if( !data_->distToolPtr_->circuitDeviceLine(
        data_->circuitBlockPtr_->circuitContext.getMILine( i ) ) )
      {
        data_->circuitBlockPtr_->handleDeviceLine(
         data_->circuitBlockPtr_->circuitContext.getMILine( i ), libSelect, libInside );
      }
    }
  }

  // Post-processing or any other work to finalize the circuit is done here.
  result = data_->handleDCSweep();
  if (!result)
    return result;

  result = data_->handleSTEPSweep();
  if (!result)
    return result;


  // call function to replace any aliased nodes on the print line with
  // actual node names
  result = substituteNodeAliases();

  if( result )
  {
    // if substituteNodeAliases returned ture, then it changed the nodes
    // on the .print line. So distribute them again.  Also need to
    // redistribute if we resolved functions on the print line.
    // Broadcast the circuit options to all procs.
    data_->distToolPtr_->updateCircuitOptions( optionsTable );
  }

  if (DEBUG_IO) {
    Xyce::dout() << "CircuitBlock::parseNetlistFilePass2 Node Alias list: " << std::endl;
    std::map<std::string,std::string>::iterator currentAliasIt_ = aliasNodeMap_.begin();
    std::map<std::string,std::string>::iterator endAliasIt_ = aliasNodeMap_.end();
    while( currentAliasIt_ != endAliasIt_ )
    {
      Xyce::dout() << "key = \"" << (*currentAliasIt_).first << "\" = \"" << (*currentAliasIt_).second << "\"" << std::endl;
      currentAliasIt_++;
    }

    print();
  }

  if (parentCircuitPtr == NULL)
  {
    data_->distToolPtr_->endDeviceLines();

    data_->distToolPtr_->checkNodeDevConflicts();

    // Here's where we call netlist copy stuff
    if (data_->commandLine_.getHangingResistor().getNetlistCopy())
    {
      data_->ssfPtr_->setLocation(data_->fileStartPosition_);
      data_->ssfPtr_->setLineNumber( data_->lineStartPosition_ );
      data_->netlistIn_.clear();
      data_->netlistIn_.seekg(0, std::ios::beg);
      data_->ssfPtr_->changeCursorLineNumber( 1 );
      data_->produceUnflattenedNetlist(netlistFileName);
    }

    if (DEBUG_IO)
      Xyce::dout() << "Done with pass 2 netlist file parsing" << std::endl;
  }

  return true; // Only get here on success.
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::parseMutualInductances
// Purpose        : Special Pass for Mutual Inductances
// Special Notes  :
// Scope          :
// Creator        : Rob Hoekstra
// Creation Date  : 08/27/04
//----------------------------------------------------------------------------
bool CircuitBlock::parseMutualInductances()
{
  if (DEBUG_IO) {
    if (parentCircuitPtr == NULL)
    {
      Xyce::dout() << "Pass MI parsing of netlist file: " << netlistFileName << std::endl;
    }
  }

  // Set the start location of the circuit or subcircuit in its
  // associated file.
  data_->ssfPtr_->setLocation(data_->fileStartPosition_);
  data_->ssfPtr_->setLineNumber( data_->lineStartPosition_ );

  // If this is the main circuit, skip over the title line and continue.
  if ( parentCircuitPtr == NULL )
  {
    std::vector<SpiceSeparatedFieldTool::StringToken> line;
    char ch_input[256];
    data_->netlistIn_.clear();
    data_->netlistIn_.seekg(0, std::ios::beg);
    readline( data_->netlistIn_, ch_input );
    data_->ssfPtr_->changeCursorLineNumber( 1 );
  }

  if( circuitContext.haveMutualInductances() )
  {
    while( data_->getLinePassMI() ) {}

    // retrieve tables and MI references from current circuit context
    std::vector<CircuitContext::MutualInductance> & MIs =
     circuitContext.getMutualInductances();
    std::vector< std::set< std::string > > & iTable = circuitContext.getSharedInductorTable();
    std::vector< std::vector< int > > & mTable = circuitContext.getAllIndexedMIs();
    std::set< std::string > & cTable = circuitContext.getAllCoupledInductors();
    std::map<std::string,double>::iterator nIter;
    std::map<std::string,double>::iterator nIter_end;
    int numMIs = MIs.size();
    int doneKey = 1;
    bool done = false;
    int imin = 0;


    mTable.push_back( std::vector< int >() );
    iTable.push_back( std::set< std::string >() );

    while (!done )
    {
      std::set<std::string> indUsed;

      mTable.push_back( std::vector< int >() );
      iTable.push_back( std::set< std::string >() );

      //Add the information for imin'th mutual inductor to all of the tables:
      MIs[imin].sharedKey = doneKey;
      mTable[doneKey].push_back(imin);
      nIter = MIs[imin].inductors.begin();
      nIter_end = MIs[imin].inductors.end();
      for( ; nIter != nIter_end; ++nIter )
      {
        indUsed.insert((*nIter).first);
        cTable.insert((*nIter).first);
        iTable[doneKey].insert((*nIter).first);
      }

      for( int i = imin + 1; i < numMIs; ++i )
      {
        if( MIs[i].model == "" && MIs[i].sharedKey == 0)
        {
          bool addit = false;

          //Search to see if the i'th mutual inductor contains inductors that
          //are already contained in the imin'th mutual inductor.  If so, add
          //these inductances to the imin'th mutual inductor
          nIter = MIs[i].inductors.begin();
          nIter_end = MIs[i].inductors.end();
          for( ; nIter != nIter_end; ++nIter )
          {
            if (indUsed.find((*nIter).first) != indUsed.end())
              addit = true;
          }

          if (addit)
          {
            MIs[i].sharedKey = doneKey;
            mTable[doneKey].push_back(i);
            nIter = MIs[i].inductors.begin();
            nIter_end = MIs[i].inductors.end();
            for( ; nIter != nIter_end; ++nIter )
            {
              indUsed.insert((*nIter).first);
              cTable.insert((*nIter).first);
              iTable[doneKey].insert((*nIter).first);
            }
            i = imin;
          }
        }
          //Whenever we add a new mutual inductor to the imin'th mutual
          //inductor, we have to start the loop all over again to make sure we
          //haven't missed anything.  That's why we set i=imin after we go
          //through the "addit" process, so that we start the loop again at
          //i = imin + 1.
      }

      done=true;

      //If there are no more mutual inductors with sharedKey tag = 0, then
      //we're done.  But, otherwise, we have to repeat this process.  imin is
      //now set to the first mutual inductor index which has a sharedKey tag of
      //0

      for (int i = 0; i < numMIs; i++)
      {
        if (MIs[i].sharedKey == 0)
        {
          imin = i;
          doneKey++;
          done=false;
          break;
        }
      }
    }

    //Now that all the MIs have been broken up into Y-devices, we have to
    //make a correction to the total number of devices.  Before we get to this
    //point, Xyce computes the total number of devices as simply being the
    //total number of device lines it encounters in parsing the netlist file,
    //but when we lump the coupled inductors and multiple K-devices into big
    //Y devices, the number of total devices decreases.  Essentially, what we
    //do here is the following:  we take the total device count, subtract the
    //total number of coupled inductors, subtract the total number of K-lines
    //found in the netlist file, and add on the number of Y-devices that we
    //formed in the above loop.

    int totalCoupledIs = 0;
    int ilength = iTable.size();

    for (int i=0; i < ilength; i++)
    {
      totalCoupledIs += iTable[i].size();
    }

    circuitContext.augmentTotalDeviceCount(numMIs,totalCoupledIs,doneKey);

    // package MIs for distribution later
    circuitContext.bundleMIs();
  }

  std::map<std::string, CircuitBlock*>::iterator itcbt = circuitBlockTable_.begin();
  for( ; itcbt != circuitBlockTable_.end(); ++itcbt )
  {
    CircuitBlock * subcircuitPtr = itcbt->second;

    // Locate the subcircuit in the netlist file. It can either be in
    // the file currently being read, or in a separate include file.
    if(subcircuitPtr->netlistFileName != netlistFileName)
    { // The subcircuit is in an include file.
      // Get SSF from Pass 1's ssf map
      if( data_->ssfMap_.count( subcircuitPtr->netlistFileName ) )
        subcircuitPtr->setSSFPtr( data_->ssfMap_[subcircuitPtr->netlistFileName].second );
      else
      {
        Report::UserError() << "Can't find include file " << subcircuitPtr->netlistFileName;
        return false;
      }
    }
    else
      subcircuitPtr->setSSFPtr( data_->ssfPtr_ );

    // Set the position of the subcircuit in its file.
    subcircuitPtr->setFilePosition(subcircuitPtr->getStartPosition());
    subcircuitPtr->setLinePosition( subcircuitPtr->getLineStartPosition() );

    // switch context
    circuitContext.setContext( subcircuitPtr->getName() );

    // Parse subckt's MIs
    subcircuitPtr->parseMutualInductances();

    // restore context
    circuitContext.restorePreviousContext();
  }

  if (DEBUG_IO) {
    print();

    if (parentCircuitPtr == NULL)
    {
      Xyce::dout() << "Done with pass MI netlist file parsing" << std::endl;
    }
  }

  return true; // Only get here on success.
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::print
// Purpose       : Print the circuit.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 09/08/2001
//--------------------------------------------------------------------------
void CircuitBlock::print()
{
  Xyce::dout() << std::endl;
  Xyce::dout() << std::endl;
  Xyce::dout() << std::endl << Xyce::section_divider << std::endl;
  Xyce::dout() << "CircuitBlock::print" << std::endl;
  Xyce::dout() << "Circuit Title: " << data_->title_ << std::endl;
  Xyce::dout() << "Circuit Name:  " << data_->name_ << std::endl;
  Xyce::dout() << std::endl;

  if ( !nodeList_.empty() )
  {
    Xyce::dout() << "Subcircuit nodes:";
    int numNodes = nodeList_.size();
    for ( int i = 0; i < numNodes; ++i )
    {
      Xyce::dout() << " " << nodeList_[i];
    }

    Xyce::dout() << std::endl;
    Xyce::dout() << std::endl;
  }

  Xyce::dout() << "Circuit Devices:" << std::endl;
  int numDevices = deviceList_.size();
  for ( int i = 0; i < numDevices; ++i )
  {
    deviceList_[i].print();
  }

  int numMutualInductorDevices = mutualInductors_.size();
  for ( int i = 0; i < numMutualInductorDevices; ++i )
  {
    mutualInductors_[i].print();
  }

  Xyce::dout() << "End Circuit Devices" << std::endl;

  Xyce::dout() << std::endl;

  //if ( !modelList.empty() )
  //{
    //Xyce::dout() << "Circuit Models:" << std::endl;
    //int numModels = modelList.size();
    //for ( int i = 0; i < numModels; ++i )
    //{
      //modelList[i].print();
    //}
    //Xyce::dout() << "End Circuit Models" << std::endl;

    //Xyce::dout() << std::endl;
  //}

  if ( !optionsTable.empty() )
  {
    Xyce::dout() << "Options: " << std::endl;
    std::list<Util::OptionBlock>::iterator optionIter = optionsTable.begin();
    std::list<Util::OptionBlock>::iterator optionIterEnd = optionsTable.end();
    for ( ; optionIter != optionIterEnd; ++optionIter )
    {
      Xyce::dout() << std::endl
                   << "Option Information" << std::endl
                   << "------------------" << std::endl
                   << std::endl
                   << "  name: " << optionIter->getName() << std::endl;

      Xyce::dout() << "  parameters: " << std::endl;
      std::list< Util::Param >::const_iterator paramIter = optionIter->getParams().begin();
      std::list< Util::Param >::const_iterator paramIterEnd = optionIter->getParams().end();
      for ( ; paramIter != paramIterEnd; ++paramIter )
      {
        Xyce::dout() << "  " << paramIter->tag() << "  ";
        Xyce::dout() << paramIter->stringValue() << std::endl;
      }
    }
    Xyce::dout() << std::endl << std::endl;
  }

  //if ( netlistParameters.getNumberOfParameters() > 0 )
  //{
    //Xyce::dout() << "Circuit Parameters:" << std::endl;
    //netlistParameters.print();
    //Xyce::dout() << std::endl;
  //}

  //if ( !params.empty() )
  //{
    //int numOptionsBlocks = params.size();
    //for ( int i = 0; i < numOptionsBlocks; ++i )
    //{
      //params[i].print();
    //}
    //Xyce::dout() << std::endl;
  //}


  //if ( !functions.empty() )
  //{
    //Xyce::dout() << "Circuit Functions:" << std::endl;
    //Xyce::dout() << "------------------" << std::endl;
    //int numFunctionBlocks = functions.size();
    //for ( int i = 0; i < numFunctionBlocks; ++i )
    //{
      //functions[i].print();
    //}
    //Xyce::dout() << "------------------" << std::endl;
    //Xyce::dout() << std::endl;
    //Xyce::dout() << std::endl;
  //}

  if ( !circuitBlockTable_.empty() )
  {
    Xyce::dout() << "Subcircuits: " << std::endl;
    std::map< std::string, CircuitBlock * >::iterator itcbt = circuitBlockTable_.begin();
    for ( ; itcbt != circuitBlockTable_.end(); ++itcbt )
    {
      itcbt->second->print();
    }
    Xyce::dout() << "End Subcircuits" << std::endl;

    Xyce::dout() << std::endl;
  }
  Xyce::dout() << std::endl << Xyce::section_divider << std::endl;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::setSSFPtr
// Purpose       :
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 09/21/2001
//--------------------------------------------------------------------------
void CircuitBlock::setSSFPtr( SpiceSeparatedFieldTool* ssfPtr )
{
  data_->ssfPtr_ = ssfPtr;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::setStartPosition
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 05/19/2003
//----------------------------------------------------------------------------
void CircuitBlock::setStartPosition()
{
  data_->fileStartPosition_ = data_->ssfPtr_->getFilePosition();
  data_->lineStartPosition_ = data_->ssfPtr_->getLineNumber();

  if (DEBUG_IO) {
    Xyce::dout() << "CircuitBlock::setStartPosition being called for file "
                 << data_->ssfPtr_->getFileName() << std::endl
                 << "  start position  = " <<  data_->fileStartPosition_ << std::endl
                 << "  start line  = " << data_->lineStartPosition_ << std::endl;
  }
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::setEndPosition
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/25/2003
//----------------------------------------------------------------------------
void CircuitBlock::setEndPosition()
{
  data_->fileEndPosition_ = data_->ssfPtr_->getFilePosition();
  data_->lineEndPosition_ = data_->ssfPtr_->getLineNumber();

  if (DEBUG_IO)
    Xyce::dout() << "CircuitBlock::setEndPosition:" << std::endl
                 << "  Setting file end position: " << data_->fileEndPosition_ << std::endl
                 << "  setting line end position: " << data_->lineEndPosition_ << std::endl;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::setFilePosition
// Purpose        : Set the location in the input file to the given position.
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/21/2003
//----------------------------------------------------------------------------
void CircuitBlock::setFilePosition(std::streampos const& position)
{

  if (DEBUG_IO)
    Xyce::dout() << "CircuitBlock::setFilePosition: Setting file position to "
                 << position <<std::endl;

  data_->ssfPtr_->setLocation(position);
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::setLinePosition
// Purpose        : Set the location of the currentLine counter in the ssfPtr
// Special Notes  :
// Scope          : public
// Creator        : Eric Rankin
// Creation Date  : 10/13/2004
//----------------------------------------------------------------------------
void CircuitBlock::setLinePosition( int const& position )
{
  data_->ssfPtr_->setLineNumber( position );
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::getStartPosition
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/25/2003
//----------------------------------------------------------------------------
const std::streampos CircuitBlock::getStartPosition() const
{
  return data_->fileStartPosition_;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::getLineStartPosition
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Eric Rankin
// Creation Date  : 10/13/2004
//----------------------------------------------------------------------------
int CircuitBlock::getLineStartPosition() const
{
  return data_->lineStartPosition_;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::getEndPosition
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/25/2003
//----------------------------------------------------------------------------
const std::streampos CircuitBlock::getEndPosition() const
{
  return data_->fileEndPosition_;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::getLineEndPosition
// Purpose        :
// Special Notes  :
// Scope          : public
// Creator        : Eric Rankin
// Creation Date  : 10/13/2004
//----------------------------------------------------------------------------
int CircuitBlock::getLineEndPosition() const
{
  return data_->lineEndPosition_;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::receiveCircuitContext
// Purpose        : Receive the circuit context (from the distribution tool).
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 07/21/2003
//----------------------------------------------------------------------------
bool CircuitBlock::receiveCircuitContext(CircuitContext & ccIn)
{
  circuitContext = ccIn;
  circuitContext.resolve( std::vector<Device::Param>() );

  return true;
}


//--------------------------------------------------------------------------
// Function      : CircuitBlock::extractSubcircuitData
// Purpose       : Extract subcircuit data from parsedLine. The bulk of
//                 the subcircuit data is stored in the circuit context.
//                 A circuit block is created to represent the subcircuit
//                 but mainly exists to help with the recursive descent
//                 through the circuit during pass 2 to instantiate devices.
//                 All that is needed for this is the subcircuit name.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 09/21/2001
//--------------------------------------------------------------------------
bool CircuitBlock::extractSubcircuitData(std::string fileName, int lineNum)
{
  int numFields = data_->parsedLine_.size();

  if ( numFields < 3 )
  {
    Report::DevelFatal0().in("CircuitBlock::extractSubcircuitData").at(netlistFileName, data_->parsedLine_[0].lineNumber_)
      << "This should have been detected earlier in parse";
    return false;
  }

  // Extract the subcircuit name.
  ExtendedString field ( data_->parsedLine_[1].string_ );
  field.toUpper();
  data_->name_ = field;

  std::map<std::string,int> dupNodes;
  for (int i=2 ; i<numFields ; ++i)
  {
    field = data_->parsedLine_[i].string_;
    field.toUpper();
    if (field == "=")
    {
      break;
    }
    if (field == "PARAMS:")
      break;
    if (field == "0")
    {
      Report::UserError0().at(fileName, lineNum) << "Ground node '0' appears in .SUBCKT line";
    }
  }

  // Once data is extracted, the info in parsedLine is no longer
  // needed, trim it now using the swap trick.
  data_->parsedLine_.clear();
  std::vector<SpiceSeparatedFieldTool::StringToken>(data_->parsedLine_).swap(data_->parsedLine_);

  return true;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::instantiateDevices
// Purpose        : Instantiate the devices related to a given (sub)circuit
//                  instance.
// Special Notes  : The context for the devices to be instantiated should
//                  have been set prior to calling this method (except in
//                  the default "main circuit" context).
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/04/2003
//----------------------------------------------------------------------------
bool CircuitBlock::instantiateDevices(std::string libSelect, std::string libInside)
{
  std::vector<SpiceSeparatedFieldTool::StringToken> line;

  while (data_->getLinePass2(line, libSelect, libInside))
  {
    // parse locally if distool does not distribute
    if( !data_->distToolPtr_->circuitDeviceLine( line ) )
    {
       handleDeviceLine( line, libSelect, libInside );
    }
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::instantiateDevice
// Purpose        : Extract data from tokenized device and place all data about
//                  the device in Topology.
// Special Notes  : This method is misnamed as it does not actually instantiate
//                  the device.  That is done after partitioning (if parallel)
//                  and is called from Topology.
//                  This actually fills the DeviceBlock which is stored in
//                  Topology.
// Scope          : private
// Creator        : Lon Waters
// Creation Date  : 07/28/2003
//----------------------------------------------------------------------------
bool CircuitBlockData::instantiateDevice(
    DeviceBlock &                          device,
    std::string &                               prefix,
    const std::map<std::string,std::string> &   nodeMap,
    const std::string &                         libSelect,
    const std::string &                         libInside)
{
  bool result;
  std::map<std::string,std::string>::const_iterator nodeMapIter;
  std::map<std::string,std::string>::const_iterator nodeMapEnd = nodeMap.end();

  // Global params were discovered by proc 0 in pass 1 and stored in root level of
  // circuitContext.  Then circuitContext was transmitted to all procs.

  // Now, before device parsing, move global params from circuitContext to deviceMgr

  // This is done here because the global parameters are needed for the parsing of
  // device lines, but are not known until after pass 1.  This point is immediately
  // before device line parsing on all processors.

  if (!globalParamsInserted_)
  {
    OptionBlock *globals = circuitBlockPtr_->circuitContext.getGlobals();
    std::list<Util::Param>::iterator opt = globals->optionData.begin();
    std::list<Util::Param>::iterator optEnd = globals->optionData.end();
    for ( ; opt != optEnd; ++opt)
    {
      devIntPtr_->addGlobalPar(*opt);
    }
    globalParamsInserted_ = true;
    circuitBlockPtr_->outputMgrPtr_->registerNodeDevNames (&circuitBlockPtr_->nodeNames_,
                                                           &circuitBlockPtr_->deviceNames_);
  }

  // Fully parse the device line and extract the relevant data.
  //result = devicePtr->extractData( &circuitBlockPtr_->circuitContext,
                                   //&circuitBlockPtr_->metadata );
  result = device.extractData( );

//DNS: this is where device name mapping info would be registered with DeviceMgr

  if (!result)
return result;

  // Map device nodes.
  if (!nodeMap.empty())
  {
    std::list<tagged_param> nodeValues;
    std::list<tagged_param>::const_iterator nodeIter = device.getNodeValues().begin();
    std::list<tagged_param>::const_iterator nodeIterEnd = device.getNodeValues().end();
    for (; nodeIter != nodeIterEnd; ++nodeIter)
    {
      nodeMapIter = nodeMap.find(nodeIter->tag);
      if (nodeMapIter != nodeMap.end())
      {
        // The device node is an external node, map it using nodeMap
        nodeValues.push_back(tagged_param(nodeMapIter->second,0));
      }
      else if (nodeIter->tag != "0" && nodeIter->tag.substr(0,2) != "$G" &&
               !circuitBlockPtr_->circuitContext.globalNode(nodeIter->tag))
      {
        // The node is internal, prepend subcircuitPrefix. Note: the
        // ground node and global nodes (0 and $G*) are unchanged.
        nodeValues.push_back(tagged_param(prefix + ":" + nodeIter->tag,0));
      }
      else
      {
        nodeValues.push_back(tagged_param(nodeIter->tag,0));
      }
    }

    device.setNodeValues(nodeValues);
  }
  // If the device is not a subcircuit instance add its data to the
  // tables, otherwise, expand the instance.
  if (!device.isSubcircuitInstance())
  {
    // Map device name.
    if (prefix != "")
      device.setName(prefix + ":" + device.getName());

    // If the device has a model, find it and instantiate it (if it has not
    // already been instantiated). Prepend the model name in the device and
    // the model with the model prefix.  Add the model to the circuit.
    if (device.getModelName() != "")
    {
      std::string modelPrefix;
      ParameterBlock* modelPtr;
      bool modelFound = circuitBlockPtr_->circuitContext.findModel(
          device.getModelName(), modelPtr, modelPrefix);
      if (modelFound)
      {
        // Add the model prefix to the device's model name.
        if (modelPrefix != "")
        {
          device.setModelName(modelPrefix + ":" + device.getModelName());
        }
      }
      else
      {
        Report::UserError0() << "Unable to find model " << device.getModelName()
                             << " for device " << device.getName();
        return false;
      }

      // ERK: tmpModel is a persistant object in the CircuitBlockData class.
      // It is cheaper to use the "=" operator than the copy constructor here,
      // as the copy constructor will require allocations and deletions every
      // time this function is called, which is a lot.
      //ParameterBlock model(*modelPtr);
      tmpModel = (*modelPtr);

      // Set the model parameter values.
      tmpModel.setParameterValues(&circuitBlockPtr_->circuitContext);

      // Add the model to the circuit.
      mainCircuitPtr_->addModel(tmpModel, modelPrefix);
    }

    // If the device is a Bsource, the nodes and sources appearing in
    // the expression must be appropriately prefixed.

    // THIS IS AN AWFUL HACK!  It means that *ANY* device that has
    // "SOLN_DEP" parameters must be hacked in here or they can't be used
    // in a subcircuit!

    if ((device.getNetlistDeviceType() == "B" || device.getNetlistDeviceType() == "S" || device.getNetlistDeviceType() == "C")&& (!nodeMap.empty()))
    {
      Device::Param* givenParameter = 0;
      if (device.getNetlistDeviceType() == "B")
      {
        Device::Param* I_parameter = device.findInstanceParameter( "I" );
        Device::Param* V_parameter = device.findInstanceParameter( "V" );
        if ( I_parameter->given() )
          givenParameter = I_parameter;
        else if ( V_parameter->given() )
          givenParameter = V_parameter;
      }
      else
      {
        Device::Param* C_parameter = 0;
        if (device.getNetlistDeviceType() == "S")
        {
            C_parameter = device.findInstanceParameter( "CONTROL" );
        }
        else if (device.getNetlistDeviceType() == "C")
        {
            C_parameter = device.findInstanceParameter( "C" );
        }

        // Note!  The capacitor could be a semiconductor capacitor, and
        // therefore might not even HAVE a C parameter.  Must guard against
        // that.
        if ( C_parameter->given() )
          givenParameter = C_parameter;
        else
          givenParameter = 0;
      }


      if (givenParameter && givenParameter->getType() == Xyce::Util::EXPR)
      {
        Util::Expression &expression = givenParameter->getValue<Util::Expression>();
        if ( (expression.get_num( XEXP_NODE ) > 0) ||
            (expression.get_num( XEXP_INSTANCE ) > 0) ||
            (expression.get_num( XEXP_LEAD ) > 0) )
        {
          // If the expression has nodes or voltage source instances, get
          // the nodes and map them appropriately or add the subcircuit
          // prefix to them.
          std::vector<std::string> names;
          std::vector<std::string> nodes;
          std::vector<std::string> instances;
          std::vector<std::string> leads;
          if ( expression.get_num( XEXP_NODE ) > 0 )
          {
            expression.get_names( XEXP_NODE, nodes );
            names.insert( names.end(), nodes.begin(), nodes.end() );
          }
          if ( expression.get_num( XEXP_INSTANCE ) > 0 )
          {
            expression.get_names( XEXP_INSTANCE, instances);
            names.insert( names.end(), instances.begin(), instances.end() );
          }
          if ( expression.get_num( XEXP_LEAD ) > 0 )
          {
            expression.get_names( XEXP_LEAD, leads);
            names.insert( names.end(), leads.begin(), leads.end() );
          }

          std::vector<std::string> actualName;
          std::string newName, tmpName;
          unsigned int i;

          nodeMapEnd = nodeMap.end();
          for (i = 0; i < names.size(); ++i )
          {
            nodeMapIter = nodeMap.find(names[i]);
            if (nodeMapIter != nodeMapEnd)
            {
              newName = nodeMapIter->second;
            }
            else if (!(i < nodes.size() && ( names[i].substr(0,2) == "$G" ||
                     circuitBlockPtr_->circuitContext.globalNode(names[i]))))
            {
              newName = prefix + ":" + names[i];
            }
            else
            {
              // which leaves the one case, where we're a global node $G...,
              // in which case we just keep the original name.
              newName = names[i];
            }

            // Replace the old node name with the new node name
            // in the expression.
            if (names[i] != newName)
            {
              if (DEBUG_IO)
                Xyce::dout() << "Replacing " << names[i] << " with " << newName << std::endl;

              actualName.push_back(newName);
              tmpName = ";" + newName;
              expression.replace_name( names[i], tmpName );
              if (DEBUG_EXPRESSION) {
                Xyce::dout() << "CircuitBlock::instantiateDevice:  After this replacement, get_expression returns "
                             << expression.get_expression() << std::endl
                             << " Parse Tree: " << std::endl;
                expression.dumpParseTree();
              }
            }
          }
          for (i=0 ; i<actualName.size() ; ++i)
          {
            if (DEBUG_IO)
              Xyce::dout() << "Replacing ;"<<actualName[i]<<" with " << actualName[i] << std::endl;

            newName = actualName[i];
            tmpName = ";" + newName;
            expression.replace_name( tmpName, newName);

            if (DEBUG_EXPRESSION) {
              Xyce::dout() << "CircuitBlock::instantiateDevice:  After this replacement get_expression returns "
                           << expression.get_expression() << std::endl;
              expression.dumpParseTree();
            }
          }

          // Reset Bsource's expression.
          if (DEBUG_IO)
            Xyce::dout() << "CircuitBlock::instantiateDevice:  After all expression handling, get_expression returns "
                         << expression.get_expression() << std::endl;
        }
      }
    }

    if (device.getNetlistDeviceType() == "L" &&
     circuitBlockPtr_->circuitContext.haveMutualInductances() )
    {
      handleMutualInductance(device);
    }

    // Add device info to the main circuit device and node tables.
    circuitBlockPtr_->addTableData(device);
  }
  else
  {
    expandSubcircuitInstance(device, libSelect, libInside);
  }

  if (DEBUG_IO)
    device.print();

  return true;
}

#ifdef Xyce_PARALLEL_MPI

//----------------------------------------------------------------------------
// Function       : CircuitBlock::getDeviceNames
// Purpose        : return a vector of the device names on this processor
// Special Notes  : returns byte count for packing
// Scope          : public
// Creator        : Dave Shirley, PSSI
// Creation Date  : 01/20/2006
//----------------------------------------------------------------------------
int CircuitBlock::getDeviceNames(std::vector<std::string> & names)
{
  int byteCount = 0;

  names.clear();

  std::map<std::string,RCP<Device::InstanceBlock> >::iterator nName;
  std::map<std::string,RCP<Device::InstanceBlock> >::iterator nNameEnd;

  std::set<std::string>::iterator iterN;

  nName = deviceNames_.begin();
  nNameEnd = deviceNames_.end();

  for ( ; nName != nNameEnd ; ++nName)
  {
    names.push_back(nName->first);
    byteCount += sizeof(int) + nName->first.size();
  }

  return byteCount;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::checkDeviceNames
// Purpose        : checks a vector of device names against the local device names
// Special Notes  :
// Scope          : public
// Creator        : Dave Shirley, PSSI
// Creation Date  : 01/20/2006
//----------------------------------------------------------------------------
void CircuitBlock::checkDeviceNames(const std::vector<std::string> & names)
{
  for (std::vector<std::string>::const_iterator nName = names.begin(); nName != names.end(); ++nName)
  {
   // find duplicates across procs
    std::map<std::string,RCP<Device::InstanceBlock> >::iterator iterDN = deviceNames_.find(*nName);

    if (iterDN != deviceNames_.end())
    {
      if (Teuchos::nonnull( iterDN->second ) )
      {
        int lastColon = nName->find_last_of( ':' );
        Report::UserError().at(deviceNames_[*nName]->netlistFileName_, deviceNames_[*nName]->lineNumber_ )
          << "Duplicate device name detected: " << nName->substr( lastColon + 1 );
      }
    }
  }
}

#endif


//-----------------------------------------------------------------------------
// Function      : CircuitBlock::addTableData
// Purpose       : Add a device to the circuit.
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 02/21/2002
//-----------------------------------------------------------------------------
void CircuitBlock::addTableData( DeviceBlock & device )
{
  if (!cmdChecked)
  {
    if (commandLine_.argExists("-syntax"))
    {
      netlistSave = false;
    }
    else
    {
      netlistSave = true;
    }
    cmdChecked = true;
  }

  if (DEBUG_IO) {
    ++devProcessedNumber;
    if (devProcessedNumber%1000 == 0)
      Xyce::dout() << devProcessedNumber << "Devices processed" << std::endl;
  }

  std::string dName;

  Teuchos::RCP<Device::InstanceBlock> iBPtr = Teuchos::rcp( new Device::InstanceBlock(device.getDeviceData().getDevBlock()) );
  dName = iBPtr->getName();

  if (netlistSave)
  {
    iBPtr->numExtVars = device.getNumberOfNodes();

    if (DEBUG_IO) {
      int lastColon = dName.find_last_of( ':' );
      Xyce::dout() << "Inserting device: " << dName
                   << " locally named " <<  dName.substr( lastColon + 1 )
                   << " from file: " << iBPtr->netlistFileName_
                   << " line: " << iBPtr->lineNumber_ << std::endl;
    }

    data_->insertionToolPtr_->insertNode(device.getDeviceData().getNodeBlock(), iBPtr);
  }
  else
  {
    data_->devIntPtr_->addDeviceToCount(device.getNetlistDeviceType());

    // std::string tmpString ( device.getNetlistDeviceType() );
    // outputMgrPtr_->addDeviceToCount(tmpString);
  }

  // save device names for syntax diagnostics
  unsigned int expectedSize = deviceNames_.size() + 1;

  deviceNames_[dName] = iBPtr;

  if ( deviceNames_.size() != expectedSize )
  {
    int lastColon = dName.find_last_of( ':' );

    Report::UserError().at(iBPtr->netlistFileName_, iBPtr->lineNumber_) << "Duplicate device name detected "<< dName.substr( lastColon + 1 );
  }

  // save node names for syntax diagnostics
  std::list<tagged_param> nodeList = device.getDeviceData().getNodeBlock().get_NodeList();
  std::list<tagged_param>::iterator ii=nodeList.begin();
  std::list<tagged_param>::iterator ii_end=nodeList.end();
  for ( ; ii != ii_end; ++ii)
  {
    if (DEBUG_IO)
      Xyce::dout() << "  Node: " << ii->tag << std::endl;

    nodeNames_.insert(ii->tag);
  }

}

//-----------------------------------------------------------------------------
// Function      : CircuitBlock::addModel
// Purpose       : Add a model to the circuit.
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 02/21/2002
//-----------------------------------------------------------------------------
void CircuitBlock::addModel( ParameterBlock & model, std::string const& modelPrefix)
{
  std::string modelName(model.getName());
  if (modelPrefix != "")
  {
    modelName = modelPrefix + ":" + modelName;
    model.setName(modelName);
  }

  typedef std::map<std::string,int> mnMap;
  std::map<std::string,int>::iterator iterMN = data_->modelNames_.find(modelName);

  if ( iterMN == data_->modelNames_.end() )
  {
    data_->modelNames_.insert(mnMap::value_type(modelName,0));
    data_->devIntPtr_->addDeviceModel(model.modelData);
  }
}

//-----------------------------------------------------------------------------
// Function      : CircuitBlock::addOptions
// Purpose       : Add a set of options corresponding to a .OPTIONS netlist
//                 line to the circuit.
// Special Notes :
// Creator       : Lon Waters, SNL
// Creation Date : 02/21/2002
//-----------------------------------------------------------------------------
void CircuitBlock::addOptions( OptionBlock const& options )
{
  if ( options.getName() == "PRINT" )
  {
    int numParameters = options.getNumberOfParameters();
    for ( int i = 0; i < numParameters; ++i )
    {
      aliasNodeMapHelper_.insert( std::string( ( options.getParameter( i ) ).uTag() ) );
    }
  }

  else if ( options.getName() == "DEVICE" )
  {
    deviceOptions = options.optionData;
  }
  else if (options.getName() == "OUTPUT-LINE")
  {
    // The parameters on a .OUTPUT line need to be added to the option block
    // with the name "OUTPUT" which was created by a prior .OPTIONS OUTPUT

    // line in the netlist. Find the option block, report an error if not
    // found.
    std::list<Util::OptionBlock>::iterator optionBlockIter;
    std::list<Util::OptionBlock>::iterator first = optionsTable.begin();
    std::list<Util::OptionBlock>::iterator last = optionsTable.end();

    for (optionBlockIter = first; optionBlockIter != last; ++optionBlockIter)
    {
      if (optionBlockIter->getName() == "OUTPUT")
        break;
    }

    if (optionBlockIter == last)
    {
      // The line number of the .OUTPUT line was stored as the 3rd parameter.
      int lineNum = (options.getParameter(2)).getImmutableValue<int>();

      // Could not find required option block, report error.
      Report::UserError0().at(netlistFileName, lineNum) << "A .OPTIONS OUTPUT line is required before any .OUTPUT line in the netlist";
    }
    else
    {
      // If we get here, all is well, add the parameters.
      optionBlockIter->getParams().push_back(options.getParameter(0));
      optionBlockIter->getParams().push_back(options.getParameter(1));
    }

    return;
  }

  if (DEBUG_IO)
    options.printDiagnostic();

  optionsTable.push_back( options.optionData );
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::handleLinePass1
// Purpose       : Determine the type of netlist line in parsedLine and
//                 handle appropriately.
//
// Special Notes : Returns false on end of file
//
// Creator       : Lon Waters
//
// Creation Date : 09/08/2001
//--------------------------------------------------------------------------
bool CircuitBlockData::handleLinePass1(
  bool & result, std::map<std::string,int> & devMap,
  std::map<std::string,int> & par, std::map<std::string,int> & fun,
  CircuitContext::ModelMap & modMap, std::map<std::string,int> & sub,
  const std::string &libSelect, std::string &libInside )
{
  result = true;

  // Get the next line of input.
  std::vector<SpiceSeparatedFieldTool::StringToken> line;
  int eof = !ssfPtr_->getLine(line,replace_ground_); // Breaks the line into fields.

  if (DEBUG_IO) {
    Xyce::dout() << " After getline, file position is " << ssfPtr_->getFilePosition() << std::endl;
    Xyce::dout() << "pass 1 read netlist line: ";
    for (unsigned int i = 0; i < line.size(); ++i)
    {
      Xyce::dout() << "\"" << line[i].string_ << "\" ";
    }
    Xyce::dout() << std::endl;
  }

  // Determine what to do with the parsed line.
  if ( line.empty() )
  {
    // Blank line.
    return !eof;
  }

  char lineType;
  bool removecomponent = false;
  ExtendedString ES1 ( line[0].string_ );
  ES1.toUpper();
  if (libSelect != libInside && libInside != "" && ES1 != ".ENDL")
  {
    return !eof;
  }
  lineType = ES1[0];
  if (lineType >= 'A' && lineType <= 'Z')
  {
    if (lineType == 'C' || lineType == 'D' || lineType == 'I' ||
        lineType == 'L' || lineType == 'R' || lineType == 'V')
    {
      if (line.size() > 2) //make sure that there are two nodes to check!
      {
        ExtendedString node1 ( line[1].string_ );
        ExtendedString node2 ( line[2].string_ );
        node1.toUpper();
        node2.toUpper();

        removecomponent = removeTwoTerminalDevice(lineType, node1, node2);
      }
    }
    else if (lineType == 'M' || lineType == 'Q')
    {
      if (line.size() > 3) //make sure that there are three nodes to check!
      {
        ExtendedString node1 ( line[1].string_ );
        ExtendedString node2 ( line[2].string_ );
        ExtendedString node3 ( line[3].string_ );
        node1.toUpper();
        node2.toUpper();
        node3.toUpper();

        removecomponent = removeThreeTerminalDevice(lineType, node1, node2,
                                                    node3);
      }
    }

    if (!removecomponent)
    {
      if (lineType != 'X' && lineType != 'K')
      {
        circuitBlockPtr_->circuitContext.incrementDeviceCount();
      }
      if (lineType == 'X')
      {
        DeviceBlock device(circuitBlockPtr_->circuitContext,
                           circuitBlockPtr_->metadata,circuitBlockPtr_->netlistFileName,
                           line );

        device.extractSubcircuitInstanceData();
        circuitBlockPtr_->circuitContext.addInstance(device.getModelName(),
                                                     circuitBlockPtr_->netlistFileName,
                                                     line[0].lineNumber_);
      }
      if (lineType == 'K')
      {
        circuitBlockPtr_->circuitContext.incrementDeviceCount();
        DeviceBlock device(circuitBlockPtr_->circuitContext,
                           circuitBlockPtr_->metadata,
                           circuitBlockPtr_->netlistFileName, line );

        // save this information for later resolution
        mainCircuitPtr_->rawMIs.insert(std::pair< CircuitContext *,
                                       DeviceBlock >(circuitBlockPtr_->circuitContext.getCurrentContextPtr(), device ) );
      }

      if (lineType == 'U')
      {
        if (line.size() < 2)
        {
          Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
            << "U device line specified with only one token: " <<  ES1 << ". Need at least two.";
          result = false;
        }
        else {
          if (line.size() < 3)
          {
            Report::UserWarning().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
              << "U device line (" << ES1 << ") specified with only two tokens. Likely need more.";
          }
          ES1 = "U" + line[0].string_ + "%" + line[1].string_;
          ES1[1] = '%';
          ES1.toUpper();
        }
      }

      if (lineType == 'Y')
      {
        if (line.size() < 2)
        {
          Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
            << "Y device line specified with only one token: " <<  ES1 << ". Need at least two.";
          result = false;
        }
        else {
          if (line.size() < 3)
          {
            Report::UserWarning().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
              << "Y device line (" << ES1 << ") specified with only two tokens. Likely need more.";
          }
          ES1 = "Y" + line[0].string_ + "%" + line[1].string_;
          ES1[1] = '%';
          ES1.toUpper();
        }
      }
    }

#ifdef Xyce_DEBUG_IO
    else
    {
      Xyce::dout() << "Netlist Parse 1:  ";
      Xyce::dout() << "removing component " << ES1 << ".  All nodes on the device";
      Xyce::dout() << " are the same."  << std::endl;
    }
#endif
  }
  else if (lineType == '.')
  {
    if( ES1 == ".DC" )
    {
      // need to store an option block for each sweep variable
      std::vector< OptionBlock > optionBlocks;
      OptionBlock options( circuitBlockPtr_->netlistFileName, line,
                           circuitBlockPtr_->metadata );
      int oBs = options.extractDCData( optionBlocks );

      // create a block for each line
      for( int i = 0; i < oBs; ++i)
      {
        circuitBlockPtr_->addOptions( optionBlocks[i] );
      }
    }
    else if ( ES1 == ".OP" || ES1 == ".OPTIONS" || ES1 == ".DCOP" ||
              ES1 == ".OUTPUT" || ES1 == ".PRINT" || ES1 == ".TRAN" || ES1 == ".TR" ||
              ES1 == ".STEP" || ES1 == ".RESULT" || ES1 == ".OBJECTIVE" ||
              ES1 == ".IC" || ES1 == ".DCVOLT" || ES1 == ".NODESET" ||
              ES1 == ".SAVE" || ES1 == ".LOAD" || ES1 == ".MPDE" ||
              ES1 == ".HB" ||  ES1 == ".AC" || ES1 == ".MEASURE" || ES1 == ".MEAS" ||
              ES1 == ".MOR" || ES1 == ".FOUR" || ES1 == ".FFT" ||
              ES1 == ".SENS" )
    {
      // Create an OptionBlock for this line, extract the parameters,
      // and add the OptionBlock to the circuit.
      OptionBlock options( circuitBlockPtr_->netlistFileName, line,
                           circuitBlockPtr_->metadata );
      result = options.extractData();
      circuitBlockPtr_->addOptions( options );
    }
    else if (ES1 == ".END")
    {
      return false;
    }
    else if (ES1 == ".ENDS")
    {
      // End the current subcircuit context.
      circuitBlockPtr_->setEndPosition();
      circuitBlockPtr_->circuitContext.endSubcircuitContext();
      return false;
    }
    else if (ES1 == ".FUNC")
    {
      // Create a FunctionBlock for the .FUNC line, extract the
      // data from line, and add the function to the circuit.
      FunctionBlock function( circuitBlockPtr_->netlistFileName, line );
      result = function.extractData();
      if (result) {
        circuitBlockPtr_->circuitContext.addFunction( function );
        ExtendedString F ( line[1].string_ );
        F.toUpper();
        if (fun[F]++ != 0)
        {
          Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
            << "Duplicate function definition detected: " <<  F;
          result = false;
        }
      }
    }

    else if (ES1 == ".INCLUDE" || ES1 == ".INC" || ES1 == ".LIB")
    {
      std::string includeFile, libSelect_child;
      libSelect_child = libSelect;
      handleIncludeLine( line, ES1, includeFile, libSelect_child, libInside );
      if (includeFile != "")
      {
        result = parseIncludeFile(includeFile, libSelect_child, devMap, par, fun, modMap, sub);
      }
    }
    else if (ES1 == ".ENDL")
    {
      handleEndlLine ( line, libSelect, libInside );
    }
    else if (ES1 == ".INITCOND" )
    {
      handleInitCond( line );
    }
    else if (ES1 == ".MODEL")
    {
      ParameterBlock* modelPtr =
        new ParameterBlock (circuitBlockPtr_->netlistFileName, line);
      result = modelPtr->extractModelData( circuitBlockPtr_->metadata );

      ExtendedString M ( line[1].string_ );
      M.toUpper();
      std::map<std::string,ParameterBlock*,LessNoCase>::iterator mp = modMap.find(M);
      if (mp == modMap.end())
      {
        // Save for potential use later, like if there are data for other temperatures
        modMap[M] = modelPtr;
        // Add the model to the circuit context. Note, the circuit context
        // will handle deleting the model.
        circuitBlockPtr_->circuitContext.addModel(modelPtr);
      }
      else
      {
        // A duplicate model name has been detected.  Hopefully this is a data
        // point at another temperature, or other independent parameter.
        ParameterBlock *pb = modMap[M];
        if (pb->getType() == modelPtr->getType() && pb->getLevel() == modelPtr->getLevel())
        {
          std::vector<Device::Param> addMP;
          addMP.push_back(Device::Param("INDEPENDENT;PARAM","TNOM"));
          pb->addParameters(addMP);
          pb->addParameters(modelPtr->getParams());
          delete modelPtr;
        }
        else
        {
          Report::UserWarning().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
            << "Duplicate model definition detected: " << M;
          result = false;
        }
      }
    }
    else if (ES1 == ".PARAM" || ES1 == ".GLOBAL_PARAM")
    {
      if ( ES1 == ".GLOBAL_PARAM" && circuitBlockPtr_->parentCircuitPtr != NULL )
      {
        Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
          << "Attempt to assign global_param inside of subcircuit";
        result = false;
      }
      else {
        OptionBlock param( circuitBlockPtr_->netlistFileName, line,
                           circuitBlockPtr_->metadata );
        result = param.extractData();
        if (result) {
          int numberOfParams = param.getNumberOfParameters();
          for (int i = 0; i < numberOfParams; ++i)
          {
            ExtendedString P ( param.getParameter(i).tag() );
            P.toUpper();
            if (!P.possibleParam())
            {
              Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
                << "Illegal parameter name: " << P;
              result = false;
            }
            else if (par[P]++ != 0)
            {
              Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
                << "Duplicate parameter definition detected: " << P;
              result = false;
            }
          }
        }
        if (ES1 == ".PARAM")
          circuitBlockPtr_->circuitContext.addParams( param );
        else
          circuitBlockPtr_->circuitContext.addGlobalParams( param );
      }
    }
    else if (ES1 == ".GLOBAL")
    {
      if (circuitBlockPtr_->parentCircuitPtr != NULL )
      {
        Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
          << "Attempt to assign global node inside of subcircuit";
      }
      if (line.size() != 2)
      {
        Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
          << "Syntax error in .global, should be .global <node>";
      }
      ExtendedString ES2(line[1].string_);
      ES2.toUpper();
      circuitBlockPtr_->circuitContext.addGlobalNode ( ES2 );
    }
    else if (ES1 == ".SAVE")
    {
      Report::UserWarning0().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
        << ".SAVE line not currently handled, statement skipped";
    }
    else if (ES1 == ".SUBCKT")
    {
      // Create a new CircuitBlock to hold the subcircuit definition.
      // Set the parentCircuitPtr of the new CircuitBlock.
      CircuitBlock* subcircuitBlockPtr =
        new CircuitBlock(
          circuitBlockPtr_->netlistFileName,
          line, commandLine_,
          circuitBlockPtr_->metadata,
          modelNames_,
          ssfMap_,
          circuitBlockPtr_->circuitContext,
          circuitBlockPtr_->outputMgrPtr_,
          circuitBlockPtr_->useCount,
          globalParamsInserted_,
          mainCircuitPtr_,
          distToolPtr_,
          insertionToolPtr_,
          devIntPtr_,
          circuitBlockPtr_->deviceNames_,
          circuitBlockPtr_->nodeNames_,
          externalNetlistParams_,
          remove_redundant_C_,
          remove_redundant_D_,
          remove_redundant_I_,
          remove_redundant_L_,
          remove_redundant_M_,
          remove_redundant_Q_,
          remove_redundant_R_,
          remove_redundant_V_,
          replace_ground_
                         );

      // Start a new subcircuit context in the circuit context.
      result = circuitBlockPtr_->circuitContext.beginSubcircuitContext(circuitBlockPtr_->netlistFileName, line);
      subcircuitBlockPtr->parentCircuitPtr = circuitBlockPtr_;
      subcircuitBlockPtr->netlistFileName = circuitBlockPtr_->netlistFileName;

      // Subcircuits must use the ssfPtr_ that was set up by the
      // main circuit for parsing the input file.
      subcircuitBlockPtr->setSSFPtr(ssfPtr_);
      subcircuitBlockPtr->setStartPosition();

      // Extract the subcircuit data from line.
      if (result) {
        result = subcircuitBlockPtr->extractSubcircuitData(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
                 && result;

        ExtendedString S ( line[1].string_ );
        S.toUpper();
        if ( circuitBlockPtr_->circuitBlockTable_.find( S ) != circuitBlockPtr_->circuitBlockTable_.end() )
        {
          Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
            << "Duplicate subcircuit definition detected: " <<  S;
          result = false;
        }

        circuitBlockPtr_->circuitBlockTable_[S] = subcircuitBlockPtr;
      }

      result = subcircuitBlockPtr->parseNetlistFilePass1(libSelect, libInside)
               && result;
    }

    else if (ES1 == ".PREPROCESS")
    {
      result=true; //ignore this line; it's job is done in the preprocess
      //phase
    }
    else
    {
      // If we get here then we have an unrecognized "." line, flag it
      // with a warning, ignore it and continue.
      Report::UserWarning0().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
        <<  "Unrecognized dot line will be ignored";
    }
  }
  else if (lineType == '*' || lineType == ' ' || lineType == '\t')
  {
    result = true;
  }
  else
  {
    Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
      << "Unrecognized line";
    result = false;
  }
  return !eof;
}


//----------------------------------------------------------------------------
// Function       : CircuitBlockData::getLinePass2
// Purpose        :
// Special Notes  :
// Scope          : private
// Creator        : Lon Waters
// Creation Date  : 08/01/2003
//----------------------------------------------------------------------------
bool CircuitBlockData::getLinePass2(
    std::vector<SpiceSeparatedFieldTool::StringToken> & line,
    const std::string &libSelect, std::string &libInside)
{
  int eof = 0;

  while (!eof)
  {
    eof = !ssfPtr_->getLine(line,replace_ground_); // Breaks the line into fields.

    if (DEBUG_IO) {
      Xyce::dout() << "pass 2 read netlist line:  "<< std::endl;
      for (unsigned int i = 0; i < line.size(); ++i)
      {
        Xyce::dout() << line[i].string_ << " ";
      }
      Xyce::dout() << std::endl;

      if (eof)
      {
        Xyce::dout() << "We are at the end of file, size of line is " << line.size()
                     <<std::endl;
      }
    }

    // better not try to do anything if getLine returned an empty line!
    if ( line.empty() )
      continue;

    // Determine what to do with the parsed line.
    ExtendedString ES1(line[0].string_);
    ES1.toUpper();
    if ( !(libSelect != libInside && libInside != "" && ES1 != ".ENDL"))
    {
      char lineType;
      lineType = ES1[0];
      bool removecomponent=false;

      if (lineType >= 'A' && lineType <= 'Z')
      {

        if (lineType == 'C' || lineType == 'D' || lineType == 'I' || lineType
          == 'L' || lineType == 'R' || lineType == 'V')
        {
          if (line.size() > 2) //make sure that there are two nodes to check!
          {
            ExtendedString node1 ( line[1].string_ );
            ExtendedString node2 ( line[2].string_ );
            node1.toUpper();
            node2.toUpper();

            removecomponent = removeTwoTerminalDevice(lineType,node1, node2);
          }
        }
        else if (lineType == 'M' || lineType == 'Q')
        {
          if (line.size() > 3) //make sure that there are three nodes to check!
          {
            ExtendedString node1 ( line[1].string_ );
            ExtendedString node2 ( line[2].string_ );
            ExtendedString node3 ( line[3].string_ );
            node1.toUpper();
            node2.toUpper();
            node3.toUpper();

            removecomponent = removeThreeTerminalDevice(lineType,node1, node2,
                    node3);
          }
        }

        if (!removecomponent)
        {
          if (lineType != 'X' && lineType != 'K' && lineType != 'L')
          {
            // check for .INITCOND values
            std::string tmpName(circuitBlockPtr_->circuitContext.getPrefix());
            if (tmpName == "")
              tmpName = ES1;
            else
              tmpName = tmpName + ":" + ES1;

            if( mainCircuitPtr_->initCondIndex.find( tmpName ) !=
              mainCircuitPtr_->initCondIndex.end() )
            {
              // append the IC=val1...valN list to device line
              line.insert( line.end(),
                ( mainCircuitPtr_->initCondIndex[tmpName] ).begin(),
                ( mainCircuitPtr_->initCondIndex[tmpName] ).end() );
            }

            return true;
          }
          else if (lineType == 'L')
          {
            std::set< std::string > & cTable =
              circuitBlockPtr_->circuitContext.getAllCoupledInductors();

            if( cTable.find( ES1 ) == cTable.end() )
            {
              // treat as a normal device line if inductor is not coupled
              return true;
            }
            continue;
          }
          else if (lineType == 'K')
          {
              // DeviceBlock device(circuitBlockPtr_->netlistFileName, line);
              // circuitBlockPtr_->mutualInductors_.push_back(device);
          }
          else if (lineType == 'X')
          {
            circuitBlockPtr_->handleDeviceLine(line, libSelect, libInside);
          }
        }
#ifdef Xyce_DEBUG_IO
        else
        {
          Xyce::dout() << "Netlist Parse 2:  ";
          Xyce::dout() << "removing component " << ES1 << ".  All nodes on the device";
          Xyce::dout() << " are the same."  << std::endl;
        }
#endif
      }
      else if (lineType == '.')
      {
        if (ES1 == ".SUBCKT")
        {
          // Jump to the end of the subcircuit.
          // Find the subcircuit corresponding to this instance.
          CircuitBlock* subcircuitPtr = circuitBlockPtr_->
            findSubcircuit(ExtendedString(line[1].string_).toUpper());

          // Set the end location of the subcircuit in its associated file.
          subcircuitPtr->setFilePosition(subcircuitPtr->getEndPosition());
          subcircuitPtr->setLinePosition( subcircuitPtr->getLineEndPosition() );

        }
        else if (ES1 == ".INCLUDE" || ES1 == ".INC" || ES1 == ".LIB")
        {
          std::string includeFile, libSelect_child;
          libSelect_child = libSelect;
          handleIncludeLine( line, ES1, includeFile, libSelect_child, libInside );
          if (includeFile != "")
          {
            int result = parseIncludeFile2(includeFile, libSelect_child);
            if (!result)
return result;
          }
        }
        else if (ES1 == ".ENDL")
        {
          handleEndlLine ( line, libSelect, libInside );
        }
        else if (ES1 == ".ENDS")
        {
          return false;
        }
        if (ES1 == ".END")
        {
          return false;
        }
      }
      else
      {
        // Do nothing for other types of lines
      }
    }
  }
  // Only get here if end of file.

  return false;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlockData::getLinePassMI
// Purpose        :
// Special Notes  :
// Scope          : private
// Creator        : Lon Waters
// Creation Date  : 08/01/2003
//----------------------------------------------------------------------------
bool CircuitBlockData::getLinePassMI()
{
  int eof = 0;

  std::vector<SpiceSeparatedFieldTool::StringToken> line;

  while (!eof)
  {

    eof = !ssfPtr_->getLine(line,replace_ground_);
  // Breaks the line into fields.
    //FLAG  Should I be using the optional boolean argument here?  I think so.
    if (DEBUG_IO) {
      Xyce::dout() << "pass MI read netlist line: ";
      for (unsigned int i = 0; i < line.size(); ++i)
      {
        Xyce::dout() << line[i].string_ << " ";
      }
      Xyce::dout() << std::endl;
    }

    // better not try to do anything if getLine returned an empty line!
    if ( !(line.empty()) )
    {
      // Determine what to do with the parsed line.
      char lineType;
      ExtendedString ES1 ( line[0].string_ );
      ES1.toUpper();
      lineType = ES1[0];

      if (lineType == 'L')
      {
        // This is an inductor, check for assoc. MI
        std::vector<CircuitContext::MutualInductance> & MIs =
          circuitBlockPtr_->circuitContext.getMutualInductances();
        int numMIs = MIs.size();
        for( int i = 0; i < numMIs; ++i )
        {
          if( MIs[i].inductors.count( ES1 ) )
          {
            device_.clear();
            device_.setParsedLine(line);
            device_.setFileName(circuitBlockPtr_->netlistFileName);
            device_.extractData( );

            int numParams = device_.getNumberOfInstanceParameters();
            for( int j = 0; j < numParams; ++j )
            {
              Device::Param param = device_.getInstanceParameter(j);
              if( param.uTag() == "L" )
              {
                MIs[i].inductors[ES1] = param.getImmutableValue<double>();

                // store terminal names associated with this inductor
                ( MIs[i].terminals[ES1] ).push_back( device_.getNodeValue( 0 ) );
                ( MIs[i].terminals[ES1] ).push_back( device_.getNodeValue( 1 ) );
              }
            }
          }
        }
        return true;
      }
      else if (lineType == '.')
      {
        // Jump to the end of the subcircuitj; parseMutualInductances will
        // handle subckt directly
        if (ES1 == ".SUBCKT")
        {
          // Find the subcircuit corresponding to this instance.
          CircuitBlock* subcircuitPtr = circuitBlockPtr_->
            findSubcircuit(ExtendedString(line[1].string_).toUpper());

          // Set the end location of the subcircuit in its associated file.
          subcircuitPtr->setFilePosition(subcircuitPtr->getEndPosition());
          subcircuitPtr->setLinePosition( subcircuitPtr->getLineEndPosition() );

        }
        else if (ES1 == ".ENDS" || ES1 == ".END")
          return false;
      }
      else
      {
        return true;
      }
    }
  }

  // Only get here if end of file.
  return false;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::setFileName
// Purpose       : Change netlist file name on slave processors
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 03/10/2006
//--------------------------------------------------------------------------
void CircuitBlock::setFileName ( std::string & fileNameIn )
{
  netlistFileName = fileNameIn;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlock::handleDeviceLine
// Purpose        : Processor zero / serial processing of devices
// Special Notes  :
// Scope          : public
// Creator        : Lon Waters
// Creation Date  : 02/21/2003
//----------------------------------------------------------------------------
bool CircuitBlock::handleDeviceLine(
    std::vector<SpiceSeparatedFieldTool::StringToken> const& deviceLine,
    const std::string &libSelect, const std::string &libInside)
{
  bool result;

  data_->device_.clear();
  data_->device_.setParsedLine(deviceLine);
  data_->device_.setFileName(netlistFileName);
  std::string prefix(circuitContext.getPrefix());
  std::map<std::string, std::string> * nodeMapPtr = circuitContext.getNodeMapPtr();
  result = data_->instantiateDevice(data_->device_, prefix, *nodeMapPtr, libSelect, libInside);

  return result;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::expandSubcircuitInstance
// Purpose       : Expand a subcircuit instance by adding the devices and
//                 device models that compose the subcircuit to the main
//                 (top level) circuit.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 12/28/2001
//--------------------------------------------------------------------------
bool CircuitBlockData::expandSubcircuitInstance(
  IO::DeviceBlock &     subcircuitInstance,
  const std::string &   libSelect,
  const std::string &   libInside)
{
  // Set subcircuitPrefix.
  std::string subcircuitPrefix("");
  if ( circuitBlockPtr_->circuitContext.getPrefix() != "" )
  {
    subcircuitPrefix = circuitBlockPtr_->circuitContext.getPrefix() +
      ":" + subcircuitInstance.getName();
  }
  else
  {
    subcircuitPrefix = subcircuitInstance.getName();
  }

  // Find the subcircuit corresponding to this instance.
  CircuitBlock* subcircuitPtr = circuitBlockPtr_->findSubcircuit(subcircuitInstance.getModelName());
  if ( subcircuitPtr == NULL )
  {
    Report::UserError0().at(circuitBlockPtr_->netlistFileName, subcircuitInstance.getParsedLine()[0].lineNumber_)
      << "Subcircuit " << subcircuitInstance.getModelName()
      << " has not been defined for instance " << subcircuitInstance.getName();
    distToolPtr_->endDeviceLines();
    return false;
  }

  // get the list of X nodes
  std::list< std::string > subcircuitInstanceNodes;
  subcircuitInstance.getAllNodeNames( subcircuitInstanceNodes );

  // Set the context for this subcircuit instance.
  std::string subcircuitName(subcircuitInstance.getModelName());
  bool cresult = circuitBlockPtr_->circuitContext.setContext(subcircuitName, subcircuitPrefix, subcircuitInstanceNodes);
  if (!cresult)
  {
    distToolPtr_->endDeviceLines();
    Report::UserError0().at(circuitBlockPtr_->netlistFileName, subcircuitInstance.getParsedLine()[0].lineNumber_)
      << "Error invoking subcircuit " << subcircuitInstance.getModelName()
      << " instance " << subcircuitInstance.getName();
    return false;
  }

  // get the list of .SUBCKT nodes
  std::vector<std::string> subcircuitNodes =
    circuitBlockPtr_->circuitContext.getNodeList();

  // Make sure the subcircuit instance and subcircuit definition agree on the
  // number of nodes.
  if ( subcircuitInstanceNodes.size() != subcircuitNodes.size() )
  {
    distToolPtr_->endDeviceLines();
    Report::UserError0() << "Number of nodes for subcircuit instance " << subcircuitInstance.getName()
                         << " does not agree with number of nodes in subcircuit "
                         << circuitBlockPtr_->circuitContext.getCurrentContextName();
    return false;
  }

  // these iterators loop over the nodes in this context (i.e. the real node
  // connected to the interface nodes)
  std::list<std::string>::iterator circuitInstanceNodeIt_ = subcircuitInstanceNodes.begin();
  std::list<std::string>::iterator endCircuitInstanceNodeIt_ = subcircuitInstanceNodes.end();

  // these iterators loop over the subcircuits interface nodes
  std::vector<std::string>::iterator subcircuitNodeInterfaceIt_ = subcircuitNodes.begin();
  std::vector<std::string>::iterator endSubcircuitNodeInterfaceIt_ = subcircuitNodes.end();

  // add the interface nodes of this subcirciut to the aliasNodeMap so that
  // we can look up the aliases later if needed
  while( ( circuitInstanceNodeIt_ != endCircuitInstanceNodeIt_ ) &&
         ( subcircuitNodeInterfaceIt_ != endSubcircuitNodeInterfaceIt_ ) )
  {
    std::string key( subcircuitPrefix + ":" + *subcircuitNodeInterfaceIt_ );

    if ( ( mainCircuitPtr_->aliasNodeMapHelper_ ).find( key ) !=
         ( mainCircuitPtr_->aliasNodeMapHelper_ ).end() )
    {
      ( mainCircuitPtr_->aliasNodeMap_ )[key] = *circuitInstanceNodeIt_;

      if (DEBUG_IO)
        Xyce::dout() << "Found node alias:  " << key << " ==> " << *circuitInstanceNodeIt_ << std::endl;
    }

    circuitInstanceNodeIt_++;
    subcircuitNodeInterfaceIt_++;
  }


  // Locate the subcircuit in the netlist file. It can either be in
  // the file currently being read, or in a separate include file.
  SpiceSeparatedFieldTool * newssf = ssfPtr_;
  if (subcircuitPtr->netlistFileName != circuitBlockPtr_->netlistFileName)
  { // The subcircuit is in an include file.
    // Get SSF from Pass 1's ssf map
    if( ssfMap_.count( subcircuitPtr->netlistFileName ) )
      newssf = ssfMap_[subcircuitPtr->netlistFileName].second;
    else
    {
      distToolPtr_->endDeviceLines();
      Report::UserError().at(circuitBlockPtr_->netlistFileName, subcircuitInstance.getParsedLine()[0].lineNumber_)
        << "Can't find include file " << subcircuitPtr->netlistFileName;
    }
  }

  subcircuitPtr->setSSFPtr( newssf );

  // Set the position of the subcircuit in its file.
  std::streampos oldLoc = newssf->getFilePosition();
  int oldLine = subcircuitInstance.getParsedLine()[0].lineNumber_;
  subcircuitPtr->setFilePosition(subcircuitPtr->getStartPosition());
  subcircuitPtr->setLinePosition( subcircuitPtr->getLineStartPosition() );

  // Resolve parameters and functions in the current context.
  bool result;
  std::vector<Device::Param> subcircuitInstanceParams;
  subcircuitInstance.getInstanceParameters(subcircuitInstanceParams);

  result = circuitBlockPtr_->circuitContext.resolve(subcircuitInstanceParams);
  if (!result)
return result;

  // Tell the distribution tool that the context has changed.
  distToolPtr_->circuitStart( subcircuitName,
                              subcircuitInstanceNodes,
                              subcircuitPrefix,
                              subcircuitInstanceParams );

  // Instantiate the devices in this subcircuit instance.
  subcircuitPtr->instantiateDevices(libSelect, libInside);
  // send MIs if present
  if( circuitBlockPtr_->circuitContext.haveMutualInductances() )
  {
    // Distribute each Y%MI?%name device line, end of circuit
    int n = circuitBlockPtr_->circuitContext.getNumMILines();
    for( int i = 0; i < n; ++i )
    {
      // parse locally if distool does not distribute; normally the
      // CircuitBlock::instantiateDevices() performs this step
      if( !distToolPtr_->circuitDeviceLine(
       circuitBlockPtr_->circuitContext.getMILine( i ) ) )
      {
        circuitBlockPtr_->handleDeviceLine(
         circuitBlockPtr_->circuitContext.getMILine( i ), libSelect, libInside );
      }
    }
  }

  // Return to previous context.
  circuitBlockPtr_->circuitContext.restorePreviousContext();
  subcircuitPtr->setFilePosition(oldLoc);
  subcircuitPtr->setLinePosition( oldLine );

  // Tell the distribution tool that the current context has ended and
  // to switch back to the previous context.
  distToolPtr_->circuitEnd();

  // Only get here on success.
  return true;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlock::findSubcircuit
// Purpose       : Search the circuitBlockTable_ of the current circuit block
//                 for the subcircuit of the given name. If it is not found,
//                 recursively search each parent subcircuit. Return a
//                 pointer to the circuit block if it is found, otherwise
//                 return NULL.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 12/28/2001
//--------------------------------------------------------------------------
CircuitBlock* CircuitBlock::findSubcircuit( std::string const& subcircuitName)
{
  // Search this circuit blocks subcircuit list.
  if ( circuitBlockTable_.find( subcircuitName ) != circuitBlockTable_.end() )
  {
      return circuitBlockTable_.find( subcircuitName )->second;
  }
  else
  {
    // The subcircuit was not found in the current circuit's subcircuit list,
    // recursively search the parent circuit's subcircuit list.
    CircuitBlock* circuitBlockPtr = NULL;
    if ( parentCircuitPtr != NULL )
    {
      return circuitBlockPtr = parentCircuitPtr->findSubcircuit( subcircuitName );
    }
    else
    {
      return NULL;
    }
  }

}


//--------------------------------------------------------------------------
// Function      : CircuitBlock::substituteNodeAliases
// Purpose       : Scans the OptionBlock data held in the class var
//                 std::list<Util::OptionBlock> optionsTable for aliased node
//                 names and then replaces those with real names
// Special Notes :
// Creator       : Rich Schiek, Electrical and Microsystems Modeling
// Creation Date : 04/28/2009
//--------------------------------------------------------------------------
bool CircuitBlock::substituteNodeAliases()
{
  bool didASubstitution = false;
  // scan through the list of options block
  std::list<Util::OptionBlock>::iterator currentOptionBlockIter = optionsTable.begin();
  std::list<Util::OptionBlock>::iterator endOptionBlockIter = optionsTable.end();

  while( currentOptionBlockIter != endOptionBlockIter )
  {
    // for now we'll just substitute aliases that are on the print line
    // we could expand this to do all of the option blocks if needed
    if (currentOptionBlockIter->getName() == "PRINT")
    {
      std::list< Util::Param >::iterator currentParamIt =  currentOptionBlockIter->getParams().begin();
      std::list< Util::Param >::iterator endParamIt =  currentOptionBlockIter->getParams().end();
      std::map<std::string,std::string>::iterator endOfAliasMapIt =  aliasNodeMap_.end();

      while( currentParamIt != endParamIt )
      {
        std::string key(currentParamIt->uTag());
        if( aliasNodeMap_.find( key ) != endOfAliasMapIt )
        {
          // replace this tag with the real node name from the alias map
          currentParamIt->setTag( aliasNodeMap_[ key ] );
          didASubstitution=true;
        }
        currentParamIt++;
      }

      break;
    }
    currentOptionBlockIter++;
  }

  return didASubstitution;
}


//--------------------------------------------------------------------------
// Function      : CircuitBlock::resolveExpressionsInOptionBlocks
// Purpose       : Scans the OptionBlock data held in the class var
//                 std::list<Util::OptionBlock> optionsTable for expressions
//                 and tries to resolve them so that a call on an
//                 expressions eval() method will return the right result.
// Special Notes :
// Creator       : Rich Schiek, Electrical Systems Modeling
// Creation Date : 02/02/2012
//--------------------------------------------------------------------------
bool CircuitBlock::resolveExpressionsInOptionBlocks()
{
  bool resolutionSuccesfull = false;
  // scan through the list of options block
  std::list<Util::OptionBlock>::iterator currentOptionBlockIter = optionsTable.begin();
  std::list<Util::OptionBlock>::iterator endOptionBlockIter = optionsTable.end();

  while( currentOptionBlockIter != endOptionBlockIter )
  {
    // for now we'll just substitute aliases that are on the print line
    // we could expand this to do all of the option blocks if needed
    if (currentOptionBlockIter->getName() == "PRINT")
    {
      std::list< Util::Param >::iterator currentParamIt =  currentOptionBlockIter->getParams().begin();
      std::list< Util::Param >::iterator endParamIt =  currentOptionBlockIter->getParams().end();
      while( currentParamIt != endParamIt )
      {
        if( Util::hasExpressionTag(*currentParamIt) )
        {
          // ok here's the odd part.  the tag in this case is the expression "{sin(v(a))}"
          // and the value has been set to zero earlier.  The code that processes Util::Param
          // objects into full expressions assumes the tag is a descriptor like "IMPORTANT_PARAM"
          // and that the value is the expression as in "{sin(v(a))}" here.  So we need to
          // make an adjustment here by copying the "tag" to the "value".
          currentParamIt->setVal( currentParamIt->tag() );
          std::vector<std::string> exceptionStrings;
          circuitContext.resolveParameter(*currentParamIt,exceptionStrings);
          // To do:  need to check if exceptionStrings comes back as non zero if it
          // couldn't fully resolve the expression.  Output a better error if that's
          // the case.
          if( exceptionStrings.size() > 0)
          {
            Xyce::dout() << "Exception strings found in processing expression in CircuitBlock::resolveExpressionsInOptionBlocks() = ";
            for( unsigned int i=0; i<exceptionStrings.size(); i++ )
              Xyce::dout() << "\"" << exceptionStrings[i] << "\" ";
            Xyce::dout() << " Ignoring for now. Please report this error." << std::endl;
          }
          resolutionSuccesfull=true;
        }
        currentParamIt++;
      }
      // can't break here as we need to set the status of the other option blocks
      // break;
    }

    currentOptionBlockIter++;
  }

  return resolutionSuccesfull;
}


//--------------------------------------------------------------------------
// Function      : CircuitBlockData::handleIncludeLine
// Purpose       : Handle a netlist .include or .lib line, add the include file
//                 to includeFiles_ for later processing.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 01/10/2001
//--------------------------------------------------------------------------
void CircuitBlockData::handleIncludeLine(
      std::vector<SpiceSeparatedFieldTool::StringToken> const& parsedLine,
      const ExtendedString & ES1, std::string& includeFile, std::string& libSelect, std::string& libInside)
{
  bool lib = ES1.substr(0,4) != ".INC";

  std::string includeFileTmp;
  
  // Get the file name indicated by the .include line.
  if ( parsedLine.size() < 2 )
  {
    Report::UserWarning0().at(circuitBlockPtr_->netlistFileName, parsedLine[0].lineNumber_)
      << ES1 << " is missing argument(s), ignoring";
    return;
  }
  else
  {
    includeFileTmp = parsedLine[1].string_;
  }

  if (!lib || (lib && (parsedLine.size() >= 3 || (ssfMap_.size() == 1 && parsedLine.size() >= 2))))
  {
    // Strip off the enclosing double quotes if they are present.
    if ( (includeFileTmp[0] == '"') &&
         (includeFileTmp[includeFileTmp.length()-1] == '"') )
    {
      includeFile = includeFileTmp.substr( 1, includeFileTmp.length()-2 );
    }
    else
    {
      includeFile = includeFileTmp;
    }
  }
  else
  {
    includeFile = "";
  }

  if ( lib )
  {
    if ( parsedLine.size() > 3)
    {
      Report::UserWarning0().at(circuitBlockPtr_->netlistFileName, parsedLine[0].lineNumber_)
        << "Extraneous data on .LIB ignored";
    }
    if ( parsedLine.size() >= 3 )
    {
      libSelect = ExtendedString(parsedLine[2].string_).toUpper();
    }
    else
    {
      libInside = ExtendedString(parsedLine[1].string_).toUpper();
    }
  }
  else
  {
    if ( parsedLine.size() >= 3)
    {
      Report::UserWarning0().at(circuitBlockPtr_->netlistFileName, parsedLine[0].lineNumber_)
        << "Extraneous data on .INCLUDE ignored";
    }
  }
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::handleStd::EndlLine
// Purpose       : Handle a netlist .endl line
// Special Notes :
// Creator       : Dave Shirley, PSSI
// Creation Date : 10/09/2009
//--------------------------------------------------------------------------
void CircuitBlockData::handleEndlLine(
      std::vector<SpiceSeparatedFieldTool::StringToken> const& parsedLine,
      const std::string& libSelect, std::string& libInside)
{
  if (libInside == "")
  {
    Report::UserError().at(circuitBlockPtr_->netlistFileName, parsedLine[0].lineNumber_)
      << ".ENDL encountered without .LIB ";
    return;
  }

  if ( parsedLine.size() >= 2 )
  {
    ExtendedString libName ( parsedLine[1].string_ );
    libName.toUpper();
    if (libName != libInside)
    {
      Report::UserError().at(circuitBlockPtr_->netlistFileName, parsedLine[0].lineNumber_)
        << ".ENDL encountered with name " << libName << ", which does not match .LIB name " << libInside;
    }
    if ( parsedLine.size() >2 )
    {
      Report::UserWarning().at(circuitBlockPtr_->netlistFileName, parsedLine[0].lineNumber_)
        << "Extraneous field(s) following library name in .ENDL";
    }
  }
  else
  {
    Report::UserError().at(circuitBlockPtr_->netlistFileName, parsedLine[0].lineNumber_)
      << ".ENDL encountered without library name, currently inside .LIB " << libInside;
  }

  libInside = "";
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::parseIncludeFile
// Purpose       : Parse each include file in includeFiles_ adding the
//                 contents to the current CircuitBlock.
// Special Notes :
// Creator       : Lon Waters
// Creation Date : 01/10/2001
//--------------------------------------------------------------------------
bool CircuitBlockData::parseIncludeFile(
  std::string const& includeFile, std::string const& libSelect,
  std::map<std::string,int> & devMap, std::map<std::string,int> & par, std::map<std::string,int> & fun,
  CircuitContext::ModelMap & modMap, std::map<std::string,int> & sub)
{
  // Save current ssfPtr_ and netlistFileName.
  SpiceSeparatedFieldTool* oldssfPtr = ssfPtr_;

  // save the old file name (parent file)
  std::string old_netlistFileName(circuitBlockPtr_->netlistFileName);


  // set the current netlist file to the name of this include file.
  circuitBlockPtr_->netlistFileName = includeFile;

  if (DEBUG_IO)
    Xyce::dout() << "CircuitBlockData::parseIncludeFile: Parsing include file: " << includeFile << std::endl;

  if( !ssfMap_.count(includeFile) )
  {
    // Create a new SpiceSeparatedFieldTool for this include file.
    std::ifstream * includeIn = new std::ifstream;
    // Using binary to avoid issues with compiler/plat
    // *fstream differences in implementation
    includeIn->open( includeFile.c_str(), std::ios::in | std::ios::binary );
    if ( !includeIn->is_open() )
    {
      Report::UserError0() << "Could not find include file " << includeFile;
      return false;
    }

    ssfPtr_ = new SpiceSeparatedFieldTool(*includeIn, includeFile, externalNetlistParams_);


    ssfMap_[includeFile] = FileSSFPair( includeIn, ssfPtr_ );
  }
  else
  {
    // we already have an ssF for this file in the map, just pull it out,
    // rewind, and proceed.
    ssfPtr_ = ssfMap_[includeFile].second;

    if (DEBUG_IO) {
      Xyce::dout() << "  CircuitBlockData::parseIncludeFile: found existing ssFT " << std::endl
                   << " \t its file name is " << ssfPtr_->getFileName() << std::endl
                   << "\t its current location is " << ssfPtr_->getFilePosition() << std::endl
                   << "\t its current line number is " << ssfPtr_->getLineNumber() << std::endl
                   << "\t Rewinding to location 0 and line 1" << std::endl;
    }

    ssfPtr_->setLocation(0);
    ssfPtr_->setLineNumber(1);
  }

  // Handle the include file lines.
  bool result = false;
  std::string libInside;
  for (;;) {
    bool line_parsed = true;

    if (handleLinePass1(line_parsed, devMap, par, fun, modMap, sub, libSelect, libInside) )
      result = result && line_parsed;
    else
      break;
  }

  // Restore old ssfPtr_ and netlistFileName.
  ssfPtr_ = oldssfPtr;

  // restore the parent file name
  circuitBlockPtr_->netlistFileName = old_netlistFileName;


  if (DEBUG_IO)
    Xyce::dout() << "CircuitBlockData::parseIncludeFile: finished with include file: " << includeFile << std::endl;

  return true;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::parseIncludeFile2
// Purpose       : Jump to include file for 2nd pass
// Special Notes :
// Creator       : Rob Hoekstra
// Creation Date : 08/02/2004
//--------------------------------------------------------------------------
bool CircuitBlockData::parseIncludeFile2(std::string const& includeFile,
             const std::string &libSelect)
{
  // Save current ssfPtr_ and netlistFileName.
  SpiceSeparatedFieldTool* oldssfPtr = ssfPtr_;
  std::string old_netlistFileName(circuitBlockPtr_->netlistFileName);
  circuitBlockPtr_->netlistFileName = includeFile;
  distToolPtr_->setFileName(includeFile);

  if (DEBUG_IO)
    Xyce::dout() << "Parsing include file Pass 2: " << includeFile << std::endl;

  // Find SSF for file
  if( !ssfMap_.count( includeFile ) )
  {
    distToolPtr_->endDeviceLines();
    Report::UserError() << "Could not find include file SSF " << includeFile;
    return false;
  }
  ssfPtr_ = ssfMap_[includeFile].second;

  //Save current position and move to beginning
  ssfPtr_->setLocation(0);
  ssfPtr_->setLineNumber( 1 );

  int result = circuitBlockPtr_->instantiateDevices(libSelect, std::string(""));
  if (!result)
    return result;

  // Restore old ssfPtr_ and netlistFileName.
  ssfPtr_ = oldssfPtr;
  circuitBlockPtr_->netlistFileName = old_netlistFileName;
  distToolPtr_->setFileName(old_netlistFileName);

  if (DEBUG_IO)
    Xyce::dout() << "Done with include file Pass 2: " << includeFile << std::endl;

  return true; // Only get here on success.
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::handleInitCond
// Purpose       : Retrieve separate IC= data from line or external file and
//               : temporarily store in CircuitBlock
// Special Notes : Validation of .initcond lines is not done here.  Semantic
//               : errors are handled during device instantiation.
// Creator       :
// Creation Date :
//--------------------------------------------------------------------------
void CircuitBlockData::handleInitCond(
 std::vector< SpiceSeparatedFieldTool::StringToken > const& parsedLine )
{
  // check for multiple .initcond lines
  if( !(mainCircuitPtr_->initCondIndex.empty()) )
  {
    Report::UserError0() << ".INITCOND line may appear only once.";
  }

  // check minimum line length
  if( parsedLine.size() < 3 )
  {
    Report::UserError0().at(circuitBlockPtr_->netlistFileName, parsedLine[0].lineNumber_ )
      << ".INITCOND line is missing information";
  }

  ExtendedString tmpType ( parsedLine[1].string_ );
  tmpType.toUpper();

  // read IC values from separate file:  .INITCOND FILE fileName|"filename"
  if( tmpType == "FILE" )
  {
    // Strip off the enclosing double quotes if they are present.
    std::string initCondFile(parsedLine[2].string_);
    if ( ( initCondFile[0] == '"' ) &&
     ( initCondFile[initCondFile.length() - 1] == '"' ) )
    {
      initCondFile = initCondFile.substr( 1, initCondFile.length() - 2 );
    }

    // open the file for reading
    std::ifstream initCondIn;
    initCondIn.open( initCondFile.c_str(), std::ios::in | std::ios::binary );
    if( !initCondIn.is_open() )
    {
      Report::UserError0() << "Could not open the .INITCOND file " << initCondFile;
      return;
    }

    // use parser to extract data from the file
    SpiceSeparatedFieldTool ssfICPtr( initCondIn, initCondFile, externalNetlistParams_ );
    std::vector< SpiceSeparatedFieldTool::StringToken > line;

    while( !initCondIn.eof() )
    {
      // tokenize lines; file ptr is advanced in getLine()
      ssfICPtr.getLine( line, replace_ground_ );

      // check for enough data on line
      if( line.size() < 4 )
      {
        Report::UserError0() << ".INITCOND file '" << initCondFile << "' is not formatted properly.";
      }
      else {

        // store tokenized line in the index
        tmpType = line[0].string_;
        tmpType.toUpper();
        mainCircuitPtr_->initCondIndex[tmpType] =
          std::vector< SpiceSeparatedFieldTool::StringToken >(
            line.begin() + 1, line.end() );
      }
    }
  }

  // read IC values from line:  .INITCOND ( fqDevName IC = val (, val)* )+
  else
  {
    std::vector< SpiceSeparatedFieldTool::StringToken >::const_iterator
      pIter, pNextIter, pEndIter;

    pIter = parsedLine.begin() + 1;
    pEndIter = parsedLine.end();

    // check for enough data on line
    if( distance( pIter, pEndIter ) < 2 )
    {
      Report::UserError0() << ".INITCOND line is not formatted properly.";
    }
    else {

      // find the next dev name
      while( pIter != pEndIter )
      {
        // point to beginning of IC=val1...valN list
        pNextIter = pIter + 3;

        // build var1..varN list
        while( pNextIter != pEndIter && (*pNextIter).string_ != "=")
        {
          ++pNextIter;
        }

        // keep checking if more data is on the line
        if( pNextIter != pEndIter )
        {
          // move back to end of list
          pNextIter -= 2;
        }

        // copy into map
        ExtendedString tmpType ( (*pIter).string_ );
        tmpType.toUpper();
        mainCircuitPtr_->initCondIndex[tmpType] =
          std::vector< SpiceSeparatedFieldTool::StringToken >(
            pIter + 1, pNextIter );

        // move to end of line
        pIter = pNextIter;
      }
    }
  }

  if (DEBUG_IO) {
    std::map< std::string, std::vector< SpiceSeparatedFieldTool::StringToken > >::const_iterator a, b;
    a=mainCircuitPtr_->initCondIndex.begin();
    b=mainCircuitPtr_->initCondIndex.end();
    Xyce::dout() << ".INITCOND line yields " << mainCircuitPtr_->initCondIndex.size() <<
      " parsed devices:  ";
    for(; a!=b;++a)
    {
      Xyce::dout() << (*a).first << " ";
      for(unsigned int i=0; i< mainCircuitPtr_->initCondIndex[(*a).first].size();++i)
        Xyce::dout() << ( (mainCircuitPtr_->initCondIndex[(*a).first])[i] ).string_;
      Xyce::dout() << " ";
    }
    Xyce::dout() << std::endl;
  }
}


//----------------------------------------------------------------------------
// Function       : CircuitBlockData::handleDCSweep
// Purpose        : Post process a DC sweep if one was specified in the
//                  netlist.
//
// Special Notes  : This function:
//
//                     (1) adds VSTART,VSTOP, and VSTEP parameters to the
//                     .PRINT options that are passed to the output mgr.
//
//                     (2) Calculates VINTERVAL for each source,  and adds
//                     this information to the device options for each
//                     source, and also to the DC options for each source.
//
//                     (3) Calculates the total number of steps (the
//                     STEPS parameter) and adds it to the .DC options
//                     block.  Note that STEPS is the final vinterval.
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/06/2002
//----------------------------------------------------------------------------

namespace {

struct GetNameEqual
{
  bool operator()(const Util::OptionBlock &op, const char *name)
  {
    return op.getName() == name;
  }
};

} // namespace <unnamed>

bool CircuitBlockData::handleDCSweep()
{
  static const char *analysisOptions_[] = {"DC", "TRAN", "TR", "MPDE", "HB", "AC", "OP", "MOR"};
  static const char *printOptions_[] = {"PRINT"};

  std::list<Util::OptionBlock>::const_iterator op_analysis_it
    = std::find_first_of(circuitBlockPtr_->optionsTable.begin(), circuitBlockPtr_->optionsTable.end(),
                         &analysisOptions_[0], &analysisOptions_[sizeof(analysisOptions_)/sizeof(analysisOptions_[0])], GetNameEqual());

  if (op_analysis_it == circuitBlockPtr_->optionsTable.end())
  {
    distToolPtr_->endDeviceLines();

    // Problem, no analysis specified.
    Report::UserError() << "No analysis specified.";
    return false;
  }

  // Add the DC sweep parameters to the PRINT options. First, find the
  // DC PRINT options.
  std::list<Util::OptionBlock>::const_iterator op_param_it
    = std::find_first_of(circuitBlockPtr_->optionsTable.begin(), circuitBlockPtr_->optionsTable.end(),
                         &printOptions_[0], &printOptions_[sizeof(printOptions_)/sizeof(printOptions_[0])], GetNameEqual());
  if (op_param_it == circuitBlockPtr_->optionsTable.end())
  {
    Report::UserWarning0() << "No print specified";
    return true;
  }

  std::list<Util::Param>::const_iterator paramIter;
  paramIter = std::find(op_param_it->getParams().begin(), op_param_it->getParams().end(), Util::Param("TYPE", "") );

  // Check for consistency between analysis type and print type.
  std::string analysisName = op_analysis_it->getName();
  std::string usVal = paramIter->usVal();
  if (analysisName == "TR")
  {
    analysisName = "TRAN"; // TR is a synonym for TRAN
  }
  if (usVal == "TR") usVal = "TRAN";

  if ( (analysisName == "TRAN" && usVal != "TRAN") ||
       (analysisName == "DC" && usVal != "DC") ||
       (analysisName == "MPDE" && usVal != "TRAN") ||
       (analysisName == "HB" && usVal != "HB") ||
       (analysisName == "AC" && usVal != "AC") ||
       (analysisName == "MOR" && usVal != "MOR") )
  {
    distToolPtr_->endDeviceLines();

    // Problem, inconsistent analysis type and print type.
    Report::UserError0() << "Analysis type " << op_analysis_it->getName() << " and print type " << paramIter->usVal() << " are inconsistent.";
    return false;
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlockData::handleSTEPSweep
// Purpose        : Post process a STEP sweep if one was specified in the
//                  netlist.
//
// Special Notes  : This function literally does  nothing.  It should
//                  probably be removed.  It was setup originally to
//                  add the "stepenable" flag to the .PRINT options, but
//                  this isn't needed anymore.
//
// Scope          :
// Creator        : Eric R. Keiter, SNL
// Creation Date  : 10/30/2003
//----------------------------------------------------------------------------
bool CircuitBlockData::handleSTEPSweep()
{

  // Find the STEP analysis options.  Combine all instances of STEP into a
  // single set of STEP options.  Unfortunately, we have to support
  // multiple STEP statements.
  std::list<Util::OptionBlock>::iterator optionBlockIter;
  std::list<Util::OptionBlock>::iterator firstSTEPIter;
  std::list<Util::OptionBlock>::iterator first = circuitBlockPtr_->optionsTable.begin();
  std::list<Util::OptionBlock>::iterator last = circuitBlockPtr_->optionsTable.end();

  // First, find the STEP options, if they exist.  Nothing is actually
  // done to the STEP options in this function.  This is only a test to
  // see if STEP options exist in this netlist.  This means that it is only
  // neccessary to find one STEP statement, even if there are several.
  std::list<Util::OptionBlock>::iterator stepOBiter = first;
  while( stepOBiter != last )
  {
    if ( stepOBiter->getName() == "STEP" )
    {
        break;
    }
    ++stepOBiter;
  }

  // if didn't find STEP, skip the rest of this function.
  if (stepOBiter == last)
  {
    return true;
  }

  return true;
}

//----------------------------------------------------------------------------
// Function       : CircuitBlockData::handleMutualInductances
// Purpose        : Post-process the mutual inductors in the current circuit,
//                  each inductor must get the list of inductors it is coupled
//                  to, the coupling coeefficient, and model information.
// Special Notes  :
// Scope          :
// Creator        : Lon Waters
// Creation Date  : 08/08/2002
//----------------------------------------------------------------------------
bool CircuitBlockData::handleMutualInductances()
{
  if ( circuitBlockPtr_->mutualInductors_.empty() )
  {
    return true; // Return with success, no mutual inductors to handle.
  }

  std::string subcircuitPrefix(circuitBlockPtr_->circuitContext.getPrefix());

  std::vector<DeviceBlock>::iterator MI_Iter;
  std::vector<DeviceBlock>::iterator firstMI =
    circuitBlockPtr_->mutualInductors_.begin();
  std::vector<DeviceBlock>::iterator lastMI =
    circuitBlockPtr_->mutualInductors_.end();

  for ( MI_Iter = firstMI; MI_Iter != lastMI; ++MI_Iter )
  {
    DeviceBlock MIDev(*MI_Iter);

    if (!MIDev.isExtracted())
      MIDev.extractData();

    // Extract the inductor names from the mutual inductor instance and
    // get the coupling coefficient. Also get a handle to each inductors
    // instance block in the main circuit's device table.
    std::vector<std::string> inductorList;
    double coupling;
    std::vector<Device::InstanceBlock*> inductorIBs;
    int numParameters = MIDev.getNumberOfInstanceParameters();
    for ( int i = 0; i < numParameters; ++i )
    {
      const Device::Param &parameter = MIDev.getInstanceParameter(i);

      if ( parameter.tag() != "COUPLING" )
      {
        std::string inductorName = parameter.uTag();
        if ( subcircuitPrefix != "" )
        {
          inductorName = subcircuitPrefix + ":" + inductorName;
        }

        inductorList.push_back( inductorName );

        if (mainCircuitPtr_->deviceTable.find(inductorName) ==
              mainCircuitPtr_->deviceTable.end())
        {
          Report::UserError0() << "Could not find inductor " << parameter.uTag()
                               << " named in the mutual inductor " << MIDev.getName() << " in the circuit";
          return false;
        }

        inductorIBs.push_back(
            &mainCircuitPtr_->deviceTable[inductorName]->getDevBlock() );
      }
      else
      {
        coupling = parameter.getImmutableValue<double>();
      }
    }

    // Add the coupling info to the inductor instance block parameter data.
    size_t numCoupledInductors = inductorList.size();
    for ( size_t i = 0; i < numCoupledInductors; ++i)
    {
      if (i == 0)
      {
        Device::Param parameter;

        parameter.setTag( "FIRSTINDUCTOR" );
        parameter.setGiven(true);
        parameter.setVal( 1 );
        inductorIBs[i]->params.push_back(parameter);
      }

      if (MIDev.getModelName() != "")
      {
        Device::Param parameter;

        parameter.setTag( "NONLINEARCOUPLING" );
        parameter.setVal( 1 );
        inductorIBs[i]->params.push_back(parameter);
        inductorIBs[i]->modelFlag = true;

        std::string modelPrefix;
        ParameterBlock* MI_modelPtr;
        bool modelFound = circuitBlockPtr_->circuitContext.
          findModel( MIDev.getModelName(), MI_modelPtr, modelPrefix );

        if (modelFound)
        {
          // Check whether the current context is dependent on subcircuit
          // parameters (via a "params:" list on the .subckt line. This
          // will affect the model prefixing.
          if (circuitBlockPtr_->circuitContext.hasSubcircuitParams() &&
            MI_modelPtr->hasExpressionValuedParams())
          {
            // If there is a subcircuit parameter dependency, we will
            // assume the worst and add the subcircuit instance prefix
            // to the model prefix. This will insure that model parameters
            // that might depend on subcircuit parameters are assigned to
            // a unique model for each subcircuit instance.

            //TVR:  No, don't assume the worst, as doing so is inconsistent
            // with instantiation of the model below!  Only do this if
            // the model actually has expression valued params.
            modelPrefix = subcircuitPrefix + ":" + modelPrefix;
          }

          // Add the model prefix to the device's model name.
          inductorIBs[i]->setModelName(MIDev.getModelName());
          if (modelPrefix != "")
          {
            inductorIBs[i]->setModelName(modelPrefix + ":" + MIDev.getModelName());
          }
        }
        else
        {
          Report::UserError0() << "Unable to find mutual inductor model named " << MIDev.getModelName();
          return false;
        }

        // ERK: tmpModel is a persistant object in the CircuitBlockData class.
        // It is cheaper to use the "=" operator than the copy constructor here,
        // as the copy constructor will require allocations and deletions every
        // time this function is called, which is a lot.
        //ParameterBlock model(*MI_modelPtr);
        ParameterBlock tmpModel = (*MI_modelPtr);

        // Set the model parameter values.
        tmpModel.setParameterValues(&circuitBlockPtr_->circuitContext);

        // Add the model to the circuit.
        mainCircuitPtr_->addModel(tmpModel, modelPrefix);
      }

      Device::Param parameter;

      parameter.setTag( "COUPLING" );
      parameter.setVal( coupling );
      parameter.setGiven(true);
      inductorIBs[i]->params.push_back(parameter);

      // For each inductor, add the set of inductors to which it is coupled
      // the those inductors inductances to its parameter list.
      for ( size_t j = 0; j < numCoupledInductors; ++j )
      {
        if (j != i)
        {
          {
            Device::Param parameter;

            parameter.setTag( "COUPLEDINDUCTOR" );
            parameter.setVal( inductorList[j] );
            parameter.setGiven (true);
            inductorIBs[i]->params.push_back(parameter);
          }

          {
            Device::Param parameter;
            parameter.setTag("COUPLEDINDUCTANCE");

            std::vector<Device::Param>::iterator paramIter;
            std::vector<Device::Param>::iterator first = inductorIBs[j]->params.begin();
            std::vector<Device::Param>::iterator last = inductorIBs[j]->params.end();
            for (paramIter = first; paramIter != last; ++paramIter)
            {
              if (paramIter->uTag() == "L")
              {
                parameter.setVal(paramIter->getImmutableValue<double>());
                parameter.setGiven (true);
                inductorIBs[i]->params.push_back(parameter);
                break;
              }
            }
          }
        }
      }
    }
  }

  return true;
}


//----------------------------------------------------------------------------
// Function       : CircuitBlockData::handleMutualInductance
// Purpose        : Post-process the mutual inductors in the current circuit,
//                  each inductor must get the list of inductors it is coupled
//                  to, the coupling coeefficient, and model information.
// Special Notes  :
// Scope          :
// Creator        : Rob Hoekstra
// Creation Date  : 08/28/04
//----------------------------------------------------------------------------
bool CircuitBlockData::handleMutualInductance( DeviceBlock & device )
{
  std::string subcircuitPrefix = circuitBlockPtr_->circuitContext.getPrefix();

  //Check for mutual inductance assoc. with this inductor
  std::vector<CircuitContext::MutualInductance> & MIs = circuitBlockPtr_->circuitContext.getMutualInductances();

  std::string name = device.getName();
  std::string::size_type pos = name.find_last_of(":");
  if (pos != std::string::npos)
    name = name.substr( pos + 1, name.length() - (pos + 1));

  for( int i = 0; i < MIs.size(); ++i)
  {
    if( MIs[i].inductors.count( name ) )
    {
      // add mutual inductance info to this inductor
      if( MIs[i].inductors.begin()->first == name )
      {
        Device::Param param;

        param.setTag( "FIRSTINDUCTOR" );
        param.setGiven(true);
        param.setVal( 1 );
        device.addInstanceParameter( param );
      }

      if( MIs[i].model != "" )
      {
        std::string modelName(MIs[i].model);

        Device::Param param;

        param.setTag( "NONLINEARCOUPLING" );
        param.setVal( 1 );
        param.setGiven(true);
        device.addInstanceParameter( param );

        std::string modelPrefix;
        ParameterBlock* MI_modelPtr;
        bool modelFound = circuitBlockPtr_->circuitContext.
              findModel( modelName, MI_modelPtr, modelPrefix );

        if(modelFound)
        {
          // Check whether the current context is dependent on subcircuit
          // parameters (via a "params:" list on the .subckt line. This
          // will affect the model prefixing.
          if (circuitBlockPtr_->circuitContext.hasSubcircuitParams() &&
                MI_modelPtr->hasExpressionValuedParams())
          {
            // If there is a subcircuit parameter dependency, we will
            // assume the worst and add the subcircuit instance prefix
            // to the model prefix. This will insure that model parameters
            // that might depend on subcircuit parameters are assigned to
            // a unique model for each subcircuit instance.
            if( subcircuitPrefix != "" )
              modelPrefix = subcircuitPrefix + ":" + modelPrefix;
          }

          // Add the model prefix to the device's model name.
          if( modelPrefix != "" ) modelName = modelPrefix + ":" + modelName;

          device.setModelName( modelName );
        }
        else
        {
          Report::UserError0() << "Unable to find mutual inductor model " << modelName;
          return false;
        }

        // ERK: tmpModel is a persistant object in the CircuitBlockData class.
        // It is cheaper to use the "=" operator than the copy constructor here,
        // as the copy constructor will require allocations and deletions every
        // time this function is called, which is a lot.
        //ParameterBlock model(*MI_modelPtr);
        ParameterBlock tmpModel = (*MI_modelPtr);

        // Set the model parameter values.
        tmpModel.setParameterValues(&circuitBlockPtr_->circuitContext);

        // Add the model to the circuit.
        mainCircuitPtr_->addModel(tmpModel, modelPrefix);
      }
      {
        Device::Param param;

        param.setTag( "COUPLING" );
        param.setVal( MIs[i].coupling );
        param.setGiven(true);
        device.addInstanceParameter( param );
      }

      // add the set of inductors to which it is coupled
      std::map<std::string,double>::iterator iterI = MIs[i].inductors.begin();
      std::map<std::string,double>::iterator  endI = MIs[i].inductors.end();
      for( ; iterI != endI; ++iterI )
      {
        std::string ci = iterI->first;

        if( name != ci )
        {
          Device::Param param;

          param.setTag( "COUPLEDINDUCTOR" );
          if( subcircuitPrefix != "" )
          {
            ci = subcircuitPrefix + ":" + ci;
          }

          param.setVal( ci );
          param.setGiven (true);
          device.addInstanceParameter( param );

          param.setTag("COUPLEDINDUCTANCE");

          if (DEBUG_IO)
            Xyce::dout() << "coupledinductance value from " << name << " to " << ci
                         << " is " << iterI->second << std::endl;

          param.setVal( iterI->second);
          param.setGiven (true);
          device.addInstanceParameter( param );
        }
      }
    }
  }

  return true;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::parsePreprocess
// Purpose       : This function introduces a preprocessing phase whereby
//                 it is determined whether Xyce should remove "redundant"
//                 devices (devices for which all of the nodes are the
//                 same).  This has to be caught before the initial parse
//                 to determine whether to add a device to the circuit or
//                 not.
//
// Special Notes : Updated 12/6/07 for additional detection:  we detect flags
//                 to create netlist files which contain resistors to ground
//                 for nodes which are connected only to one device terminal,
//                 and/or resistors to ground for nodes that have no DC path
//                 to ground.
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/05/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::parsePreprocess(const std::string & netlistFileName)
{
  if (DEBUG_IO)
    Xyce::dout() << "Preprocess of netlist file: " << netlistFileName << std::endl;

  //Get the first line of input.
  std::vector<SpiceSeparatedFieldTool::StringToken> line;
  int eof = !ssfPtr_->getLine(line); //Breaks the line into fields.
  int removecounter = 0;
  int replacecounter = 0;
  int onetermcounter = 0;
  int nodcpathcounter = 0;

  while (!eof)
  {
    if (DEBUG_IO) {
      Xyce::dout() << " After getline, file position is " << ssfPtr_->getFilePosition() << std::endl;
      Xyce::dout() << "preprocess of netlist line: ";
      for (unsigned int i = 0; i < line.size(); ++i)
      {
      Xyce::dout() << line[i].string_ << " ";
      }
      Xyce::dout() << std::endl;
    }

    if ( !line.empty() )
    {
      ExtendedString ES1 ( line[0].string_ );
      ES1.toUpper();

      if ( ES1 != ".PREPROCESS" )
      {
        //do nothing
      }
      else if (line.size() < 3)
      {
        Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
          << "Too few parameters specified in .PREPROCESS statement.";
      }
      else
      {
        ExtendedString preprocarg ( line[1].string_ );
        preprocarg.toUpper();

        if (preprocarg == "REMOVEUNUSED")
        {
          if (removecounter != 0)
          {
            Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
              << "Multiple .PREPROCESS REMOVEUNUSED statements.";
          }
          else
          {
            removecounter++;
            ExtendedString removeparam ( ES1 );
            bool anyparamsremoved = false;

            for (unsigned int i = 2; i < line.size(); ++i)
            {
              removeparam=line[i].string_;
              removeparam.toUpper();
              if (removeparam == "C")
              {
                remove_redundant_C_ = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "D")
              {
                remove_redundant_D_ = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "I")
              {
                remove_redundant_I_ = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "L")
              {
                remove_redundant_L_ = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "M")
              {
                remove_redundant_M_ = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "Q")
              {
                remove_redundant_Q_ = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "R")
              {
                remove_redundant_R_ = true;
                anyparamsremoved = true;
              }
              else if (removeparam == "V")
              {
                remove_redundant_V_ = true;
                anyparamsremoved = true;
              }
              else if (removeparam == ",")
              {
                //skip commas
              }
              else
              {
                Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
                  << "Unknown argument type " << removeparam << " in .PREPROCESS REMOVEUNUSED statement.";
              }
            }
            if (anyparamsremoved == false)
            {
              //didn't find any parameters on the line
              Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
                << "No remove parameters specified in .PREPROCESS REMOVEUNUSED statement.";
            }
          }
        }
        else if (preprocarg == "REPLACEGROUND")
        {
          if (replacecounter != 0)
          {
            Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
              << "Multiple .PREPROCESS REPLACEGROUND statements.";
          }
          else
          {
            replacecounter++;
            if (line.size() > 3)
            {
              Report::UserWarning().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
                << "Additional parameters in .PREPROCESS REPLACEGROUND statement.  Ignoring.";
            }

            ExtendedString replaceparam(line[2].string_);
            replaceparam.toUpper();

            if (replaceparam == "TRUE")
            {
              replace_ground_=true;
            }
            else if (replaceparam != "FALSE")
            {
              Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
                << "Unknown argument " << replaceparam << " in .PREPROCESS REPLACEGROUND statement.";
            }
          }
        }
        else if (preprocarg == "ADDRESISTORS")
        {
          commandLine_.getHangingResistor().setNetlistCopy(true);

          if (line.size() > 4)
          {
            Report::UserWarning().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
              << "Additional parameters in .PREPROCESS ADDRESISTORS statement.  Ignoring.";
          }
          else if (line.size() < 4)
          {
            Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
              << "Missing resistance value in .PREPROCESS ADDRESISTORS statement.";
          }

          ExtendedString netlistparam(line[2].string_);
          netlistparam.toUpper();
          ExtendedString resistanceparam(line[3].string_);
          resistanceparam.toUpper();

          if (netlistparam == "ONETERMINAL")
          {
            if (onetermcounter != 0)
            {
              Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
                << "Multiple .PREPROCESS ADDRESISTORS ONETERMINAL statements.";
            }
            else
            {
              onetermcounter++;
              commandLine_.getHangingResistor().setOneTerm(true);
              commandLine_.getHangingResistor().setOneTermRes(resistanceparam);
            }
          }
          else if (netlistparam == "NODCPATH")
          {

            if (nodcpathcounter != 0)
            {
              Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
                << "Multiple .PREPROCESS ADDRESISTORS NODCPATH statements.";
            }
            else
            {
              nodcpathcounter++;
              commandLine_.getHangingResistor().setNoDCPath(true);
              commandLine_.getHangingResistor().setNoDCPathRes(resistanceparam);
            }
          }
          else
          {
            Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
              <<  "Unknown argument " << netlistparam << " in .PREPROCESS ADDRESISTORS statement.";
          }
        }
        else
        {
          Report::UserError().at(circuitBlockPtr_->netlistFileName, line[0].lineNumber_)
            << "Unknown keyword " <<  preprocarg << " specified in .PREPROCESS statement.";
        }
      }
    }
    eof = !ssfPtr_->getLine(line); //get the next line
  }

  if (DEBUG_IO) {
    Xyce::dout() << std::endl << "Unused components to be removed.  (1 means remove "
                 << "redundancies,"  << std::endl <<" 0 means do not remove redundancies): "
                 << std::endl << std::endl

                 << "Remove Unused C:  " << remove_redundant_C_ << std::endl
                 << "Remove Unused D:  " << remove_redundant_D_ << std::endl
                 << "Remove Unused I:  " << remove_redundant_I_ << std::endl
                 << "Remove Unused L:  " << remove_redundant_L_ << std::endl
                 << "Remove Unused M:  " << remove_redundant_M_ << std::endl
                 << "Remove Unused Q:  " << remove_redundant_Q_ << std::endl
                 << "Remove Unused R:  " << remove_redundant_R_ << std::endl
                 << "Remove Unused V:  " << remove_redundant_V_ << std::endl << std::endl

                 << "Replace Ground Flag set to:  " << replace_ground_ << std::endl << std::endl
                 << "Netlist copy Flag set to:  "
                 << commandLine_.getHangingResistor().getNetlistCopy() << std::endl << std::endl
                 << "One terminal Flag set to:  "
                 << commandLine_.getHangingResistor().getOneTerm() << std::endl << std::endl
                 << "No DC Path Flag set to:  "
                 << commandLine_.getHangingResistor().getNoDCPath() << std::endl << std::endl
                 << "One terminal resistance:  "
                 << commandLine_.getHangingResistor().getOneTermRes() << std::endl << std::endl
                 << "No DC path resistance:  "
                 << commandLine_.getHangingResistor().getNoDCPathRes() << std::endl << std::endl

                 << "Done with preprocess netlist file parsing." << std::endl << std::endl;
  }

  return true;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::removeTwoTerminalDevice
// Purpose       : Given a two terminal device, this function checks to see
//                 if both nodes on the device are the same and, if so,
//                 decides whether or not the device should be removed
//                 from the circuit.  The decision to remove is based upon
//                 whether the specific device type was specified during the
//                 preprocessing phase as a device for which redundancies
//                 should be removed.  E.g., if the lines
//
//                 .PREPROCESS REMOVEUNUSED C
//                 C1 1 1 1
//                 R1 2 2 1
//
//                 appear in the netlist, the capacitor C1 will be removed
//                 from the netlist, whereas the resistor R1 will not.
//
// Special Notes : This function only determines *whether* a specific device
//                 should be removed; the process of removing the device from
//                 the netlist takes place in the handleLinePass1 and
//                 getLinePass2 functions.
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/08/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::removeTwoTerminalDevice(const char linetype, const
ExtendedString & node1, const ExtendedString & node2)
{
  bool result=false;

  if (node1 == node2)
  {
    if (remove_redundant_C_  && linetype == 'C')
    {
      result=true;
    }
    else if (remove_redundant_D_  && linetype == 'D')
    {
      result=true;
    }
    else if (remove_redundant_I_  && linetype == 'I')
    {
      result=true;
    }
    else if (remove_redundant_L_  && linetype == 'L')
    {
      result=true;
    }
    else if (remove_redundant_R_  && linetype == 'R')
    {
      result=true;
    }
    else if (remove_redundant_V_  && linetype == 'V')
    {
      result=true;
    }
  }
  return result;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::removeThreeTerminalDevice
// Purpose       : Given a three terminal device, this function checks to see
//                 if all three nodes on the device are the same and, if so,
//                 decides whether or not the device should be removed
//                 from the circuit.  The decision to remove is based upon
//                 whether the specific device type was specified during the
//                 preprocessing phase as a device for which redundancies
//                 should be removed.  E.g., if the lines
//
//                 .PREPROCESS REMOVEUNUSED M
//                 M1 1 1 1 4 NMOS
//                 Q1 2 2 2 5 NPN
//
//                 appear in the netlist, the MOSFET M1 will be removed
//                 from the netlist, whereas the BJT Q1 will not.
//
// Special Notes : This function only determines *whether* a specific device
//                 should be removed; the process of removing the device from
//                 the netlist takes place in the handleLinePass1 and
//                 getLinePass2 functions.
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/10/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::removeThreeTerminalDevice(const char linetype,
  const ExtendedString & node1,
  const ExtendedString & node2,
  const ExtendedString & node3)
{
  bool result=false;

  if (node1 == node2 && node2 == node3)
  {
    if (remove_redundant_M_ == true  && linetype == 'M')
    {
      result=true;
    }
    else if (remove_redundant_Q_ == true && linetype == 'Q')
    {
      result=true;
    }
  }
  return result;
}
/*
//--------------------------------------------------------------------------
// Function      : CircuitBlockData::removeC
// Purpose       : Accessor function for boolean variable remove_redundant_C
//
// Special Notes :
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/10/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::removeC()
{
  return remove_redundant_C_;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::removeD
// Purpose       : Accessor function for boolean variable remove_redundant_D
//
// Special Notes :
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/10/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::removeD()
{
  return remove_redundant_D_;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::removeI
// Purpose       : Accessor function for boolean variable remove_redundant_I
//
// Special Notes :
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/10/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::removeI()
{
  return remove_redundant_I_;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::removeL
// Purpose       : Accessor function for boolean variable remove_redundant_L
//
// Special Notes :
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/10/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::removeL()
{
  return remove_redundant_L_;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::removeM
// Purpose       : Accessor function for boolean variable remove_redundant_M
//
// Special Notes :
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/10/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::removeM()
{
  return remove_redundant_M_;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::removeQ
// Purpose       : Accessor function for boolean variable remove_redundant_Q
//
// Special Notes :
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/10/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::removeQ()
{
  return remove_redundant_Q_;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::removeR
// Purpose       : Accessor function for boolean variable remove_redundant_R
//
// Special Notes :
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/10/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::removeR()
{
  return remove_redundant_R_;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::removeV
// Purpose       : Accessor function for boolean variable remove_redundant_V
//
// Special Notes :
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/10/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::removeV()
{
  return remove_redundant_V_;
}

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::replgnd
// Purpose       : Accessor function for boolean variable replace_ground_
//
// Special Notes :
//
// Creator       : Keith Santarelli
//
// Creation Date : 10/11/2007
//--------------------------------------------------------------------------
bool CircuitBlockData::replgnd()
{
  return replace_ground_;
}
*/

//--------------------------------------------------------------------------
// Function      : CircuitBlockData::produceUnflattenedNetlist
// Purpose       : Generates a copy of the current netlist in the file
//                 netlistFilename_copy.cir.  This is used to create netlist
//                 files that contain resistors that connect ground to
//                 "dangling" nodes (nodes which either don't have a dc path
//                 to ground, or which are only connected to one device
//                 terminal).
//
// Special Notes :
//
// Creator       : Keith Santarelli, Electrical and Microsystems Modeling
//
// Creation Date : 12/5/07
//--------------------------------------------------------------------------
void CircuitBlockData::produceUnflattenedNetlist(const std::string & netlistFileName)
{

  if (DEBUG_IO)
    Xyce::dout() << "Producing copy of netlist file: " << netlistFileName << std::endl;

  // Create the output file stream
  std::string netlistCopy(netlistFileName);
  netlistCopy += "_xyce.cir";
  std::ofstream copyFile(netlistCopy.c_str());

  // Some error checking in case we can't open the file.
  if(copyFile.fail())
  {
    Report::UserError0() << ".PREPROCESS NETLISTCOPY cannot open output file " << netlistCopy;
    return;
  }


  //Create date/time stamp
  const time_t now = time( NULL );
  char timeDate[ 40 ];
  strftime( timeDate, 40, "TIME='%I:%M:%S %p' DATE='%b %d, %Y' ",
    localtime( &now ) );

  //Create title line (not the same as title line of original netlist file!)
  std::string header("XYCE-generated Netlist file copy:  ");
  header += timeDate;
  copyFile << header << std::endl;

  //Add the original title:
  copyFile << "*Original Netlist Title:  " << std::endl << std::endl;
  copyFile << "*";

  std::vector<SpiceSeparatedFieldTool::StringToken> separatedLine;
  int eof = !ssfPtr_->getLineWithComments(separatedLine);

  //Add the original title text:
  for (unsigned int i = 0; i < separatedLine.size(); ++i)
  {
    copyFile << separatedLine[i].string_;
  }

  copyFile << std::endl;

  eof = !ssfPtr_->getLineWithComments(separatedLine);

  bool endflag=false;
  bool addresistbool=false;
  ExtendedString firstarg("");
  ExtendedString addresistarg("");

  if ( !(separatedLine.empty()) )
  {
    firstarg=separatedLine[0].string_;
    firstarg.toUpper();
  }

  while (!eof && !endflag)
  {
    if (firstarg != ".PREPROCESS")
    {
      for (unsigned int i = 0; i < separatedLine.size(); ++i)
      {
        copyFile << separatedLine[i].string_;
      }
    }
    else
    {
      for (unsigned int i = 0; i < separatedLine.size(); ++i)
      {
        addresistarg = separatedLine[i].string_;
        addresistarg.toUpper();
        if (addresistarg == "ADDRESISTORS")
        {
          addresistbool = true;
        }
      }

      if (!addresistbool)
      {
        for (unsigned int i = 0; i < separatedLine.size(); ++i)
        {
          copyFile << separatedLine[i].string_;
        }
      }
      else
      {
        copyFile << "*";
        for (unsigned int i = 0; i < separatedLine.size(); ++i)
        {
          copyFile << separatedLine[i].string_;
        }
        copyFile << "*Xyce:  \".PREPROCESS ADDRESISTORS\" statement";
        copyFile << " automatically commented out in netlist copy.";
        copyFile << std::endl;
      }
    }

    //Get the next line of input
    eof = !ssfPtr_->getLineWithComments(separatedLine);

    if ( !(separatedLine.empty()) )
    {
      firstarg = separatedLine[0].string_;
    }
    else
    {
      firstarg="";
    }
    firstarg.toUpper();

    //We don't reproduce anything after a .END statement!
    if ( !(separatedLine.empty()) && firstarg == ".END")
    {
      endflag = true;
    }
  }

 copyFile.close();
}

void CircuitBlock::fixupYDeviceNames()
{
  // add Y device name aliases
  std::map<std::string, Teuchos::RefCountPtr<Device::InstanceBlock> >::iterator it_DN = deviceNames_.begin();
  std::map<std::string, Teuchos::RefCountPtr<Device::InstanceBlock> >::const_iterator ite_DN = deviceNames_.end();
  int firstP, lastP;
  for ( ; it_DN != ite_DN; ++it_DN)
  {
    if ((*it_DN).first.find_last_of('%') != std::string::npos)
    {
      firstP =(*it_DN).first.find_first_of('%');
      lastP =(*it_DN).first.find_last_of('%');
      deviceNames_[(*it_DN).first.substr(lastP+1)] = Teuchos::RCP< Device::InstanceBlock >();
      if (firstP != lastP)
      {
       deviceNames_[(*it_DN).first.substr(0, firstP) +
                    (*it_DN).first.substr(firstP+1, lastP-firstP-1) +
                     " " +
                    (*it_DN).first.substr(lastP+1)] = Teuchos::RCP< Device::InstanceBlock >();
      }
    }
  }
}

} // namespace IO
} // namespace Xyce
