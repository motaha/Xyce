#
# A simple python program to convert bugzilla's XML
# to a format Merlin and MS Project can read
#
# usage:
#
# python bugzillaToMerlin.py bugzillaXmlFile
#
# the output is in a file called bugzillaXmlFile_merlin.xml
#
# How do you get bugzilla xml?  Do a search for some bugs in bugzilla and 
# then click on the "xml" button.  In Firefox, you can then save the 
# resulting xml in a file.  In Safari, you must "View Source" and then 
# save that as a file.
# 
# This script takes the latest change date on a bug as the project start
# date to approximate the idea that you're looking at a bunch of issues 
# and trying to see how long it will take to get a critical one done.
#
# A bug's estimated time and actual time are used to guess at how long an
# issue will take to resolve.  If no value was given in the bug report,
# then 100 hours is assumed just to avoid having lots of bugs with zero
# time to fix them in place.
#
# Finally, this script tries to link all bugs that depend on each other in
# the input xml file.  If some of the input bugs depend on bugs that were
# not in the input file, new bugs with the correct bug number but the
# summary title of "unknown" are added to the output.  This way you can
# see that an issue depends on something and go and look that up in bugzilla
# if you want to know more. 

import sys
import string

def makeTextNode( document, parent, tag, text ):
  # makes a text node called "tag" with the data "text"
  # and attaches it to the xml parent "parent"
  element = document.createElement(tag)
  elementText = document.createTextNode(text)
  element.appendChild(elementText)
  parent.appendChild(element)  
  
def insertCalendarXml( parent ):
  # this ugly stuff is a required calendar element that the project
  # files need.  We just parse it and store it in our xml output
  calendarXmlString="""<Calendars>
        <Calendar>
            <UID>1</UID>
            <WeekDays>
                <WeekDay>
                    <DayWorking>0</DayWorking>
                    <DayType>1</DayType>
                </WeekDay>
                <WeekDay>
                    <DayWorking>1</DayWorking>
                    <WorkingTimes>
                        <WorkingTime>
                            <FromTime>08:00:00</FromTime>
                            <ToTime>12:00:00</ToTime>
                        </WorkingTime>
                        <WorkingTime>
                            <FromTime>13:00:00</FromTime>
                            <ToTime>17:00:00</ToTime>
                        </WorkingTime>
                    </WorkingTimes>
                    <DayType>2</DayType>
                </WeekDay>
                <WeekDay>
                    <DayWorking>1</DayWorking>
                    <WorkingTimes>
                        <WorkingTime>
                            <FromTime>08:00:00</FromTime>
                            <ToTime>12:00:00</ToTime>
                        </WorkingTime>
                        <WorkingTime>
                            <FromTime>13:00:00</FromTime>
                            <ToTime>17:00:00</ToTime>
                        </WorkingTime>
                    </WorkingTimes>
                    <DayType>3</DayType>
                </WeekDay>
                <WeekDay>
                    <DayWorking>1</DayWorking>
                    <WorkingTimes>
                        <WorkingTime>
                            <FromTime>08:00:00</FromTime>
                            <ToTime>12:00:00</ToTime>
                        </WorkingTime>
                        <WorkingTime>
                            <FromTime>13:00:00</FromTime>
                            <ToTime>17:00:00</ToTime>
                        </WorkingTime>
                    </WorkingTimes>
                    <DayType>4</DayType>
                </WeekDay>
                <WeekDay>
                    <DayWorking>1</DayWorking>
                    <WorkingTimes>
                        <WorkingTime>
                            <FromTime>08:00:00</FromTime>
                            <ToTime>12:00:00</ToTime>
                        </WorkingTime>
                        <WorkingTime>
                            <FromTime>13:00:00</FromTime>
                            <ToTime>17:00:00</ToTime>
                        </WorkingTime>
                    </WorkingTimes>
                    <DayType>5</DayType>
                </WeekDay>
                <WeekDay>
                    <DayWorking>1</DayWorking>
                    <WorkingTimes>
                        <WorkingTime>
                            <FromTime>08:00:00</FromTime>
                            <ToTime>12:00:00</ToTime>
                        </WorkingTime>
                        <WorkingTime>
                            <FromTime>13:00:00</FromTime>
                            <ToTime>17:00:00</ToTime>
                        </WorkingTime>
                    </WorkingTimes>
                    <DayType>6</DayType>
                </WeekDay>
                <WeekDay>
                    <DayWorking>0</DayWorking>
                    <DayType>7</DayType>
                </WeekDay>
            </WeekDays>
            <Name>Standard</Name>
            <IsBaseCalendar>1</IsBaseCalendar>
            <BaseCalendarUID>-1</BaseCalendarUID>
        </Calendar>
    </Calendars>"""
  calDom = xml.dom.minidom.parseString(calendarXmlString)
  parent.appendChild(calDom.firstChild)
  calDom.unlink()
  
  
if __name__ == '__main__':
  import sys
  import os.path
  import xml.dom.minidom
  
  # sys.argv[0] is the python program's file name
  inputFile = sys.argv[1]
  (baseName,extention) = os.path.splitext(inputFile)
  outputFileName = baseName + "_merlin.xml" 
  outputFile = open( outputFileName, 'w' )
  
  # try to read the input xml file
  inputDom = xml.dom.minidom.parse(inputFile)
  
  # to simplify some of the later processing, search the input xml
  # for a bug with the latest delta time stamp.  Use this bug's timestamp
  # as the start date for the project.
  
  # Also, generate a list of Bug ID's in the input XML file.  When we add
  # the depends_on information to the output XML we must exclude any
  # bugs that are not part of the input, or generate new bugs for these
  # unknow ones.  Otherwise Merlin (and probably Project)
  # will not make the remaing connections

  # get all of the bug elements
  bugList = inputDom.getElementsByTagName("bug")
  projectStartDate=""
  startingBug=""
  inputBugIDList=[]
  for aBug in bugList:
    for bugChild in aBug.childNodes:
      # look for the id and creation time stamp
      
      if bugChild.nodeName == "bug_id":
        bugChild.normalize()
        bugID = bugChild.firstChild.data
        inputBugIDList.append( bugID )
      
      if bugChild.nodeName == "delta_ts":
        bugChild.normalize()
        creationTimeStamp = bugChild.firstChild.data[0:10] + "T" + bugChild.firstChild.data[12:]
        if creationTimeStamp.count(":") < 2:
          creationTimeStamp = creationTimeStamp + ":00"
    
    if (projectStartDate=="") or (projectStartDate < creationTimeStamp):
      projectStartDate = creationTimeStamp
      startingBug = bugID
  
  # ok, now startingBug and projectStartDate are set.
  # print "First bug is ", startingBug, " with time stamp ", projectStartDate
        
  # here's where we'll store the output
  domImp = xml.dom.minidom.getDOMImplementation()
  outputDom = domImp.createDocument(None, None, None)
  outputDomProject = outputDom.createElement("Project")
  outputDomProject.setAttribute("xmlns", "http://schemas.microsoft.com/project")
  outputDom.appendChild( outputDomProject )

  makeTextNode( outputDom, outputDomProject, "Title", "Xyce Bugzilla Items")
  makeTextNode( outputDom, outputDomProject, "ScheduleFromStart", startingBug )
  makeTextNode( outputDom, outputDomProject, "StartDate", projectStartDate )
  makeTextNode( outputDom, outputDomProject, "CalendarUID", "1" )
  insertCalendarXml( outputDomProject )
  
  # all bugs go into a container of Tasks
  outputTaskList = outputDom.createElement("Tasks")
  outputDomProject.appendChild(outputTaskList)
  
  # the first Task is a summary of the whole project
  outputBug = outputDom.createElement("Task")
  outputTaskList.appendChild(outputBug)
  
  makeTextNode( outputDom, outputBug, "UID", "0")   
  makeTextNode( outputDom, outputBug, "ID", "0")
  makeTextNode( outputDom, outputBug, "Name", "Xyce Bugzilla Items")
  makeTextNode( outputDom, outputBug, "OutlineLevel", "0")
  makeTextNode( outputDom, outputBug, "CalendarUID", "1")

  bugCount = 0
  bugsToBeCreated = []
  for aBug in bugList:
    # zero out the data we will extract from this bug xml object
    dependsOnList = []
    bugCount += 1
    estimatedTime = 0
    actualTime = 0
    for bugChild in aBug.childNodes:
      # look for the elements of a bug we want to pass
      # to the output document
      
      if bugChild.nodeName == "bug_id":
        bugChild.normalize()
        bugID = bugChild.firstChild.data
      
      if bugChild.nodeName == "short_desc":
        bugChild.normalize()
        shortDescription = bugChild.firstChild.data
      
      if bugChild.nodeName == "creation_ts":
        bugChild.normalize()
        creationTimeStamp = bugChild.firstChild.data[0:10] + "T" + bugChild.firstChild.data[12:]
        if creationTimeStamp.count(":") < 2:
          creationTimeStamp = creationTimeStamp + ":00"
      
      if bugChild.nodeName == "delta_ts":
        bugChild.normalize()
        deltaTimeStamp = bugChild.firstChild.data[0:10] + "T" + bugChild.firstChild.data[12:]
        if deltaTimeStamp.count(":") < 2:
          deltaTimeStamp = deltaTimeStamp + ":00"
      
      if bugChild.nodeName == "reporter":
        reporter = bugChild.getAttribute("name")
        
      if bugChild.nodeName == "assigned_to":
        assignedTo = bugChild.getAttribute("name")
        
      if bugChild.nodeName == "dependson":
        bugChild.normalize()
        dep = bugChild.firstChild.data
        dependsOnList.append( dep )
        if( inputBugIDList.count( dep ) == 0 ):
          # the dependency is not in the input list, so create it later
          bugsToBeCreated.append( dep )       
        
      if bugChild.nodeName == "estimated_time":
        bugChild.normalize()
        estimatedTime = bugChild.firstChild.data
      
      if bugChild.nodeName == "remaining_time":
        bugChild.normalize()
        remainingTime = bugChild.firstChild.data
      
      if bugChild.nodeName == "actual_time":
        bugChild.normalize()
        actualTime = bugChild.firstChild.data
      
    # ok we have the bugzilla info, now save it in the
    # output xml data structure
    outputBug = outputDom.createElement("Task")
    outputTaskList.appendChild(outputBug)
    
    # convert the bug count to a string value
    bugCountAsString = '%s' % bugCount
    
    # I tried to use the creation time and delta timestamp to estimate
    # actual time taken, but this doesn't work well when times get stacked up
    
    # use the difference in time from creationTimeStamp and deltaTimeStamp
    # to get duration estimate
    #tsCreation = datetime(*(strptime(creationTimeStamp, "%Y-%m-%dT%H:%M:%S")[0:6]))
    #tsDelta = datetime(*(strptime(deltaTimeStamp, "%Y-%m-%dT%H:%M:%S")[0:6]))
    #
    # in Python 2.5 we can just do this
    #tsCreation = datetime.strptime(creationTimeStamp, "%Y-%m-%dT%H:%M:%S" )
    #tsDelta = datetime.strptime(deltaTimeStamp, "%Y-%m-%dT%H:%M:%S" )
    #
    #duration = tsDelta - tsCreation
    #hours = int(duration.seconds / (60*60))
    #minutes = int( (duration.seconds - 60*60*hours)/60 )
    #seconds = duration.seconds - 60*60*hours - 60*minutes
    #hours = hours + 24 * duration.days
    #durationString = 'PT%sH%sM%sS' % ( hours, minutes, seconds )
    
    duration = max( [ float(estimatedTime), float(actualTime) ] ) 
    if duration <= 0:
      duration = 100
    durationString = 'PT%sH%sM%sS' % ( duration, 0, 0 )
    
    makeTextNode( outputDom, outputBug, "UID", bugID )     
    makeTextNode( outputDom, outputBug, "ID", bugCountAsString )
    makeTextNode( outputDom, outputBug, "Name", "Bug " + bugID + ":  " + shortDescription )
    makeTextNode( outputDom, outputBug, "OutlineLevel", "1" )
    makeTextNode( outputDom, outputBug, "CalendarUID", "1" )
    makeTextNode( outputDom, outputBug, "Start", deltaTimeStamp )
    makeTextNode( outputDom, outputBug, "Duration", durationString )
    
    # attached depends on data if it is there
    if len(dependsOnList) > 0:  
      for dep in dependsOnList:
        outputPredElement = outputDom.createElement("PredecessorLink")
        outputBug.appendChild(outputPredElement)
        makeTextNode( outputDom, outputPredElement, "PredecessorUID", dep )
        makeTextNode( outputDom, outputPredElement, "Type", "1" )
  
  # cycle through bugsToBeCreated list and make any bugs needed for 
  # dependency info
  for aBug in bugsToBeCreated:
    bugCount += 1
    outputBug = outputDom.createElement("Task")
    outputTaskList.appendChild(outputBug)
    
    # convert the bug count to a string value
    bugCountAsString = '%s' % bugCount
    
    makeTextNode( outputDom, outputBug, "UID", aBug )      
    
    makeTextNode( outputDom, outputBug, "UID", bugID )     
    makeTextNode( outputDom, outputBug, "ID", bugCountAsString )
    makeTextNode( outputDom, outputBug, "Name", "Bug " + aBug + ":  unknown" )
    makeTextNode( outputDom, outputBug, "OutlineLevel", "1" )
    makeTextNode( outputDom, outputBug, "CalendarUID", "1" )
    makeTextNode( outputDom, outputBug, "Start", projectStartDate )
    makeTextNode( outputDom, outputBug, "Duration", "PT40H0M0S" )
  
  #outputFile.write( outputDom.toprettyxml("  ") )
  outputFile.write( outputDom.toxml() )
  outputFile.close()
  
  # need to clean up dom objects when done
  inputDom.unlink()
  outputDom.unlink()