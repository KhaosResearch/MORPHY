//
// File: ApplicationTools.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 21 16:19 2005
// From old file created on: Sun Dec 14 09:36:26 2003
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

This software is a computer program whose purpose is to provide basal and 
utilitary classes. This file belongs to the Bio++ Project.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "ApplicationTools.h"

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/******************************************************************************/

OutputStream* ApplicationTools::error   = new StdErr();
OutputStream* ApplicationTools::message = new StdOut();
OutputStream* ApplicationTools::warning = new StdOut();
time_t ApplicationTools::startTime;
size_t ApplicationTools::terminalWidth = 80;
float ApplicationTools::terminalSplit = 0.5;
bool ApplicationTools::interactive = true;

/******************************************************************************/

bool ApplicationTools::parameterExists(
  const std::string & parameterName,
  std::map<std::string, std::string>& params)
{
  return (params.find(parameterName) != params.end() && !TextTools::isEmpty(params[parameterName]));
}
  
/******************************************************************************/

std::string ApplicationTools::getAFilePath(
  const std::string& parameter,
  std::map<std::string, std::string>& params,
  bool isRequired,
  bool mustExist,
  const std::string& suffix,
  bool suffixIsOptional) throw (Exception)
{
  string filePath = getStringParameter(parameter, params, "none", suffix, suffixIsOptional, false);
  if (filePath == "") filePath = "none";
  if (filePath == "none" && isRequired)
  {
    throw Exception("You must specify a file for this parameter: " + parameter + (suffixIsOptional ? "" : suffix));
  }
  if(filePath == "none") return filePath;
  if(mustExist && !FileTools::fileExists(filePath))
  {
    throw Exception("File does not exists: " + filePath);
  }
  return filePath;
}

/******************************************************************************/

double ApplicationTools::getDoubleParameter(
  const std::string & parameterName,
  std::map<std::string, std::string> & params,
  double defaultValue,
  const std::string & suffix,
  bool suffixIsOptional,
  bool warn)
{
  double dParam = defaultValue;
  if (parameterExists(parameterName + suffix, params))
  {
    dParam = TextTools::toDouble(params[parameterName + suffix]);
  }
  else if (suffixIsOptional && parameterExists(parameterName, params))
  {
    dParam = TextTools::toDouble(params[parameterName]);
  }
  else if(warn)
  {
    displayWarning("Parameter " + parameterName + suffix + " not specified. Default used instead: " + TextTools::toString(defaultValue));
  }
  return dParam;
}

/******************************************************************************/

int ApplicationTools::getIntParameter(
  const std::string & parameterName,
  std::map<std::string, std::string> & params,
  int defaultValue,
  const std::string & suffix,
  bool suffixIsOptional,
  bool warn)
{
  int iParam = defaultValue;
  if (parameterExists(parameterName + suffix, params)) {
    iParam = TextTools::toInt(params[parameterName + suffix]);
  } else if(suffixIsOptional && parameterExists(parameterName, params)) {
    iParam = TextTools::toInt(params[parameterName]);
  } else if (warn) {
    displayWarning("Parameter " + parameterName + suffix + " not specified. Default used instead: " + TextTools::toString(defaultValue));
  }
  return iParam;
}

/******************************************************************************/

std::string ApplicationTools::getStringParameter(
  const std::string& parameterName,
  std::map<std::string, std::string>& params,
  const std::string& defaultValue,
  const std::string& suffix,
  bool suffixIsOptional,
  bool warn)
{
  string sParam = defaultValue;
  if (parameterExists(parameterName + suffix, params)) {
    sParam = params[parameterName + suffix];
  } else if (suffixIsOptional && parameterExists(parameterName, params)) {
    sParam = params[parameterName];
  } else if (warn) {
    displayWarning("Parameter " + parameterName + " not specified. Default used instead: " + defaultValue);
  }
  return sParam;
}

/******************************************************************************/

bool ApplicationTools::getBooleanParameter(
  const std::string& parameterName,
  std::map<std::string, std::string>& params,
  bool defaultValue,
  const std::string& suffix,
  bool suffixIsOptional,
  bool warn)
{
  string sParam;
  bool bParam = defaultValue;
  if (parameterExists(parameterName + suffix, params))
  {
    sParam = params[parameterName + suffix];
  }
  else if (suffixIsOptional && parameterExists(parameterName, params))
  {
    sParam = params[parameterName];
  }
  else {
    if (warn)
    {
      displayWarning("Parameter " + parameterName + " not specified. Default used instead: " + TextTools::toString(defaultValue));
    }
    return bParam;
  }
  if ((sParam == "true") 
   || (sParam == "TRUE")
   || (sParam == "t")
   || (sParam == "T")
   || (sParam == "yes")
   || (sParam == "YES")
   || (sParam == "y")
   || (sParam == "Y")
   || (sParam == "1"))
    bParam = true;
  else if ((sParam == "false") 
        || (sParam == "FALSE")
        || (sParam == "f")
        || (sParam == "F")
        || (sParam == "no")
        || (sParam == "NO")
        || (sParam == "n")
        || (sParam == "N")
        || (sParam == "0"))
    bParam = false;
  else throw Exception("ApplicationTools::getBooleanParameter. Wrong description:" + sParam);
 
  return bParam;
}

/******************************************************************************/

void ApplicationTools::displayMessage(const std::string& text) { if(message) (*message << text).endLine(); }
    
void ApplicationTools::displayError(const std::string& text) { if(error) (*error << "ERROR!!! " << text).endLine(); }
    
void ApplicationTools::displayWarning(const std::string& text) { if(warning) (*warning << "WARNING!!! " << text).endLine(); }

void ApplicationTools::displayTask(const std::string& text, bool eof)
{
  if (message)
  {
    *message << TextTools::resizeRight(text, static_cast<size_t>(static_cast<float>(terminalWidth) * terminalSplit - 1), '.') << ": ";
    if (eof) message->endLine();
    else     message->flush();
  }
}
    
void ApplicationTools::displayTaskDone() { if(message) (*message << "Done.").endLine(); }

/******************************************************************************/

void ApplicationTools::displayGauge(size_t iter, size_t total, char symbol, const std::string& mes)
{
  if (!message) return;
  if (total == 0) return;//We do not display anything in that case.
  size_t width = static_cast<size_t>(static_cast<float>(terminalWidth) * terminalSplit - 2);
  string gauge = string(static_cast<size_t>((1. * static_cast<double>(iter) / static_cast<double>(total)) * static_cast<double>(width)), symbol);
  string info = mes;
  if (interactive)
  {
    string fill = string(width - gauge.length(), ' ');
    gauge = "[" + gauge + fill + "] " + TextTools::resizeLeft(TextTools::toString(100 * iter / total), 3) + "%";
    if (mes.size() > terminalWidth - gauge.size())
      info = TextTools::resizeRight(mes, terminalWidth - gauge.size());
    *message << '\r' + info + gauge;
    message->flush();
  }
  else
  {
    if (iter == 0)
    {
      *message << "[";
      message->flush();
      return;
    }
    size_t step = static_cast<size_t>(ceil(1. * static_cast<double>(total) / static_cast<double>(width)));
    if (iter >= total)
    {
      size_t fill = static_cast<size_t>(static_cast<float>(terminalWidth) * terminalSplit) - (total - 1) / step - 1;
      *message << TextTools::resizeLeft("]", fill, symbol);
      message->flush();
      return;
    }
    size_t x = iter % step;
    if (x == 0) { *message << symbol; message->flush(); }
  }
}

/******************************************************************************/

void ApplicationTools::displayUnlimitedGauge(size_t iter, const std::string& mes)
{
  if (!message) return;
  string chars = "-/-\\";
  string info = mes;
  if (interactive)
  {
    unsigned int i = iter % 4;
    *message << '\r' << info << chars[i] << " " << TextTools::toString(iter);
    message->flush();
  }
  else
  {
    if (iter == 0)
      *message << info;
    *message << "*";
    message->flush();
    return;
  }
}

/******************************************************************************/
  
void ApplicationTools::displayTime(const std::string& msg)
{
  time_t endTime;
  time(&endTime);
  if (message)
  {
    double nsec = difftime(endTime, startTime);
    double nmin = floor(nsec/60.);
    double nhou = floor(nmin/60.);
    double nday = floor(nhou/24.);
    nhou = nhou - nday * 24;
    nmin = nmin - (nday * 24 + nhou) * 60;
    nsec = nsec - ((nday * 24 + nhou) * 60 + nmin) * 60;
    *message << msg << " ";
    *message << nday << "d, ";
    *message << nhou << "h, ";
    *message << nmin << "m, ";
    *message << nsec << "s.";
    message->endLine();
  }
}

/******************************************************************************/
  
double ApplicationTools::getTime()
{
  time_t endTime;
  time(&endTime);
  return difftime(endTime, startTime);
}

/******************************************************************************/

