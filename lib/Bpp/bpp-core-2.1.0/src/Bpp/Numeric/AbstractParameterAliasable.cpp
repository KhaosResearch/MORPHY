//
// File: AbstractParameterAliasable.cpp
// Created by: Julien Dutheil
// Created on: Sun Mar 29 09:10 2009
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 19, 2004)

   This software is a computer program whose purpose is to provide classes
   for numerical calculus.

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

#include "AbstractParameterAliasable.h"
#include "VectorTools.h"

using namespace bpp;
using namespace std;

AbstractParameterAliasable::AbstractParameterAliasable(const AbstractParameterAliasable& ap) :
  AbstractParametrizable(ap),
  independentParameters_(ap.independentParameters_),
  aliasListenersRegister_()
{
  // Actualize the register with adequate pointers:
  for (map<string, AliasParameterListener*>::const_iterator it = ap.aliasListenersRegister_.begin();
       it != ap.aliasListenersRegister_.end();
       it++)
  {
    AliasParameterListener* listener = it->second->clone();
    listener->setParameterList(&getParameters_());
    aliasListenersRegister_[it->first] = listener;
    //Now correct parameters with appropriate pointers:
    for (unsigned int i = 0; i < getNumberOfParameters(); ++i) {
      if (getParameters_()[i].hasParameterListener(it->first)) {
        getParameters_()[i].removeParameterListener(it->first);
        getParameters_()[i].addParameterListener(listener, false);
      }
    }
  }
}

AbstractParameterAliasable& AbstractParameterAliasable::operator=(const AbstractParameterAliasable& ap)
{
  independentParameters_ = ap.independentParameters_;

  // Actualize the register with adequate pointers:
  for (map<string, AliasParameterListener*>::const_iterator it = ap.aliasListenersRegister_.begin();
       it != ap.aliasListenersRegister_.end();
       it++)
  {
    AliasParameterListener* listener = it->second->clone();
    listener->setParameterList(&getParameters_());
    aliasListenersRegister_[it->first] = listener;
    //Now correct parameters with appropriate pointers:
    for (unsigned int i = 0; i < getNumberOfParameters(); ++i) {
      if (getParameters_()[i].hasParameterListener(it->first)) {
        getParameters_()[i].removeParameterListener(it->first);
        getParameters_()[i].addParameterListener(listener, false);
      }
    }
  }
  return *this;
}

AbstractParameterAliasable::~AbstractParameterAliasable()
{
  // Delete the registry content:
  for (map<string, AliasParameterListener*>::iterator it = aliasListenersRegister_.begin();
       it != aliasListenersRegister_.end();
       it++)
  {
    delete it->second;
  }
}

void AbstractParameterAliasable::aliasParameters(const std::string& p1, const std::string& p2)
throw (ParameterNotFoundException, Exception)
{
  // In case this is the first time we call this method:
  if (getNumberOfParameters() > 0 && independentParameters_.size() == 0)
    independentParameters_ = getParameters();

  if (!hasParameter(p1))
    throw ParameterNotFoundException("AbstractParameterAliasable::aliasParameters", p1);
  if (!hasParameter(p2))
    throw ParameterNotFoundException("AbstractParameterAliasable::aliasParameters", p2);
  if (!independentParameters_.hasParameter(getNamespace() + p2))
    throw Exception("AbstractParameterAliasable::aliasParameters. Parameter " + p2 + " is already aliased to a parameter and can't be aliased twice.");

  string id = "__alias_" + p2 + "_to_" + p1;
  string idCheck = "__alias_" + p1 + "_to_" + p2;
  if (aliasListenersRegister_.find(idCheck) != aliasListenersRegister_.end())
    throw Exception("AbstractParameterAliasable::aliasParameters. Trying to alias parameter " + p2 + " to " + p1 + ", but parameter " + p1 + " is already aliased to parameter " + p2 + ".");
  Parameter* param1 = &getParameter_(p1);
  Parameter* param2 = &getParameter_(p2);
  if (!param1->hasConstraint())
  {
    if (param2->hasConstraint())
      throw Exception("AbstractParameterAliasable::aliasParameters. Cannot alias parameter " + p2 + " to " + p1 + ", because the constraints attached to these two parameters are different.");
  }
  else
  // We use a small trick here, we test the constraints on the basis of their string description (C++ does not provide a default operator==() :( ).
  if (param2->hasConstraint() && (param1->getConstraint()->getDescription() != param2->getConstraint()->getDescription()))
    throw Exception("AbstractParameterAliasable::aliasParameters. Cannot alias parameter " + p2 + " to " + p1 + ", because the constraints attached to these two parameters are different.");

  // Every thing seems ok, let's create the listener and register it:
  AliasParameterListener* aliasListener = new AliasParameterListener(id, getParameters().whichParameterHasName(getNamespace() + p2), &getParameters_());
  aliasListenersRegister_[id] = aliasListener;
  // Now we add it to the appropriate parameter, that is p1.
  // The parameter will not own the listener, the bookkeeping being achieved by the register:
  param1->addParameterListener(aliasListener, false);
  // Finally we remove p2 from the list of independent parameters:
  independentParameters_.deleteParameter(getNamespace() + p2);
}

void AbstractParameterAliasable::unaliasParameters(const std::string& p1, const std::string& p2)
throw (ParameterNotFoundException, Exception)
{
  if (!hasParameter(p1))
    throw ParameterNotFoundException("AbstractParameterAliasable::unaliasParameters", p1);
  if (!hasParameter(p2))
    throw ParameterNotFoundException("AbstractParameterAliasable::unaliasParameters", p2);

  string id = "__alias_" + p2 + "_to_" + p1;
  map<string, AliasParameterListener*>::iterator it = aliasListenersRegister_.find(id);
  if (it == aliasListenersRegister_.end())
    throw Exception("AbstractParameterAliasable::unaliasParameters. Parameter " + p2 + " is not aliased to parameter " + p1 + ".");
  // Remove the listener:
  getParameter_(p1).removeParameterListener(id);
  delete it->second;
  aliasListenersRegister_.erase(it);
  // Finally we re-add p2 to the list of independent parameters:
  independentParameters_.addParameter(getParameter(p2));
}

void AbstractParameterAliasable::setNamespace(const std::string& prefix)
{
  string currentName;
  // First we correct the independent parameter list
  for (unsigned int i = 0; i < independentParameters_.size(); i++)
  {
    currentName = independentParameters_[i].getName();
    if (TextTools::startsWith(currentName, getNamespace()))
      independentParameters_[i].setName(prefix + currentName.substr(getNamespace().size()));
    else
      independentParameters_[i].setName(prefix + currentName);
  }
  // Then we modify all the listeners
  for (map<string, AliasParameterListener*>::iterator it = aliasListenersRegister_.begin();
       it != aliasListenersRegister_.end();
       it++)
  {
    currentName = it->second->getName();
    if (TextTools::startsWith(currentName, getNamespace()))
      it->second->rename(prefix + currentName.substr(getNamespace().size()));
    else
      it->second->rename(prefix + currentName);
  }
  // Finally we notify the mother class:
  AbstractParametrizable::setNamespace(prefix);
}

vector<string> AbstractParameterAliasable::getAlias(const string& name) const
{
  vector<string> aliases;
  for (map<string, AliasParameterListener*>::const_iterator it = aliasListenersRegister_.begin();
       it != aliasListenersRegister_.end();
       it++)
  {
    if (it->second->getName() == name)
    {
      string alias = it->second->getAlias();
      aliases.push_back(alias);
      vector<string> chainAliases = getAlias(alias);
      VectorTools::append(aliases, chainAliases);
    }
  }
  return aliases;
}

