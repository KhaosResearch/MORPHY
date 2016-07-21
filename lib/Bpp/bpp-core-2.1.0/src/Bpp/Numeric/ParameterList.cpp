//
// File: ParameterList.cpp
// Created by: Julien Dutheil
// Created on: Wed Oct 15 18:17:29 2003
//

/*
   Copyright or © or Copr. Julien Dutheil, (November 19, 2004)

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

#include "ParameterList.h"

using namespace bpp;

#include <iostream>
#include <algorithm>
using namespace std;

/** Copy constructor: *********************************************************/

ParameterList::ParameterList(const ParameterList& pl) :
  parameters_(pl.size())
{
  // Now copy all parameters:
  for (unsigned int i = 0; i < size(); i++)
  {
    parameters_[i] = dynamic_cast<Parameter*>(pl.parameters_[i]->clone());
  }
}

/** Assignation operator: *****************************************************/

ParameterList& ParameterList::operator=(const ParameterList& pl)
{
  // First delete all parameters:
  reset();

  // Then resize the vector:
  parameters_.resize(pl.size());

  // Now copy all parameters:
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    parameters_[i] = dynamic_cast<Parameter*>(pl.parameters_[i]->clone());
  }

  return *this;
}

/** Destructor: ***************************************************************/

ParameterList::~ParameterList()
{
  // Delete all parameter objects.
  reset();
}

/******************************************************************************/
const Parameter& ParameterList::getParameter(const std::string& name) const throw (ParameterNotFoundException)
{
  for (unsigned int i = 0; i < size(); i++)
  {
    const Parameter* p = parameters_[i];
    if (p->getName() == name) return *p;
  }
  throw ParameterNotFoundException("ParameterList::getParameter('name').", name);
}

/******************************************************************************/
double ParameterList::getParameterValue(const std::string& name) const throw (ParameterNotFoundException)
{
  for (unsigned int i = 0; i < size(); i++)
  {
    const Parameter* p = parameters_[i];
    if (p->getName() == name) return p->getValue();
  }
  throw ParameterNotFoundException("ParameterList::getParameterValue('name').", name);
}

/******************************************************************************/
Parameter& ParameterList::getParameter(const std::string& name) throw (ParameterNotFoundException)
{
  for (unsigned int i = 0; i < size(); i++)
  {
    Parameter* p = parameters_[i];
    if (p->getName() == name) return *p;
  }
  throw ParameterNotFoundException("ParameterList::getParameter('name').", name);
}

/******************************************************************************/
ParameterList ParameterList::subList(const std::vector<std::string>& names) const throw (ParameterNotFoundException)
{
  ParameterList pl;
  for (unsigned int i = 0; i < names.size(); i++)
  {
    Parameter param = getParameter(names[i]);
    pl.addParameter(param);
  }
  return pl;
}

/******************************************************************************/
ParameterList ParameterList::subList(const std::string& name) const throw (ParameterNotFoundException)
{
  ParameterList pl;
  Parameter param = getParameter(name);
  pl.addParameter(param);
  return pl;
}

/******************************************************************************/
ParameterList ParameterList::subList(const std::vector<size_t>& parameters) const
{
  ParameterList pl;
  for (unsigned int i = 0; i < parameters.size(); i++)
  {
    if (parameters[i] < size()) pl.parameters_.push_back(dynamic_cast<Parameter*>(parameters_[parameters[i]]->clone()));
  }
  return pl;
}

/******************************************************************************/
ParameterList ParameterList::subList(size_t parameter) const
{
  ParameterList pl;
  if (parameter < size()) pl.parameters_.push_back(dynamic_cast<Parameter*>(parameters_[parameter]->clone()));
  return pl;
}

/******************************************************************************/
ParameterList ParameterList::getCommonParametersWith(const ParameterList& params) const
{
  ParameterList pl;
  for (unsigned int i = 0; i < params.size(); i++)
  {
    const Parameter& p = params[i];
    if (hasParameter(p.getName()))
      pl.parameters_.push_back(dynamic_cast<Parameter*>(p.clone()));                                                                        
    // We use push_back instead of addParameter because we are sure the name is not duplicated.
  }

  return pl;
}

/******************************************************************************/

std::vector<std::string> ParameterList::getParameterNames() const
{
  vector<string> pNames(size());
  for (unsigned int i = 0; i < size(); i++)
  {
    pNames[i] = parameters_[i]->getName();
  }
  return pNames;
}

/******************************************************************************/

void ParameterList::addParameter(const Parameter& param) throw (ParameterException)
{
  if (hasParameter(param.getName()))
    throw ParameterException("ParameterList::addParameter. Parameter with name '" + param.getName() + "' already exists.", &param);
  parameters_.push_back(dynamic_cast<Parameter*>(param.clone()));
}

/******************************************************************************/

void ParameterList::addParameter(Parameter* param) throw (ParameterException)
{
  if (hasParameter(param->getName()))
    throw ParameterException("ParameterList::addParameter. Parameter with name '" + param->getName() + "' already exists.", param);
  parameters_.push_back(param);
}

/******************************************************************************/

void ParameterList::setParameter(size_t index, const Parameter& param) throw (IndexOutOfBoundsException)
{
  if (index >= size()) throw IndexOutOfBoundsException("ParameterList::setParameter.", index, 0, size());
  delete parameters_[index];
  parameters_[index] = dynamic_cast<Parameter*>(param.clone());
}

/******************************************************************************/

void ParameterList::includeParameters(const ParameterList& params)
{
  for (unsigned int i = 0; i < params.size(); i++)
  {
    if (hasParameter(params[i].getName()))
      setParameterValue(params[i].getName(), params[i].getValue());
    else
      parameters_.push_back(dynamic_cast<Parameter*>(params[i].clone()));
  }
}

/******************************************************************************/

void ParameterList::addParameters(const ParameterList& params)
throw (ParameterException)
{
  for (unsigned int i = 0; i < params.size(); i++)
  {
    addParameter(params[i]);
  }
}

/******************************************************************************/

void ParameterList::setParameterValue(const string& name, double value)
throw (ParameterNotFoundException, ConstraintException)
{
  Parameter* p = &getParameter(name);
  p->setValue(value);
}

/******************************************************************************/

void ParameterList::setAllParametersValues(const ParameterList& params)
throw (ParameterNotFoundException, ConstraintException)
{
  // First we check if all values are correct:
  for (vector<Parameter*>::iterator it = parameters_.begin(); it < parameters_.end(); it++)
  {
    const Parameter* p = &params.getParameter((*it)->getName());
    if ((*it)->hasConstraint() && !(*it)->getConstraint()->isCorrect(p->getValue()))
      throw ConstraintException("ParameterList::setParametersValues()", *it, p->getValue());
  }

  // If all values are ok, we set them:
  for (vector<Parameter*>::iterator it = parameters_.begin(); it < parameters_.end(); it++)
  {
    const Parameter* p = &params.getParameter((*it)->getName());
    (*it)->setValue(p->getValue());
  }
}

/******************************************************************************/

void ParameterList::setParametersValues(const ParameterList& params)
{
  // First we check if all values are correct:
  for (vector<Parameter*>::const_iterator it = params.parameters_.begin(); it < params.parameters_.end(); it++)
  {
    if (hasParameter((*it)->getName()))
    {
      Parameter* p = &getParameter((*it)->getName());
      if (p->hasConstraint() && !p->getConstraint()->isCorrect((*it)->getValue()))
        throw ConstraintException("ParameterList::setParametersValues()", p, (*it)->getValue());
    }
  }

  // If all values are ok, we set them:
  {
    for (vector<Parameter*>::const_iterator it = params.parameters_.begin(); it < params.parameters_.end(); it++)
    {
      if (hasParameter((*it)->getName()))
      {
        Parameter* p = &getParameter((*it)->getName());
        p->setValue((*it)->getValue());
      }
    }
  }
}

/******************************************************************************/

bool ParameterList::testParametersValues(const ParameterList& params) const
{
  // First we check if all values are correct:
  for (vector<Parameter*>::const_iterator it = params.parameters_.begin(); it < params.parameters_.end(); it++)
  {
    if (hasParameter((*it)->getName()))
    {
      const Parameter* p = &getParameter((*it)->getName());
      if (p->hasConstraint() && !p->getConstraint()->isCorrect((*it)->getValue()))
        throw ConstraintException("ParameterList::matchParametersValues()", p, (*it)->getValue());
    }
  }

  // If all values are ok, we test them:
  bool ch = 0;

  for (vector<Parameter*>::const_iterator it = params.parameters_.begin(); it < params.parameters_.end(); it++)
  {
    if (hasParameter((*it)->getName()))
    {
      const Parameter* p = &getParameter((*it)->getName());
      if (p->getValue() != (*it)->getValue())
        ch |= 1;
    }
  }
  return ch;
}

/******************************************************************************/

bool ParameterList::matchParametersValues(const ParameterList& params)
throw (ConstraintException)
{
  // First we check if all values are correct:
  for (vector<Parameter*>::const_iterator it = params.parameters_.begin(); it < params.parameters_.end(); it++)
  {
    if (hasParameter((*it)->getName()))
    {
      Parameter* p = &getParameter((*it)->getName());
      if (p->hasConstraint() && !p->getConstraint()->isCorrect((*it)->getValue()))
        throw ConstraintException("ParameterList::matchParametersValues()", p, (*it)->getValue());
    }
  }

  // If all values are ok, we set them:
  bool ch = 0;

  for (vector<Parameter*>::const_iterator it = params.parameters_.begin(); it < params.parameters_.end(); it++)
  {
    if (hasParameter((*it)->getName()))
    {
      Parameter* p = &getParameter((*it)->getName());
      if (p->getValue() != (*it)->getValue())
        ch |= 1;
      p->setValue((*it)->getValue());
    }
  }
  return ch;
}

/******************************************************************************/
void ParameterList::setAllParameters(const ParameterList& params)
throw (ParameterNotFoundException)
{
  for (vector<Parameter*>::iterator it = parameters_.begin(); it < parameters_.end(); it++)
  {
    const Parameter* p = &params.getParameter((*it)->getName());
    **it = *p;
  }
}

/******************************************************************************/
void ParameterList::setParameters(const ParameterList& params)
throw (ParameterNotFoundException)
{
  for (vector<Parameter*>::const_iterator it = params.parameters_.begin(); it < params.parameters_.end(); it++)
  {
    Parameter* p = &getParameter((*it)->getName());
    *p = **it;
  }
}

/******************************************************************************/
bool ParameterList::hasParameter(const std::string& name) const
{
  for (unsigned int i = 0; i < size(); i++)
  {
    const Parameter* p = parameters_[i];
    if (p->getName() == name)
      return true;
  }
  return false;
}

/******************************************************************************/
void ParameterList::matchParameters(const ParameterList& params)
{
  for (vector<Parameter*>::const_iterator it = params.parameters_.begin(); it < params.parameters_.end(); it++)
  {
    if (hasParameter((*it)->getName()))
    {
      Parameter* p = &getParameter((*it)->getName());
      *p = **it;
    }
  }
}

/******************************************************************************/
void ParameterList::deleteParameter(const std::string& name) throw (ParameterNotFoundException)
{
  for (unsigned int i = 0; i < size(); i++)
  {
    Parameter* p = parameters_[i];
    if (p->getName() == name)
    {
      delete p;
      parameters_.erase(parameters_.begin() + i);
      return;
    }
  }
  throw ParameterNotFoundException("ParameterList::deleteParameter", name);
}

/******************************************************************************/
void ParameterList::deleteParameters(const std::vector<std::string>& names) throw (ParameterNotFoundException)
{
  for (unsigned int i = 0; i < names.size(); i++)
  {
    deleteParameter(names[i]);
  }
}

/******************************************************************************/
void ParameterList::deleteParameter(size_t index) throw (IndexOutOfBoundsException)
{
  if (index >= size()) throw IndexOutOfBoundsException("ParameterList::deleteParameter.", index, 0, size());
  Parameter* p = parameters_[index];
  delete p;
  parameters_.erase(parameters_.begin() + index);
}

/******************************************************************************/
void ParameterList::deleteParameters(const std::vector<size_t>& indices) throw (IndexOutOfBoundsException)
{
  vector<size_t> tmp(indices);
  sort(tmp.begin(), tmp.end());
  for (vector<size_t>::reverse_iterator i = tmp.rbegin(); i != tmp.rend(); i++)
  {
    size_t index = *i;
    if (index >= size()) throw IndexOutOfBoundsException("ParameterList::deleteParameter.", index, 0, size());
    Parameter* p = parameters_[index];
    delete p;
    parameters_.erase(parameters_.begin() + index);
  }
}

/******************************************************************************/
size_t ParameterList::whichParameterHasName(const std::string& name) const throw (ParameterNotFoundException)
{
  for (size_t i = 0; i < size(); i++)
  {
    if (parameters_[i]->getName() == name) return i;
  }
  throw ParameterNotFoundException("ParameterList::whichParameterHasName.", name);
}

/******************************************************************************/
void ParameterList::printParameters(OutputStream& out) const
{
  (out << "Name:\tValue:\tConstraint:").endLine();
  (out << "_________________________________________________").endLine();
  for (unsigned int i = 0; i < size(); i++)
  {
    out << parameters_[i]->getName();
    out << "\t" << parameters_[i]->getValue();
    out << (parameters_[i]->hasConstraint() ? "\t" + parameters_[i]->getConstraint()->getDescription() : string(""));
    out.endLine();
  }
}

/******************************************************************************/
void ParameterList::reset()
{
  for (unsigned int i = 0; i < size(); i++)
  {
    delete parameters_[i];
  }
  parameters_.resize(0);
}

/******************************************************************************/

