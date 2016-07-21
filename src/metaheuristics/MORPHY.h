//  MORPHY.h
//
//  Author:
//       Cristian Zambrano-Vega <czambrano@uteq.edu.ec>
//
//  Copyright (c) 2016 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __MORPHY__
#define __MORPHY__

#include <Algorithm.h>
#include <Problem.h>
#include <SolutionSet.h>
#include <Distance.h>
#include <Ranking.h>
#include <CrowdingComparator.h>

//#include <QualityIndicator.h>

/**
  * @class MORPHY
  * @brief This class implements the MORPHY algorithm
**/

class MORPHY : public Algorithm {

public:
  MORPHY(Problem * problem);
  double RFDistance(Solution * sol1, Solution * sol2, bool normalized, int NumberOfTaxas);
  SolutionSet * execute();

};


#endif /* __MORPHY__ */
