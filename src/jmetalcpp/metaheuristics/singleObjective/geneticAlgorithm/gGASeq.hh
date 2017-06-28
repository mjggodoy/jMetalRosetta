//  gGA.h
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
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

#ifndef __GGASEQ_H__
#define __GGASEQ_H__

#include <jmetalcpp/core/Algorithm.hh>
#include <jmetalcpp/core/Problem.hh>
#include <jmetalcpp/core/SolutionSet.hh>
#include <jmetalcpp/util/comparators/ObjectiveComparator.hh>

class gGASeq : public Algorithm {
//private:
  //int populationSize_; Maria 12/5/2017 These variables are never used
  //int maxEvaluations_;
public:
  gGASeq(Problem * problem);
  SolutionSet * execute();
};

#endif /* GGA_H_ */
