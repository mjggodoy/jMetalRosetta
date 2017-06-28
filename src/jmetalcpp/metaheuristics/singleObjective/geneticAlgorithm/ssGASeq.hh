//  ssGA.h
//
//  Author:
//       María Jesús García Godoy <mjgarciag@lcc.uma.es>
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

#ifndef __SSGASEQ_H__
#define __SSGASEQ_H__

#include <jmetalcpp/core/Algorithm.hh>
#include <jmetalcpp/core/Problem.hh>
#include <jmetalcpp/core/SolutionSet.hh>
#include <jmetalcpp/util/comparators/ObjectiveComparator.hh>
#include <jmetalcpp/operators/selection/WorstSolutionSelection.hh>

/**
 * Class implementing a steady-state genetic algorithm
 */
class ssGASeq : public Algorithm {

public:
  ssGASeq(Problem * problem);
  SolutionSet * execute();

};

#endif /* __SSGA_H__ */
