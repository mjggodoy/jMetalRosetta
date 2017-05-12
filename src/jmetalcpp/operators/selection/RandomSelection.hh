//  RandomSelection.h
//
//  Author:
//       Cristian Zambrano V. <cristian_uteq@hotmail.com>
//
//  Copyright (c) 2013 Antonio J. Nebro, Juan J. Durillo
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

#ifndef __RANDOM_SELECTION__
#define __RANDOM_SELECTION__

#include <jmetalcpp/core/SolutionSet.hh>
#include <jmetalcpp/operators/selection/Selection.hh>
#include <jmetalcpp/util/comparators/Comparator.hh>
#include <jmetalcpp/util/comparators/DominanceComparator.hh>
#include <jmetalcpp/util/PseudoRandom.hh>

/**
 * This class implements an binary tournament selection operator
 */
class RandomSelection : public Selection {

//private:
 // Comparator * comparator_; //Maria 12/5/2017: This private variable is never used.

public:
  RandomSelection(map<string, void *> parameters);
  ~RandomSelection();
  void *execute(void *);

};

#endif
