//  Real.cpp
//
//  Author:
//       Juan J. Durillo <durillo@lcc.uma.es>
//       Esteban López-Camacho <esteban@lcc.uma.es>
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


#include "RealJM.hh"


/** 
  * Empty constructor. 
  * It will only initialize all the variables.
 **/
RealJM::RealJM() {
  value_ = 0.0;
} // Real


/**
 * Lower and Upper bounds based constructor.
 * It will initialize the upper and lower bounds of the 
 * variable and it will initialize its value to a random
 * value between those upper and lower bound values.
 **/
RealJM::RealJM(double lowerBound, double upperBound) {
  value_ = PseudoRandom::randDouble(lowerBound, upperBound);
  lowerBound_ = lowerBound;
  upperBound_ = upperBound;
} // Real


/**
 * Constructor
 */
RealJM::RealJM(Variable * variable) {
  lowerBound_ = variable->getLowerBound();
  upperBound_ = variable->getUpperBound();
  value_      = variable->getValue();
} // Real


/**
 * Destructor
 */
RealJM::~RealJM() { /* do nothing */ }


/**
 * Gets the value of the <code>Real</code> variable.
 * @return the value.
 */
double RealJM::getValue() {
  return value_;
} // getValue


/**
 * Sets the value of the variable.
 * @param value The value.
 */
void RealJM::setValue(double value) {
  value_ = value;
} // setValue


/**
 * Returns a exact copy of the <code>Real</code> variable
 * @return the copy
 */
Variable *RealJM::deepCopy() {
  return new RealJM(this);
} // deepCopy


/**
 * Gets the lower bound of the variable.
 * @return the lower bound.
 */
double RealJM::getLowerBound() {
  return lowerBound_;
} // getLowerBound


/**
 * Gets the upper bound of the variable.
 * @return the upper bound.
 */
double RealJM::getUpperBound() {
  return upperBound_;
} // getUpperBound


/**
 * Sets the lower bound of the variable.
 * @param lowerBound The lower bound.
 */
void RealJM::setLowerBound(double bound) {
  lowerBound_ = bound;
} // setLowerBound


/**
 * Sets the upper bound of the variable.
 * @param upperBound The upper bound.
 */
void RealJM::setUpperBound(double bound) {
  upperBound_ = bound;
} // setUpperBound


/**
 * Returns a string representing the object
 * @return The string
 */
string RealJM::toString(){
  std::ostringstream stringStream;
  stringStream << value_ ;
  string aux = stringStream.str() + " ";

  return aux ;
} // toString
