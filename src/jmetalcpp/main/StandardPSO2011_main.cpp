//  StandardPSO2011_main.cpp
//
//  Author:
//       Esteban López-Camacho <esteban@lcc.uma.es>
//       Antonio J. Nebro <antonio@lcc.uma.es>
//
//  Copyright (c) 2014 Antonio J. Nebro
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

#include <Problem.h>
#include <Algorithm.h>
#include <Solution.h>
#include <iostream>
#include <StandardPSO2011.h>
#include <time.h>
#include <PolynomialMutation.h>
#include <ProblemFactory.h>

int main(int argc, char ** argv) {

  clock_t t_ini, t_fin;

  Problem   * problem   ; // The problem to solve
  Algorithm * algorithm ; // The algorithm to use
  //Operator  * mutation  ; // "Turbulence" operator

  if (argc>=2) {
    problem = ProblemFactory::getProblem(argc, argv);
    cout << "Selected problem: " << problem->getName() << endl;
  } else {
    cout << "No problem selected." << endl;
    cout << "Default problem will be used: Sphere" << endl;
    problem = ProblemFactory::getProblem(const_cast<char *>("Sphere"));
  }

  algorithm = new StandardPSO2011(problem);

  // Algorithm parameters
  int swarmSize = 100;
  int maxIterations = 8000;
  int numberOfParticlesToInform = 3;
  algorithm->setInputParameter("swarmSize",&swarmSize);
  algorithm->setInputParameter("maxIterations",&maxIterations);
  algorithm->setInputParameter("numberOfParticlesToInform", &numberOfParticlesToInform);

  // Execute the Algorithm
  t_ini = clock();
  SolutionSet * population = algorithm->execute();
  t_fin = clock();
  double secs = (double) (t_fin - t_ini);
  secs = secs / CLOCKS_PER_SEC;

  // Result messages
  cout << "Total execution time: " << secs << "s" << endl;
  cout << "Variables values have been written to file VAR" << endl;
  population->printVariablesToFile("VAR");
  cout << "Objectives values have been written to file FUN" << endl;
  population->printObjectivesToFile("FUN");

  delete population;
  delete algorithm;

} // main
