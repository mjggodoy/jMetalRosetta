//  DE.cpp
//
//  Author:
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


#include "DESeq.hh"


/*
 * This class implements a differential evolution algorithm.
 */


/**
 * Constructor
 * @param problem Problem to solve
 */
DESeq::DESeq(Problem *problem) : Algorithm(problem) {
} // DE


/*
 * Runs of the DE algorithm.
 * @return a <code>SolutionSet</code> that is a set of non dominated solutions
 * as a result of the algorithm execution
 */
SolutionSet * DESeq::execute() {

  int populationSize;
  int maxEvaluations;
  int evaluations;

  SolutionSet * population;
  SolutionSet * offspringPopulation;

  Operator * selectionOperator;
  Operator * crossoverOperator;


  Comparator * comparator;
  comparator = new ObjectiveComparator(0); // Single objective comparator

  // Differential evolution parameters


  Solution ** parent;

  //Read the parameters
  populationSize = *(int *) getInputParameter("populationSize");
  maxEvaluations = *(int *) getInputParameter("maxEvaluations");
  // Maria: Use as an input the population:

  //if (getInputParameter("population") == NULL) cout << "EXITO" << endl;
  //cout << getInputParameter("population") << endl;

  population = (SolutionSet *)getInputParameter("population");

  //if (population==NULL) cout << "FRACASO (EXITO?)" << endl;

  //cout << population << endl;

  selectionOperator = operators_["selection"];
  crossoverOperator = operators_["crossover"];

  //Initialize the variables
  
  evaluations = 0;

  if (population) {

    cout << "Testing if the population is getting correctly" << endl;
    cout << "Population's size: " << population->size() << endl;
    cout << population->get(0) << endl;
    cout << "Maria: DESeq: Population does exist. We evaluate it:" << endl;
    // Evaluate the input population
    for (int i = 0; i < populationSize; i++) {
      Solution * solution = population->get(i);
      problem_->evaluate(solution);
      problem_->evaluateConstraints(solution);
      evaluations++;
    } //for
    
  } else {

    population  = new SolutionSet(populationSize);
    cout << "Maria: DESeq: Population does not exist. We evaluate it:" << endl;
    // Create the initial solutionSet
    //cout << "AQUI" << endl;
    Solution * newSolution;
    for (int i = 0; i < populationSize; i++) {
      //cout << i << endl;
      newSolution = new Solution(problem_);
      //cout << "AQUI "<< endl;
      problem_->evaluate(newSolution);
      //cout << "AQUI "<< endl;
      problem_->evaluateConstraints(newSolution);
      //cout << "AQUI "<< endl;
      evaluations++;
      //cout << "AQUI "<< endl;
      population->add(newSolution);
      //cout << "AQUI "<< endl;

    } //for
    
  }

  // Create the initial solutionSet
  //Solution * newSolution;
  

  // Generations ...
  population->sort(comparator) ;
  while (evaluations < maxEvaluations) {

   // Create the offSpring solutionSet
   offspringPopulation  = new SolutionSet(populationSize);

   //offspringPopulation.add(new Solution(population.get(0))) ;

   for (int i = 0; i < populationSize; i++) {
     // Obtain parents. Two parameters are required: the population and the
     //                 index of the current individual
     void ** object1 = new void*[2];
     object1[0] = population;
     object1[1] = &i;
     parent = (Solution **) (selectionOperator->execute(object1));
     delete[] object1;

     Solution * child;

     // Crossover. Two parameters are required: the current individual and the
     //            array of parents
     void ** object2 = new void*[2];
     object2[0] = population->get(i);
     object2[1] = parent;
     child = (Solution *) (crossoverOperator->execute(object2));
     delete[] object2;
     delete[] parent;

     problem_->evaluate(child);

     evaluations++;
     cout << "EVALUATIONS = " << evaluations << endl;

     if (comparator->compare(population->get(i), child) < 0) {
       offspringPopulation->add(new Solution(population->get(i)));
       delete child;
     } else {
       offspringPopulation->add(child);
     }
   } // for

   // The offspring population becomes the new current population
   for (int i = 0; i < populationSize; i++) {
     delete population->get(i);
   }
   population->clear();
   for (int i = 0; i < populationSize; i++) {
     population->add(offspringPopulation->get(i));
   }
   offspringPopulation->clear();
   delete offspringPopulation;
   population->sort(comparator);

   cout << "DE: Best objective from this iteration -> " << population->get(0)->getObjective(0) << 
    " (EVALUATIONS = " << evaluations << ")" << endl;
  } // while

  delete comparator;

  // Return a population with the best individual
  //SolutionSet * resultPopulation = new SolutionSet(1);
  //resultPopulation->add(new Solution(population->get(0)));
  //delete population;
  
  return population;

} // execute
