//  AbInitio.hh
//
//  Authors:
//       Maria Jesus Garcia Godoy <mjgarciag@lcc.uma.es>
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

#ifndef __ABINITIO__
#define __ABINITIO__

#include <jmetalcpp/core/Problem.hh>
#include <math.h>
#include <jmetalcpp/encodings/solutionType/RealSolutionType.hh>
#include <jmetalcpp/core/Solution.hh>


class AbInitio : public Problem {
   
public:
	/// @brief Constructor 1 of the Abinito problem
	AbInitio(string solutionType, int numberOfVariables = 10);
    
    /// @brief Destructor
	~AbInitio();
	void evaluate(Solution *solution);
};

private: 

    void configureEvaluation();

private:

	int evals; 						//Number of Evaluations
	int rma_stage_sample;           // Number of Rosetta's stage
	int RMA_ITERATIONS;				// Iterations at each stage of RMA
    int STAGE1_ITERATIONS;			// Iterations for Rosetta stage1 
	int STAGE2_ITERATIONS;			// Iterations for Rosetta stage2
	int STAGE3_ITERATIONS;			// Iterations for Rosetta stage3 
	int STAGE4_ITERATIONS;			// Iterations for Rosetta stage4
    int MAX_EVALUATIONS_STAGE1;		// Max allowed evaluations in stage1		
	int MAX_EVALUATIONS_STAGE2;		// Max allowed evaluations in stage2		
	int MAX_EVALUATIONS_STAGE3;		// Max allowed evaluations in stage3		
	int MAX_EVALUATIONS_STAGE4;		// Max allowed evaluations in stage4


#endif

/**
  * @class AbInitio
  * @brief Class representing problem Sphere
 **/