//  AbInitio.cc
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


// Maria: Dependences from Rosetta's code

#include "AbInitio.hh"
#include <jmetalcpp/core/Solution.hh>
#include <protocols/abinitio/Protocol.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/Sequence.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>


// Maria: Dependences from C++ code

#include <string>
#include <utility/io/util.hh>

using core::Size;
using namespace core;



/**
 * Constructor.
 * Creates a new instance of the AbInitio problem.
 * @param solutionType The solution type must "Real", "BinaryReal, and "ArrayReal".
 * @param numberOfVariables Number of variables of the problem
 */

 AbInitio::AbInitio(string solutionType, ProtocolOP ab, std::string const& sequence, int numberOfVariables) {
    rosetta_abinitio = ab;
	numberOfVariables_   = numberOfVariables;
	numberOfObjectives_  = 1; // monoobjective problem so the number of objectives is just one.
	numberOfConstraints_ = 0;
	problemName_= "AbInitio";
	
    lowerLimit_ = new double[numberOfVariables_];
	if (lowerLimit_ == NULL) {
		cout << "Abinitio::Abinitio. Error reserving memory for storing the array of lower limits" << endl;
		exit(-1) ;
	}	
	
	upperLimit_ = new double[numberOfVariables_];
	if (upperLimit_ == NULL) {
		cout << "Abinitio::Abinitio. Error reserving memory for storing the array of upper limits" << endl;
		exit(-1) ;
	}

    //Maria: Range of values for the problem
    
    for (int i=0; i<numberOfVariables_; i++) {
        
        lowerLimit_[i] = -180;
        upperLimit_[i] = 180;    
    }
	
    // TODO: Solution type initialization
    solutionType_ = new RealSolutionType(this);	
    
    //Inizialite number of evaluations:
    rma_stage_sample = 1;
    RMA_ITERATIONS = 100;
    STAGE1_ITERATIONS = 200;
    STAGE2_ITERATIONS = 200;
    STAGE3_ITERATIONS = 200;
    STAGE4_ITERATIONS = 400;
    MAX_EVALUATIONS_STAGE1 = STAGE1_ITERATIONS * 500 * RMA_ITERATIONS;
    MAX_EVALUATIONS_STAGE2 = MAX_EVALUATIONS_STAGE1 + STAGE2_ITERATIONS * 500 * RMA_ITERATIONS;
    MAX_EVALUATIONS_STAGE3 = MAX_EVALUATIONS_STAGE2 + STAGE3_ITERATIONS * 500 * RMA_ITERATIONS;
    MAX_EVALUATIONS_STAGE4 = MAX_EVALUATIONS_STAGE3 + STAGE4_ITERATIONS * 500 * RMA_ITERATIONS;

} // AbInitio


/**
 * Destructor
 */
AbInitio::~AbInitio() {
  delete [] lowerLimit_ ;
  delete [] upperLimit_ ;
  delete solutionType_ ;
} // ~AbInitio


void AbInitio::evaluate(Solution *solution) {

    iterations = STAGE1_ITERATIONS * RMA_ITERATIONS;
    variable_temp = ( temp_strategy == "VT" ) ? true : false;


    //Maria 30-5-17: Create pose from sequence and initialize the angles in -150, 150 and 180.
    core::pose::PoseOP pose = createPose(sequence);
    
    Variable **variables = solution->getDecisionVariables();
    
    if((int)pose->size()*3 != (int)solution->getNumberOfVariables()){
 
        cout << "The number of variables do not equal to the size of the problem: " << "\n" << "Number of angles per aminoacids in protein" << pose->size()*3 
        << " Number of variables in jMetal solution: " << solution->getNumberOfVariables() <<endl;
        exit(-1);
    }


    //Maria 31-5-17: Getting the Rosetta's solution.

    for ( int pos = 0; pos <= numberOfVariables_; pos++ ) {

        pose->set_phi(pos, variables[pos*3]->getValue());
        pose->set_psi(pos, variables[pos*3+1]->getValue());
        pose->set_omega(pos, variables[pos*3+2]->getValue());

    }


    if(rma_stage_sample==1){

        std::cout << "Maria:  Evaluation in Stage 1" << evals << std::endl;

        // Maria: Evaluation of pose (Stage1)
        rosetta_abinitio->mgf_apply_STAGE1(*pose, iterations, do_recover, variable_temp );
        

    }else if(rma_stage_sample==2){

        std::cout << "Maria:  Evaluation in Stage 2" << evals << std::endl;

         // Maria: Evaluation of pose (Stage2)
        rosetta_abinitio->mgf_apply_STAGE2(*pose, iterations, do_recover, variable_temp );


    }else if(rma_stage_sample==3){

        
         std::cout << "Maria:  Evaluation in Stage 3" << evals << std::endl;

        // Maria: Evaluation of pose (Stage3)
        rosetta_abinitio->mgf_apply_STAGE3(*pose, iterations, do_recover, variable_temp );


    }else if(rma_stage_sample==4){

        std::cout << "Maria:  Evaluation in Stage 4" << evals << std::endl;

        // Maria: Evaluation of pose (Stage4)
        rosetta_abinitio->mgf_apply_STAGE4(*pose, iterations, do_recover, variable_temp );

    } else {
        // ERROR
        exit(-1);
    }


     // Maria: Retrieving the energy from pose

    for (int i = 0; i < numberOfObjectives_; i++) {
           
            solution->setObjective(i,pose->energies().total_energy());
    }

    //Maria: get all angles after evaluation and write them in jMetal array of real values.

    for ( int pos = 0; pos <= numberOfVariables_; pos++ ) {

        variables[pos*3]->setValue(pose->phi(pos*3));
        variables[pos*3+1]->setValue(pose->psi(pos*3+1));
        variables[pos*3+2]->setValue(pose->omega(pos*3+2));

    }

    evals++;
    configureEvaluation();

}


void AbInitio::configureEvaluation(){

    if(evals<=MAX_EVALUATIONS_STAGE1) {

        rma_stage_sample=1;

    }else if(evals <= MAX_EVALUATIONS_STAGE2 ){

        rma_stage_sample=2;

    }else if(evals<= MAX_EVALUATIONS_STAGE3){

        rma_stage_sample=3;

    }else if(evals<= MAX_EVALUATIONS_STAGE4){

        rma_stage_sample=4;

    }
}

core::pose::PoseOP AbInitio::createPose(std::string const& sequence){

    core::pose::PoseOP pose = core::pose::PoseOP( new core::pose::Pose );
    core::pose::make_pose_from_sequence(*pose, sequence, *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )));

    for ( Size pos = 1; pos <= pose->size(); pos++ ) {
		
		if ( ! pose->residue(pos).is_protein() ) continue;
		pose->set_phi( pos, -150 );
		pose->set_psi( pos, 150);
		pose->set_omega( pos, 180 );
	}

    return pose;
}







