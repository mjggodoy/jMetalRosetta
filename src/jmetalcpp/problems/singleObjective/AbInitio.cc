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


// Maria: Dependences from C++ language

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

 AbInitio::AbInitio(string solutionType, ProtocolOP ab,  Pose & fold_pose, int numberOfVariables, 
	string const strategy, int rma_stage_sample) {
    stage = rma_stage_sample;
    pose = fold_pose;
    rosetta_abinitio = ab;
	numberOfVariables_   = numberOfVariables;
	numberOfObjectives_  = 1; // monoobjective problem so the number of objectives is just one.
	numberOfConstraints_ = 0;
	problemName_= "AbInitio";

    if(stage==1){

        stage_iterations=200;
    }

    if(stage==2){

        stage_iterations=10;

    }

    if(stage==3){

        stage_iterations=10;

    }

    if(stage==4){

        stage_iterations=40;

    }

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

     
    
    //MAX_EVALUATIONS_STAGE1 = JMETAL_ITERATIONS_STAGE1 * population_size; //500*100=50000
    //MAX_EVALUATIONS_STAGE3 = JMETAL_ITERATIONS_STAGE3 * population_size; // 150000
    //MAX_EVALUATIONS_STAGE4 = JMETAL_ITERATIONS_STAGE4 * population_size; // 200000

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


//cout << stage_iterations << endl;

if( temp_strategy != "FT" && temp_strategy != "VT"){
    
        std::cout << "\n\tERROR: Unrecognised format for temp_strategy " << temp_strategy << std::endl << std::endl;
		exit(1);
    
    }else{
        
        variable_temp = ( temp_strategy == "VT" ) ? true : false;
    }
  
    Variable **variables = solution->getDecisionVariables();
    //std::cout << "fold_pose's size: " << pose.size() << std::endl;
    //std::cout << "number of variables: " << solution->getNumberOfVariables() << std::endl;
    //std::cout << "ESTOY AQUI 0" << std::endl;

    if((int)pose.size()*3 != (int)solution->getNumberOfVariables()){
 
        cout << "The number of variables does not equal to the size of the problem: " << "\n" << "Number of angles per aminoacids in protein: " << pose.size()*3 
        << " Number of variables in jMetal solution: " << solution->getNumberOfVariables() <<endl;
        exit(-1);
    }

    //std::cout << "ESTOY AQUI 1" << std::endl;
    //Maria 31-5-17: Getting the Rosetta's solution.

    for ( int pos = 0; pos < (int)pose.size(); pos++ ) {

        //std::cout << "pos = " << pos << std::endl;
        //Sstd::cout << 

        pose.set_phi(pos+1, variables[pos*3]->getValue());
        pose.set_psi(pos+1, variables[pos*3+1]->getValue());
        pose.set_omega(pos+1, variables[pos*3+2]->getValue());

    }

    pose.fraglength = 9;
        //std::cout << "ESTOY AQUI 2" << std::endl;

    if(stage==1){

        //stage_iterations=200;

        std::cout << "Maria:  Evaluation in Stage 1: " << stage << " stage_iterations " << stage_iterations << std::endl;
        
        // Maria: Evaluation of pose (Stage1)
        rosetta_abinitio->mgf_apply_STAGE1(pose, stage_iterations, do_recover, variable_temp ); 
        

    }else if(stage==2){

        //stage_iterations=200;
        
        std::cout << "Maria:  Evaluation in Stage 2: " << stage <<  " stage_iterations " << stage_iterations << std::endl;

         // Maria: Evaluation of pose (Stage2)
        rosetta_abinitio->mgf_apply_STAGE2(pose, stage_iterations, do_recover, variable_temp ); 


    }else if(stage==3){

        //stage_iterations=200;

        std::cout << "Maria:  Evaluation in Stage 3: "  << stage << " stage_iterations " << stage_iterations << std::endl;

        // Maria: Evaluation of pose (Stage3)
        rosetta_abinitio->mgf_apply_STAGE3(pose, stage_iterations, do_recover, variable_temp ); 


    }else if(stage==4){

        //stage_iterations=400;

        std::cout << "Maria:  Evaluation in Stage 4: " << stage << " stage_iterations " << stage_iterations << std::endl;

        // Maria: Evaluation of pose (Stage4)
        rosetta_abinitio->mgf_apply_STAGE4(pose, stage_iterations, do_recover, variable_temp );

    } else {

        std::cout << " Strategy not found. Please provide an strategy: " << stage << std::endl;

        // ERROR
        exit(-1);
    }

     // Maria: Retrieving the energy from pose
     //std::cout << "ESTOY AQUI 3" << std::endl;

    for (int i = 0; i < numberOfObjectives_; i++) {
           
            solution->setObjective(i,pose.energies().total_energy());
    }

    //Maria: get all angles after evaluation and write them in jMetal array of real values.
   // std::cout << "ESTOY AQUI 4" << std::endl;

    for ( int pos = 0; pos < (int)pose.size(); pos++ ) {

        variables[pos*3]->setValue(pose.phi(pos+1));
        variables[pos*3+1]->setValue(pose.psi(pos+1));
        variables[pos*3+2]->setValue(pose.omega(pos+1));

    }

    //Maria: evals
   // evals++;
    
    //Maria: Recaluating 
    //revaluate(evals, solution);
     //Maria: Configuring evaluations 
    //configureEvaluation();
}

 



//core::pose::PoseOP AbInitio::createPose(std::string const& sequence){

    //core::pose::PoseOP pose = core::pose::PoseOP( new core::pose::Pose );
    //core::pose::make_pose_from_sequence(*pose, sequence, *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )));

    //for ( Size pos = 1; pos <= pose->size(); pos++ ) {
		
		//if ( ! pose->residue(pos).is_protein() ) continue;
		//pose->set_phi( pos, -150 );
		//pose->set_psi( pos, 150);
		//pose->set_omega( pos, 180 );
	//}

    //return pose;
//}







