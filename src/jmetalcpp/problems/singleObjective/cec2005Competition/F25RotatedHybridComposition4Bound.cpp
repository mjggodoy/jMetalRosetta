//  F25RotatedHybridComposition4Bound.cpp
//
//  Authors:
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


#include <F25RotatedHybridComposition4Bound.h>

// Fixed (class) parameters
const string F25RotatedHybridComposition4Bound::FUNCTION_NAME = "Rotated Hybrid Composition Function 4 without bounds";
// TODO: Cambiar ruta
const string F25RotatedHybridComposition4Bound::DEFAULT_FILE_DATA = "../../data/cec2005CompetitionResources/supportData/hybrid_func4_data.txt";
const string F25RotatedHybridComposition4Bound::DEFAULT_FILE_MX_PREFIX = "../../data/cec2005CompetitionResources/supportData/hybrid_func4_M_D";
const string F25RotatedHybridComposition4Bound::DEFAULT_FILE_MX_SUFFIX = ".txt";

const double F25RotatedHybridComposition4Bound::m_sigma[NUM_FUNC] = {
  2.0,  2.0,  2.0,  2.0,  2.0,
  2.0,  2.0,  2.0,  2.0,  2.0
};

const double F25RotatedHybridComposition4Bound::m_lambda[NUM_FUNC] = {
    10.0,       5.0/20.0,   1.0,   5.0/32.0,   1.0,
    5.0/100.0,  5.0/50.0,   1.0,   5.0/100.0,  5.0/100.0
};

const double F25RotatedHybridComposition4Bound::m_func_biases[NUM_FUNC] = {
  0.0,    100.0,  200.0,  300.0,  400.0,
  500.0,  600.0,  700.0,  800.0,  900.0
};


/**
 * Constructor.
 */
F25RotatedHybridComposition4Bound::F25RotatedHybridComposition4Bound(int dimension, double bias)
    : F25RotatedHybridComposition4Bound(dimension, bias, DEFAULT_FILE_DATA, getFileMxName(DEFAULT_FILE_MX_PREFIX, dimension, DEFAULT_FILE_MX_SUFFIX)) {
} // F25RotatedHybridComposition4Bound


/**
 * Constructor
 */
F25RotatedHybridComposition4Bound::F25RotatedHybridComposition4Bound(int dimension, double bias, string file_data, string file_m)
    : TestFunc(dimension, bias, FUNCTION_NAME) {

  // Note: dimension starts from 0
  m_o = new double*[NUM_FUNC];
  m_M = new double**[NUM_FUNC];

  m_testPoint = new double[m_dimension];
  m_testPointM = new double[m_dimension];
  m_fmax = new double[NUM_FUNC];

  m_w  = new double[NUM_FUNC];
  m_z  = new double*[NUM_FUNC];
  m_zM = new double*[NUM_FUNC];

  for (int i=0; i<NUM_FUNC; i++) {
    m_o[i]  = new double[m_dimension];
    m_z[i]  = new double[m_dimension];
    m_zM[i] = new double[m_dimension];
    m_M[i] = new double*[m_dimension];
    for (int j=0; j<m_dimension; j++) {
      m_M[i][j] = new double[m_dimension];
    }
  }

  // Load the shifted global optimum
  Benchmark::loadMatrixFromFile(file_data, NUM_FUNC, m_dimension, m_o);
  // Load the matrix
  Benchmark::loadNMatrixFromFile(file_m, NUM_FUNC, m_dimension, m_dimension, m_M);

  // Initialize the hybrid composition job object
  theJob = new MyHCJob(NUM_FUNC);
  theJob->num_func = NUM_FUNC;
  theJob->num_dim = m_dimension;
  theJob->C = 2000.0;
  for (int i=0; i<NUM_FUNC; i++) {
    theJob->sigma[i]  = m_sigma[i];
    theJob->biases[i] = m_func_biases[i];
    theJob->lambda[i] = m_lambda[i];
  }
  theJob->o = m_o;
  theJob->M = m_M;
  theJob->w = m_w;
  theJob->z = m_z;
  theJob->zM = m_zM;
  // Calculate/estimate the fmax for all the functions involved
  for (int i=0; i<NUM_FUNC; i++) {
    for (int j=0; j<m_dimension; j++) {
      m_testPoint[j] = (5.0 / m_lambda[i]);
    }
    Benchmark::rotate(m_testPointM, m_testPoint, m_M[i], m_dimension);
    m_fmax[i] = fabs(theJob->basic_func(i, m_testPointM, m_dimension));
  }
  theJob->fmax = m_fmax;

} // F25RotatedHybridComposition4Bound


/**
 * Destructor
 */
F25RotatedHybridComposition4Bound::~F25RotatedHybridComposition4Bound() {

  for (int i=0; i<NUM_FUNC; i++) {
    delete [] m_o[i];
    delete [] m_z[i];
    delete [] m_zM[i];
  }
  for (int i=0; i<NUM_FUNC; i++) {
    for (int j=0; j<m_dimension; j++) {
      delete [] m_M[i][j];
    }
    delete [] m_M[i];
  }

  delete [] m_o;
  delete [] m_M;

  delete [] m_testPoint;
  delete [] m_testPointM;
  delete [] m_fmax;

  delete [] m_w;
  delete [] m_z;
  delete [] m_zM;

  delete theJob;

} // ~F25RotatedHybridComposition4Bound

F25RotatedHybridComposition4Bound::MyHCJob::MyHCJob(int numFunc)
    : HCJob(numFunc) { }

double F25RotatedHybridComposition4Bound::MyHCJob::basic_func(int func_no, double* x, int length) {
  double result = 0.0;
  switch(func_no) {
    case 0:
      result = Benchmark::weierstrass(x,length);
      break;
    case 1:
      result = Benchmark::EScafferF6(x,length);
      break;
    case 2:
      result = Benchmark::F8F2(x,length);
      break;
    case 3:
      result = Benchmark::ackley(x,length);
      break;
    case 4:
      result = Benchmark::rastrigin(x,length);
      break;
    case 5:
      result = Benchmark::griewank(x,length);
      break;
    case 6:
      result = Benchmark::EScafferF6NonCont(x,length);
      break;
    case 7:
      result = Benchmark::rastriginNonCont(x,length);
      break;
    case 8:
      result = Benchmark::elliptic(x,length);
      break;
    case 9:
      result = Benchmark::sphere_noise(x,length);
      break;
    default:
      cerr << "func_no is out of range." << endl;
      exit(-1);
  }
  return (result);
}


/**
 * Function body
 */
double F25RotatedHybridComposition4Bound::f(double * x) {
  double result = 0.0;

  result = Benchmark::hybrid_composition(x, theJob, m_dimension);

  result += m_bias;

  return result;
}


string F25RotatedHybridComposition4Bound::getFileMxName(string prefix, int dimension, string suffix) {
  std::stringstream sstm;
  sstm << prefix << dimension << suffix;
  return sstm.str();
}
