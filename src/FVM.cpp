#include<iostream>
#include<fstream>

#include <vector>
#include <string>

#include "coord.h"
#include "functions.h"
#include "initiateFields.h"
#include "calculateGeom.h"
#include "calculateOperators.h"

#include "outputToTecplot.h"

int main(void) {

  std::fstream inputData("inputData.dat");

  std::string meshFileName; //Name of file, containing mesh
  bool linear;              //Using linear or non-linear functions for initialisation of P and V fields
  std::string gradMethod;   //Gradient calculation methods: Gauss, iterGauss, OLS
  std::string scheme;   //Scheme for calculating of Pressure on a Face: Central, Upwind, SOU
  std::string outputFileName; //Name of output file
  bool correction;
  size_t nIter;   //Number of iterations for iteration Gauss method
  size_t order;   //Order of derivatives' approximation on boarders in laplacian calculation

  //********************** Read Input Data **************************************

  if (inputData.is_open()) {

    inputData >> meshFileName;
    inputData >> linear;
    if ((linear != 0) && (linear != 1)) {
      std::cout << "Invalid type of initialization functions";
      return 0;
    }

    inputData >> gradMethod;
    if ((gradMethod != "Gauss") && (gradMethod != "iterGauss") && (gradMethod != "OLS")) {
      std::cout << "Invalid name of gradient calculation method";
      return 0;
    }

    inputData >> nIter;

    inputData >> scheme;
    if ((scheme != "Central") && (scheme != "Upwind") && (scheme != "SOU")) {
      std::cout << "Invalid name of scalar interpolation scheme";
      return 0;
    }
    inputData >> order;
    if ((order != 1) && (order != 2) ) {
      std::cout << "Invalid order of derivatives approximation";
      return 0;
    }
    inputData >> correction; 
    inputData >> outputFileName;

    inputData.close();
  }
  else {
    std::cout << "Invalid input file" << std::endl;
    return 0;
  }
  //********************** Read Mesh File **************************************
  std::fstream inputMesh(meshFileName);

  size_t NI, NJ; //Number of nodes in I and J directions

  matrixVec mesh;   //Vector array containing mesh nodes 

  std::cout << "Mesh will be read from " << meshFileName << std::endl;

  if (inputMesh.is_open()) {

    std::cout << "Reading nodes number..." << std::endl;

    inputMesh >> NI >> NJ;
    std::cout << "nI = " << NI << " nJ = " << NJ << std::endl;

    std::cout << "Reading nodes coordinates..." << std::endl;
   
    mesh.resize(NI);

    for (size_t i = 0; i != NI; i++)
      mesh[i].resize(NJ);

    for (size_t j = 0; j != NJ; j++)
      for (size_t i = 0; i != NI; i++)
        inputMesh >> mesh[i][j].x >> mesh[i][j].y;

    inputMesh.close();
  }
  else {
    std::cout << "Invalid Mesh file";
    return 0;
  }

  //********************** Allocate All Arrays  **********************************

  std::cout << "Allocate arrays..." << std::endl;

  matrixScalar P(NI + 1, std::vector<double>(NJ + 1));  //Pressure
  matrixVec gradP(NI + 1, std::vector<coord>(NJ + 1));  //Pressure gradient
  matrixVec gradPExact(NI + 1, std::vector<coord>(NJ + 1));  //Exact pressure gradient
  matrixScalar gradError(NI + 1, std::vector<double>(NJ + 1));  //Pressure gradient error

  matrixScalar lapP(NI + 1, std::vector<double>(NJ + 1));  //Pressure laplasian
  matrixScalar lapPExact(NI + 1, std::vector<double>(NJ + 1));  //Exact pressure laplasian
  matrixScalar lapPError(NI + 1, std::vector<double>(NJ + 1));  //Pressure laplasian error

  matrixVec V(NI + 1, std::vector<coord>(NJ + 1));      //Velocity
  matrixScalar divV(NI + 1, std::vector<double>(NJ + 1));  //Velocity divergence
  matrixScalar divVExact(NI + 1, std::vector<double>(NJ + 1));  //Exact velocity divergence
  matrixScalar divVError(NI + 1, std::vector<double>(NJ + 1));  //Velocity divergence error

  matrixVec PV(NI + 1, std::vector<coord>(NJ + 1));      //Velocity x Pressure
  matrixScalar divPV(NI + 1, std::vector<double>(NJ + 1));  //Velocity x Pressure divergence
  matrixScalar divPVExact(NI + 1, std::vector<double>(NJ + 1));  //Exact Velocity x Pressure divergence
  matrixScalar divPVError(NI + 1, std::vector<double>(NJ + 1));  //Velocity x Pressure divergence error

  matrixScalar cellVolumes(NI - 1, std::vector<double>(NJ - 1));  //Cell Volumes
  matrixVec cellCenter(NI + 1, std::vector<coord>(NJ + 1)); //Cell Centers
  matrixVec iFaceCenter(NI, std::vector<coord>(NJ - 1));  //Face Centers for I-faces
  matrixVec iFaceVector(NI, std::vector<coord>(NJ - 1));  //Face Vectors for I-faces
  matrixVec jFaceCenter(NI - 1, std::vector<coord>(NJ));  //Face Centers for J-faces
  matrixVec jFaceVector(NI - 1, std::vector<coord>(NJ));  //Face Vectors for J-faces

//********************** Calculate Metric **************************************

  std::cout << "Calculate metric..." << std::endl;

  calculateMetric(mesh, cellVolumes, cellCenter, iFaceCenter, iFaceVector, jFaceCenter, jFaceVector);

  //********************** Initiate Fields ***************************************

  std::cout << "Initiate fields..." << std::endl;

  InitiatePressureField(P, cellCenter, linear);
  InitiateVelocityField(V, cellCenter, linear);
  InitiatePVField(PV, P, V);

  //********************** Calculate Pressure Gradient ******************************

  std::cout << "Calculate gradient..." << std::endl;
  std::cout << gradMethod << " method" << std::endl;
  CalculateExactGradientField(gradPExact, cellCenter, linear);

  if (gradMethod == "Gauss") {
    calculateGradientGauss(P, gradP, cellVolumes, cellCenter, iFaceCenter, iFaceVector, jFaceCenter, jFaceVector);
    CalculateErrorField(gradP, gradPExact, gradError);
    std::cout << matrixMax(gradError) << std::endl;
  }

  else if (gradMethod == "iterGauss") {
    std::cout << "Iter" << " maxErr" << std::endl;

    for (size_t k = 0; k != nIter; k++) {
      calculateGradientIterGauss(P, gradP, cellVolumes, cellCenter, iFaceCenter, iFaceVector, jFaceCenter, jFaceVector);
      CalculateErrorField(gradP, gradPExact, gradError);
      double m = matrixMax(gradError);

      std::cout << k << " " << matrixMax(gradError) << std::endl;
    }
  }

  else if (gradMethod == "OLS") {
    calculateGradientOLS(P, gradP, cellVolumes, cellCenter, iFaceCenter, iFaceVector, jFaceCenter, jFaceVector);
    CalculateErrorField(gradP, gradPExact, gradError);
    std::cout << matrixMax(gradError) << std::endl;
  }

  //********************** Calculate Velocity Divergence ****************************

  std::cout << "Calculate divergence..." << std::endl;
  CalculateExactDivergenceField(divVExact, cellCenter, linear);
  calculateDivergence(V, divV, cellVolumes, cellCenter, iFaceCenter, iFaceVector, jFaceCenter, jFaceVector);
  CalculateErrorField(divV, divVExact, divVError);

  std::cout << matrixMax(divVError) << std::endl;

  //********************** Calculate Velocity x Pressure  Divergence ******************

  std::cout << "Calculate divergence of pV..." << std::endl;
  std::cout << "Interpolation of the scalar on faces will be performed according to " << scheme << " scheme" << std::endl;
  CalculateExactPVDivergenceField(divPVExact, cellCenter, linear);

  calculatePVDivergence(scheme, P, V, divPV, cellVolumes, \
      cellCenter, iFaceCenter, iFaceVector, jFaceCenter, jFaceVector, gradP);

  CalculateErrorField(divPV, divPVExact, divPVError);
  
  std::cout << matrixMax(divPVError) << std::endl;

  //********************** Calculate Velocity Divergence ****************************

  std::cout << "Calculate laplacian..." << std::endl;

  std::cout << "Approximation of the derivatives on boarders will be with " << order << " order" << std::endl;

  if (correction) {
    std::cout << "Skew correction will be used" << std::endl;

  }

  CalculateExactLaplacianField(lapPExact, cellCenter, linear);
  calculateLaplacian(P, lapP, gradP, cellVolumes, cellCenter, iFaceCenter, iFaceVector, jFaceCenter, jFaceVector, correction, order);
  CalculateErrorField(lapP, lapPExact, lapPError);
  std::cout << matrixMax(lapPError) << std::endl;

  //********************** Output Fields *****************************************
  std::cout << "Output fields in file: " << outputFileName << std::endl;

  ouputToTecplot(outputFileName, mesh, P, V, PV, \
    gradP, divV, divPV, lapP, \
    gradPExact, divVExact, divPVExact, lapPExact, \
    gradError, divVError, divPVError, lapPError);

  return 0;
}
