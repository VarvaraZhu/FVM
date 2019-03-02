#include "initiateFields.h"
#include <iostream>
#include "functions.h"

//Initialization of Pressure Field
void InitiatePressureField(matrixScalar &P, const matrixVec& cellCenter, const bool& linear){
	if((P.size() != cellCenter.size()) || (P[0].size() != cellCenter[0].size()))
		std::cout << "Array sizes does not match" << std::endl;
  	else{
    	for (size_t i = 0; i != P.size(); i++){
      		for(size_t j = 0; j != P[0].size(); j++)
        		P[i][j] = Pressure(cellCenter[i][j], linear);
    		}
  	}
  	return;
}

//Initialization of Velocity Field
void InitiateVelocityField(matrixVec &V, const matrixVec& cellCenter, const bool& linear) {
  if ((V.size() != cellCenter.size()) || (V[0].size() != cellCenter[0].size()))
    std::cout << "Array sizes does not match" << std::endl;
  else {
    for (size_t i = 0; i != V.size(); i++) {
      for (size_t j = 0; j != V[0].size(); j++)
        V[i][j] = Velocity(cellCenter[i][j], linear);
    }
  }
  return;
}

//Initialization of Pressure x Velocity Field
void InitiatePVField(matrixVec &PV, const matrixScalar &P, const matrixVec& V) {
  for (size_t i = 0; i != P.size(); i++) {
    for (size_t j = 0; j != P[0].size(); j++)
      PV[i][j] = V[i][j] * P[i][j];
  }
  return;
}

//Initialization of Exact Pressure Gradient Field
void CalculateExactGradientField(matrixVec &gradPExact, const matrixVec& cellCenter, const bool& linear) {
  if ((gradPExact.size() != cellCenter.size()) || (gradPExact[0].size() != cellCenter[0].size()))
    std::cout << "Array sizes does not match" << std::endl;
  else {
    for (size_t i = 0; i != gradPExact.size(); i++) {
      for (size_t j = 0; j != gradPExact[0].size(); j++)
        gradPExact[i][j] = gradientExact(cellCenter[i][j], linear);
    }
  }
  return;
}

//Initialization of Exact Pressure Laplacian Field
void CalculateExactLaplacianField(matrixScalar &lapPExact, const matrixVec& cellCenter, const bool& linear){
  if ((lapPExact.size() != cellCenter.size()) || (lapPExact[0].size() != cellCenter[0].size()))
    std::cout << "Array sizes does not match" << std::endl;
  else {
    for (size_t i = 0; i != lapPExact.size(); i++) {
      for (size_t j = 0; j != lapPExact[0].size(); j++)
        lapPExact[i][j] = laplacianExact(cellCenter[i][j], linear);
    }
  }
  return;
}

//Initialization of Velocity Divergence Field
void CalculateExactDivergenceField(matrixScalar &divVExact, const matrixVec& cellCenter, const bool& linear) {
  if ((divVExact.size() != cellCenter.size()) || (divVExact[0].size() != cellCenter[0].size()))
    std::cout << "Array sizes does not match" << std::endl;
  else {
    for (size_t i = 0; i != divVExact.size(); i++) {
      for (size_t j = 0; j != divVExact[0].size(); j++)
        divVExact[i][j] = divergenceExact(cellCenter[i][j], linear);
    }
  }
  return;
}

//Initialization of Pressure x Velocity Divergence Field
void CalculateExactPVDivergenceField(matrixScalar &divPVExact, const matrixVec& cellCenter, const bool& linear) {

  if ((divPVExact.size() != cellCenter.size()) || (divPVExact[0].size() != cellCenter[0].size()))
    std::cout << "Array sizes does not match" << std::endl;
  else {
    for (size_t i = 0; i != divPVExact.size(); i++) {
      for (size_t j = 0; j != divPVExact[0].size(); j++)
        divPVExact[i][j] = divergencePVExact(cellCenter[i][j], linear);
    }
  }
  return;
}

//Calculation of Scalar Error
void CalculateErrorField(const matrixVec &gradP, const matrixVec &gradPExact, matrixScalar &Error) {
  for (size_t i = 0; i != gradPExact.size(); i++) {
    for (size_t j = 0; j != gradPExact[0].size(); j++)
      Error[i][j] = error(gradP[i][j], gradPExact[i][j]);
  }
  return;
}

//Calculation of vector Error
void CalculateErrorField(const matrixScalar &A, const matrixScalar &AExact, matrixScalar &Error) {
  for (size_t i = 0; i != AExact.size(); i++) {
    for (size_t j = 0; j != AExact[0].size(); j++)
      Error[i][j] = error(A[i][j], AExact[i][j]);
  }
  return;
}