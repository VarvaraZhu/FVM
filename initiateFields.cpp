#include "initiateFields.h"
#include <iostream>
#include "functions.h"

void InitiatePressureField(matrixScalar &P, const matrixVec& cellCenter){
	if((P.size() != cellCenter.size()) || (P[0].size() != cellCenter[0].size()))
		std::cout << "Array sizes does not match" << std::endl;
  	else{
    	for (size_t i = 0; i != P.size(); i++){
      		for(size_t j = 0; j != P[0].size(); j++)
        		P[i][j] = Pressure(cellCenter[i][j]);
    		}
  	}
  	return;
}
void InitiateVelocityField(matrixVec &V, const matrixVec& mesh){ return;}