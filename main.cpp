#include<iostream>
#include<fstream>

#include <vector>
#include <string>

#include "functions.h"
#include "calculateGeom.h"
#include "initiateFields.h"
#include "coord.h"
#include "calculateOperators.h"

int main(int argc, char* argv[]){

//  Mesh file's name is read from cmd and contained in argv[1]

   /* if (argc > 1)
        std::cout << "Mesh will be read from " << argv[1] << std::endl;

    else {
        std::cout << "Please, set the mesh file" << std::endl;
        return 0;
    }*/

    std::string outputFileName = "output.plt";
    std::string inputFileName = "base.msh";

    std::fstream inputMesh(inputFileName); // name of file with computational mesh

    size_t NI, NJ;

  //********************** Read Mesh File **************************************
  matrixVec mesh;

  if (inputMesh.is_open()){

    std::cout << "Reading nodes number..." << std::endl;

    inputMesh >> NI >> NJ;
    std::cout << "nI = " << NI << " nJ = " << NJ << std::endl;

    std::cout << "Reading Mesh file..." << std::endl;
    //Vector array containing mesh nodes
    mesh.resize(NI);

    for (size_t i = 0; i != NI; i++){
      for(size_t j = 0; j != NJ; j++) {
		      mesh[i].resize(NJ);
             inputMesh >> mesh[i][j].x >> mesh[i][j].y;
           }
    }
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

  matrixScalar cellVolumes(NI - 1, std::vector<double>(NJ - 1));  //Cell Volumes
  matrixVec cellCenter(NI + 1, std::vector<coord>(NJ + 1)); //Cell Centers
  matrixVec iFaceCenter(NI, std::vector<coord>(NJ - 1));  //Face Centers for I-faces
  matrixVec iFaceVector(NI, std::vector<coord>(NJ - 1));  //Face Vectors for I-faces
  matrixVec jFaceCenter(NI - 1, std::vector<coord>(NJ));  //Face Centers for J-faces
  matrixVec jFaceVector(NI - 1, std::vector<coord>(NJ));  //Face Vectors for J-faces

//********************** Calculate Metric **************************************

  std::cout << "Calculate metric..." << std::endl;

  calculateGeom(mesh, cellVolumes, cellCenter, iFaceCenter, iFaceVector, jFaceCenter, jFaceVector);

//********************** Initiate Fields ***************************************

std::cout << "Initiate fields..." << std::endl;

InitiatePressureField(P, cellCenter);

//********************** Calculate Gradient ************************************
std::cout << "Calculate gradient..." << std::endl;

calculateGradient(P, gradP, cellVolumes, cellCenter, iFaceCenter, iFaceVector, jFaceCenter, jFaceVector);

//********************** Output Fields *****************************************
std::cout << "Output fields in file: " << outputFileName << std::endl;

   std::ofstream outputData;
   outputData.open(outputFileName);
   if (outputData.is_open())
   {
     outputData << "VARIABLES =  \"X\",\"Y\",\"P\",\"gradPx\",\"gradPy\"" << std::endl;
     outputData << "ZONE T=\"ZONE 001\"" << std::endl;

     outputData << "I=" << NI << " J=" << NJ <<", DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)"<< std::endl;

     //outputData.setf(ios::std::scientific);

     for (int j = 0; j != NJ; j++){
       for (int i = 0; i != NI; i++) {
         outputData << mesh[i][j].x << " ";
       }
     }
     outputData << std::endl;

     for (int j = 0; j != NJ; j++){
       for (int i = 0; i != NI; i++) {
         outputData << mesh[i][j].y << " ";
       }
     }
     outputData << std::endl;

     for (int j = 1; j != NJ; j++){
       for (int i = 1; i != NI; i++) {
         outputData << P[i][j] << " ";
       }
     }

     for (int j = 1; j != NJ; j++){
       for (int i = 1; i != NI; i++) {
         outputData << gradP[i][j].x << " ";
       }
     }

     for (int j = 1; j != NJ; j++){
       for (int i = 1; i != NI; i++) {
         outputData << gradP[i][j].y << " ";
       }
     }

     outputData << std::endl;

   }
   else {
     std::cout << "Cannot open output file" << std::endl;
     return 0;
   }

       return 0;
}
