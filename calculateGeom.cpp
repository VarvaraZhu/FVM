#include "calculateGeom.h"
#include<iostream>

void calculateGeom(const matrixVec& mesh, \
  matrixScalar& cellVolumes, matrixVec& cellCenter, matrixVec& iFaceCenter, \
  matrixVec& iFaceVector, matrixVec& jFaceCenter, matrixVec& jFaceVector){

    coord r;

//****************** Face centers and face vectors *****************************
//i - direction
//NI - 1, NJ -2 -> NI - 1, NJ - 1------>все ок
  std::cout << "Face centers and face vectors : i" << std::endl;
     for(size_t i = 0; i != iFaceVector.size(); i++){
        for(size_t j = 0; j != iFaceVector[0].size(); j++){
            r = mesh[i][j + 1] - mesh[i][j];  // r - vector from one node to another
            iFaceVector[i][j].x = r.y;        //iFaceVector = r rotated on 90 degree
            iFaceVector[i][j].y = - r.x;      //iFaceVector directer to increasing i

            iFaceCenter[i][j] = (mesh[i][j] + mesh[i][j + 1]) * 0.5;
         }
     }

//j - direction
//NI - 2, NJ - 1 -> NI - 1, NJ - 1------>все ок
std::cout << "Face centers and face vectors : j" << std::endl;

    for(size_t i = 0; i != jFaceVector.size(); i++){
        for(size_t j = 0; j != jFaceVector[0].size(); j++){
            r = mesh[i + 1][j] - mesh[i][j];  // r - vector from one node to another
            jFaceVector[i][j].x = - r.y;      //jFaceVector = r rotated on -90 degree
            jFaceVector[i][j].y = r.x;        //jFaceVector directer to increasing j

            jFaceCenter[i][j] = (mesh[i][j] + mesh[i + 1][j]) * 0.5;
        }
    }

//********************** Cell Volumes  *****************************************
//NI - 2, NJ -2 -> NI - 1, NJ - 1------>все ок
std::cout << "Cell Volumes" << std::endl;

    for(size_t i = 0; i != cellVolumes.size(); i++){
        for(size_t j = 0; j != cellVolumes[0].size(); j++){
          r = mesh[i + 1][j + 1] - mesh[i][j];
          cellVolumes[i][j] = 0.5 * (iFaceVector[i][j] * r) + 0.5 * (jFaceVector[i][j] * r);
        }
    }

//********************** Cell Centers  *****************************************
//For inner cells: center of contour
std::cout << "Cell Centers for inner cells " << std::endl;

//NI - 2, NJ - 2 -> NI - 1, NJ - 1------>все ок

    for(size_t i = 1; i != cellCenter.size() - 1; i ++){
      for(size_t j = 1; j != cellCenter[0].size() - 1; j ++){

         cellCenter[i][j] = (iFaceCenter[i - 1][j - 1] * norm2(iFaceVector[i - 1][j - 1]) +\
                          iFaceCenter[i][j - 1] * norm2(iFaceVector[i][j - 1]) +\
                            jFaceCenter[i - 1][j - 1] * norm2(jFaceVector[i - 1][j - 1]) +\
                              jFaceCenter[i - 1][j] * norm2(jFaceVector[i - 1][j])) * \
                              (1.0 /(norm2(iFaceVector[i - 1][j - 1]) + norm2(iFaceVector[i][j - 1])+\
                                + norm2(jFaceVector[i - 1][j - 1]) + norm2(jFaceVector[i - 1][j])));
      }
    }

//For dummy cells on boundaries: cell center = face centers
//i - boundaries
std::cout << "Cell Centers for boundaries: i " << std::endl;

//NI - 2, NJ - 2 -> NI - 1, NJ - 1------>все ок
    for(size_t j = 1; j != cellCenter[0].size() - 1; j++){
        cellCenter[0][j] = iFaceCenter[0][j - 1];
        cellCenter[cellCenter.size() - 1][j] = iFaceCenter[iFaceCenter.size() - 1][j - 1];
    }

//j - boundaries
std::cout << "Cell Centers for boundaries: j " << std::endl;

    for(size_t i = 1; i != cellCenter.size() - 1; i++){
        cellCenter[i][0] = jFaceCenter[i - 1][0];
        cellCenter[i][cellCenter[0].size() - 1] = jFaceCenter[i - 1][jFaceCenter[0].size() - 1];
    }

    return;
  }
