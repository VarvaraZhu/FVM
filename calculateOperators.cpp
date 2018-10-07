#include "calculateOperators.h"
#include <iostream>
#include "functions.h"

void calculateGradient(const matrixScalar& scalarField, matrixVec& gradient,\
    const matrixScalar& cellVolumes, const matrixVec& cellCenter,\
        const matrixVec& iFaceCenter, const matrixVec& iFaceVector,\
            const matrixVec& jFaceCenter, const matrixVec& jFaceVector){

    std::vector<coord> SF(4), RF(4);
    std::vector<coord> Ncell(4);  //индексы соседних ячеек

    if((scalarField.size() != gradient.size()) || (scalarField[0].size() != gradient[0].size()))
        std::cout << "Array sizes does not match" << std::endl;
    else{
        for(size_t i = 1; i != gradient.size() - 1; i++)
            for(size_t j = 1; j != gradient[0].size() - 1; j++){
              std::cout << i <<" "<< j << std::endl;

                SF[0] = -iFaceVector[i - 1][j - 1];
                SF[1] = iFaceVector[i][j - 1];
                SF[2] = -jFaceVector[i - 1][j - 1];
                SF[3] = jFaceVector[i - 1][j];

                RF[0] = iFaceCenter[i - 1][j - 1];
                RF[1] = iFaceCenter[i][j - 1];
                RF[2] = jFaceCenter[i - 1][j - 1];
                RF[3] = jFaceCenter[i - 1][j];

                Ncell[0] = coord(i - 1, j);
                Ncell[1] = coord(i + 1, j);
                Ncell[2] = coord(i, j - 1);
                Ncell[3] = coord(i, j + 1);

                for(size_t k = 0; k != 4; k++ ){
                  std::cout << k<< std::endl;

                  size_t in = Ncell[k].x;
                  size_t jn = Ncell[k].y;
                  //double PF = 0.5 * (scalarField[i][j] + scalarField[in][jn]);
                   double DC = norm2((RF[k] - cellCenter[i][j]));
                   double DN = norm2((RF[k] - cellCenter[in][jn]));

                   double PF = linearInterpolation(DC, DN, scalarField[i][j], scalarField[in][jn]);
                   gradient[i][j] = gradient[i][j] + SF[k] * PF;
                }
                std::cout << cellVolumes[i][j] << std::endl;

                gradient[i][j] = gradient[i][j] * (1 / cellVolumes[i - 1][j - 1]);
              }
        }
    return;
}
