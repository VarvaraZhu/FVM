#include <iostream>

#include "functions.h"
#include "calculateOperators.h"

#include <string>
//Calculate gradien using Gauss method

void calculateGradientGauss(const matrixScalar& Pressure, matrixVec& gradient, \
                        const matrixScalar& cellVolumes, const matrixVec& cellCenter,\
                          const matrixVec& iFaceCenter, const matrixVec& iFaceVector,\
                            const matrixVec& jFaceCenter, const matrixVec& jFaceVector){

  std::vector<coord> SF(4), RF(4);    //RF - Center of current Face
                                       //SF - normal to current Face
  std::vector<coord> Ncell(4);          //indices of neigbouhood cells

 
    for(size_t i = 1; i != gradient.size() - 1; i++)
      for(size_t j = 1; j != gradient[0].size() - 1; j++){

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

        for (size_t k = 0; k != 4; k++) {

          size_t in = Ncell[k].x;
          size_t jn = Ncell[k].y;

          //PF - Pressure value on current Face

          double DC = norm2((RF[k] - cellCenter[i][j]));
          double DN = norm2((RF[k] - cellCenter[in][jn]));

          double PF = linearInterpolation(DC, DN, Pressure[i][j], Pressure[in][jn]);

          gradient[i][j] = gradient[i][j] + SF[k] * PF;
        }

        gradient[i][j] = gradient[i][j] * (1 / cellVolumes[i - 1][j - 1]);
      }
  
  return;
}

//Calculate gradien using Gauss method with iterations

void calculateGradientIterGauss(const matrixScalar& Pressure, matrixVec& gradient, \
                                  const matrixScalar& cellVolumes, const matrixVec& cellCenter, \
                                    const matrixVec& iFaceCenter, const matrixVec& iFaceVector, \
                                      const matrixVec& jFaceCenter, const matrixVec& jFaceVector) {

  std::vector<coord> SF(4), RF(4);    //RF - Center of current Face
                                       //SF - normal to current Face
  std::vector<coord> Ncell(4);          //indices of neigbouhood cells

    for (size_t i = 1; i != gradient.size() - 1; i++)
      for (size_t j = 1; j != gradient[0].size() - 1; j++) {

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

        coord GP = coord(0, 0);

        for (size_t Face = 0; Face != 4; Face++) {

          size_t in = Ncell[Face].x;
          size_t jn = Ncell[Face].y;

          double DC = norm2((RF[Face] - cellCenter[i][j]));
          double DN = norm2((RF[Face] - cellCenter[in][jn]));

          double PFM = linearInterpolation(DC, DN, Pressure[i][j], Pressure[in][jn]);
          coord GradM = linearInterpolation(DC, DN, gradient[i][j], gradient[in][jn]);
          coord RM = linearInterpolation(DC, DN, cellCenter[i][j], cellCenter[in][jn]);

          double PF = PFM + (RF[Face] - RM) * GradM;

          GP = GP + SF[Face] * PF;
        }

        gradient[i][j] = GP * (1 / cellVolumes[i - 1][j - 1]);

      }
 
  return;
}

//Calculate gradien using Ordinary Least Squared method

void calculateGradientOLS(const matrixScalar& Pressure, matrixVec& gradient, \
                            const matrixScalar& cellVolumes, const matrixVec& cellCenter, \
                              const matrixVec& iFaceCenter, const matrixVec& iFaceVector, \
                                const matrixVec& jFaceCenter, const matrixVec& jFaceVector) {

    std::vector<coord> Ncell(4);          //indices of neigbouhood cells
    std::vector<double> H(2);
    std::vector<std::vector<double>> G(2, std::vector<double>(2));

    for (size_t i = 1; i != gradient.size() - 1; i++)
      for (size_t j = 1; j != gradient[0].size() - 1; j++) {

        Ncell[0] = coord(i - 1, j);
        Ncell[1] = coord(i + 1, j);
        Ncell[2] = coord(i, j - 1);
        Ncell[3] = coord(i, j + 1);

        H[0] = 0;
        H[1] = 0;

        G[0][0] = 0;
        G[0][1] = 0;
        G[1][0] = 0;
        G[1][1] = 0;

        for (size_t k = 0; k != 4; k++) {
          
          size_t iN = Ncell[k].x;
          size_t jN = Ncell[k].y;

          coord PN = cellCenter[iN][jN] - cellCenter[i][j];

          double PN2 = PN * PN;

          H[0] += (Pressure[iN][jN] - Pressure[i][j]) * PN.x / PN2;
          H[1] += (Pressure[iN][jN] - Pressure[i][j]) * PN.y / PN2;

          G[0][0] += PN.x * PN.x / PN2;
          G[0][1] += PN.x * PN.y / PN2;
          G[1][0] += PN.y * PN.x / PN2;
          G[1][1] += PN.y * PN.y / PN2;
        }

        double det = G[0][0] * G[1][1] - G[1][0] * G[0][1];

        gradient[i][j].x = (H[0] * G[1][1] - H[1] * G[0][1]) / det;
        gradient[i][j].y = (H[1] * G[0][0] - H[0] * G[1][0]) / det;

      }
  
  return;
}


void calculateDivergence(const matrixVec& Velocity, matrixScalar& divergence, \
                           const matrixScalar& cellVolumes, const matrixVec& cellCenter, \
                             const matrixVec& iFaceCenter, const matrixVec& iFaceVector, \
                                const matrixVec& jFaceCenter, const matrixVec& jFaceVector) {

  std::vector<coord> SF(4), RF(4);    //RF - Center of current Face
                                        //SF - normal to current Face
  std::vector<coord> Ncell(4);          //indices of neigbouhood cells


    for (size_t i = 1; i != divergence.size() - 1; i++)
      for (size_t j = 1; j != divergence[0].size() - 1; j++) {

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

        for (size_t k = 0; k != 4; k++) {

          size_t in = Ncell[k].x;
          size_t jn = Ncell[k].y;

          //VF - Velocity on current Face

          double DC = norm2((RF[k] - cellCenter[i][j]));
          double DN = norm2((RF[k] - cellCenter[in][jn]));

          coord VF = linearInterpolation(DC, DN, Velocity[i][j], Velocity[in][jn]);

          divergence[i][j] = divergence[i][j] + SF[k] * VF;
        }

        divergence[i][j] = divergence[i][j] * (1 / cellVolumes[i - 1][j - 1]);
      }
  
  return;
}

void calculatePVDivergence(const std::string& scheme, const matrixScalar& Pressure, const matrixVec& Velocity, matrixScalar& divergence, \
                            const matrixScalar& cellVolumes, const matrixVec& cellCenter, \
                              const matrixVec& iFaceCenter, const matrixVec& iFaceVector, \
                                const matrixVec& jFaceCenter, const matrixVec& jFaceVector, \
                                  const matrixVec& gradP) {

  std::vector<coord> SF(4), RF(4);    //RF - Center of current Face
                                        //SF - normal to current Face
  std::vector<coord> Ncell(4);          //indices of neigbouhood cells

    for (size_t i = 1; i != divergence.size() - 1; i++)
      for (size_t j = 1; j != divergence[0].size() - 1; j++) {

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

        for (size_t k = 0; k != 4; k++) {

          size_t in = Ncell[k].x;
          size_t jn = Ncell[k].y;

          //VF - Velocity on current Face

          double DC = norm2((RF[k] - cellCenter[i][j]));
          double DN = norm2((RF[k] - cellCenter[in][jn]));

          coord VF = linearInterpolation(DC, DN, Velocity[i][j], Velocity[in][jn]);
          double PF;

          if (scheme == "Central")
            PF = linearInterpolation(DC, DN, Pressure[i][j], Pressure[in][jn]);

          else if (scheme == "Upwind") {
            double GF = VF * SF[k];
            if (GF > 0)
              PF = Pressure[i][j];
            else {
              if (DN < 1e-7)
                PF = 2 * Pressure[in][jn] - Pressure[i][j];
              else
                PF = Pressure[in][jn];
            }
          }

          else if (scheme == "SOU") {
            double GF = VF * SF[k];
            if (GF > 0) {
              coord r = RF[k] - cellCenter[i][j];
              PF = Pressure[i][j] + gradP[i][j] * r;
            }
            else {
              if (DN < 1e-7){
                double PO = 2 * Pressure[in][jn] - Pressure[i][j];
                double G1 = gradP[i][j] * (cellCenter[i][j] - RF[k]);
                double GF = Pressure[i][j] - Pressure[in][jn];
                PF = PO + (4 * GF - 3 * G1);
              }
                
              else {
                coord r = RF[k] - cellCenter[in][jn];
                PF = Pressure[in][jn] + gradP[in][jn] * r;
              }
            }
          }

          divergence[i][j] = divergence[i][j] + SF[k] * (VF * PF);
        }

        divergence[i][j] = divergence[i][j] * (1 / cellVolumes[i - 1][j - 1]);
      }
  
  return;
}


void calculateLaplacian(const matrixScalar& Pressure, matrixScalar& laplacian,
 const matrixVec& gradP, \
                          const matrixScalar& cellVolumes,
 const matrixVec& cellCenter, \
                            const matrixVec& iFaceCenter,
 const matrixVec& iFaceVector, \
                              const matrixVec& jFaceCenter, const matrixVec& jFaceVector, bool &correction, size_t order) {

  std::vector<coord> SF(4), RF(4);    //RF - Center of current Face
                                        //SF - normal to current Face
  std::vector<coord> Ncell(4), Ocell(4);          //indices of neigbouhood cells


  for (size_t i = 1; i != laplacian.size() - 1; i++)
    for (size_t j = 1; j != laplacian[0].size() - 1; j++) {

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

      Ocell[0] = coord(i + 1, j);
      Ocell[1] = coord(i - 1, j);
      Ocell[2] = coord(i, j + 1);
      Ocell[3] = coord(i, j - 1);

      laplacian[i][j] = 0;

      for (size_t k = 0; k != 4; k++) {

        size_t in = Ncell[k].x;
        size_t jn = Ncell[k].y;

        size_t io = Ocell[k].x;
        size_t jo = Ocell[k].y;

        double DNC = norm2(cellCenter[in][jn] - cellCenter[i][j]);
        coord RNC = (cellCenter[in][jn] - cellCenter[i][j]) * (1 / DNC);
        double DC = norm2((RF[k] - cellCenter[i][j]));
        double DN = norm2((RF[k] - cellCenter[in][jn]));

        coord NF = SF[k] * (1 / norm2(SF[k]));
        double dpdn;
        coord GF;
     
        if (DN < 1e-7) {
          //Первый порядок точности аппроксимации производной на границах
          if (order == 1)
            dpdn = (Pressure[in][jn] - Pressure[i][j]) / DC;
          // аппроксимация со вторым порядком точности
          //dpdn = (9 * (Pressure[in][jn] - Pressure[i][j]) - (Pressure[in][jn] - Pressure[io][jo])) / (6 * DNC); 
          else if (order == 2) {
            double dpdn_center = gradP[i][j] * NF;
            dpdn = 5. / 3. * (Pressure[in][jn] - Pressure[i][j]) / DNC - 2. / 3. * dpdn_center;
          }

          GF = gradP[i][j];
        }
        else {
          dpdn = (Pressure[in][jn] - Pressure[i][j]) / DNC;
          GF = linearInterpolation(DC, DN, gradP[i][j], gradP[in][jn]);
        }

        if (correction == 1)
          dpdn = dpdn + (NF - RNC) * GF;

        laplacian[i][j] = laplacian[i][j] + dpdn * norm2(SF[k]);
      }

      laplacian[i][j] = laplacian[i][j] * (1 / cellVolumes[i - 1][j - 1]);
    }

  return;
}