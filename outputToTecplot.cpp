#include "tecioFunctions.h"

#include "outputToTecplot.h"

void ouputToTecplot(std::string outputFileName, const matrixVec& mesh, matrixScalar& P, const matrixVec& V, const matrixVec& PV, \
  const matrixVec& gradP, const matrixScalar& divV, const matrixScalar& divPV, const matrixScalar& lapP, \
  const matrixVec& gradPExact, const matrixScalar& divVExact, const matrixScalar& divPVExact, const matrixScalar& lapPExact, \
  const matrixScalar& gradError, const matrixScalar& divVError, const matrixScalar& divPVError, const matrixScalar& lapPError) {



  std::vector<std::string> variables = { "X", "Y", "P", "Vx", "Vy", "PVx", "PVy",\
                                        "gradPx","gradPy","divV", "divPV", "lapP",\
                                          "exactGradPx","exactGradPy", "exactDivV", "exactDivPV", "exactLapP",\
                                            "GradError", "DivError" , "DivPVError", "lapPError" }; // List of variables' names

  size_t NI = mesh.size(), NJ = mesh[0].size();

  openBinStream(outputFileName, variables, NI, NJ);
  createZone(NI, NJ);

  matrixScalar outX(NI, std::vector<double>(NJ)), outY(NI, std::vector<double>(NJ)), outP(NI - 1, std::vector<double>(NJ - 1)), outVx(NI - 1, std::vector<double>(NJ - 1)), outVy(NI - 1, std::vector<double>(NJ - 1)), outPVx(NI - 1, std::vector<double>(NJ - 1)), outPVy(NI - 1, std::vector<double>(NJ - 1)), \
    outgradPx(NI - 1, std::vector<double>(NJ - 1)), outgradPy(NI - 1, std::vector<double>(NJ - 1)), outdivV(NI - 1, std::vector<double>(NJ - 1)), outdivPV(NI - 1, std::vector<double>(NJ - 1)), outlapP(NI - 1, std::vector<double>(NJ - 1)), \
    outexactGradPx(NI - 1, std::vector<double>(NJ - 1)), outexactGradPy(NI - 1, std::vector<double>(NJ - 1)), outexactDivV(NI - 1, std::vector<double>(NJ - 1)), outexactDivPV(NI - 1, std::vector<double>(NJ - 1)), outexactLapP(NI - 1, std::vector<double>(NJ - 1)), \
    outGradError(NI - 1, std::vector<double>(NJ - 1)), outDivError(NI - 1, std::vector<double>(NJ - 1)), outDivPVError(NI - 1, std::vector<double>(NJ - 1)), outlapPError(NI - 1, std::vector<double>(NJ - 1));


  for (int i = 0; i != NI; i++)
    for (int j = 0; j != NJ; j++) {
      outX[i][j] = mesh[i][j].x;
      outY[i][j] = mesh[i][j].y;
    }

  for (int i = 1; i != NI; i++)
    for (int j = 1; j != NJ; j++) {
      outP[i - 1][j - 1] = P[i][j];
      outVx[i - 1][j - 1] = V[i][j].x;
      outVy[i - 1][j - 1] = V[i][j].y;
      outPVx[i - 1][j - 1] = PV[i][j].x;
      outPVy[i - 1][j - 1] = PV[i][j].y;

      outgradPx[i - 1][j - 1] = gradP[i][j].x;
      outgradPy[i - 1][j - 1] = gradP[i][j].y;
      outdivV[i - 1][j - 1] = divV[i][j];
      outdivPV[i - 1][j - 1] = divPV[i][j];
      outlapP[i - 1][j - 1] = lapP[i][j];

      outexactGradPx[i - 1][j - 1] = gradPExact[i][j].x;
      outexactGradPy[i - 1][j - 1] = gradPExact[i][j].y;
      outexactDivV[i - 1][j - 1] = divVExact[i][j];
      outexactDivPV[i - 1][j - 1] = divPVExact[i][j];
      outexactLapP[i - 1][j - 1] = lapPExact[i][j];

      outGradError[i - 1][j - 1] = gradError[i][j];
      outDivError[i - 1][j - 1] = divVError[i][j];
      outDivPVError[i - 1][j - 1] = divPVError[i][j];
      outlapPError[i - 1][j - 1] = lapPError[i][j];
    }

  writeToBinStream(outX);
  writeToBinStream(outY);
  writeToBinStream(outP);
  writeToBinStream(outVx);
  writeToBinStream(outVy);
  writeToBinStream(outPVx);
  writeToBinStream(outPVy);

  writeToBinStream(outgradPx);
  writeToBinStream(outgradPy);
  writeToBinStream(outdivV);
  writeToBinStream(outdivPV);
  writeToBinStream(outlapP);

  writeToBinStream(outexactGradPx);
  writeToBinStream(outexactGradPy);
  writeToBinStream(outexactDivV);
  writeToBinStream(outexactDivPV);
  writeToBinStream(outexactLapP);
  writeToBinStream(outGradError);
  writeToBinStream(outDivError);
  writeToBinStream(outDivPVError);
  writeToBinStream(outlapPError);


  closeBinStream();

  return;
}