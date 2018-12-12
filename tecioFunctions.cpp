#include "tecioFunctions.h"

#include <TECIO.h>

#include <fstream>
#include <iostream>

#include <cassert>

double *holder;

bool openBinStream(const std::string& fileName,
                     const std::vector<std::string>& vars,
                        size_t NJ, size_t NI) {

  INTEGER4 Debug = 0;
  INTEGER4 VIsDouble = 1;
  INTEGER4 FileFormat = 0;
  INTEGER4 I = 0;

  std::string varsStr = "";

  for (int i = 0; i < vars.size(); ++i) {

    if (i != 0) {
      varsStr += " ";
    }

    varsStr += vars[i];
  }

  I = TECINI112((char*)"IJ Ordered Zones",
                (char*)varsStr.c_str(),
                (char*)fileName.c_str(),
                (char*)".",
                &FileFormat,
                &Debug,
                &VIsDouble);

  if (I != 0) {
    return false;
  }

  holder = new double[NJ * NI];
  return true;
}

bool createZone(size_t NI, size_t NJ) {
  
  INTEGER4 I = 0;

  INTEGER4 ZoneType = 0; // Ordered Zone
  
  INTEGER4 IMax = NI;
  INTEGER4 JMax = NJ;
  INTEGER4 KMax = 1;

  INTEGER4 ICellMax = 0;
  INTEGER4 JCellMax = 0;
  INTEGER4 KCellMax = 0;

  double SolutionTime = 0;
  INTEGER4 StrandID = 0;
  INTEGER4 ParentZone = 0;

  INTEGER4 IsBlock = 1;
  INTEGER4 NumFaceConnections = 0;
  INTEGER4 FaceNeighborMode = 0;

  INTEGER4 TotalNumFaceNodes = 0;
  INTEGER4 NumConnectedBoundaryFaces = 0;
  INTEGER4 TotalNumBndryConnections = 0;

  INTEGER4 DIsDouble = 1;

  INTEGER4 TotalNumBndryFaces = 0;
  INTEGER4 ShrConn = 0;
  int ValueLocation[21] = { 1, 1, 0, 0, 0,\
                            0, 0, 0, 0, 0,\
                            0, 0, 0, 0, 0,\
                            0, 0, 0, 0, 0, 0};

  //INTEGER4 ValueLocation = 1;
  I = TECZNE112((char*)"Zone 1",
                  &ZoneType,
                  &IMax,
                  &JMax,
                  &KMax,
                  &ICellMax,
                  &JCellMax,
                  &KCellMax,
                  &SolutionTime,
                  &StrandID,
                  &ParentZone,
                  &IsBlock,
                  &NumFaceConnections,
                  &FaceNeighborMode,
                  &TotalNumFaceNodes,
                  &NumConnectedBoundaryFaces,
                  &TotalNumBndryConnections,
                  NULL,
                  ValueLocation,
                  NULL,
                  &ShrConn);

  if (I != 0) {

    return false;

  }

  return true;

}

bool writeToBinStream(const std::vector<std::vector<double>>& fieldData) {

  assert(!fieldData.empty());

  size_t NI = fieldData.size(),
            NJ = fieldData[0].size();

  INTEGER4 size = NI * NJ;
  INTEGER4 DIsDouble = 1;
  INTEGER4 I = 0;

  for (size_t i = 0; i < NI; i++) {
    for (size_t j = 0; j < NJ; j++) {
      holder[j * NI + i] = fieldData[i][j];
    }
  }

  I = TECDAT112(&size, holder, &DIsDouble);

  if (I != 0) {
    return false;
  }

  return true;
}

void closeBinStream() {

  TECEND112();

  delete holder;
}
