#pragma once
#include <vector>

bool openBinStream(const std::string& fileName,
                      const std::vector<std::string>& vars,
                          size_t NI, size_t NJ);

bool createZone(size_t NI, size_t NJ);
bool writeToBinStream(const std::vector<std::vector<double>>& fieldData);

void closeBinStream();