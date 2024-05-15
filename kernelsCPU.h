#ifndef KERNELSCPU_H
#define KERNELSCPU_H
#include <vector>
#include <complex>
#include "utility.h"
void generateRotationMatrices(std::vector<float>, std::vector<float> &, int,int);
void transposeBlockMatrices(std::vector<float> &, std::vector<float>, int,int,int);
void compTransferFunc(std::vector<float>, std::vector<std::complex<float>> &, int, paramContainer);
void scale_columns(std::vector<std::complex<float>> &, int, int, int);
void shrinkArrays(std::vector<std::complex<float>>, std::vector<std::complex<float>> &,int, int, int);
void prodEigs(std::vector<float>, std::vector<float> &,int, int,int);
void det2GC(std::vector<float>, std::vector<float>, std::vector<float> &, int,int);
#endif
