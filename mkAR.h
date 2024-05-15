#ifndef MKAR_H
#define MKAR_H
#include "dataContainers.h"
#include <vector>

void mkAR(dataList, std::vector<int>, ARmodel &, dataList &);
void orthonormalizeR(dataList, dataList &, matrix &);
void rotate_model(ARmodel &, matrix L);


#endif
