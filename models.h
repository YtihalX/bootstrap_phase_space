#ifndef MODELS_H
#define MODELS_H

#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <stdint.h>
using Eigen::SelfAdjointEigenSolver;
using std::max;

bool double_well(double E, double xsq, uint64_t size);
bool double_well_single(double E, double xsq, uint64_t size);
bool double_well_fixed(double x1, double x3, uint64_t size);
bool toda(double E, double ex, uint64_t size);
bool toda_single(double E, double ex, uint64_t size);
bool harmonics(double E, uint64_t size);
bool harmonics_single(double E, uint64_t size);
bool coulomb(double E, double sp, uint64_t size);
bool coulomb_single(double E, uint64_t size);

#endif // !MODELS_H
