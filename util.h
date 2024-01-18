#ifndef UTIL_H
#define UTIL_H

#include <stdint.h>
#include "routines.h"
#include <eigen3/Eigen/Dense>
#include <vector>

using std::vector;
using std::array;

void linspace(double *x, double start, double end, uint64_t len);
Eigen::MatrixXd mat_single(double *x, uint64_t k);
Eigen::MatrixXd mat_double(double **xp, uint64_t k);
uint64_t bootstrap_single(vector<double> *x_axis, vector<double> *y_axis, bool (*model)(double, uint64_t));
uint64_t bootstrap_double(vector<double> *x_axis, vector<double> *y_axis, bool (*model)(double, double, uint64_t));

#endif // !UTIL_H

