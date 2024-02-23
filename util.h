#ifndef UTIL_H
#define UTIL_H

#include "routines.h"
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MPRealSupport>
#include <stdint.h>
#include <vector>

using std::array;
using std::vector;
typedef Eigen::Matrix<mpfr::mpreal, Eigen::Dynamic, Eigen::Dynamic> MatrixXmp;

extern int64_t num_x, num_y, kk;
extern double plot_range[2][2];
void linspace(double *x, double start, double end, uint64_t len);
MatrixXmp mat_single(double *x, uint64_t k);
MatrixXmp mat_double(double **xp, uint64_t k);
uint64_t bootstrap_single(vector<double> *x_axis, vector<double> *y_axis,
			  bool (*model)(double, uint64_t));
uint64_t bootstrap_double(vector<double> *x_axis, vector<double> *y_axis,
			  bool (*model)(double, double, uint64_t));

#endif // !UTIL_H
