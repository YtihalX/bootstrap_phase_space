#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <stdint.h>
#include <stdio.h>
bool anharmonics(double E, double xsq, uint64_t size) {
  uint64_t xlen = std::max(2 * size - 1, 8ul);
  double **xp = (double **)malloc(xlen * sizeof(double *));
  for (uint64_t i = 0; i < xlen; i++)
    xp[i] = (double *)malloc((2 * size - 1) * sizeof(double));
  xp[0][0] = 1;
  xp[1][0] = 0;
  xp[2][0] = xsq;
  xp[3][0] = 0;
  for (uint64_t i = 4; i < xlen; i++) {
    xp[i][0] = ((i - 2) * xp[i - 2][0] + (i - 3) * xp[i - 4][0] * E) / (i - 1);
  }
  for (uint64_t i = 0; i < 2 * size - 1; i++) {
    xp[i][1] = 0;
  }
  for (uint64_t j = 2; j < 2 * size - 1; j++) {
    for (uint64_t i = 0; i < 4; i++) {
      xp[i][j] = E * xp[i][j - 2] + xp[i + 2][j - 2] - xp[i + 4][j - 2];
    }
    for (uint64_t i = 4; i < 2 * size - 1; i++) {
      xp[i][j] = ((i + j - 2) * xp[i - 2][j] + E * (i - 3) * xp[i - 4][j]) /
		 (i + 2 * j - 1);
    }
  }
  //    for (uint64_t i = 0; i < xlen; i++) {
  //        for (uint64_t j = 0; j < 2*size - 1; j++) {
  //            printf("%f ", xp[i][j]);
  //        }
  //        puts("");
  //    }
  Eigen::MatrixXd mat(size * size, size * size);
  for (uint64_t i = 0; i < size * size; i++) {
    for (uint64_t j = i; j < size * size; j++) {
      uint64_t a = i / size;
      uint64_t b = i % size;
      uint64_t c = j / size;
      uint64_t d = j % size;
      mat(i, j) = xp[a + c][b + d];
      mat(j, i) = mat(i, j);
    }
  }
  // std::cout<<mat<<'\n';
  Eigen::EigenSolver<Eigen::MatrixXd> solver(mat);
  return solver.eigenvalues().real().minCoeff() >= 0;
}

bool toda(double E, double ex, uint64_t size) {
  uint64_t xlen = std::max(2 * size - 1, 4ul);
  double **xp = (double **)malloc(xlen * sizeof(double *));
  for (uint64_t i = 0; i < xlen; i++)
    xp[i] = (double *)malloc((2 * size - 1) * sizeof(double));
  xp[0][0] = 1;
  xp[1][0] = ex;
  for (int64_t i = 2; i < xlen; i++) {
    xp[i][0] = (2 * E * (i - 1) * xp[i - 1][0] + (-2 * i + 3) * xp[i - 2][0]) /
	       (2 * i - 1);
    printf("what? %f\n", xp[i][0]);
  }
  for (int64_t i = 0; i < xlen; i++)
    xp[i][1] = 0;
  for (int64_t j = 2; j < 2 * size - 1; j++) {
    for (int64_t i = 1; i < 3; i++)
      xp[i][j] = E * xp[i][j - 2] - xp[i + 1][j - 2] - xp[i - 1][j - 2];
    xp[0][j] = ((3 + j) * xp[2][j] - 2 * E * xp[1][j]) / (j - 1);
    for (int64_t i = 3; i < xlen; i++)
      xp[i][j] =
	  (2 * E * (i - 1) * xp[i - 1][j] + (j - 2 * i + 3) * xp[i - 2][j]) /
	  (2 * i + j - 1);
  }
  for (int64_t i = 0; i < xlen; i++) {
    for (int64_t j = 0; j < 2 * size - 1; j++) {
      printf("%.4f ", xp[i][j]);
    }
    puts("");
  }
  Eigen::MatrixXd mat(size * size, size * size);
  for (uint64_t i = 0; i < size * size; i++) {
    for (uint64_t j = i; j < size * size; j++) {
      uint64_t a = i / size;
      uint64_t b = i % size;
      uint64_t c = j / size;
      uint64_t d = j % size;
      mat(i, j) = xp[a + c][b + d];
      mat(j, i) = mat(i, j);
    }
  }
  std::cout << mat << '\n';
  for (uint64_t i = 0; i < xlen; i++)
    free(xp[i]);
  free(xp);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
  return solver.eigenvalues().minCoeff() >= 0;
}

int main() { std::cout << toda(0, 0.3, 4) << '\n'; }
