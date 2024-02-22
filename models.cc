#include "models.h"
#include "util.h"

bool double_well(double E, double xsq, uint64_t size) {
  uint64_t xlen = max(2 * size - 1, 8ul);
  double **xp = (double **)malloc(xlen * sizeof(double *));
  for (uint64_t i = 0; i < xlen; i++)
    xp[i] = (double *)malloc((2 * size - 1) * sizeof(double));
  xp[0][0] = 1;
  xp[1][0] = 0.;
  xp[2][0] = xsq;
  xp[3][0] = 0.;
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
    for (uint64_t i = 4; i < xlen; i++) {
      xp[i][j] = ((i + j - 2) * xp[i - 2][j] + E * (i - 3) * xp[i - 4][j]) /
		 (i + 2 * j - 1);
    }
  }
  auto mat = mat_double(xp, size);
  for (uint64_t i = 0; i < xlen; i++)
    free(xp[i]);
  free(xp);
  SelfAdjointEigenSolver<MatrixXd> solver(mat);
  return solver.eigenvalues().minCoeff() >= 0;
}

bool double_well_fixed(double x1, double x3, uint64_t size) {
  double E = -0.1;
  uint64_t xlen = max(2 * size - 1, 8ul);
  double **xp = (double **)malloc(xlen * sizeof(double *));
  for (uint64_t i = 0; i < xlen; i++)
    xp[i] = (double *)malloc((2 * size - 1) * sizeof(double));
  xp[0][0] = 1;
  xp[1][0] = x1;
  xp[2][0] = 0.4052;
  xp[3][0] = x3;
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
    for (uint64_t i = 4; i < xlen; i++) {
      xp[i][j] = ((i + j - 2) * xp[i - 2][j] + E * (i - 3) * xp[i - 4][j]) /
		 (i + 2 * j - 1);
    }
  }
  auto mat = mat_double(xp, size);
  for (uint64_t i = 0; i < xlen; i++)
    free(xp[i]);
  free(xp);
  SelfAdjointEigenSolver<MatrixXd> solver(mat);
  return solver.eigenvalues().minCoeff() >= 0;
}

bool double_well_single(double E, double xsq, uint64_t size) {
  double *x = (double *)malloc((2 * size - 1) * sizeof(double));
  x[0] = 1;
  x[1] = 0;
  x[2] = xsq;
  x[3] = 0;
  for (uint64_t i = 4; i < 2 * size - 1; i++) {
    x[i] = ((i - 2) * x[i - 2] + (i - 3) * x[i - 4] * E) / (i - 1);
  }
  auto mat = mat_single(x, size);
  free(x);
  SelfAdjointEigenSolver<MatrixXd> solver(mat);
  return solver.eigenvalues().minCoeff() >= 0;
}

bool harmonics_single(double E, uint64_t size) {
  double *x = (double *)malloc((2 * size - 1) * sizeof(double));
  x[0] = 1;
  x[1] = 0;
  for (int64_t i = 2; i < 2 * size - 1; i++) {
    x[i] = E * (i - 1) / i * x[i - 2];
  }
  auto mat = mat_single(x, size);
  free(x);
  SelfAdjointEigenSolver<MatrixXd> solver(mat);
  return solver.eigenvalues().minCoeff() >= 0;
}

bool harmonics(double E, uint64_t size) {
  uint64_t xlen = max(2 * size - 1, 4ul);
  double **xp = (double **)malloc(xlen * sizeof(double *));
  for (uint64_t i = 0; i < xlen; i++)
    xp[i] = (double *)malloc((2 * size - 1) * sizeof(double));
  xp[0][0] = 1;
  xp[1][0] = 0;
  for (int64_t i = 2; i < xlen; i++)
    xp[i][0] = E * (i - 1) / i * xp[i - 2][0];
  for (int64_t i = 0; i < xlen; i++)
    xp[i][1] = 0;
  for (int64_t j = 2; j < 2 * size - 1; j++) {
    for (int64_t i = 0; i < 2; i++)
      xp[i][j] = E * xp[i][j - 2] - xp[i + 2][j - 2];
    for (int64_t i = 2; i < xlen; i++)
      xp[i][j] = E * (i - 1) / (i + j) * xp[i - 2][j];
  }
  auto mat = mat_double(xp, size);
  for (uint64_t i = 0; i < xlen; i++)
    free(xp[i]);
  free(xp);
  SelfAdjointEigenSolver<MatrixXd> solver(mat);
  return solver.eigenvalues().minCoeff() >= 0;
}

bool coulomb_single(double E, uint64_t size) {
  double *x = (double *)malloc((2 * size - 1) * sizeof(double));
  x[0] = 1;
  x[1] = -1 - 3. / 4 / E;
  for (int64_t i = 2; i < 2 * size - 1; i++)
    x[i] = (2 * i * x[i - 2] - (1 + 2 * i) * x[i - 1]) / 2 / (i + 1) / E;
  auto mat = mat_single(x, size);
  free(x);
  SelfAdjointEigenSolver<MatrixXd> solver(mat);
  return solver.eigenvalues().minCoeff() >= 0;
}

bool coulomb(double E, uint64_t size) {
  uint64_t xlen = max((2 * size - 1), 4lu);
  double **xp = (double **)malloc(xlen * sizeof(double *));
  for (uint64_t i = 0; i < xlen; i++)
    xp[i] = (double *)malloc((2 * size - 1) * sizeof(double));
  xp[0][0] = 1;
  xp[1][0] = -1 - 3. / 4 / E;
  for (int64_t i = 2; i < xlen; i++)
    xp[i][0] =
	(2 * i * xp[i - 2][0] - (1 + 2 * i) * xp[i - 1][0]) / 2 / (i + 1) / E;
  for (uint64_t i = 0; i < xlen; i++)
    xp[i][1] = 0;
  for (int64_t j = 2; j < 2 * size - 1; j++) {
    for (int64_t i = 2; i < xlen; i++)
      xp[i][j] = E * xp[i][j - 2] + xp[i - 1][j - 2] - xp[i - 2][j - 2];
    for (int64_t i = 1; i >= 0; i--)
      xp[i][j] =
	  (2 * (i + 3) * E * xp[i + 2][j] - (j - 2 * i - 5) * xp[i + 1][j]) /
	  2 / (i - j + 2);
  }
  auto mat = mat_double(xp, size);
  for (uint64_t i = 0; i < xlen; i++)
    free(xp[i]);
  free(xp);
  SelfAdjointEigenSolver<MatrixXd> solver(mat);
  return solver.eigenvalues().minCoeff() >= 0;
}

bool toda(double E, double ex, uint64_t size) {
  uint64_t xlen = max(2 * size - 1, 4ul);
  double **xp = (double **)malloc(xlen * sizeof(double *));
  for (uint64_t i = 0; i < xlen; i++)
    xp[i] = (double *)malloc((2 * size - 1) * sizeof(double));
  xp[0][0] = 1;
  xp[1][0] = ex;
  for (int64_t i = 2; i < xlen; i++)
    xp[i][0] = (2 * E * (i - 1) * xp[i - 1][0] + (-2 * i + 3) * xp[i - 2][0]) /
	       (2 * i - 1);
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
  auto mat = mat_double(xp, size);
  for (uint64_t i = 0; i < xlen; i++)
    free(xp[i]);
  free(xp);
  SelfAdjointEigenSolver<MatrixXd> solver(mat);
  return solver.eigenvalues().minCoeff() >= 0;
}

bool toda_single(double E, double ex, uint64_t size) {
  double *x = (double *)malloc((2 * size - 1) * sizeof(double));
  x[0] = 1.;
  x[1] = ex;
  for (int64_t i = 2; i < 2 * size - 1; i++) {
    x[i] = (2 * E * (i - 1) * x[i - 1] + (-2 * i + 3) * x[i - 2]) / (2 * i - 1);
  }
  auto mat = mat_single(x, size);
  free(x);
  SelfAdjointEigenSolver<MatrixXd> solver(mat);
  return solver.eigenvalues().minCoeff() >= 0;
}

bool trig_single(double E, double ex, uint64_t size) {
  double *x = (double *)malloc((2 * size - 1) * sizeof(double));
  int64_t offset = size - 1;
  x[offset] = 1;
  x[1 + offset] = ex;
  for (int64_t i = -1; i > -size; i--) {
    x[i + offset] =
	((2 * i + 3) * x[i + 2 + offset] - 2 * (i + 1) * x[i + 1 + offset]) /
	(-1 - 2 * i);
  }
  for (int64_t i = 2; i < size; i++) {
    x[i + offset] = (2 * (i - 1) * E * x[i - 1 + offset] +
		     (3 - 2 * i) * x[i - 2 + offset]) /
		    (2 * i - 1);
  }
  MatrixXd mat(size, size);
  for (uint64_t i = 0; i < size; i++) {
    for (uint64_t j = 0; j < size; j++) {
      mat(i, j) = x[-i + j + offset];
    }
  }
  free(x);
  SelfAdjointEigenSolver<MatrixXd> solver(mat);
  return solver.eigenvalues().minCoeff() >= 0;
}

// bool trig(double E, double ex, uint64_t size) {
//     double **xp = (double **)malloc((2*size - 2)*sizeof(double *));
//     for (uint64_t i = 0; i < 2*size - 1; i++) xp[i] = (double*)malloc((2*size
//     - 1)*sizeof(double)); int64_t offset = size - 1; xp[offset][0] = 1; xp[1
//     + offset][0] = ex; for (int64_t i = -1; i > -size; i--) xp[i + offset][0]
//     = ((2*i + 3)*xp[i + 2 + offset][0] - 2*(i + 1)*xp[i + 1 + offset][0])/(-
//     1 - 2*i); for (int64_t i = 2; i < size; i++) xp[i + offset][0] = (2*(i -
//     1)*E*xp[i - 1 + offset][0] + (3 - 2*i)*xp[i - 2 + offset][0])/(2*i - 1);
//     for (int64_t i = 0; i < 2*size - 1; i++) xp[i][1] = 0;
//     for (int64_t j = 2; j < 2*size - 1; j++) {
//         for (int64_t i = 0; i < 2; i++) {
//             xp[i + offset][j] = E*xp[i + offset][j - 2] - xp[i + offset +
//             1][j - 2] - xp[i + offset - 1][j - 2];
//         }
//         for (int64_t i = -1; i > -size; i--) xp[i + offset][j] = ((2*i + j +
//         3)*xp[i + 2 + offset][j] - 2*(i + 1)*xp[i + 1 + offset][j])/(j - 2*i
//         - 1); for (int64_t i = 2; i < size; i++) xp[i + offset][j] = (2*(i -
//         1)*E*xp[i - 1 + offset][j] + (j - 2*i + 3)*xp[i - 2 +
//         offset][j])/(2*i + j - 1);
//     }
//     MatrixXd mat(size*size, size*size);
//     for (int64_t i = 0; i < size*size; i++) {
//         for (int64_t j = i; j < size*size; j++) {
//             int64_t a = i/size;
//             int64_t b = i%size;
//             int64_t c = j/size;
//             int64_t d = j%size;
//         }
//     }
// }
