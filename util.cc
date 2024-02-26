#include "util.h"
void linspace(double *x, double start, double end, uint64_t len) {
  double step = (end - start) / len;
  x[0] = start;
  for (int i = 1; i < len; i++) {
    x[i] = x[i - 1] + step;
  }
  return;
}

Eigen::MatrixXd mat_single(double *x, uint64_t k) {
  Eigen::MatrixXd mat(k, k);
  for (uint64_t i = 0; i < k; i++) {
    for (uint64_t j = i; j < k; j++) {
      mat(i, j) = x[i + j];
      mat(j, i) = x[i + j];
    }
  }
  return mat;
}

Eigen::MatrixXd mat_double(double **xp, uint64_t k) {
  Eigen::MatrixXd mat(k * k, k * k);
  for (uint64_t i = 0; i < k * k; i++) {
    for (uint64_t j = i; j < k * k; j++) {
      uint64_t a = i / k;
      uint64_t b = i % k;
      uint64_t c = j / k;
      uint64_t d = j % k;
      mat(i, j) = xp[a + c][b + d];
      mat(j, i) = mat(i, j);
    }
  }
  return mat;
}

MatrixXmp mat_double_mpfr(double **xp, uint64_t k) {
  MatrixXmp mat(k * k, k * k);
  for (uint64_t i = 0; i < k * k; i++) {
    for (uint64_t j = i; j < k * k; j++) {
      uint64_t a = i / k;
      uint64_t b = i % k;
      uint64_t c = j / k;
      uint64_t d = j % k;
      mat(i, j) = xp[a + c][b + d];
      mat(j, i) = mat(i, j);
    }
  }
  return mat;
}

uint64_t bootstrap_single(vector<double> *x_axis, vector<double> *y_axis,
			  bool (*model)(double, uint64_t)) {

  double *x_parameter = (double *)malloc(num_x * sizeof(double));
  linspace(x_parameter, plot_range[0][0], plot_range[0][1], num_x);
  //  for (uint64_t i = 0; i < num_x; i++) {
  //    fprintf(stderr, "%f ", x_parameter[i]);
  //  }
  //  fprintf(stderr, "\n");
  // printf("%f, %f\n", x_range[124], x_range.back());

  vector<double> x_param_allowed[NUM_THREADS];

  pthread_t threads[NUM_THREADS];
  struct ThreadArg pargs[NUM_THREADS];

  for (int i = 0; i < NUM_THREADS; i++) {
    x_param_allowed[i].reserve(num_x);
    pargs[i].x_range = x_parameter;
    pargs[i].x_len = num_x;
    pargs[i].y_range = NULL;
    pargs[i].y_len = 0;
    pargs[i].x_param_allowed = x_param_allowed + i;
    pargs[i].y_param_allowed = NULL;
    pargs[i].thread_id = i;
    pargs[i].model = (void *)model;

    if (pthread_create(threads + i, NULL, routine_single,
		       (void *)(pargs + i)) != 0) {
      perror("pthread creation error");
      exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; i < NUM_THREADS; i++) {
    if (pthread_join(threads[i], NULL) != 0) {
      perror("pthread join error");
      exit(EXIT_FAILURE);
    }
  }

  uint64_t len = num_x;
  uint64_t res = len % NUM_THREADS;
  for (uint64_t i = len - res; i < len; i++) {
    if (model(x_parameter[i], kk)) {
      x_param_allowed[0].push_back(x_parameter[i]);
    }
  }

  uint64_t length = 0;
  for (uint64_t i = 0; i < NUM_THREADS; i++) {
    length += x_param_allowed[i].size();
  }
  // printf("length: %lu\n", length);
  x_axis->reserve(length);
  y_axis->assign(length, 0);
  for (int i = 0; i < NUM_THREADS; i++) {
    x_axis->insert(x_axis->end(), x_param_allowed[i].begin(),
		   x_param_allowed[i].end());
  }
  free(x_parameter);

  return length;
}

uint64_t bootstrap_double(vector<double> *x_axis, vector<double> *y_axis,
			  bool (*model)(double, double, uint64_t)) {

  double *x_parameter = (double *)malloc(num_x * sizeof(double));
  double *y_parameter = (double *)malloc(num_y * sizeof(double));
  linspace(x_parameter, plot_range[0][0], plot_range[0][1], num_x);
  linspace(y_parameter, plot_range[1][0], plot_range[1][1], num_y);
  // linspace(x_parameter, -6, 6, num_x);
  // linspace(y_parameter, -3, 3, num_y);
  // printf("%f, %f\n", x_range[124], x_range.back());

  vector<double> x_param_allowed[NUM_THREADS];
  vector<double> y_param_allowed[NUM_THREADS];

  pthread_t threads[NUM_THREADS];
  struct ThreadArg pargs[NUM_THREADS];

  for (int i = 0; i < NUM_THREADS; i++) {
    y_param_allowed[i].reserve(num_y);
    x_param_allowed[i].reserve(num_x);
    pargs[i].x_range = x_parameter;
    pargs[i].x_len = num_x;
    pargs[i].y_range = y_parameter;
    pargs[i].y_len = num_y;
    pargs[i].x_param_allowed = x_param_allowed + i;
    pargs[i].y_param_allowed = y_param_allowed + i;
    pargs[i].thread_id = i;
    pargs[i].model = (void *)model;

    if (pthread_create(threads + i, NULL, routine, (void *)(pargs + i)) != 0) {
      perror("pthread creation error");
      exit(EXIT_FAILURE);
    }
  }

  for (int i = 0; i < NUM_THREADS; i++) {
    if (pthread_join(threads[i], NULL) != 0) {
      perror("pthread join error");
      exit(EXIT_FAILURE);
    }
  }

  uint64_t res = num_y * num_x % NUM_THREADS;
  for (uint64_t j = num_y - res; j < num_y; j++) {
    uint64_t i = num_x - 1;
    if (model(x_parameter[i], y_parameter[j], kk)) {
      y_param_allowed[0].push_back(y_parameter[j]);
      x_param_allowed[0].push_back(x_parameter[i]);
    }
  }

  uint64_t length = 0;
  for (uint64_t i = 0; i < NUM_THREADS; i++) {
    length += y_param_allowed[i].size();
  }
  // printf("length: %lu\n", length);
  x_axis->reserve(length);
  y_axis->reserve(length);
  for (int i = 0; i < NUM_THREADS; i++) {
    x_axis->insert(x_axis->end(), x_param_allowed[i].begin(),
		   x_param_allowed[i].end());
    y_axis->insert(y_axis->end(), y_param_allowed[i].begin(),
		   y_param_allowed[i].end());
  }

  free(x_parameter);
  free(y_parameter);

  return length;
}
