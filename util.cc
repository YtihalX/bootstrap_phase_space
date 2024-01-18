#include "util.h"
void linspace(double *x, double start, double end, uint64_t len) {
    double step = (end - start)/len;
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
    Eigen::MatrixXd mat(k*k, k*k);
    for (uint64_t i = 0; i < k*k; i++) {
        for (uint64_t j = i; j < k*k; j++) {
            uint64_t a = i/k;
            uint64_t b = i%k;
            uint64_t c = j/k;
            uint64_t d = j%k;
            mat(i, j) = xp[a + c][b + d];
            mat(j, i) = mat(i, j);
        }
    }
    return mat;
}

uint64_t bootstrap_single(vector<double> *x_axis, vector<double> *y_axis, bool (*model)(double, uint64_t)) {

    array<double, NUM_X> x_parameter;
    linspace(x_parameter.data(), 0., 10., x_parameter.size());
    // printf("%f, %f\n", x_range[124], x_range.back());

    vector<double> x_param_allowed[NUM_THREADS];

    pthread_t threads[NUM_THREADS];
    struct ThreadArg pargs[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; i++) {
        x_param_allowed[i].reserve(x_parameter.size());
        pargs[i].x_range = x_parameter.data();
        pargs[i].x_len = x_parameter.size();
        pargs[i].y_range = NULL;
        pargs[i].y_len = 0;
        pargs[i].x_param_allowed = x_param_allowed + i;
        pargs[i].y_param_allowed = NULL;
        pargs[i].thread_id = i;
        pargs[i].model = (void *)model;

        if (pthread_create(threads + i, NULL, routine_single, (void *)(pargs + i)) != 0) {
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

    uint64_t len = x_parameter.size();
    uint64_t res = len%NUM_THREADS;
    for (uint64_t i = len - res; i < len; i++) {
        if (model(x_parameter[i], SIZE_MAT)) {
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
        x_axis->insert(x_axis->end(), x_param_allowed[i].begin(), x_param_allowed[i].end());
    }

    return length;

}

uint64_t bootstrap_double(vector<double> *x_axis, vector<double> *y_axis, bool (*model)(double, double, uint64_t)) {

    array<double, NUM_Y> y_parameter;
    array<double, NUM_X> x_parameter;
    linspace(x_parameter.data(), -0.27, 2., x_parameter.size());
    linspace(y_parameter.data(), 0., 0.9, y_parameter.size());
    // printf("%f, %f\n", x_range[124], x_range.back());

    vector<double> x_param_allowed[NUM_THREADS];
    vector<double> y_param_allowed[NUM_THREADS];

    pthread_t threads[NUM_THREADS];
    struct ThreadArg pargs[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; i++) {
        y_param_allowed[i].reserve(y_parameter.size());
        x_param_allowed[i].reserve(x_parameter.size());
        pargs[i].x_range = x_parameter.data();
        pargs[i].x_len = x_parameter.size();
        pargs[i].y_range = y_parameter.data();
        pargs[i].y_len = y_parameter.size();
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

    uint64_t res = y_parameter.size()*x_parameter.size()%NUM_THREADS;
    for (uint64_t j = y_parameter.size() - res; j < y_parameter.size(); j++) {
        uint64_t i = x_parameter.size() - 1;
        if (model(x_parameter[i], y_parameter[j], SIZE_MAT)) {
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
        x_axis->insert(x_axis->end(), x_param_allowed[i].begin(), x_param_allowed[i].end());
        y_axis->insert(y_axis->end(), y_param_allowed[i].begin(), y_param_allowed[i].end());
    }

    return length;
}
