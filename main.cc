#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <vector>
#include <array>
#include <stdint.h>
#include <eigen3/Eigen/Dense>
#include <algorithm>

using namespace ::std;

#define NUM_THREADS 12
#define NUM_X 10000
#define NUM_E 10000
#define SIZE_MAT 9
#define DEBUG_THREAD 1

void *routine(void *arg);
void *toda_routine(void *arg);
void linspace(double *x, double start, double end, uint64_t len);
bool double_well(double E, double xsq, uint64_t size);
bool toda(double E, double ex, uint64_t size);
bool double_well_single(double E, double xsq, uint64_t size);
bool double_well_fixed(double x1, double x3, uint64_t size);

struct ThreadArg {
    bool (*model)(double, double, uint64_t);
    double *x_range;
    double *E_range;
    uint64_t x_len;
    uint64_t E_len;
    vector<double> *y_param_allowed;
    vector<double> *x_param_allowed;
    int thread_id;
};

int main() {
    array<double, NUM_X> y_parameter;
    array<double, NUM_E> x_parameter;
    linspace(y_parameter.data(), 0., 0.5, y_parameter.size());
    linspace(x_parameter.data(), -0.28, 0.3, x_parameter.size());
    // printf("%f, %f\n", x_range[124], x_range.back());

    vector<double> y_param_allowed[NUM_THREADS];
    vector<double> x_param_allowed[NUM_THREADS];

    pthread_t threads[NUM_THREADS];
    struct ThreadArg pargs[NUM_THREADS];

    for (int i = 0; i < NUM_THREADS; i++) {
        y_param_allowed[i].reserve(y_parameter.size());
        x_param_allowed[i].reserve(x_parameter.size());
        pargs[i].E_range = x_parameter.data();
        pargs[i].E_len = x_parameter.size();
        pargs[i].x_range = y_parameter.data();
        pargs[i].x_len = y_parameter.size();
        pargs[i].x_param_allowed = x_param_allowed + i;
        pargs[i].y_param_allowed = y_param_allowed + i;
        pargs[i].thread_id = i;
        pargs[i].model = &double_well_single;

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
        if (double_well(x_parameter[i], y_parameter[j], SIZE_MAT)) {
            y_param_allowed[0].push_back(y_parameter[j]);
            x_param_allowed[0].push_back(x_parameter[i]);
        }
    }
    
    uint64_t length = 0;
    for (uint64_t i = 0; i < NUM_THREADS; i++) {
        length += y_param_allowed[i].size();
    }
    // printf("length: %lu\n", length);
    vector<double> y_axis;
    vector<double> x_axis;
    y_axis.reserve(length);
    x_axis.reserve(length);
    for (int i = 0; i < NUM_THREADS; i++) {
        y_axis.insert(y_axis.end(), y_param_allowed[i].begin(), y_param_allowed[i].end());
        x_axis.insert(x_axis.end(), x_param_allowed[i].begin(), x_param_allowed[i].end());
    }

    FILE *data = fopen("data.csv", "w");
    if (data != NULL) {
        fprintf(data, "E,xsq\n");
        for (uint64_t i = 0; i < length; i++) {
            fprintf(data, "%f,%f\n", x_axis[i], y_axis[i]);
        }
        fclose(data);
    }
    return 0;
}

void linspace(double *x, double start, double end, uint64_t len) {
    double step = (end - start)/len;
    x[0] = start;
    for (int i = 1; i < len; i++) {
        x[i] = x[i - 1] + step;
    }
    return;
}

void *routine(void *arg) {
    ThreadArg parg = *(ThreadArg *)arg;
    uint64_t step = parg.x_len*parg.E_len/NUM_THREADS;
    uint64_t start = parg.thread_id * step;
    uint64_t end = start + step;
    uint64_t i = start/parg.x_len;
    uint64_t j = start%parg.x_len;
    uint64_t ist = i;
    uint64_t ied = end/parg.x_len;
    uint64_t istep = (ied - ist)/10;
    fprintf(stderr, "id: %u, start: %lu, end: %lu, i: %lu, j: %lu\n\n", parg.thread_id, start, end, i, j);
    while (j < parg.x_len) {
        // if (parg.thread_id == DEBUG_THREAD) printf("i: %lu, j: %lu, E: %f, x: %f\n", i, j, parg.E_range[i], parg.x_range[j]);
        if (parg.model(parg.E_range[i], parg.x_range[j], SIZE_MAT)) {
            // if (parg.thread_id == 0) printf("are you ok?\n");
            parg.y_param_allowed->push_back(parg.x_range[j]);
            parg.x_param_allowed->push_back(parg.E_range[i]);
        }
        j++;
    }
    while (i < end/parg.x_len) {
        for (uint64_t idx = 1; idx < 10; idx++) {
            if (i == ist + istep*idx) printf("thread %u: %lu%%\n\n", parg.thread_id, idx*10);
        }
        for (j = 0; j < parg.x_len; j++) {
        // if (parg.thread_id == DEBUG_THREAD) printf("i: %lu, j: %lu, E: %f, x: %f\n", i, j, parg.E_range[i], parg.x_range[j]);
            if (parg.model(parg.E_range[i], parg.x_range[j], SIZE_MAT)) {
            // if (parg.thread_id == 0) printf("are you ok?\n");
                parg.y_param_allowed->push_back(parg.x_range[j]);
                parg.x_param_allowed->push_back(parg.E_range[i]);
            }
        }
        i++;
    }
    for (j = 0; j < end%parg.x_len; j++) {
        // if (parg.thread_id == DEBUG_THREAD) printf("i: %lu, j: %lu, E: %f, x: %f\n", i, j, parg.E_range[i], parg.x_range[j]);
        if (parg.model(parg.E_range[i], parg.x_range[j], SIZE_MAT)) {
            // if (parg.thread_id == 0) printf("are you ok?\n");
            parg.y_param_allowed->push_back(parg.x_range[j]);
            parg.x_param_allowed->push_back(parg.E_range[i]);
        }
    }
    fprintf(stderr, "id: %u, len: %lu\n", parg.thread_id, parg.y_param_allowed->size());
    return NULL;
}

bool double_well(double E, double xsq, uint64_t size) {
    uint64_t xlen = max(2*size - 1, 8ul);
    double **xp = (double **)malloc(xlen*sizeof(double *));
    for (uint64_t i = 0; i < xlen; i++) xp[i] = (double*)malloc((2*size - 1)*sizeof(double));
    xp[0][0] = 1;
    xp[1][0] = 0.;
    xp[2][0] = xsq;
    xp[3][0] = 0.;
    for (uint64_t i = 4; i < xlen; i++) {
        xp[i][0] = ((i - 2)*xp[i - 2][0] + (i - 3)*xp[i - 4][0]*E)/(i - 1);
    }
    for (uint64_t i = 0; i < 2*size - 1; i++) {
        xp[i][1] = 0;
    }
    for (uint64_t j = 2; j < 2*size - 1; j++) {
        for (uint64_t i = 0; i < 4; i++) {
            xp[i][j] = E*xp[i][j - 2] + xp[i + 2][j - 2] - xp[i + 4][j - 2];
        }
        for (uint64_t i = 4; i < xlen; i++) {
            xp[i][j] = ((i + j - 2)*xp[i - 2][j] + E*(i - 3)*xp[i - 4][j])/(i + 2*j - 1);
        }
    }
    Eigen::MatrixXd mat(size*size, size*size);
    for (uint64_t i = 0; i < size*size; i++) {
        for (uint64_t j = i; j < size*size; j++) {
            uint64_t a = i/size;
            uint64_t b = i%size;
            uint64_t c = j/size;
            uint64_t d = j%size;
            mat(i, j) = xp[a + c][b + d];
            mat(j, i) = mat(i, j);
        }
    }
    for (uint64_t i = 0; i < xlen; i++) free(xp[i]);
    free(xp);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
    return solver.eigenvalues().minCoeff() >= 0;
}

bool double_well_fixed(double x1, double x3, uint64_t size) {
    double E = -0.1;
    uint64_t xlen = max(2*size - 1, 8ul);
    double **xp = (double **)malloc(xlen*sizeof(double *));
    for (uint64_t i = 0; i < xlen; i++) xp[i] = (double*)malloc((2*size - 1)*sizeof(double));
    xp[0][0] = 1;
    xp[1][0] = x1;
    xp[2][0] = 0.4052;
    xp[3][0] = x3;
    for (uint64_t i = 4; i < xlen; i++) {
        xp[i][0] = ((i - 2)*xp[i - 2][0] + (i - 3)*xp[i - 4][0]*E)/(i - 1);
    }
    for (uint64_t i = 0; i < 2*size - 1; i++) {
        xp[i][1] = 0;
    }
    for (uint64_t j = 2; j < 2*size - 1; j++) {
        for (uint64_t i = 0; i < 4; i++) {
            xp[i][j] = E*xp[i][j - 2] + xp[i + 2][j - 2] - xp[i + 4][j - 2];
        }
        for (uint64_t i = 4; i < xlen; i++) {
            xp[i][j] = ((i + j - 2)*xp[i - 2][j] + E*(i - 3)*xp[i - 4][j])/(i + 2*j - 1);
        }
    }
    Eigen::MatrixXd mat(size*size, size*size);
    for (uint64_t i = 0; i < size*size; i++) {
        for (uint64_t j = i; j < size*size; j++) {
            uint64_t a = i/size;
            uint64_t b = i%size;
            uint64_t c = j/size;
            uint64_t d = j%size;
            mat(i, j) = xp[a + c][b + d];
            mat(j, i) = mat(i, j);
        }
    }
    for (uint64_t i = 0; i < xlen; i++) free(xp[i]);
    free(xp);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
    return solver.eigenvalues().minCoeff() >= 0;
}

bool double_well_single(double E, double xsq, uint64_t size) {
    double *x = (double *)malloc((2*size - 1)*sizeof(double));
    x[0] = 1;
    x[1] = 0;
    x[2] = xsq;
    x[3] = 0;
    for (uint64_t i = 4; i < 2*size - 1; i++) {
        x[i] = ((i - 2)*x[i - 2] + (i - 3)*x[i - 4]*E)/(i - 1);
    }
    Eigen::MatrixXd mat(size, size);
    for (uint64_t i = 0; i < size; i++) {
        for (uint64_t j = i; j < size; j++) {
            mat(i, j) = x[i + j];
            mat(j, i) = x[i + j];
        }
    }
    free(x);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
    return solver.eigenvalues().minCoeff() >= 0;
}

bool toda(double E, double ex, uint64_t size) {
    uint64_t xlen = max(2*size - 1, 4ul);
    double **xp = (double **)malloc(xlen*sizeof(double *));
    for (uint64_t i = 0; i < xlen; i++) xp[i] = (double*)malloc((2*size - 1)*sizeof(double));
    xp[0][0] = 1;
    xp[1][0] = ex;
    for (int64_t i = 2; i < xlen; i++) xp[i][0] = (2*E*(i - 1)*xp[i - 1][0] + (- 2*i + 3)*xp[i - 2][0])/(2*i - 1);
    for (int64_t i = 0; i < xlen; i++) xp[i][1] = 0;
    for (int64_t j = 2; j < 2*size - 1; j++) {
        for (int64_t i = 1; i < 3; i++) xp[i][j] = E*xp[i][j - 2] - xp[i + 1][j - 2] - xp[i - 1][j - 2];
        xp[0][j] = ((3 + j)*xp[2][j] - 2*E*xp[1][j])/(j - 1);
        for (int64_t i = 3; i < xlen; i++) xp[i][j] = (2*E*(i - 1)*xp[i - 1][j] + (j - 2*i + 3)*xp[i - 2][j])/(2*i + j - 1);
    }
    Eigen::MatrixXd mat(size*size, size*size);
    for (uint64_t i = 0; i < size*size; i++) {
        for (uint64_t j = i; j < size*size; j++) {
            uint64_t a = i/size;
            uint64_t b = i%size;
            uint64_t c = j/size;
            uint64_t d = j%size;
            mat(i, j) = xp[a + c][b + d];
            mat(j, i) = mat(i, j);
        }
    }
    for (uint64_t i = 0; i < xlen; i++) free(xp[i]);
    free(xp);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
    return solver.eigenvalues().minCoeff() >= 0;
}

bool trig_single(double E, double ex, uint64_t size) {
    double *x = (double *)malloc((2*size - 1)*sizeof(double));
    int64_t offset = size - 1;
    x[offset] = 1;
    x[1 + offset] = ex;
    for (int64_t i = -1; i > -size; i--) {
        x[i + offset] = ((2*i + 3)*x[i + 2 + offset] - 2*(i + 1)*x[i + 1 + offset])/(- 1 - 2*i);
    }
    for (int64_t i = 2; i < size; i++) {
        x[i + offset] = (2*(i - 1)*E*x[i - 1 + offset] + (3 - 2*i)*x[i - 2 + offset])/(2*i - 1);
    }
    Eigen::MatrixXd mat(size, size);
    for (uint64_t i = 0; i < size; i++) {
        for (uint64_t j = 0; j < size; j++) {
            mat(i, j) = x[- i + j + offset];
        }
    }
    free(x);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(mat);
    return solver.eigenvalues().minCoeff() >= 0;
}

//bool trig(double E, double ex, uint64_t size) {
//    double **xp = (double **)malloc((2*size - 2)*sizeof(double *));
//    for (uint64_t i = 0; i < 2*size - 1; i++) xp[i] = (double*)malloc((2*size - 1)*sizeof(double));
//    int64_t offset = size - 1;
//    xp[offset][0] = 1;
//    xp[1 + offset][0] = ex;
//    for (int64_t i = -1; i > -size; i--) xp[i + offset][0] = ((2*i + 3)*xp[i + 2 + offset][0] - 2*(i + 1)*xp[i + 1 + offset][0])/(- 1 - 2*i);
//    for (int64_t i = 2; i < size; i++) xp[i + offset][0] = (2*(i - 1)*E*xp[i - 1 + offset][0] + (3 - 2*i)*xp[i - 2 + offset][0])/(2*i - 1);
//    for (int64_t i = 0; i < 2*size - 1; i++) xp[i][1] = 0;
//    for (int64_t j = 2; j < 2*size - 1; j++) {
//        for (int64_t i = 0; i < 2; i++) {
//            xp[i + offset][j] = E*xp[i + offset][j - 2] - xp[i + offset + 1][j - 2] - xp[i + offset - 1][j - 2];
//        }
//        for (int64_t i = -1; i > -size; i--) xp[i + offset][j] = ((2*i + j + 3)*xp[i + 2 + offset][j] - 2*(i + 1)*xp[i + 1 + offset][j])/(j - 2*i - 1);
//        for (int64_t i = 2; i < size; i++) xp[i + offset][j] = (2*(i - 1)*E*xp[i - 1 + offset][j] + (j - 2*i + 3)*xp[i - 2 + offset][j])/(2*i + j - 1);
//    }
//    Eigen::MatrixXd mat(size*size, size*size);
//    for (int64_t i = 0; i < size*size; i++) {
//        for (int64_t j = i; j < size*size; j++) {
//            int64_t a = i/size;
//            int64_t b = i%size;
//            int64_t c = j/size;
//            int64_t d = j%size;
//        }
//    }
//}

