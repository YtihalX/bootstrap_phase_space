#ifndef ROUTINES_H
#define ROUTINES_H

#define NUM_THREADS 12
#define SIZE_MAT 4
#define NUM_Y 10000lu
#define NUM_X 10000lu
#define DEBUG_THREAD 1

#include <pthread.h>
#include <vector>
#include <stdint.h>
#include <stdio.h>

void *routine(void *arg);
void *routine_single(void *arg);

struct ThreadArg {
    void *model;
    double *y_range;
    double *x_range;
    uint64_t y_len;
    uint64_t x_len;
    std::vector<double> *y_param_allowed;
    std::vector<double> *x_param_allowed;
    int thread_id;
};

#endif // !ROUTINES_H
