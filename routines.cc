#include "routines.h"

void *routine(void *arg) {
  ThreadArg parg = *(ThreadArg *)arg;
  uint64_t step = parg.y_len * parg.x_len / NUM_THREADS;
  uint64_t start = parg.thread_id * step;
  uint64_t end = start + step;
  uint64_t i = start / parg.y_len;
  uint64_t j = start % parg.y_len;
  uint64_t ist = i;
  uint64_t ied = end / parg.y_len;
  uint64_t istep = (ied - ist) / 10;
  bool (*model)(double, double, uint64_t) =
      (bool (*)(double, double, uint64_t))parg.model;
  fprintf(stderr, "id: %u, start: %lu, end: %lu, i: %lu, j: %lu\n\n",
	  parg.thread_id, start, end, i, j);
  while (j < parg.y_len) {
    // if (parg.thread_id == DEBUG_THREAD) printf("i: %lu, j: %lu, E: %f, x:
    // %f\n", i, j, parg.E_range[i], parg.x_range[j]);
    if (model(parg.x_range[i], parg.y_range[j], kk)) {
      // if (parg.thread_id == 0) printf("are you ok?\n");
      parg.x_param_allowed->push_back(parg.x_range[i]);
      parg.y_param_allowed->push_back(parg.y_range[j]);
    }
    j++;
  }
  while (i < end / parg.y_len) {
    for (uint64_t idx = 1; idx < 10; idx++) {
      if (i == ist + istep * idx)
	fprintf(stderr, "thread %u: %lu%%\n\n", parg.thread_id, idx * 10);
    }
    for (j = 0; j < parg.y_len; j++) {
      // if (parg.thread_id == DEBUG_THREAD) printf("i: %lu, j: %lu, E: %f, x:
      // %f\n", i, j, parg.E_range[i], parg.x_range[j]);
      if (model(parg.x_range[i], parg.y_range[j], kk)) {
	// if (parg.thread_id == 0) printf("are you ok?\n");
	parg.x_param_allowed->push_back(parg.x_range[i]);
	parg.y_param_allowed->push_back(parg.y_range[j]);
      }
    }
    i++;
  }
  for (j = 0; j < end % parg.y_len; j++) {
    // if (parg.thread_id == DEBUG_THREAD) printf("i: %lu, j: %lu, E: %f, x:
    // %f\n", i, j, parg.E_range[i], parg.x_range[j]);
    if (model(parg.x_range[i], parg.y_range[j], kk)) {
      // if (parg.thread_id == 0) printf("are you ok?\n");
      parg.x_param_allowed->push_back(parg.x_range[i]);
      parg.y_param_allowed->push_back(parg.y_range[j]);
    }
  }
  fprintf(stderr, "id: %u, len: %lu\n", parg.thread_id,
	  parg.y_param_allowed->size());
  return NULL;
}

void *routine_single(void *arg) {
  ThreadArg parg = *(ThreadArg *)arg;
  uint64_t step = parg.x_len / NUM_THREADS;
  uint64_t start = parg.thread_id * step;
  uint64_t end = start + step;
  uint64_t ist = start;
  uint64_t istep = step / 10;
  bool (*model)(double, uint64_t) = (bool (*)(double, uint64_t))parg.model;
  fprintf(stderr, "id: %u, start: %lu, end: %lu\n\n", parg.thread_id, start,
	  end);
  while (start < end) {
    for (uint64_t i = 1; i < 10; i++) {
      if (start == ist + istep * i)
	fprintf(stderr, "thread %u: %lu%%\n\n", parg.thread_id, i * 10);
    }
    if (model(parg.x_range[start], kk))
      parg.x_param_allowed->push_back(parg.x_range[start]);
    start++;
  }
  fprintf(stderr, "id: %u, len: %lu\n", parg.thread_id,
	  parg.x_param_allowed->size());
  return NULL;
}
