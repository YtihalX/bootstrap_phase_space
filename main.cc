#include "models.h"
#include "util.h"

#define MODEL 1

using namespace ::std;

int64_t num_x, num_y, kk;
double plot_range[2][2];

int main(int argc, char *argv[]) { // =Double_or_Single= =model= =num_x= =num_y=
				   // =k= =range_x range_x= =range_y range_y=
  if (argc != 10) {
    fprintf(stderr, "invalid usage, supposed to be 9 numbers\n");
    return 1;
  }
  int64_t args[5];
  for (uint64_t i = 0; i < 5; i++) {
    char *end;
    int64_t res = strtol(argv[i + 1], &end, 10);
    if (*end != '\0') {
      fprintf(stderr, "invalid usage, supposed to be 9 numbers\n");
      return 2;
    }
    args[i] = res;
  }
  num_x = args[2];
  num_y = args[3];
  kk = args[4];
  for (uint64_t i = 0; i < 4; i++) {
    char *end;
    double res = strtod(argv[i + 6], &end);
    if (*end != '\0') {
      fprintf(stderr, "invalid usage, supposed to be 9 numbers\n");
      return 3;
    }
    plot_range[i / 2][i % 2] = res;
    fprintf(stderr, "%f\n", res);
  }
  vector<double> x_axis;
  vector<double> y_axis;
  void *models[8] = {
      (void *)double_well_single, (void *)toda_single, (void *)harmonics_single,
      (void *)coulomb_single,	  (void *)double_well, (void *)toda,
      (void *)harmonics,	  (void *)coulomb,
  };
  uint64_t length;
  if (args[0] == 0) {
    length = bootstrap_single(&x_axis, &y_axis,
			      (bool (*)(double, uint64_t))models[args[1]]);
  } else {
    length = bootstrap_double(
	&x_axis, &y_axis, (bool (*)(double, double, uint64_t))models[args[1]]);
  }

  printf("E,xsq\n");
  for (uint64_t i = 0; i < length; i++) {
    printf("%.*f,%.*f\n", std::numeric_limits<double>::digits10, x_axis[i],
	   std::numeric_limits<double>::digits10, y_axis[i]);
  }
  return 0;
}
