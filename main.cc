#include "util.h"
#include "models.h"

#define MODEL 0

using namespace ::std;

int main() {
    vector<double> x_axis;
    vector<double> y_axis;
    #if MODEL == 0
    auto length = bootstrap_double(&x_axis, &y_axis, toda);
    #elif MODEL == 1
    auto length = bootstrap_single(&x_axis, &y_axis, coulomb);
    #endif

    FILE *data = fopen("data.csv", "w");
    if (data != NULL) {
        fprintf(data, "E,xsq\n");
        for (uint64_t i = 0; i < length; i++) {
            fprintf(data, "%.*f,%.*f\n", std::numeric_limits<double>::digits10, x_axis[i], std::numeric_limits<double>::digits10, y_axis[i]);
        }
        fclose(data);
    }
    return 0;
}

