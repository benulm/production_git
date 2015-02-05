#ifndef DRIVER_
#define DRIVER_
#include <string>

const int print = 100000;
const double tol = 3e-5;
const std::string outputattatch = "obstac";
std::string outputdir;
const int n_steps_between_avgs = 100;

#define TMS_USED_
#define ROUGHNESS_

#define PRINT_STATE_
//#define TIME_AVERAGE_

//#define ENTROPIC_COLLIDE_
#define KBC_COLLIDE_

#endif // DRIVER_
