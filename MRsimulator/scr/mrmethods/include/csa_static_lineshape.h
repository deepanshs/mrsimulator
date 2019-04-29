

#include "c_array.h"
#include <time.h>
#include "powder_setup.h"
#include "OCPowderScheme.h"
#include "core.h"

extern void lineshape_csa_static(
        double * spec,
        double * cpu_time_,                  
        int m,
        int nt,
        double fstart,
        double fwidth,
        double iso,
        double aniso,
        double eta,
        int octa,
        int npros);


