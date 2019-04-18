
#include "OCPowderScheme.h"
#include "sphere_lebedev_rule.h"
#include "lebedev131.h"
#include "c_array.h"
#include "leboct31.h"

OCEulerAngle *OCEulerAngleCreatePowderAngles(char *filename, uint64_t *Nangles, double **weights)
{
    /* Read in the angles for the powder average */
    FILE *fp;
    if ((fp = fopen(filename, "r")) == NULL) {
        printf("Cannot open powder average angles file %s\n",filename);
        exit(1);
    }
    /* Next find out how many lines are in file */
    uint64_t Nangle=1;
    int c;
    while((c=fgetc(fp))!=EOF) {
        char cc = c;
        if((cc=='\n')||(cc=='\r')) (Nangle)++;
    }
    rewind(fp);
    
    OCEulerAngle *powderAngles = malloc(sizeof(OCEulerAngle)*Nangle);
    *weights = malloc(sizeof(double)*Nangle);
    
    
    char buffer[256];
    fgets(buffer, 256, fp);
    for(int i=0;i<=255;i++)
        if((buffer[i]=='\r')||(buffer[i]=='\n')) buffer[i] = 0;
    
    double alpha, beta, gamma, weight;
    
    int angleset = sscanf(buffer, "%lf %lf %lf %lf",
                      &alpha,
                      &beta,
                      &gamma,
                      &weight);
    
    rewind(fp);
    angleset -=1;
    char string[120];
    sprintf(string,"\nReading in %llu Powder Average %d Angle Set from file %s...\n",Nangle,angleset,filename);
//    Log(WRITE_LOG, string);
    
    if(angleset==3) {
        for(uint64_t index=0;index<Nangle;index++) {
            fscanf(fp, "%lg %lg %lg %lg",&alpha,&beta,&gamma,&weight);
            OCEulerAngle angle = {alpha, beta, gamma};
            powderAngles[index] = angle;
            *weights[index] = weight;
        }
    }
    else {
        for(uint64_t index=0;index<Nangle;index++) {
            gamma = 0.0;
            fscanf(fp, "%lg %lg %lg",&alpha,&beta,&weight);
            OCEulerAngle angle = {alpha, beta, gamma};
            powderAngles[index] = angle;
            *weights[index] = weight;
        }
    }
    fclose(fp);
    return powderAngles;
}


OCPowderScheme OCCreatePowderScheme(int scheme, int size)
{
    double alpha, beta;
    double pi = 3.14159265358979323846;
    double deg_to_rad = pi/180.;
    
    if (scheme == 1){
        double *x = createDouble1DArray(size);
        double *y = createDouble1DArray(size);
        double *z = createDouble1DArray(size);
        double *w = createDouble1DArray(size);

        OCEulerAngle *powderAngles = malloc(sizeof(OCEulerAngle)*size);
        ld_by_order(size, x, y, z, w);
        for(int i=0; i<size; i++){
            xyz_to_tp(x[i], y[i], z[i], &beta, &alpha);
            OCEulerAngle angle = {alpha*deg_to_rad, beta*deg_to_rad, 0.0};
            powderAngles[i] = angle;
        }
        OCPowderScheme powder_scheme = {
            powderAngles,
            w,
            size
        };
        // powder_scheme.name = "1";
        // powder_scheme.size = size;
        // powder_scheme.angles = powderAngles;
        // powder_scheme.weights = w;

        destroyDouble1DArray(x);
        destroyDouble1DArray(y);
        destroyDouble1DArray(z);
        return powder_scheme;
    }

    if (scheme == 2){
        int N = Lebedev_Laikov_npoint(size);
        printf("Lebedev_Laikov_npoint %i", N);
        double *alpha = createDouble1DArray(N);
        double *beta = createDouble1DArray(N);
        double *weights = createDouble1DArray(N);

        
        Lebedev_quad_arrays(N, alpha, beta, weights);

        OCEulerAngle *powderAngles = malloc(sizeof(OCEulerAngle)*N);
        for(int i=0; i<N; i++){
            OCEulerAngle angle = {alpha[i], beta[i], 0.0};
            powderAngles[i] = angle;
        }
        OCPowderScheme powder_scheme = {
            powderAngles,
            weights,
            N
        };
        // powder_scheme.name = "1";
        // powder_scheme.size = size;
        // powder_scheme.angles = powderAngles;
        // powder_scheme.weights = w;

        destroyDouble1DArray(alpha);
        destroyDouble1DArray(beta);
        return powder_scheme;
    }
    // lebedev octahedral l=31 //
    if (scheme == 3){
        // double *alpha = createDouble1DArray(31);
        // double *beta = createDouble1DArray(31);
        // double *weights = createDouble1DArray(31);


        // OCEulerAngle *powderAngles = malloc(sizeof(OCEulerAngle)*31);
        // for(int i=0; i<31; i++){
        //     OCEulerAngle angle = {alpha[i]*deg_to_rad, beta[i]*deg_to_rad, 0.0};
        //     powderAngles[i] = angle;
        // }
        OCPowderScheme powder_scheme = leboct31();
        // powder_scheme.name = "1";
        // powder_scheme.size = size;
        // powder_scheme.angles = powderAngles;
        // powder_scheme.weights = w;

        // destroyDouble1DArray(alpha);
        // destroyDouble1DArray(beta);
        return powder_scheme;
    }
}