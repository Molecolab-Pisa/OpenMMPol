#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int countLines(char *fin){
    FILE *fp = fopen(fin, "r");
    
    char c;
    int lines = 1;

    if(fp == NULL) return 0;

    do{
        c = fgetc(fp);
        if(c == '\n') lines++;
    }while(c != EOF);

    fclose(fp);
  
    return lines - 1;
}

double **read_ef(char *fin){
    double **ef;
    
    int pol_atoms = countLines(fin);

    ef = (double **) malloc(sizeof(double *) * 3 * pol_atoms);

    FILE *fp = fopen(fin, "r");
    for(int i =0; i < pol_atoms; i++){
        ef[i] = (double *) malloc(sizeof(double) * 3);
        fscanf(fp, "%lf %lf %lf", &(ef[i][0]), &(ef[i][1]),  &(ef[i][2]));
    }
    fclose(fp);

    return ef;
}

int main(int argc, char **argv){
    if(argc != 4 && argc != 5){
        printf("Syntax expected\n");
        printf("    $ test_init.exe <INPUT FILE> <OUTPUT FILE> <SOLVER> [<ELECTRIC FIELD FILE>]\n");
        return 1;
    }
    
    int n_ipd, pol_atoms, solver;
    double *_ipd;
    double ***ipd, *electric_field;
    double **external_ef;
    int32_t *polar_mm;

    OMMP_SYSTEM_PRT my_system = ommp_init_mmp(argv[1]);
    // ommp_set_verbose(OMMP_VERBOSE_DEBUG);

    n_ipd = ommp_get_n_ipd(my_system);
    pol_atoms = ommp_get_pol_atoms(my_system);

    if(strcmp(argv[3], "inversion") == 0){
        printf("Solving with matrix inversion\n");
        solver = OMMP_SOLVER_INVERSION;
    }
    else if(strcmp(argv[3], "cg") == 0){
        printf("Solving with conjugate gradients\n");
        solver = OMMP_SOLVER_CG;
    }
    else if(strcmp(argv[3], "diis") == 0){
        printf("Solving with Jacobi-DIIS\n");
        solver = OMMP_SOLVER_DIIS;
    }
    else{
        printf("Unrecognized solver (\"%s\"). Exiting.\n", argv[3]);
        return 1;
    }
    
    electric_field = (double *) malloc(sizeof(double) * 3 * pol_atoms);
    polar_mm = (int32_t *) ommp_get_polar_mm(my_system);
   
    if(argc == 5)
        external_ef = read_ef(argv[4]);

    for(int j = 0; j < pol_atoms; j++)
        for(int k = 0; k < 3; k++)
            if(argc == 5)
                electric_field[j*3+k] = external_ef[polar_mm[j]-1][k];
            else
                electric_field[j*3+k] = 0.0;

    ommp_set_external_field(my_system, electric_field, solver);

    // Get induced point dipoles
    _ipd = (double *) ommp_get_ipd(my_system);
    ipd = (double ***) malloc(sizeof(double **) * n_ipd );

    for(int i = 0; i < n_ipd; i++){
        ipd[i] = (double **) malloc(sizeof(double *) * pol_atoms);

        for(int j = 0; j < pol_atoms; j++){
            ipd[i][j] = &(_ipd[i*pol_atoms*3 + j*3]);
        }
    }

    FILE *fp = fopen(argv[2], "w+");

    for(int k = 0; k < n_ipd; k++){
        // fprintf(fp, "IPD SET %d\n", k+1);
        
        for(int i = 0; i < pol_atoms; i++){
            for(int j = 0; j < 3; j++)
                fprintf(fp, "%20.12e", ipd[k][i][j]);
            fprintf(fp, "\n");
        }
    }
    
    fclose(fp);
    free(electric_field);
    
    for(int i = 0; i < n_ipd; i++)
        free(ipd[i]);
    free(ipd);
    
    if(argc == 5)
        free(external_ef);
    ommp_terminate(my_system);
    
    return 0;
}
