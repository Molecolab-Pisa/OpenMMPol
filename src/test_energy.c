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
    if(argc != 4 && argc != 3){
        printf("Syntax expected\n");
        printf("    $ test_init.exe <INPUT FILE> <OUTPUT FILE> [<ELECTRIC FIELD FILE>]\n");
        return 0;
    }
    
    int n_ipd, pol_atoms;
    double E_MMMM, E_MMPOL;
    double *electric_field;
    double **external_ef;
    int32_t *polar_mm;

    mmpol_init_mmp(argv[1]);
    set_verbose(OMMP_VERBOSE_NONE);
    n_ipd = get_n_ipd();
    pol_atoms = get_pol_atoms();
    
    electric_field = (double *) malloc(sizeof(double) * n_ipd * 3 * pol_atoms);
    polar_mm = (int32_t *) get_polar_mm();
    
    if(argc == 4)
        external_ef = read_ef(argv[3]);

    for(int i = 0; i < n_ipd; i++)
        for(int j = 0; j < pol_atoms; j++)
            for(int k = 0; k < 3; k++)
                if(argc == 4)
                    electric_field[i*pol_atoms*3+j*3+k] = external_ef[polar_mm[j]-1][k];
                else
                    electric_field[i*pol_atoms*3+j*3+k] = 0.0;

    do_qmmm(electric_field, OMMP_SOLVER_DEFAULT);

    get_energy(&E_MMMM, &E_MMPOL); 

    FILE *fp = fopen(argv[2], "w+");

    fprintf(fp, "%20.12e\n", E_MMMM);
    fprintf(fp, "%20.12e\n", E_MMPOL);
    
    fclose(fp);
    return 0;
}
