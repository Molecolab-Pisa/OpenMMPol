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
    
    int pol_atoms;
    double E_MMMM, E_MMPOL;
    double *electric_field;
    double **external_ef;
    int32_t *polar_mm;

    OMMP_SYSTEM_PRT my_system = ommp_init_mmp(argv[1]);
    pol_atoms = ommp_get_pol_atoms(my_system);
    
    electric_field = (double *) malloc(sizeof(double) * 3 * pol_atoms);
    polar_mm = (int32_t *) ommp_get_polar_mm(my_system);
   
    if(argc == 4)
        external_ef = read_ef(argv[3]);

    for(int j = 0; j < pol_atoms; j++)
        for(int k = 0; k < 3; k++)
            if(argc == 4)
                electric_field[j*3+k] = external_ef[polar_mm[j]][k];
            else
                electric_field[j*3+k] = 0.0;
    
    E_MMMM = ommp_get_fixedelec_energy(my_system);
    ommp_set_external_field(my_system, electric_field, OMMP_SOLVER_DEFAULT);
    E_MMPOL = ommp_get_polelec_energy(my_system);

    FILE *fp = fopen(argv[2], "w+");

    fprintf(fp, "%20.12e\n", E_MMMM);
    fprintf(fp, "%20.12e\n", E_MMPOL);
    
    fclose(fp);
    free(electric_field);
    ommp_terminate(my_system);
    
    return 0;
}
