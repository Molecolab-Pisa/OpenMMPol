#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"
#define NELEMENTS 11

int symbol2z(char *s){
    const char elements[NELEMENTS][2] = {"X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"};
    for(int i=0; i < NELEMENTS; i++)
        if(strcmp(elements[i], s) == 0) return i;
    return -1;
}

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

void read_qm(char *fin, int32_t *nqm, double ***coord, int32_t **atype, 
             int32_t **z, double **q){
   
    *nqm = countLines(fin);
    char element[2];

    *coord = (double **) malloc(sizeof(double *) * 3 * (*nqm));
    double *_coord = (double *) malloc(sizeof(double) * 3 * (*nqm));
    *q = (double *) malloc(sizeof(double) * (*nqm));
    *z = (int32_t *) malloc(sizeof(int32_t) * (*nqm));
    *atype = (int32_t *) malloc(sizeof(int32_t) * (*nqm));

    FILE *fp = fopen(fin, "r");
    for(int i=0; i < (*nqm); i++){
        (*coord)[i] = &(_coord[i*3]);
        fscanf(fp, "%s %lf %lf %lf %d", element, &((*coord)[i][0]), 
                                        &((*coord)[i][1]),&((*coord)[i][2]),
                                        &((*atype)[i]));
        (*z)[i] = symbol2z(element);
        (*q)[i] = (double)(*z)[i];
        printf("%d\n", (*z)[i]);
    }
}

int main(int argc, char **argv){
    if(argc != 5){
        printf("Syntax expected\n");
        printf("    $ test_geomgrad_xyz.exe <XYZ-MM FILE> <XYZ-QM FILE> <PRM FILE> <OUTPUT FILE>\n");
        return 0;
    }
    
    ommp_set_verbose(OMMP_VERBOSE_DEBUG);
    OMMP_SYSTEM_PRT my_system = ommp_init_xyz(argv[1], argv[3]);

    double **cqm, *qqm;
    int32_t *zqm, *qmatype, nqm;
    read_qm(argv[2], &nqm, &cqm, &qmatype, &zqm, &qqm);
    for(int i=0; i <nqm; i++)
        printf("%d [%d %3.1f] %12.6f %12.6f %12.6f\n", i, zqm[i], qqm[i], cqm[i][0], cqm[i][1], cqm[i][2]);
    OMMP_QM_HELPER_PRT my_qmh = ommp_init_qm_helper(nqm, &(cqm[0]), qqm, zqm);

    ommp_qm_helper_set_attype(my_qmh, qmatype);
    ommp_qm_helper_init_vdw_prm(my_qmh, "/home/mattia/LibEnv/pyscf-openmmpol/examples/qmmm/amoeba09.prm");

    pol_atoms = ommp_get_pol_atoms(my_system);
    
    electric_field = (double *) malloc(sizeof(double) * 3 * pol_atoms);
    polar_mm = (int32_t *) ommp_get_polar_mm(my_system);
   
    for(int j = 0; j < pol_atoms; j++)
        for(int k = 0; k < 3; k++)
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
