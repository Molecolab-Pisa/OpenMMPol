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
    if(argc != 5 && argc != 4){
        printf("Syntax expected\n");
        printf("    $ test_init_xyz.exe <XYZ FILE> <PRM FILE> <OUTPUT FILE> [<ELECTRIC FIELD>]\n");
        return 0;
    }
    
    int n_ipd, pol_atoms;
    double E_MMMM, E_MMPOL, E_VDW, E_BND;
    double *electric_field;
    double **external_ef;
    int32_t *polar_mm;

    set_verbose(OMMP_VERBOSE_DEBUG);
    mmpol_init_xyz(argv[1], argv[2]);
    
    n_ipd = get_n_ipd();
    pol_atoms = get_pol_atoms();
    
    electric_field = (double *) malloc(sizeof(double) * n_ipd * 3 * pol_atoms);
    polar_mm = (int32_t *) get_polar_mm();
    
    if(argc == 5)
        external_ef = read_ef(argv[4]);

    for(int i = 0; i < n_ipd; i++)
        for(int j = 0; j < pol_atoms; j++)
            for(int k = 0; k < 3; k++)
                if(argc == 5)
                    electric_field[i*pol_atoms*3+j*3+k] = external_ef[polar_mm[j]-1][k];
                else
                    electric_field[i*pol_atoms*3+j*3+k] = 0.0;

    do_qmmm(electric_field, OMMP_SOLVER_DEFAULT);

    get_energy(&E_MMMM, &E_MMPOL); 
    get_vdw_energy(&E_VDW);
    get_bond_energy(&E_BND);

    FILE *fp = fopen(argv[3], "w+");

    E_MMMM *= AU2KCALMOL;
    E_MMPOL *= AU2KCALMOL;
    E_VDW *= AU2KCALMOL;
    E_BND *= AU2KCALMOL;

    fprintf(fp, "MM-MM:  %20.12e\n", E_MMMM);
    fprintf(fp, "MM-Pol: %20.12e\n", E_MMPOL);
    fprintf(fp, "VDW:    %20.12e\n", E_VDW);
    fprintf(fp, "BOND:   %20.12e\n", E_BND);
    
    fclose(fp);
    
    return 0;
}
