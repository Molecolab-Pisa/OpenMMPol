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
    if(argc != 4){
        printf("Syntax expected\n");
        printf("    $ test_init_xyz.exe <XYZ FILE> <PRM FILE> <OUTPUT FILE>\n");
        return 0;
    }
    
    int n_ipd, pol_atoms;
    double eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, ebt, eat, 
           ett, ev, er, edsp, ec, ecd, ed, em, ep, ect, erxf, es, elf, eg, ex;
    
    double *electric_field;
    double **external_ef;

    set_verbose(OMMP_VERBOSE_DEBUG);
    mmpol_init_xyz(argv[1], argv[2]);
    
    n_ipd = get_n_ipd();
    pol_atoms = get_pol_atoms();
    
    electric_field = (double *) malloc(sizeof(double) * n_ipd * 3 * pol_atoms);
    
    for(int i = 0; i < n_ipd; i++)
        for(int j = 0; j < pol_atoms; j++)
            for(int k = 0; k < 3; k++)
                electric_field[i*pol_atoms*3+j*3+k] = 0.0;

    do_qmmm(electric_field, OMMP_SOLVER_DEFAULT);

    get_energy(&em, &ep); 
    get_vdw_energy(&ev);
    get_bond_energy(&eb);
    get_angle_energy(&ea);
    get_urey_energy(&eub);

    FILE *fp = fopen(argv[3], "w+");

    em *= AU2KCALMOL;
    ep *= AU2KCALMOL;
    ev *= AU2KCALMOL;
    eb *= AU2KCALMOL;
    ea *= AU2KCALMOL;
    eub *= AU2KCALMOL;

    fprintf(fp, "EM      %20.12e\n", em);
    fprintf(fp, "EP      %20.12e\n", ep);
    fprintf(fp, "EV      %20.12e\n", ev);
    fprintf(fp, "EB      %20.12e\n", eb);
    fprintf(fp, "EA      %20.12e\n", ea);
    fprintf(fp, "EUB     %20.12e\n", eub);
    
    fclose(fp);
    
    return 0;
}
