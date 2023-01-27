#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

void numerical_geomgrad(OMMP_SYSTEM_PRT s, double (*ene_f)(OMMP_SYSTEM_PRT), double **g){
    double *_new_c, **new_c, *cmm, tmp, dd;
    int mm_atoms;

    dd = 1e-5;
    mm_atoms = ommp_get_mm_atoms(s);
    cmm = ommp_get_cmm(s);

    _new_c = (double *) malloc(sizeof(double) * 3 * mm_atoms);
    new_c = (double **) malloc(sizeof(double *) * mm_atoms);
    for(int i=0; i < mm_atoms; i++){
        new_c[i] = &(_new_c[i*3]);
        for(int j=0; j < 3; j++) 
            new_c[i][j] = cmm[i*3+j];
    }

    for(int i=0; i < mm_atoms; i++){
        for(int j=0; j < 3; j++){
            new_c[i][j] += dd;
            ommp_update_coordinates(s, _new_c);
            tmp = ene_f(s);

            new_c[i][j] -= 2*dd;
            ommp_update_coordinates(s, _new_c);
            tmp -= ene_f(s);
            g[i][j] = tmp / (2*dd);
            
            new_c[i][j] += dd;
            ommp_update_coordinates(s, _new_c);
        }
    }
}

int main(int argc, char **argv){
    if(argc != 4 && argc != 3){
        printf("Syntax expected\n");
        printf("    $ test_init.exe <INPUT FILE> <OUTPUT FILE> [<ELECTRIC FIELD FILE>]\n");
        return 0;
    }
    
    int pol_atoms, mm_atoms, retcode=0;
    double *electric_field;
    double **external_ef, **grad_num, *_grad_num, **grad_ana, *_grad_ana;
    double delta, Mdelta;
    int32_t *polar_mm;

    OMMP_SYSTEM_PRT my_system = ommp_init_mmp(argv[1]);
    ommp_set_verbose(OMMP_VERBOSE_NONE);
    pol_atoms = ommp_get_pol_atoms(my_system);
    mm_atoms = ommp_get_mm_atoms(my_system);
    
    electric_field = (double *) malloc(sizeof(double) * 3 * pol_atoms);
    polar_mm = (int32_t *) ommp_get_polar_mm(my_system);

    _grad_num = (double *) malloc(sizeof(double) * 3 * mm_atoms);
    grad_num = (double **) malloc(sizeof(double *) * mm_atoms);
    _grad_ana = (double *) malloc(sizeof(double) * 3 * mm_atoms);
    grad_ana = (double **) malloc(sizeof(double *) * mm_atoms);
    for(int i = 0; i < mm_atoms; i++){
        grad_ana[i] = &(_grad_ana[i*3]);
        grad_num[i] = &(_grad_num[i*3]);
        for(int j=0; j < 3; j++)
            grad_ana[i][j] = grad_num[i][j] = 0.0;
    }
   
    if(argc == 4){
        printf("Currently unsupported");
        return 0;
        external_ef = read_ef(argv[3]);
    }

    for(int j = 0; j < pol_atoms; j++)
        for(int k = 0; k < 3; k++)
            if(argc == 4)
                electric_field[j*3+k] = external_ef[polar_mm[j]][k];
            else
                electric_field[j*3+k] = 0.0;
    
    FILE *fp = fopen(argv[2], "w+");
    
    numerical_geomgrad(my_system, ommp_get_fixedelec_energy, grad_num);
    ommp_fixedelec_geomgrad(my_system, _grad_ana);
    
    fprintf(fp, "DELTA NUM - ANA FIXEDELEC\n");
    Mdelta = 0.0;
    for(int i = 0; i < mm_atoms; i++){
        for(int j=0; j < 3; j++){
            delta = grad_num[i][j] - grad_ana[i][j];
            if(fabs(delta) > Mdelta) Mdelta = fabs(delta);
            fprintf(fp, "%+12.8g ", delta);
        }
        fprintf(fp, "\n");
    }

    if(Mdelta > 1e-11){
        fprintf(fp, "Numerical-Analytical gradients difference is too large (fixed).");
        retcode = retcode + 1;
    }
    
    for(int i = 0; i < mm_atoms; i++)
        for(int j=0; j < 3; j++)
            grad_ana[i][j] = grad_num[i][j] = 0.0;
    
    numerical_geomgrad(my_system, ommp_get_polelec_energy, grad_num);
    ommp_polelec_geomgrad(my_system, _grad_ana);
   
    fprintf(fp, "DELTA NUM - ANA POLELEC\n");
    Mdelta = 0.0;
    for(int i = 0; i < mm_atoms; i++){
        for(int j=0; j < 3; j++){
            delta = grad_num[i][j] - grad_ana[i][j];
            if(fabs(delta) > Mdelta) Mdelta = fabs(delta);
            fprintf(fp, "%+12.8g ", delta);
        }
        fprintf(fp, "\n");
    }

    if(Mdelta > 1e-8){
        fprintf(fp, "Numerical-Analytical gradients difference is too large (pol).");
        retcode = retcode + 2;
    }

    fclose(fp);
    free(electric_field);
    ommp_terminate(my_system);
    
    return retcode;
}

