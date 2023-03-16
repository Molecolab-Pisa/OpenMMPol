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

    dd = 1e-6;
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

int num_ana_compare(OMMP_SYSTEM_PRT sys,
                    double (*ene_f)(OMMP_SYSTEM_PRT), 
                    void (*grad_f)(OMMP_SYSTEM_PRT, double *), 
                    FILE *fp, const char *name, double del){
    
    int mm_atoms = ommp_get_mm_atoms(sys);
    
    double *_grad_num = (double *) malloc(sizeof(double) * 3 * mm_atoms);
    double **grad_num = (double **) malloc(sizeof(double *) * mm_atoms);
    double *_grad_ana = (double *) malloc(sizeof(double) * 3 * mm_atoms);
    double **grad_ana = (double **) malloc(sizeof(double *) * mm_atoms);

    for(int i = 0; i < mm_atoms; i++){
        grad_ana[i] = &(_grad_ana[i*3]);
        grad_num[i] = &(_grad_num[i*3]);
        for(int j=0; j < 3; j++)
            grad_ana[i][j] = grad_num[i][j] = 0.0;
    }

    numerical_geomgrad(sys, ene_f, grad_num);
    grad_f(sys, _grad_ana);
   
    fprintf(fp, "DELTA NUM - ANA %s\n", name);
    double Mdelta = 0.0, delta;
    
    for(int i = 0; i < mm_atoms; i++){
        fprintf(fp, "(A) ");
        for(int j=0; j < 3; j++){
            fprintf(fp, "%+12.8g ", grad_ana[i][j]);
        }
        fprintf(fp, "\n");
        fprintf(fp, "(N) ");
        for(int j=0; j < 3; j++){
            fprintf(fp, "%+12.8g ", grad_num[i][j]);
        }
        fprintf(fp, "\n");
        fprintf(fp, "(D) ");
        for(int j=0; j < 3; j++){
            delta = grad_num[i][j] - grad_ana[i][j];
            if(fabs(delta) > Mdelta) Mdelta = fabs(delta);
            fprintf(fp, "%+12.8g ", delta);
        }
        fprintf(fp, "\n\n");
    }


    if(Mdelta > del){
        fprintf(fp, "WARNING delta (%.3g) > max_delta (%.3g)\n", Mdelta, del);
        return 1;
    }else
        return 0;
}

int main(int argc, char **argv){
    if(argc != 4){
        printf("Syntax expected\n");
        printf("    $ test_geomgrad_xyz.exe <XYZ FILE> <PRM FILE> <OUTPUT FILE>\n");
        return 0;
    }
    
    int pol_atoms, mm_atoms, rc=0;
    double **grad_num, *_grad_num, **grad_ana, *_grad_ana;
    double delta, Mdelta;

    OMMP_SYSTEM_PRT my_system = ommp_init_xyz(argv[1], argv[2]);
    ommp_set_verbose(OMMP_VERBOSE_NONE);
    pol_atoms = ommp_get_pol_atoms(my_system);
    mm_atoms = ommp_get_mm_atoms(my_system);
    
    //electric_field = (double *) malloc(sizeof(double) * 3 * pol_atoms);
    //polar_mm = (int32_t *) ommp_get_polar_mm(my_system);

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
   
    FILE *fp = fopen(argv[3], "w+");
   
    /*rc += num_ana_compare(my_system, ommp_get_fixedelec_energy, ommp_fixedelec_geomgrad,
                          fp, "fixedelec", 1e-11);
    
    rc += num_ana_compare(my_system, ommp_get_polelec_energy, ommp_polelec_geomgrad,
                          fp, "polelec", 1e-8);
    
    rc += num_ana_compare(my_system, ommp_get_bond_energy, ommp_bond_geomgrad,
                          fp, "bond", 1e-11);

    rc += num_ana_compare(my_system, ommp_get_angle_energy, ommp_angle_geomgrad,
                          fp, "angle", 1e-11);

    rc += num_ana_compare(my_system, ommp_get_urey_energy, ommp_urey_geomgrad,
                          fp, "urey", 1e-11);
    
    rc += num_ana_compare(my_system, ommp_get_torsion_energy, ommp_torsion_geomgrad,
                          fp, "torsion", 1e-11);

    rc += num_ana_compare(my_system, ommp_get_imptorsion_energy, ommp_imptorsion_geomgrad,
                          fp, "imptorsion", 1e-11);
    
    rc += num_ana_compare(my_system, ommp_get_strbnd_energy, ommp_strbnd_geomgrad,
                          fp, "strbnd", 1e-11);
   
    rc += num_ana_compare(my_system, ommp_get_angtor_energy, ommp_angtor_geomgrad,
                          fp, "angtor", 1e-11);
    
    rc += num_ana_compare(my_system, ommp_get_opb_energy, ommp_opb_geomgrad,
                          fp, "opb", 1e-11);*/
    
    rc += num_ana_compare(my_system, ommp_get_strtor_energy, ommp_strtor_geomgrad,
                          fp, "strtor", 1e-10);
    fclose(fp);
    ommp_terminate(my_system);
     
    return rc;
}

