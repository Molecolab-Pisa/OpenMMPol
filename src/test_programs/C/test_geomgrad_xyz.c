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

int ana_grd_print(OMMP_SYSTEM_PRT sys,
                  void (*grad_f)(OMMP_SYSTEM_PRT, double *), 
                  FILE *fp, const char *name){
    
    int mm_atoms = ommp_get_mm_atoms(sys);
    
    double *_grad_ana = (double *) malloc(sizeof(double) * 3 * mm_atoms);
    double **grad_ana = (double **) malloc(sizeof(double *) * mm_atoms);

    for(int i = 0; i < mm_atoms; i++){
        grad_ana[i] = &(_grad_ana[i*3]);
        for(int j=0; j < 3; j++)
            grad_ana[i][j] = 0.0;
    }

    grad_f(sys, _grad_ana);
   
    fprintf(fp, "Grad %s\n", name);
    
    for(int i = 0; i < mm_atoms; i++){
        for(int j=0; j < 3; j++){
            fprintf(fp, "%+12.8g ", grad_ana[i][j]*AU2KCALMOL*ANG2AU);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
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
   
    //ana_grd_print(my_system, ommp_full_geomgrad, fp, "ETOT");
    
    ana_grd_print(my_system, ommp_vdw_geomgrad, fp, "EV");
    ana_grd_print(my_system, ommp_fixedelec_geomgrad, fp, "EM");
    ana_grd_print(my_system, ommp_polelec_geomgrad, fp, "EP");
    ana_grd_print(my_system, ommp_bond_geomgrad, fp, "EB");
    ana_grd_print(my_system, ommp_angle_geomgrad, fp, "EA");
    ana_grd_print(my_system, ommp_urey_geomgrad, fp, "EUB");
    ana_grd_print(my_system, ommp_torsion_geomgrad, fp, "ET");
    ana_grd_print(my_system, ommp_imptorsion_geomgrad, fp, "EIT");
    ana_grd_print(my_system, ommp_strbnd_geomgrad, fp, "EBA");
    ana_grd_print(my_system, ommp_angtor_geomgrad, fp, "EAT");
    ana_grd_print(my_system, ommp_opb_geomgrad, fp, "EOPB");
    ana_grd_print(my_system, ommp_strtor_geomgrad, fp, "EBT");
    ana_grd_print(my_system, ommp_tortor_geomgrad, fp, "ETT");
    ana_grd_print(my_system, ommp_pitors_geomgrad, fp, "EPT");
    
    fclose(fp);
    ommp_terminate(my_system);
     
    return rc;
}

