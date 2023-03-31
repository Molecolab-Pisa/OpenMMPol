#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if(argc != 3){
        printf("Syntax expected\n");
        printf("    $ test_geomgrad.exe <INPUT FILE> <OUTPUT FILE>\n");
        return 0;
    }
    
    int pol_atoms, mm_atoms, retcode=0;
    double **grad_ana, *_grad_ana;

    OMMP_SYSTEM_PRT my_system = ommp_init_mmp(argv[1]);
    ommp_set_verbose(OMMP_VERBOSE_NONE);
    
    pol_atoms = ommp_get_pol_atoms(my_system);
    mm_atoms = ommp_get_mm_atoms(my_system);
    
    _grad_ana = (double *) malloc(sizeof(double) * 3 * mm_atoms);
    grad_ana = (double **) malloc(sizeof(double *) * mm_atoms);
    for(int i = 0; i < mm_atoms; i++){
        grad_ana[i] = &(_grad_ana[i*3]);
        for(int j=0; j < 3; j++)
            grad_ana[i][j] = 0.0;
    }
   
    FILE *fp = fopen(argv[2], "w+");
    
    ommp_fixedelec_geomgrad(my_system, _grad_ana);
    
    fprintf(fp, "Grad EM\n");
    for(int i = 0; i < mm_atoms; i++){
        for(int j=0; j < 3; j++)
            fprintf(fp, "%+12.8g ", grad_ana[i][j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    for(int i = 0; i < mm_atoms; i++)
        for(int j=0; j < 3; j++)
            grad_ana[i][j] = 0.0;
    
    ommp_polelec_geomgrad(my_system, _grad_ana);
   
    fprintf(fp, "Grad EP\n");
    for(int i = 0; i < mm_atoms; i++){
        for(int j=0; j < 3; j++)
            fprintf(fp, "%+12.8g ", grad_ana[i][j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    fclose(fp);
    ommp_terminate(my_system);
    
    return retcode;
}

