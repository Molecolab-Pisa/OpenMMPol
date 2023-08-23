#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

void ana_grd_print(OMMP_SYSTEM_PRT sys,
                  void (*grad_f)(OMMP_SYSTEM_PRT, double *), 
                  const char *name){
    
    int mm_atoms = ommp_get_mm_atoms(sys);
    
    double *_grad_ana = (double *) malloc(sizeof(double) * 3 * mm_atoms);
    double **grad_ana = (double **) malloc(sizeof(double *) * mm_atoms);
    char msg[OMMP_STR_CHAR_MAX];

    for(int i = 0; i < mm_atoms; i++){
        grad_ana[i] = &(_grad_ana[i*3]);
        for(int j=0; j < 3; j++)
            grad_ana[i][j] = 0.0;
    }

    grad_f(sys, _grad_ana);
   
    sprintf(msg, "Grad %s", name);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    
    for(int i = 0; i < mm_atoms; i++){
        sprintf(msg, "%+20.8g %+20.8g %+20.8g", 
                grad_ana[i][0]*OMMP_AU2KCALMOL*OMMP_ANG2AU,
                grad_ana[i][1]*OMMP_AU2KCALMOL*OMMP_ANG2AU,
                grad_ana[i][2]*OMMP_AU2KCALMOL*OMMP_ANG2AU);
        ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    }
    ommp_message("", OMMP_VERBOSE_NONE, "TEST-OUT");
}

int main(int argc, char **argv){
    if(argc != 2){
        printf("Syntax expected\n");
        printf("    $ test_SI_geomgrad.exe <JSON FILE>\n");
        return 0;
    }
    
    OMMP_SYSTEM_PRT my_system;
    OMMP_QM_HELPER_PRT my_qmh;
    ommp_smartinput(argv[1], &my_system, &my_qmh);
    
    ana_grd_print(my_system, ommp_full_geomgrad, "ETOT");
    ana_grd_print(my_system, ommp_vdw_geomgrad, "EV");
    ana_grd_print(my_system, ommp_fixedelec_geomgrad, "EM");
    ana_grd_print(my_system, ommp_polelec_geomgrad, "EP");
    ana_grd_print(my_system, ommp_bond_geomgrad, "EB");
    ana_grd_print(my_system, ommp_angle_geomgrad, "EA");
    ana_grd_print(my_system, ommp_urey_geomgrad, "EUB");
    ana_grd_print(my_system, ommp_torsion_geomgrad, "ET");
    ana_grd_print(my_system, ommp_imptorsion_geomgrad, "EIT");
    ana_grd_print(my_system, ommp_strbnd_geomgrad, "EBA");
    ana_grd_print(my_system, ommp_angtor_geomgrad, "EAT");
    ana_grd_print(my_system, ommp_opb_geomgrad, "EOPB");
    ana_grd_print(my_system, ommp_strtor_geomgrad, "EBT");
    ana_grd_print(my_system, ommp_tortor_geomgrad, "ETT");
    ana_grd_print(my_system, ommp_pitors_geomgrad, "EPT");
    
    if(my_qmh != NULL) ommp_terminate_qm_helper(my_qmh);
    ommp_terminate(my_system);
    
    return 0;
}
