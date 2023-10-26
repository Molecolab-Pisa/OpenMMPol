#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
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
    if(fp == NULL){
        free(ef);
        return NULL;
    }
    for(int i =0; i < pol_atoms; i++){
        ef[i] = (double *) malloc(sizeof(double) * 3);
        fscanf(fp, "%lf %lf %lf", &(ef[i][0]), &(ef[i][1]),  &(ef[i][2]));
    }
    fclose(fp);

    return ef;
}

int main(int argc, char **argv){
    if(argc != 3 && argc != 4){
        printf("Syntax expected\n");
        printf("    $ test_SI_potential.exe <JSON FILE> <OUTPUT FILE> [<EF_FILE>]\n");
        return 1;
    }
    
    double eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, ebt, eat, etot,
           ett, ev, er, edsp, ec, ecd, ed, em, ep, ect, erxf, es, elf, eg, ex, 
           evqmmm, etotmm, eqm;

    char msg[OMMP_STR_CHAR_MAX];
    
    ommp_time_push();
    ommp_time_push();

    OMMP_SYSTEM_PRT my_system, fake_qm;
    OMMP_QM_HELPER_PRT my_qmh;
    ommp_smartinput(argv[1], &my_system, &my_qmh);
    // override output file set in JSON file
    ommp_set_outputfile(argv[2]);
    
    bool use_external_ef;
    use_external_ef = false;
    
    if(my_qmh != NULL){
        // A QM part is present!
        printf("Timing test only works with MM part");
    }

    // If an external field is present (this is mutally exclusive with a QM part)
    int pol_atoms, mm_atoms;
    double *electric_field, *qm_ef, **external_ef;

    int32_t *polar_mm = (int32_t *) ommp_get_polar_mm(my_system);
    pol_atoms = ommp_get_pol_atoms(my_system);
    mm_atoms = ommp_get_mm_atoms(my_system);
    
    electric_field = (double *) malloc(sizeof(double) * 3 * pol_atoms);
    double *grd = (double *) malloc(sizeof(double) * 3 * mm_atoms);

    if(use_external_ef){
        external_ef = read_ef(argv[3]);
        if(external_ef == NULL) return 1;
    }
    
    for(int j = 0; j < pol_atoms; j++)
        for(int k = 0; k < 3; k++){
            electric_field[j*3+k] = 0.0;
            if(use_external_ef)
                electric_field[j*3+k] += external_ef[polar_mm[j]][k];
        }
    
    ommp_time_pull("Initialization");
    em = ommp_get_fixedelec_energy(my_system);
    ommp_set_external_field(my_system, electric_field, OMMP_SOLVER_NONE, OMMP_MATV_NONE);
    ep = ommp_get_polelec_energy(my_system);
    

    ev = ommp_get_vdw_energy(my_system);
    ommp_time_push();
    eb = ommp_get_bond_energy(my_system);
    ea = ommp_get_angle_energy(my_system);
    eba = ommp_get_strbnd_energy(my_system);
    eub = ommp_get_urey_energy(my_system);
    eopb = ommp_get_opb_energy(my_system);
    ept = ommp_get_pitors_energy(my_system);
    et = ommp_get_torsion_energy(my_system);
    ett = ommp_get_tortor_energy(my_system);
    eat = ommp_get_angtor_energy(my_system);
    ebt = ommp_get_strtor_energy(my_system);
    eit = ommp_get_imptorsion_energy(my_system);
    ommp_time_pull("Total energy bonded");

    ommp_fixedelec_geomgrad(my_system, grd);
    ommp_polelec_geomgrad(my_system, grd);

    ommp_vdw_geomgrad(my_system, grd);
    ommp_time_push();
    ommp_bond_geomgrad(my_system, grd);
    ommp_angle_geomgrad(my_system, grd);
    ommp_strbnd_geomgrad(my_system, grd);
    ommp_urey_geomgrad(my_system, grd);
    ommp_torsion_geomgrad(my_system, grd);
    ommp_imptorsion_geomgrad(my_system, grd);
    ommp_angtor_geomgrad(my_system, grd);
    ommp_opb_geomgrad(my_system, grd);
    ommp_strtor_geomgrad(my_system, grd);
    ommp_tortor_geomgrad(my_system, grd);
    ommp_pitors_geomgrad(my_system, grd);

    ommp_time_pull("Total grad bonded");

    free(electric_field);
    if(my_qmh != NULL) ommp_terminate_qm_helper(my_qmh);
    ommp_terminate(my_system);
    
    ommp_time_pull("Total exec");
    return 0;
}
