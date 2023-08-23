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
    if(argc != 2 && argc != 3){
        printf("Syntax expected\n");
        printf("    $ test_init_xyz.exe <JSON FILE> [<EF_FILE>]\n");
        return 0;
    }
    
    OMMP_SYSTEM_PRT my_system;
    OMMP_QM_HELPER_PRT my_qmh;
    ommp_smartinput(argv[1], &my_system, &my_qmh);
    
    int pol_atoms;
    double eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, ebt, eat, etot,
           ett, ev, er, edsp, ec, ecd, ed, em, ep, ect, erxf, es, elf, eg, ex, evqmmm;
    
    double *electric_field, **external_ef;
    char msg[OMMP_STR_CHAR_MAX];

    int32_t *polar_mm = (int32_t *) ommp_get_polar_mm(my_system);
    pol_atoms = ommp_get_pol_atoms(my_system);
    
    electric_field = (double *) malloc(sizeof(double) * 3 * pol_atoms);
    
    if(argc == 3)
        external_ef = read_ef(argv[2]);
    
    for(int j = 0; j < pol_atoms; j++)
        for(int k = 0; k < 3; k++)
            if(argc == 3)
                electric_field[j*3+k] = external_ef[polar_mm[j]][k];
            else
                electric_field[j*3+k] = 0.0;
    
    em = ommp_get_fixedelec_energy(my_system);
    ommp_set_external_field(my_system, electric_field, OMMP_SOLVER_NONE, OMMP_MATV_NONE);
    ep = ommp_get_polelec_energy(my_system);
    
    ev = ommp_get_vdw_energy(my_system);
    
    if(my_qmh != NULL)
        evqmmm = ommp_qm_helper_vdw_energy(my_qmh, my_system);
    else
        evqmmm = 0.0;
    
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
    etot = ommp_get_full_energy(my_system);

    eaa = 0.0;
    eopd = 0.0;
    eid = 0.0;  
    er = 0.0;
    edsp = 0.0;
    ec = 0.0;
    ecd = 0.0;
    ed = 0.0;
    ect = 0.0;
    erxf = 0.0;
    es = 0.0;
    elf = 0.0;
    eg = 0.0;
    ex = 0.0;
    
    em *= OMMP_AU2KCALMOL;
    ep *= OMMP_AU2KCALMOL;
    ev *= OMMP_AU2KCALMOL;
    eb *= OMMP_AU2KCALMOL;
    ea *= OMMP_AU2KCALMOL;
    eba *= OMMP_AU2KCALMOL;
    eub *= OMMP_AU2KCALMOL;
    eopb *= OMMP_AU2KCALMOL;
    ept *= OMMP_AU2KCALMOL;
    et *= OMMP_AU2KCALMOL;
    ett *= OMMP_AU2KCALMOL;
    eat *= OMMP_AU2KCALMOL;
    ebt *= OMMP_AU2KCALMOL;
    eit *= OMMP_AU2KCALMOL;
    etot *= OMMP_AU2KCALMOL;

    sprintf(msg, "EM      %20.12e", em);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EP      %20.12e", ep);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EV      %20.12e", ev);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EVQMMM  %20.12e", evqmmm);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EB      %20.12e", eb);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EA      %20.12e", ea);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EBA     %20.12e", eba);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EUB     %20.12e", eub);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EOPB    %20.12e", eopb);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EPT     %20.12e", ept);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "ET      %20.12e", et);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "ETT     %20.12e", ett);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EAA     %20.12e", eaa); 
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EOPD    %20.12e", eopd);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EID     %20.12e", eid); 
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EIT     %20.12e", eit); 
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EBT     %20.12e", ebt); 
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EAT     %20.12e", eat); 
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "ER      %20.12e", er);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EDSP    %20.12e", edsp);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EC      %20.12e", ec);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "ECD     %20.12e", ecd);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "ED      %20.12e", ed);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "ECT     %20.12e", ect);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "ERXF    %20.12e", erxf);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "ES      %20.12e", es);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "ELF     %20.12e", elf);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EG      %20.12e", eg);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "EX      %20.12e", ex);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    sprintf(msg, "ETOT      %20.12e", etot);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    
    // Print also induced point dipoles
    int n_ipd = ommp_get_n_ipd(my_system);
    double *_ipd = (double *) ommp_get_ipd(my_system);
    double ***ipd = (double ***) malloc(sizeof(double **) * n_ipd );

    for(int i = 0; i < n_ipd; i++){
        ipd[i] = (double **) malloc(sizeof(double *) * pol_atoms);

        for(int j = 0; j < pol_atoms; j++){
            ipd[i][j] = &(_ipd[i*pol_atoms*3 + j*3]);
        }
    }

    for(int k = 0; k < n_ipd; k++){
        for(int i = 0; i < pol_atoms; i++){
            sprintf(msg, "%20.12e %20.12e %20.12e",
                    ipd[k][i][0], ipd[k][i][1], ipd[k][i][2]);
            ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-IPD");
        }
    }
    
    for(int i = 0; i < n_ipd; i++)
        free(ipd[i]);
    free(ipd);
    free(electric_field);
    if(my_qmh != NULL) ommp_terminate_qm_helper(my_qmh);
    ommp_terminate(my_system);
    
    return 0;
}
