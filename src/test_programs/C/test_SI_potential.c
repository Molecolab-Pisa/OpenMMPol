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
    
    OMMP_SYSTEM_PRT my_system, fake_qm;
    OMMP_QM_HELPER_PRT my_qmh;
    ommp_smartinput(argv[1], &my_system, &my_qmh);
    // override output file set in JSON file
    ommp_set_outputfile(argv[2]);
    
    bool use_qm = false, use_fake_qm = false, use_external_ef;
    use_external_ef = (argc == 4);
    
    if(my_qmh != NULL){
        // A QM part is present!
        use_qm = true;

        if(ommp_use_linkatoms(my_system)){
            // Since LA are used, I have to create a 
            // fake system for QM in order to have the
            // complete energy; parameters should be set 
            // anyway.
            use_fake_qm = true;
        }
    }

    if(use_fake_qm){
        // Cherrypick the parameter file use for QM from
        // the smart input
        char *prm_file, addr[] = "qm/prm_file/path";
        ommp_smartinput_cpstr(argv[1], addr, &prm_file);
        //Create the fake qm system
        fake_qm = ommp_system_from_qm_helper(my_qmh, prm_file);
        // In order to make a safe interaction, remove all the polarizabilities
        // on the fake qm system
        int natm = ommp_get_pol_atoms(fake_qm);
        int32_t *nopol = malloc(sizeof(int32_t) * natm);
        for(int j=0; j < natm; j++)
            nopol[j] = j;
        ommp_turn_pol_off(fake_qm, natm, nopol);
        free(nopol);
    }
    
    // If an external field is present (this is mutally exclusive with a QM part)
    int pol_atoms;
    double *electric_field, *qm_ef, **external_ef;

    int32_t *polar_mm = (int32_t *) ommp_get_polar_mm(my_system);
    pol_atoms = ommp_get_pol_atoms(my_system);
    
    electric_field = (double *) malloc(sizeof(double) * 3 * pol_atoms);
    if(use_fake_qm){
        qm_ef = (double *) malloc(sizeof(double) * 3 * pol_atoms);
        for(int j=0; j < 3 * pol_atoms; j++)
            qm_ef[j] = 0.0;
        ommp_field_mm2ext(fake_qm, pol_atoms, ommp_get_cpol(my_system), qm_ef);
    }

    if(use_external_ef){
        external_ef = read_ef(argv[3]);
        if(external_ef == NULL) return 1;
    }
    
    for(int j = 0; j < pol_atoms; j++)
        for(int k = 0; k < 3; k++){
            electric_field[j*3+k] = 0.0;
            if(use_external_ef)
                electric_field[j*3+k] += external_ef[polar_mm[j]][k];
            if(use_fake_qm)
                electric_field[j*3+k] += qm_ef[j*3+k];
        }
    
    em = ommp_get_fixedelec_energy(my_system);
    ommp_set_external_field(my_system, electric_field, OMMP_SOLVER_NONE, OMMP_MATV_NONE);
    ep = ommp_get_polelec_energy(my_system);
    
    ev = ommp_get_vdw_energy(my_system);
    if(use_qm)
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
    etotmm = ommp_get_full_energy(my_system);
    if(use_fake_qm)
        eqm = ommp_get_full_energy(fake_qm);
    else
        eqm = 0.0;
    etot = etotmm + eqm + evqmmm;

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
    
    sprintf(msg, "EM      %20.12e", em);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EP      %20.12e", ep);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EV      %20.12e", ev);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EVQMMM  %20.12e", evqmmm);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EB      %20.12e", eb);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EA      %20.12e", ea);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EBA     %20.12e", eba);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EUB     %20.12e", eub);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EOPB    %20.12e", eopb);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EPT     %20.12e", ept);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "ET      %20.12e", et);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "ETT     %20.12e", ett);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EAA     %20.12e", eaa); 
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EOPD    %20.12e", eopd);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EID     %20.12e", eid); 
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EIT     %20.12e", eit); 
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EBT     %20.12e", ebt); 
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EAT     %20.12e", eat); 
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "ER      %20.12e", er);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EDSP    %20.12e", edsp);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EC      %20.12e", ec);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "ECD     %20.12e", ecd);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "ED      %20.12e", ed);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "ECT     %20.12e", ect);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "ERXF    %20.12e", erxf);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "ES      %20.12e", es);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "ELF     %20.12e", elf);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EG      %20.12e", eg);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EX      %20.12e", ex);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "ETOTMM    %20.12e", etotmm);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "EQM       %20.12e", eqm);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    sprintf(msg, "ETOT      %20.12e", etot);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-ENE");
    
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
