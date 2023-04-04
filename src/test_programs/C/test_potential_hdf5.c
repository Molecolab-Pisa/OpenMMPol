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
    if(argc != 3){
        printf("Syntax expected\n");
        printf("    $ test_init_xyz.exe <HDF5 FILE> <OUTPUT FILE>\n");
        return 0;
    }
    
    int pol_atoms;
    double eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, ebt, eat, etot, 
           ett, ev, er, edsp, ec, ecd, ed, em, ep, ect, erxf, es, elf, eg, ex;
    
    double *electric_field;

    OMMP_SYSTEM_PRT my_system = ommp_init_hdf5(argv[1], "system");
    
    pol_atoms = ommp_get_pol_atoms(my_system);
    
    electric_field = (double *) malloc(sizeof(double) * 3 * pol_atoms);
    
    for(int j = 0; j < pol_atoms; j++)
        for(int k = 0; k < 3; k++)
            electric_field[j*3+k] = 0.0;

    em = ommp_get_fixedelec_energy(my_system);
    ommp_set_external_field(my_system, electric_field, OMMP_SOLVER_DEFAULT);
    ep = ommp_get_polelec_energy(my_system);
    
    ev =   ommp_get_vdw_energy(my_system);
    eb =   ommp_get_bond_energy(my_system);
    ea =   ommp_get_angle_energy(my_system);
    eba =  ommp_get_strbnd_energy(my_system);
    eub =  ommp_get_urey_energy(my_system);
    eopb = ommp_get_opb_energy(my_system);
    ept =  ommp_get_pitors_energy(my_system);
    et =   ommp_get_torsion_energy(my_system);
    ett =  ommp_get_tortor_energy(my_system);
    eat =  ommp_get_angtor_energy(my_system);
    ebt =  ommp_get_strtor_energy(my_system);
    eit = ommp_get_imptorsion_energy(my_system);
    etot = ommp_get_full_energy(my_system);

    FILE *fp = fopen(argv[2], "w+");

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
    
    em *= AU2KCALMOL;
    ep *= AU2KCALMOL;
    ev *= AU2KCALMOL;
    eb *= AU2KCALMOL;
    ea *= AU2KCALMOL;
    eba *= AU2KCALMOL;
    eub *= AU2KCALMOL;
    eopb *= AU2KCALMOL;
    ept *= AU2KCALMOL;
    et *= AU2KCALMOL;
    ett *= AU2KCALMOL;
    eat *= AU2KCALMOL;
    ebt *= AU2KCALMOL;
    eit *= AU2KCALMOL;
    etot *= AU2KCALMOL;

    fprintf(fp, "EM      %20.12e\n", em);
    fprintf(fp, "EP      %20.12e\n", ep);
    fprintf(fp, "EV      %20.12e\n", ev);
    fprintf(fp, "EB      %20.12e\n", eb);
    fprintf(fp, "EA      %20.12e\n", ea);
    fprintf(fp, "EBA     %20.12e\n", eba);
    fprintf(fp, "EUB     %20.12e\n", eub);
    fprintf(fp, "EOPB    %20.12e\n", eopb);
    fprintf(fp, "EPT     %20.12e\n", ept);
    fprintf(fp, "ET      %20.12e\n", et);
    fprintf(fp, "ETT     %20.12e\n", ett);

    fprintf(fp, "EAA     %20.12e\n", eaa); 
    fprintf(fp, "EOPD    %20.12e\n", eopd);
    fprintf(fp, "EID     %20.12e\n", eid); 
    fprintf(fp, "EIT     %20.12e\n", eit); 
    fprintf(fp, "EBT     %20.12e\n", ebt); 
    fprintf(fp, "EAT     %20.12e\n", eat); 
    fprintf(fp, "ER      %20.12e\n", er);
    fprintf(fp, "EDSP    %20.12e\n", edsp);
    fprintf(fp, "EC      %20.12e\n", ec);
    fprintf(fp, "ECD     %20.12e\n", ecd);
    fprintf(fp, "ED      %20.12e\n", ed);
    fprintf(fp, "ECT     %20.12e\n", ect);
    fprintf(fp, "ERXF    %20.12e\n", erxf);
    fprintf(fp, "ES      %20.12e\n", es);
    fprintf(fp, "ELF     %20.12e\n", elf);
    fprintf(fp, "EG      %20.12e\n", eg);
    fprintf(fp, "EX      %20.12e\n", ex);
    fprintf(fp, "ETOT      %20.12e\n", etot);
    
    fclose(fp);
    free(electric_field);
    ommp_terminate(my_system);
    
    return 0;
}
