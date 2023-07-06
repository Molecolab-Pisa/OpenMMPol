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

double full_qmmm_energy(OMMP_SYSTEM_PRT qm_sys, OMMP_QM_HELPER_PRT qmh, OMMP_SYSTEM_PRT mm_sys){
    double ene = 0.0;

    ene += ommp_get_full_energy(mm_sys);
    ene += ommp_get_full_energy(qm_sys);
    ene += ommp_qm_helper_vdw_energy(qmh, mm_sys);

    return ene;
}

void update_qmmm_coordinates(OMMP_SYSTEM_PRT qm_sys, OMMP_QM_HELPER_PRT qmh, OMMP_SYSTEM_PRT mm_sys, double *new_mm_c, double *new_qm_c){
    ommp_update_coordinates(mm_sys, new_mm_c);
    ommp_qm_helper_update_coord(qmh, new_qm_c);
    ommp_update_link_atoms_position(qmh, mm_sys);
    ommp_update_coordinates(qm_sys, ommp_qm_helper_get_cqm(qmh));
}

void numerical_geomgrad(OMMP_SYSTEM_PRT qm_sys, OMMP_QM_HELPER_PRT qmh, OMMP_SYSTEM_PRT mm_sys, double *gqm, double *gmm){

    double *new_c_mm, *new_c_qm, *cmm, *cqm, tmp, dd;
    int nmm, nqm;

    dd = 1e-6;
    nmm = ommp_get_mm_atoms(mm_sys);
    nqm = ommp_get_mm_atoms(qm_sys);
    cmm = ommp_get_cmm(mm_sys);
    cqm = ommp_get_cmm(qm_sys);

    new_c_mm = (double *) malloc(sizeof(double) * 3 * nmm);
    for(int i=0; i < nmm; i++)
        for(int j=0; j < 3; j++) 
            new_c_mm[i*3+j] = cmm[i*3+j];
    
    new_c_qm = (double *) malloc(sizeof(double) * 3 * nqm);
    for(int i=0; i < nqm; i++)
        for(int j=0; j < 3; j++) 
            new_c_qm[i*3+j] = cqm[i*3+j];

    for(int i=0; i < nmm; i++){
        for(int j=0; j < 3; j++){
            new_c_mm[i*3+j] += dd;
            update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
            tmp = full_qmmm_energy(qm_sys, qmh, mm_sys);

            new_c_mm[i*3+j] -= 2*dd;
            update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
            tmp -= full_qmmm_energy(qm_sys, qmh, mm_sys);
            gmm[i*3+j] = tmp / (2*dd);
            
            new_c_mm[i*3+j] += dd;
            update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
        }
    }
 
    for(int i=0; i < nqm; i++){
        for(int j=0; j < 3; j++){
            new_c_qm[i*3+j] += dd;
            update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
            tmp = full_qmmm_energy(qm_sys, qmh, mm_sys);

            new_c_qm[i*3+j] -= 2*dd;
            update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
            tmp -= full_qmmm_energy(qm_sys, qmh, mm_sys);
            gqm[i*3+j] = tmp / (2*dd);
            
            new_c_qm[i*3+j] += dd;
            update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
        }
    }
}

int main(int argc, char **argv){
    if(argc != 6){
        printf("Syntax expected\n");
        printf("    $ test_geomgrad_xyz.exe <XYZ-MM FILE> <XYZ-QM FILE> <PRM FILE> <LA-FILE> <OUTPUT FILE>\n");
        return 0;
    }
    
    double *qqm, *electric_field;
    int32_t *zqm, *qmatype, nqm, nmm, pol_atoms;
    double eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, ebt, eat, etot,
           ett, ev, er, edsp, ec, ecd, ed, em, ep, ect, erxf, es, elf, eg, ex,
           evqmmm;

    // Create OMMP System for MM subsystem
    ommp_set_verbose(OMMP_VERBOSE_DEBUG);
    OMMP_SYSTEM_PRT my_system = ommp_init_xyz(argv[1], argv[3]);
    nmm = ommp_get_mm_atoms(my_system);

    // Now create a OMMP system for QM subsystem to simulate QM gradients
    OMMP_SYSTEM_PRT qm_sys = ommp_init_xyz(argv[2], argv[3]);

    nqm = ommp_get_mm_atoms(qm_sys);
    zqm = ommp_get_zmm(qm_sys);
    qmatype = ommp_get_attypemm(qm_sys);
    qqm = (double *) malloc(sizeof(double) * nqm);
    for(int i=0; i<nqm; i++)
        qqm[i] = (double) zqm[i];

    OMMP_QM_HELPER_PRT my_qmh = ommp_init_qm_helper(nqm,
                                                    ommp_get_cmm(qm_sys),
                                                    qqm,
                                                    zqm);

    ommp_qm_helper_set_attype(my_qmh, qmatype);
    ommp_qm_helper_init_vdw_prm(my_qmh, argv[3]);
    
    int nla = countLines(argv[4]);
    FILE *fp = fopen(argv[4], "r");
    int32_t imm, ila, iqm;
    for(int i = 0; i < nla; i++){
        fscanf(fp, "%d %d %d", &imm, &iqm, &ila);

        ommp_create_link_atom(my_qmh, my_system, imm, iqm, ila, 
                            argv[3], OMMP_DEFAULT_LA_DIST, 
                            OMMP_DEFAULT_LA_N_EEL_REMOVE);
    }

    fclose(fp);
    ommp_update_coordinates(qm_sys, ommp_qm_helper_get_cqm(my_qmh));

    pol_atoms = ommp_get_pol_atoms(my_system);
    
    electric_field = (double *) malloc(sizeof(double) * 3 * pol_atoms);
   
    for(int j = 0; j < pol_atoms; j++)
        for(int k = 0; k < 3; k++)
            electric_field[j*3+k] = 0.0;
    
    ommp_set_external_field(my_system, electric_field, OMMP_SOLVER_DEFAULT);
    
    // Fake QM gradients, they are actually the gradients of QM subsystem
    // trated at MM level, but in this context we do not actually care
    // about where the gradients do come from.
    double *gradqm = malloc(sizeof(double) * 3 * nqm);
    ommp_full_geomgrad(qm_sys, gradqm);

    double *gradmm = malloc(sizeof(double) * 3 * nmm);
    ommp_full_geomgrad(my_system, gradmm);

    double *tmpmm = malloc(sizeof(double) * 3 * nmm);
    double *tmpqm = malloc(sizeof(double) * 3 * nqm);

    ommp_qm_helper_vdw_geomgrad(my_qmh, my_system, tmpqm, tmpmm);
    for(int i = 0; i < 3 * nmm; i++)
        gradmm[i] += tmpmm[i];
    for(int i = 0; i < 3 * nqm; i++)
        gradqm[i] += tmpqm[i];
    
    ommp_qm_helper_link_atom_geomgrad(my_qmh, my_system, tmpqm, tmpmm, gradqm);
    for(int i = 0; i < 3 * nmm; i++)
        gradmm[i] += tmpmm[i];
    for(int i = 0; i < 3 * nqm; i++)
        gradqm[i] += tmpqm[i];

    // Compute numerical gradients
    double *gradqm_num = malloc(sizeof(double) * 3 * nqm);
    double *gradmm_num = malloc(sizeof(double) * 3 * nmm);
    numerical_geomgrad(qm_sys, my_qmh, my_system, gradqm_num, gradmm_num);

    double Mdelta = 0.0, delta, del = 1e-7;

    fp = fopen(argv[5], "w+");

    fprintf(fp, "DELTA ANA - NUM (MM) \n");
    for(int i = 0; i < nmm; i++){
        fprintf(fp, "[%5d] (A) ", i+1);
        for(int j=0; j < 3; j++){
            fprintf(fp, "%+12.8g ", gradmm[i*3+j]);
        }
        fprintf(fp, "\n");
        fprintf(fp, "        (N) ");
        for(int j=0; j < 3; j++){
            fprintf(fp, "%+12.8g ", gradmm_num[i*3+j]);
        }
        fprintf(fp, "\n");
        fprintf(fp, "        (D) ");
        for(int j=0; j < 3; j++){
            delta = gradmm_num[i*3+j] - gradmm[i*3+j];
            if(fabs(delta) > Mdelta) Mdelta = fabs(delta);
            fprintf(fp, "%+12.8g ", delta);
        }
        fprintf(fp, "\n\n");
    }
    
    fprintf(fp, "DELTA ANA - NUM (QM) \n");
    for(int i = 0; i < nqm; i++){
        fprintf(fp, "[%5d] (A) ", i+1);
        for(int j=0; j < 3; j++){
            fprintf(fp, "%+12.8g ", gradqm[i*3+j]);
        }
        fprintf(fp, "\n");
        fprintf(fp, "        (N) ");
        for(int j=0; j < 3; j++){
            fprintf(fp, "%+12.8g ", gradqm_num[i*3+j]);
        }
        fprintf(fp, "\n");
        fprintf(fp, "        (D) ");
        for(int j=0; j < 3; j++){
            delta = gradqm_num[i*3+j] - gradqm[i*3+j];
            if(fabs(delta) > Mdelta) Mdelta = fabs(delta);
            fprintf(fp, "%+12.8g ", delta);
        }
        fprintf(fp, "\n\n");
    }

    fclose(fp);

    free(electric_field);
    free(gradqm_num);
    free(gradmm_num);
    free(gradqm);
    free(gradmm);
    free(tmpmm);
    free(tmpqm);
    free(qqm);
    ommp_terminate_qm_helper(my_qmh);
    ommp_terminate(qm_sys);
    ommp_terminate(my_system);
    
    if(Mdelta > del){
        fprintf(fp, "WARNING delta (%.3g) > max_delta (%.3g)\n", Mdelta, del);
        return 1;
    }else
        return 0;
}
