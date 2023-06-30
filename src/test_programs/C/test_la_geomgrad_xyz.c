#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

int main(int argc, char **argv){
    if(argc != 6){
        printf("Syntax expected\n");
        printf("    $ test_geomgrad_xyz.exe <XYZ-MM FILE> <XYZ-QM FILE> <PRM FILE> <LA-FILE> <OUTPUT FILE>\n");
        return 0;
    }
    
    ommp_set_verbose(OMMP_VERBOSE_DEBUG);
    OMMP_SYSTEM_PRT my_system = ommp_init_xyz(argv[1], argv[3]);

    double *qqm, *electric_field;
    int32_t *zqm, *qmatype, nqm, pol_atoms;
    double eb, ea, eba, eub, eaa, eopb, eopd, eid, eit, et, ept, ebt, eat, etot,
           ett, ev, er, edsp, ec, ecd, ed, em, ep, ect, erxf, es, elf, eg, ex,
           evqmmm;
    
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
    
    FILE *fp = fopen(argv[4], "r");
    int32_t imm, ila, iqm;
    fscanf(fp, "%d %d %d", &imm, &iqm, &ila);
    fclose(fp);

    ommp_create_link_atom(my_qmh, my_system, imm, iqm, ila, 
                          argv[3], OMMP_DEFAULT_LA_DIST, 
                          OMMP_DEFAULT_LA_N_EEL_REMOVE);

    ommp_update_coordinates(qm_sys, ommp_qm_helper_get_cqm(my_qmh));

    pol_atoms = ommp_get_pol_atoms(my_system);
    
    electric_field = (double *) malloc(sizeof(double) * 3 * pol_atoms);
   
    for(int j = 0; j < pol_atoms; j++)
        for(int k = 0; k < 3; k++)
            electric_field[j*3+k] = 0.0;
    
    ommp_set_external_field(my_system, electric_field, OMMP_SOLVER_DEFAULT);
    
    free(electric_field);
    ommp_terminate_qm_helper(my_qmh);
    ommp_terminate(qm_sys);
    ommp_terminate(my_system);
    
    return 0;
}
