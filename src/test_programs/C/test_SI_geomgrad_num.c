#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "openmmpol.h"

void print_qmmm_grad(char *name, int32_t mm_atoms, int32_t qm_atoms, 
                     double *grad_mm, double *grad_qm){
    char msg[OMMP_STR_CHAR_MAX];
    
    sprintf(msg, "Grad %s", name);
    ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    
    for(int i = 0; i < mm_atoms; i++){
        sprintf(msg, "MM:%-8d %+20.8g %+20.8g %+20.8g", i,
                grad_mm[i*3+0]*OMMP_AU2KCALMOL*OMMP_ANG2AU,
                grad_mm[i*3+1]*OMMP_AU2KCALMOL*OMMP_ANG2AU,
                grad_mm[i*3+2]*OMMP_AU2KCALMOL*OMMP_ANG2AU);
        ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
    }
    if(grad_qm != NULL && qm_atoms > 0){
        for(int i = 0; i < qm_atoms; i++){
            sprintf(msg, "QM:%-8d %+20.8g %+20.8g %+20.8g", i,
                    grad_qm[i*3+0]*OMMP_AU2KCALMOL*OMMP_ANG2AU,
                    grad_qm[i*3+1]*OMMP_AU2KCALMOL*OMMP_ANG2AU,
                    grad_qm[i*3+2]*OMMP_AU2KCALMOL*OMMP_ANG2AU);
            ommp_message(msg, OMMP_VERBOSE_NONE, "TEST-OUT");
        }
    }
    ommp_message("", OMMP_VERBOSE_NONE, "TEST-OUT");
}

void update_qmmm_coordinates(OMMP_SYSTEM_PRT qm_sys, 
                             OMMP_QM_HELPER_PRT qmh, 
                             OMMP_SYSTEM_PRT mm_sys, 
                             double *new_mm_c, double *new_qm_c){
    // Complete procedure to update coordinates also in presence of a QM system represented
    // by a fake MM with link atoms to the main system.
    ommp_update_coordinates(mm_sys, new_mm_c);
    if(qmh != NULL){
        ommp_qm_helper_update_coord(qmh, new_qm_c);
        if(ommp_use_linkatoms(mm_sys)){
            ommp_update_link_atoms_position(qmh, mm_sys);
        }
        if(qm_sys != NULL){
            ommp_update_coordinates(qm_sys, ommp_qm_helper_get_cqm(qmh));
        }
    }
}

void num_geomgrad(OMMP_SYSTEM_PRT qm_sys,
                  OMMP_QM_HELPER_PRT qmh,
                  OMMP_SYSTEM_PRT mm_sys,
                  double (*ene_f)(OMMP_SYSTEM_PRT,
                                  OMMP_QM_HELPER_PRT,
                                  OMMP_SYSTEM_PRT),
                  double *gqm, double *gmm, double dd){

    double *new_c_mm, *new_c_qm, *cmm, *cqm, tmp;
    int nmm, nqm;
    bool use_qm, use_fake_qm;

    use_qm = (qmh != NULL);
    use_fake_qm = (qm_sys != NULL);

    nmm = ommp_get_mm_atoms(mm_sys);
    if(use_qm) nqm = ommp_qm_helper_get_qm_atoms(qmh);
    cmm = ommp_get_cmm(mm_sys);
    if(use_qm) cqm = ommp_qm_helper_get_cqm(qmh);

    new_c_mm = (double *) malloc(sizeof(double) * 3 * nmm);
    for(int i=0; i < nmm; i++)
        for(int j=0; j < 3; j++) 
            new_c_mm[i*3+j] = cmm[i*3+j];
    
    if(use_qm){
        new_c_qm = (double *) malloc(sizeof(double) * 3 * nqm);
        for(int i=0; i < nqm; i++)
            for(int j=0; j < 3; j++) 
                new_c_qm[i*3+j] = cqm[i*3+j];
    }

    for(int i=0; i < nmm; i++){
        for(int j=0; j < 3; j++){
            new_c_mm[i*3+j] += dd;
            update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
            tmp = ene_f(qm_sys, qmh, mm_sys);

            new_c_mm[i*3+j] -= 2*dd;
            update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
            tmp -= ene_f(qm_sys, qmh, mm_sys);
            gmm[i*3+j] = tmp / (2*dd);
            
            new_c_mm[i*3+j] += dd;
            update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
        }
    }
    
    if(use_qm){
        for(int i=0; i < nqm; i++){
            for(int j=0; j < 3; j++){
                new_c_qm[i*3+j] += dd;
                update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
                tmp = ene_f(qm_sys, qmh, mm_sys);

                new_c_qm[i*3+j] -= 2*dd;
                update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
                tmp -= ene_f(qm_sys, qmh, mm_sys);
                gqm[i*3+j] = tmp / (2*dd);
                
                new_c_qm[i*3+j] += dd;
                update_qmmm_coordinates(qm_sys, qmh, mm_sys, new_c_mm, new_c_qm);
            }
        }
    }
}

void num_grd_print(OMMP_SYSTEM_PRT qm_sys,
                   OMMP_QM_HELPER_PRT qmh,
                   OMMP_SYSTEM_PRT mm_sys,
                   double (*ene_f)(OMMP_SYSTEM_PRT,
                                   OMMP_QM_HELPER_PRT,
                                   OMMP_SYSTEM_PRT),
                   char *name){

    double *gmm, *gqm = NULL;
    int nmm, nqm;
    bool use_qm, use_fake_qm;

    use_qm = (qmh != NULL);
    use_fake_qm = (qm_sys != NULL);

    nmm = ommp_get_mm_atoms(mm_sys);
    if(use_qm) nqm = ommp_qm_helper_get_qm_atoms(qmh);

    gmm = (double *) malloc(sizeof(double) * 3 * nmm);
    if(use_qm) gqm = (double *) malloc(sizeof(double) * 3 * nqm);

    num_geomgrad(qm_sys, qmh, mm_sys, ene_f, gqm, gmm, 1e-6);
    print_qmmm_grad(name, nmm, nqm, gmm, gqm);
}

double ommptest_full_ene(OMMP_SYSTEM_PRT fakeqm, OMMP_QM_HELPER_PRT qmh, OMMP_SYSTEM_PRT sys){
    double ene = 0.0;

    ene += ommp_get_full_energy(sys);
    ene += ommp_get_full_energy(fakeqm);
    ene += ommp_qm_helper_vdw_energy(qmh, sys);

    return ene;
}

double ommptest_em_ene(OMMP_SYSTEM_PRT fakeqm, OMMP_QM_HELPER_PRT qmh, OMMP_SYSTEM_PRT sys){
    // Just for signature consistency!
    return ommp_get_fixedelec_energy(sys);
}

double ommptest_ep_ene(OMMP_SYSTEM_PRT fakeqm, OMMP_QM_HELPER_PRT qmh, OMMP_SYSTEM_PRT sys){
    // Just for signature consistency!
    return ommp_get_polelec_energy(sys);
}

int main(int argc, char **argv){
    if(argc != 2){
        printf("Given a JSON SmartInput file, it computes the geoemtrical derivatives \n");
        printf("with respect to the MM atoms coordinates using numerical differentiation.\n");
        printf("If a QM system is defined in the JSON file and atom types are provided, \n");
        printf("it is converted in a non-pol MM system interacting with the main one only \n");
        printf("generating an electric field that is used for induced dipole calculation \n");
        printf("(but not in electrostatic interaction!), with VdW interaction and, if \n");
        printf("defined, by link-atoms.\n\n");
        printf("Syntax expected\n");
        printf("    $ test_SI_geomgrad_num.exe <JSON FILE>\n");
        return 0;
    }
    
    OMMP_SYSTEM_PRT my_system, fake_qm = NULL;
    OMMP_QM_HELPER_PRT my_qmh = NULL;
    ommp_smartinput(argv[1], &my_system, &my_qmh);
    
    bool use_qm, use_fake_qm;
    
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
        ommp_smartimput_cpstr(argv[1], addr, &prm_file);
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

    num_grd_print(fake_qm, my_qmh, my_system , ommptest_em_ene, "EM");
    num_grd_print(fake_qm, my_qmh, my_system , ommptest_ep_ene, "EP");
    
    if(my_qmh != NULL) ommp_terminate_qm_helper(my_qmh);
    ommp_terminate(my_system);
    
    return 0;
}
