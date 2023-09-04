#ifndef _OPENMMPOL
#define _OPENMMPOL

#include <stdbool.h>
#include <stdint.h>
//#include <math.h>
#include <openmmpol_const.h>

typedef void* OMMP_SYSTEM_PRT;
typedef void* OMMP_QM_HELPER_PRT;

#ifdef __cplusplus
extern "C" {
#endif

extern OMMP_SYSTEM_PRT ommp_init_mmp(const char *);
extern OMMP_SYSTEM_PRT ommp_init_xyz(const char *, const char *);
extern void ommp_set_default_solver(OMMP_SYSTEM_PRT, int32_t);
extern void ommp_set_default_matv(OMMP_SYSTEM_PRT, int32_t);
extern void ommp_save_mmp(OMMP_SYSTEM_PRT, const char *, int32_t);
extern void ommp_set_frozen_atoms(OMMP_SYSTEM_PRT, int32_t, const int32_t *);
extern void ommp_turn_pol_off(OMMP_SYSTEM_PRT, int32_t, const int32_t *);
extern void ommp_terminate(OMMP_SYSTEM_PRT);
#ifdef USE_HDF5
extern void ommp_save_as_hdf5(OMMP_SYSTEM_PRT, const char *, const char *);
extern void ommp_checkpoint(OMMP_SYSTEM_PRT, const char *, const char *);
extern OMMP_SYSTEM_PRT ommp_init_hdf5(const char *, const char *);
#endif
extern void ommp_set_verbose(int32_t);
extern void ommp_set_outputfile(const char *);
extern void ommp_close_outputfile(void);
extern void ommp_message(const char *, int32_t, const char *);
extern void ommp_fatal(const char *);
extern void ommp_print_summary(OMMP_SYSTEM_PRT);
extern void ommp_print_summary_to_file(OMMP_SYSTEM_PRT, const char *);

extern double ommp_get_polelec_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_fixedelec_energy(OMMP_SYSTEM_PRT);
extern void ommp_set_external_field(OMMP_SYSTEM_PRT, const double *, int32_t, int32_t);
extern void ommp_set_external_field_nomm(OMMP_SYSTEM_PRT, const double *, int32_t, int32_t);

extern void ommp_potential_mmpol2ext(OMMP_SYSTEM_PRT, int32_t, const double *, double *);
extern void ommp_potential_mm2ext(OMMP_SYSTEM_PRT, int32_t, const double *, double *);
extern void ommp_potential_pol2ext(OMMP_SYSTEM_PRT, int32_t, const double *, double *);
extern void ommp_field_mmpol2ext(OMMP_SYSTEM_PRT, int32_t, const double *, double *);
extern void ommp_field_mm2ext(OMMP_SYSTEM_PRT, int32_t, const double *, double *);
extern void ommp_field_pol2ext(OMMP_SYSTEM_PRT, int32_t, const double *, double *);

extern double ommp_get_vdw_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_bond_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_angle_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_angtor_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_strtor_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_strbnd_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_opb_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_pitors_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_torsion_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_imptorsion_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_tortor_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_urey_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_full_bnd_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_full_ele_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_full_energy(OMMP_SYSTEM_PRT);

extern bool ommp_ff_is_amoeba(OMMP_SYSTEM_PRT);

extern int32_t ommp_get_n_ipd(OMMP_SYSTEM_PRT);
extern int32_t ommp_get_ld_cart(OMMP_SYSTEM_PRT);
extern int32_t ommp_get_mm_atoms(OMMP_SYSTEM_PRT);
extern int32_t ommp_get_pol_atoms(OMMP_SYSTEM_PRT);

extern double *ommp_get_cmm(OMMP_SYSTEM_PRT);
extern int32_t *ommp_get_zmm(OMMP_SYSTEM_PRT);
extern int32_t *ommp_get_attypemm(OMMP_SYSTEM_PRT);
extern double *ommp_get_cpol(OMMP_SYSTEM_PRT);
extern double *ommp_get_q(OMMP_SYSTEM_PRT);
extern double *ommp_get_ipd(OMMP_SYSTEM_PRT);
extern int32_t *ommp_get_polar_mm(OMMP_SYSTEM_PRT);

extern bool ommp_use_frozen(OMMP_SYSTEM_PRT);
extern bool *ommp_get_frozen(OMMP_SYSTEM_PRT);

extern bool ommp_use_linkatoms(OMMP_SYSTEM_PRT);

extern void ommp_update_coordinates(OMMP_SYSTEM_PRT, const double *);

extern void ommp_full_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_full_bnd_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_fixedelec_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_polelec_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_rotation_geomgrad(OMMP_SYSTEM_PRT, const double *, const double *, double *);
extern void ommp_vdw_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_bond_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_angle_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_strbnd_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_urey_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_torsion_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_imptorsion_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_angtor_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_opb_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_strtor_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_tortor_geomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_pitors_geomgrad(OMMP_SYSTEM_PRT, double *);

extern OMMP_QM_HELPER_PRT ommp_init_qm_helper(int32_t, const double *, const double *, const int32_t *);
extern void ommp_terminate_qm_helper(OMMP_QM_HELPER_PRT);
extern void ommp_qm_helper_set_frozen_atoms(OMMP_QM_HELPER_PRT, int32_t, const int32_t *);
extern void ommp_qm_helper_update_coord(OMMP_QM_HELPER_PRT, const double *);

extern void ommp_prepare_qm_ele_ene(OMMP_SYSTEM_PRT, OMMP_QM_HELPER_PRT);
extern void ommp_prepare_qm_ele_grd(OMMP_SYSTEM_PRT, OMMP_QM_HELPER_PRT);
extern double *ommp_qm_helper_get_V_m2n(OMMP_QM_HELPER_PRT);
extern double *ommp_qm_helper_get_V_p2n(OMMP_QM_HELPER_PRT);
extern double *ommp_qm_helper_get_E_m2n(OMMP_QM_HELPER_PRT);
extern double *ommp_qm_helper_get_E_p2n(OMMP_QM_HELPER_PRT);
extern double *ommp_qm_helper_get_E_n2p(OMMP_QM_HELPER_PRT);
extern double *ommp_qm_helper_get_G_n2p(OMMP_QM_HELPER_PRT);
extern double *ommp_qm_helper_get_E_n2m(OMMP_QM_HELPER_PRT);
extern double *ommp_qm_helper_get_G_n2m(OMMP_QM_HELPER_PRT);
extern double *ommp_qm_helper_get_H_n2m(OMMP_QM_HELPER_PRT);
extern double *ommp_qm_helper_get_cqm(OMMP_QM_HELPER_PRT);
extern int32_t ommp_qm_helper_get_npol(OMMP_QM_HELPER_PRT);
extern int32_t ommp_qm_helper_get_nmm(OMMP_QM_HELPER_PRT);
extern int32_t ommp_qm_helper_get_qm_atoms(OMMP_QM_HELPER_PRT);
extern bool *ommp_qm_helper_get_frozen(OMMP_QM_HELPER_PRT);
extern bool ommp_qm_helper_use_frozen(OMMP_QM_HELPER_PRT);
extern void ommp_qm_helper_init_vdw_prm(OMMP_QM_HELPER_PRT, const char *);
extern void ommp_qm_helper_set_attype(OMMP_QM_HELPER_PRT, const int32_t *);
extern void ommp_qm_helper_init_vdw(OMMP_QM_HELPER_PRT, const double *, const double *,
                                    const double *, const char *, const char *, 
                                    const char *, const char *, const char *);
extern double ommp_qm_helper_vdw_energy(OMMP_QM_HELPER_PRT, OMMP_SYSTEM_PRT);  
extern void ommp_qm_helper_vdw_geomgrad(OMMP_QM_HELPER_PRT, OMMP_SYSTEM_PRT, 
                                          double *, double *);
extern void  ommp_qm_helper_link_atom_geomgrad(OMMP_QM_HELPER_PRT, OMMP_SYSTEM_PRT, 
                                               double *, double *, const double *);
extern bool ommp_qm_helper_use_nonbonded(OMMP_QM_HELPER_PRT);

extern void ommp_qm_helper_init_link_atom(OMMP_QM_HELPER_PRT, OMMP_SYSTEM_PRT);
extern int32_t ommp_create_link_atom(OMMP_QM_HELPER_PRT, OMMP_SYSTEM_PRT,
                                     int32_t, int32_t, int32_t, const char *, 
                                     double, int32_t);
extern void ommp_get_link_atom_coordinates(OMMP_SYSTEM_PRT, int32_t, double *);
extern void ommp_update_link_atoms_position(OMMP_QM_HELPER_PRT, OMMP_SYSTEM_PRT);

extern void ommp_smartinput(const char *, OMMP_SYSTEM_PRT, OMMP_QM_HELPER_PRT);
extern void ommp_smartinput_cpstr(const char *, char *, char **);
extern OMMP_SYSTEM_PRT ommp_system_from_qm_helper(OMMP_QM_HELPER_PRT, const char *);
#ifdef __cplusplus
}
#endif
#endif
