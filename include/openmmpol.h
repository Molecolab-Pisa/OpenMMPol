#include <stdbool.h>
#include <stdint.h>

#define OMMP_VERBOSE_DEBUG 3
#define OMMP_VERBOSE_HIGH 2 
#define OMMP_VERBOSE_LOW 1 
#define OMMP_VERBOSE_NONE 0

#define OMMP_FF_AMOEBA 1
#define OMMP_FF_WANG_AL 0 
#define OMMP_FF_WANG_DL 0 

#define OMMP_SOLVER_CG 1
#define OMMP_SOLVER_DIIS 2
#define OMMP_SOLVER_INVERSION 3
#define OMMP_SOLVER_DEFAULT OMMP_SOLVER_CG

#define OMMP_MATV_INCORE 1
#define OMMP_MATV_DIRECT 2
#define OMMP_MATV_DEFAULT OMMP_MATV_DIRECT

#define OMMP_AMOEBA_D 1
#define OMMP_AMOEBA_P 2

#define AU2KCALMOL 627.5096080306

typedef void* OMMP_SYSTEM_PRT;

#ifdef __cplusplus
extern "C" {
#endif

extern OMMP_SYSTEM_PRT ommp_init_mmp(const char *);
extern OMMP_SYSTEM_PRT ommp_init_xyz(const char *, const char *);
extern void ommp_save_mmp(OMMP_SYSTEM_PRT, const char *, int32_t);
extern void ommp_terminate(OMMP_SYSTEM_PRT);
#ifdef USE_HDF5
extern void ommp_save_as_hdf5(OMMP_SYSTEM_PRT, const char *, const char *);
extern void ommp_checkpoint(OMMP_SYSTEM_PRT, const char *, const char *);
extern OMMP_SYSTEM_PRT ommp_init_hdf5(const char *, const char *);
#endif
extern void ommp_set_verbose(int32_t);
extern void ommp_print_summary(OMMP_SYSTEM_PRT);
extern void ommp_print_summary_to_file(OMMP_SYSTEM_PRT, const char *);

extern double ommp_get_polelec_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_fixedelec_energy(OMMP_SYSTEM_PRT);
extern void ommp_set_external_field(OMMP_SYSTEM_PRT, const double *, int32_t);
extern void ommp_set_external_field_nomm(OMMP_SYSTEM_PRT, const double *, int32_t);
//extern void ommp_potential_mm2ext(int32_t, double *, double *);
//extern void ommp_potential_mmpol2ext(int32_t, double *, double *);
//extern void ommp_potential_pol2ext(int32_t, double *, double *);

extern double ommp_get_vdw_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_bond_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_angle_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_angtor_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_strtor_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_strbnd_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_opb_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_pitors_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_torsion_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_tortor_energy(OMMP_SYSTEM_PRT);
extern double ommp_get_urey_energy(OMMP_SYSTEM_PRT);

extern bool ommp_ff_is_amoeba(OMMP_SYSTEM_PRT);

extern int32_t ommp_get_n_ipd(OMMP_SYSTEM_PRT);
extern int32_t ommp_get_ld_cart(OMMP_SYSTEM_PRT);
extern int32_t ommp_get_mm_atoms(OMMP_SYSTEM_PRT);
extern int32_t ommp_get_pol_atoms(OMMP_SYSTEM_PRT);

extern double *ommp_get_cmm(OMMP_SYSTEM_PRT);
extern double *ommp_get_cpol(OMMP_SYSTEM_PRT);
extern double *ommp_get_q(OMMP_SYSTEM_PRT);
extern double *ommp_get_ipd(OMMP_SYSTEM_PRT);
extern int32_t *ommp_get_polar_mm(OMMP_SYSTEM_PRT);

extern void ommp_update_coordinates(OMMP_SYSTEM_PRT, const double *);

extern void ommp_fixedelec_numgeomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_fixedelec_anageomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_polelec_numgeomgrad(OMMP_SYSTEM_PRT, double *);
extern void ommp_polelec_anageomgrad(OMMP_SYSTEM_PRT, double *);

#ifdef __cplusplus
}
#endif
