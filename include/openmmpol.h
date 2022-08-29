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

extern void ommp_init_mmp(char *);
extern void ommp_init_xyz(char *, char *);
extern void ommp_terminate(void);
#ifdef USE_HDF5
extern void ommp_write_hdf5(char *);
extern void ommp_init_hdf5(char *);
#endif
extern void ommp_set_verbose(int32_t);
extern void ommp_print_summary(void);
extern void ommp_print_summary_to_file(char *);

extern void ommp_get_polelec_energy(double *);
extern void ommp_get_fixedelec_energy(double *);
extern void ommp_set_external_field(double *, int32_t);

extern void ommp_get_vdw_energy(double *);
extern void ommp_get_bond_energy(double *);
extern void ommp_get_angle_energy(double *);
extern void ommp_get_angtor_energy(double *);
extern void ommp_get_strtor_energy(double *);
extern void ommp_get_strbnd_energy(double *);
extern void ommp_get_opb_energy(double *);
extern void ommp_get_pitors_energy(double *);
extern void ommp_get_torsion_energy(double *);
extern void ommp_get_tortor_energy(double *);
extern void ommp_get_urey_energy(double *);

extern int32_t get_n_ipd(void);
extern int32_t ommp_get_ld_cart(void);
extern int32_t ommp_get_mm_atoms(void);
extern int32_t ommp_get_pol_atoms(void);
extern void *ommp_get_cmm(void);
extern void *ommp_get_cpol(void);
extern void *ommp_get_q(void);
extern void *ommp_get_ipd(void);
extern void *ommp_get_polar_mm(void);

extern bool ommp_ff_is_amoeba(void);
