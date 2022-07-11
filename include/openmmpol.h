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

extern void mmpol_init_mmp(char *);
extern void mmpol_init_xyz(char *, char *);
extern void do_mm(void);
extern void do_qmmm(double *, int32_t);
extern void get_energy(double *, double *);
extern void get_vdw_energy(double *);
extern void get_bond_energy(double *);
extern void get_urey_energy(double *);
extern void write_hdf5(char *);
extern void print_summary(void);
extern void print_summary_to_file(char *);

extern int32_t get_n_ipd(void);
extern int32_t get_ld_cart(void);
extern int32_t get_mm_atoms(void);
extern int32_t get_pol_atoms(void);

extern void *get_cmm(void);
extern void *get_cpol(void);
extern void *get_q(void);
extern void *get_ipd(void);
extern void *get_polar_mm(void);

extern bool is_amoeba(void);

extern void set_verbose(int32_t);
