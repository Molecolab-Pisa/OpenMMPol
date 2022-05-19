#include <stdbool.h>
#include <stdint.h>

extern void w_mmpol_init(char *);
extern void do_mm(void);
extern void do_qmmm(double *, double *, int32_t, int32_t, int32_t, int32_t);
extern void restart(void);
extern void get_energy(double, double);
extern void write_hdf5(char *);

extern int32_t get_n_ipd(void);
extern int32_t get_ld_cart(void);
extern int32_t get_mm_atoms(void);
extern int32_t get_pol_atoms(void);

extern void *get_cmm(void);
extern void *get_cpol(void);
extern void *get_q(void);
extern void *get_ipd(void);

extern bool is_amoeba(void);

