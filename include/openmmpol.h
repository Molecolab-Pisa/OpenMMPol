#include <stdbool.h>
#include <stdint.h>

extern void w_mmpol_init(char *);
extern void do_mm(void);
extern void do_qmmm(void);
extern void restart(void);
extern void get_energy(double, double);

extern int32_t get_n_ipd(void);
extern int32_t get_ld_cart(void);
extern int32_t get_mm_atoms(void);
extern int32_t get_pol_atoms(void);

//! cmm cpol q ipd

extern bool is_amoeba(void);
