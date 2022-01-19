#ifndef QLG_H
#define QLG_H

#include <complex.h>

typedef struct qubit {
  _Complex double zero, one;
} qubit;

extern void qlg_x(struct qubit q);
extern void qlg_y(struct qubit q);
extern void qlg_z(struct qubit q);
extern void qlg_h(struct qubit q);
extern void qlg_cnot(struct qubit q, qubit c);

extern struct qubit qlg_zero();
extern struct qubit qlg_one();
extern struct qubit qlg_plus();
extern struct qubit qlg_minus();
extern struct qubit qlg_r();
extern struct qubit qlg_l();

extern void qlg_set_zero(struct qubit q);
extern void qlg_set_one(struct qubit q);
extern void qlg_set_plus(struct qubit q);
extern void qlg_set_minus(struct qubit q);
extern void qlg_set_r(struct qubit q);
extern void qlg_set_l(struct qubit q);

extern unsigned short qlg_is_zero(struct qubit q);
extern unsigned short qlg_is_one(struct qubit q);
extern unsigned short qlg_is_plus(struct qubit q);
extern unsigned short qlg_is_minus(struct qubit q);
extern unsigned short qlg_is_r(struct qubit q);
extern unsigned short qlg_is_l(struct qubit q);
extern unsigned short qlg_is_x_basis(struct qubit q);
extern unsigned short qlg_is_y_basis(struct qubit q);
extern unsigned short qlg_is_z_basis(struct qubit q);

#endif
