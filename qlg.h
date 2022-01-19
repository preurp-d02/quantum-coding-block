#ifndef QLG_H
#define QLG_H

#include <complex.h>

#define QLG_I 0
#define QLG_C 1
#define QLG_X 2
#define QLG_Y 3
#define QLG_Z 4
#define QLG_H 5
#define QLG_CX 6
#define QLG_CY 7
#define QLG_CZ 8
#define QLG_CH 9

typedef struct qlg_qc {
  _Complex double state[2][2][2];
} qlg_qc;

typedef struct qlg_qubit {
  _Complex double state[2];
} qlg_qubit;

extern struct qlg_qc init(qlg_qubit *x);
extern double prob(struct qlg_qc qc, struct qlg_qc phi);
extern struct qlg_qc measure(struct qlg_qc qc);

extern struct qlg_qc qlg_apply(struct qlg_qc qc, int op1, int op2, int op3);
extern struct qlg_qc qlg_apply_arr(struct qlg_qc qc, int *op);
extern struct qlg_qc qlg_i(struct qlg_qc qc, int x);
extern struct qlg_qc qlg_x(struct qlg_qc qc, int x);
extern struct qlg_qc qlg_y(struct qlg_qc qc, int x);
extern struct qlg_qc qlg_z(struct qlg_qc qc, int x);
extern struct qlg_qc qlg_h(struct qlg_qc qc, int x);
extern struct qlg_qc qlg_cx(struct qlg_qc qc, int c, int x);
extern struct qlg_qc qlg_cy(struct qlg_qc qc, int c, int x);
extern struct qlg_qc qlg_cz(struct qlg_qc qc, int c, int x);
extern struct qlg_qc qlg_ch(struct qlg_qc qc, int c, int x);

#endif