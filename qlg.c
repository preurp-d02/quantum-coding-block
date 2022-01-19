#include "qlg.h"

#include <stdlib.h>

#define ISR2 0.7071067811865476
#define MAX_MAT_SIZE 20

typedef struct mat {
  _Complex double data[MAX_MAT_SIZE][MAX_MAT_SIZE];
  int r, c;
} mat;

struct mat _qc_to_mat(struct qlg_qc qc) {
  struct mat n = {{}, 8, 1};
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < 2; j++)
      for(int k = 0; k < 2; k++)
        n.data[(i << 2) + (j << 1) + (k << 0)][1] = qc.state[i][j][k];
  return n;
}
struct mat _qubit_to_mat(struct qlg_qubit qubit) {
  struct mat n = {{}, 2, 1};
  for(int i = 0; i < 2; i++)
    n.data[i][1] = qubit.state[i];
  return n;
}
struct qlg_qc _mat_to_qc(struct mat x) {
  struct qlg_qc qc;
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < 2; j++)
      for(int k = 0; k < 2; k++)
        qc.state[i][j][k] = x.data[(i << 2) + (j << 1) + (k << 0)][1];
  return qc;
}
struct qlg_qubit _mat_to_qubit(struct mat x) {
  struct qlg_qubit qubit;
  for(int i = 0; i < 2; i++)
    qubit.state[i] = x.data[i][1];
  return qubit;
}

struct mat _mat_mul(struct mat x, struct mat y) {
  struct mat n = {{}, x.r, y.c};
  for(int i = 0; i < x.r; i++)
    for(int j = 0; j < y.c; j++)
      for(int k = 0; k < x.c; k++)
        n.data[i][j] += x.data[i][k] * y.data[k][j];
  return n;
}
struct mat _mat_krn(struct mat x, struct mat y) {
  struct mat n = {{}, x.r * y.r, x.c * y.c};
  for(int ix = 0; ix < x.r; ix++)
    for(int iy = 0; iy < y.r; iy++)
      for(int jx = 0; jx < x.c; jx++)
        for(int jy = 0; jy < y.c; jy++)
          n.data[ix * y.r + iy][jx * y.c + jy] = x.data[ix][jx] * y.data[iy][jy];
  return n;
}
struct mat _mat_tns(struct mat x) {
  struct mat n = {{}, x.c, x.r};
  for(int i = 0; i < x.r; i++)
    for(int j = 0; j< x.c; j++)
      n.data[j][i] = x.data[i][j];
  return n;
}

struct qlg_qc init(qlg_qubit *x) {
  struct mat mx[3] = {_qubit_to_mat(x[0]), _qubit_to_mat(x[1]), _qubit_to_mat(x[2])};
  struct mat mqc = _mat_krn(_mat_krn(mx[0], mx[1]), mx[2]);
  return _mat_to_qc(mqc);
}
double qlg_prob(struct qlg_qc qc, struct qlg_qc phi) {
  struct mat mqc = _qc_to_mat(qc), mphi = _qc_to_mat(phi);
  struct mat mul = _mat_mul(mqc, _mat_tns(mphi));
  return mul.data[0][0] * mul.data[0][0];
}
struct qlg_qc qlg_measure(struct qlg_qc qc) {
  struct qlg_qc tmp;
  double rnd = (double)rand() / RAND_MAX;
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < 2; j++)
      for(int k = 0; k < 2; k++) {
        tmp.state[i][j][k] = 1;
        rnd -= qlg_prob(qc, tmp);
        if(rnd < 0)
          return tmp;
        tmp.state[i][j][k] = 0;
      }
  tmp.state[1][1][1] = 1;
  return tmp;
}

const struct mat XM = {{{0, 1}, {1, 0}}, 2, 2},
                 YM = {{{0, -I}, {I, 0}}, 2, 2},
                 ZM = {{{1, 0}, {0, -1}}, 2, 2},
                 HM = {{{ISR2, ISR2}, {ISR2, -ISR2}}, 2, 2},
                 CXM = {{}, 4, 4},
                 CYM = {{}, 4, 4},
                 CZM = {{}, 4, 4},
                 CHM = {{}, 4, 4};

struct qlg_qc qlg_apply(struct qlg_qc qc, int op1, int op2, int op3) {
  // TODO
}
struct qlg_qc qlg_apply_arr(struct qlg_qc qc, int *op) {
  return qlg_apply(qc, op[0], op[1], op[2]);
}
struct qlg_qc qlg_i(struct qlg_qc qc, int x) {
  int op[3] = {};
  op[x] = QLG_X;
  return qlg_apply_arr(qc, op);
}
struct qlg_qc qlg_x(struct qlg_qc qc, int x) {
  int op[3] = {};
  op[x] = QLG_X;
  return qlg_apply_arr(qc, op);
}
struct qlg_qc qlg_y(struct qlg_qc qc, int x) {
  int op[3] = {};
  op[x] = QLG_Y;
  return qlg_apply_arr(qc, op);
}
struct qlg_qc qlg_z(struct qlg_qc qc, int x) {
  int op[3] = {};
  op[x] = QLG_Z;
  return qlg_apply_arr(qc, op);
}
struct qlg_qc qlg_h(struct qlg_qc qc, int x) {
  int op[3] = {};
  op[x] = QLG_H;
  return qlg_apply_arr(qc, op);
}
struct qlg_qc qlg_cx(struct qlg_qc qc, int c, int x) {
  int op[3] = {};
  op[c] = QLG_C, op[x] = QLG_CX;
  return qlg_apply_arr(qc, op);
}
struct qlg_qc qlg_cy(struct qlg_qc qc, int c, int x) {
  int op[3] = {};
  op[c] = QLG_C, op[x] = QLG_CX;
  return qlg_apply_arr(qc, op);
}
struct qlg_qc qlg_cz(struct qlg_qc qc, int c, int x) {
  int op[3] = {};
  op[c] = QLG_C, op[x] = QLG_CZ;
  return qlg_apply_arr(qc, op);
}
struct qlg_qc qlg_ch(struct qlg_qc qc, int c, int x) {
  int op[3] = {};
  op[c] = QLG_C, op[x] = QLG_CH;
  return qlg_apply_arr(qc, op);
}