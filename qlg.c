#include "qlg.h"

#define ERRSQ 0.000001
#define ISR2 0.7071067811865476

#define abssq(x) (creal(x) * creal(x) + cimag(x) * cimag(x))

void qlg_x(struct qubit q) {
  _Complex double tmp = q.zero;
  q.zero = q.one;
  q.one = tmp;
  return;
}
void qlg_y(struct qubit q) {
  _Complex double tmp = q.zero;
  q.zero = -q.one * _Complex_I;
  q.one = tmp * _Complex_I;
  return;
}
void qlg_z(struct qubit q) {
  q.one = -q.one;
  return;
}
void qlg_h(struct qubit q) {
  _Complex double tmp = q.zero;
  q.zero = (q.zero + q.one) * ISR2;
  q.one = (tmp - q.one) * ISR2;
  return;
}
void qlg_cnot(struct qubit q, qubit c) {
  // TODO
}

void qlg_set_zero(struct qubit q) { q.zero = 1, q.one = 0; }
void qlg_set_one(struct qubit q) { q.zero = 0, q.one = 1; }
void qlg_set_plus(struct qubit q) { q.zero = ISR2, q.one = ISR2; }
void qlg_set_minus(struct qubit q) { q.zero = ISR2, q.one = -ISR2; }
void qlg_set_r(struct qubit q) { q.zero = ISR2, q.one = ISR2 * _Complex_I; }
void qlg_set_l(struct qubit q) { q.zero = ISR2, q.one = -ISR2 * _Complex_I; }

unsigned short qlg_is_zero(struct qubit q) {
  return abssq(q.zero - 1) <= ERRSQ && abssq(q.one - 0) <= ERRSQ;
}
unsigned short qlg_is_one(struct qubit q) {
  return abssq(q.zero - 0) <= ERRSQ && abssq(q.one - 1) <= ERRSQ;
}
unsigned short qlg_is_plus(struct qubit q) {
  return abssq(q.zero - ISR2) <= ERRSQ && abssq(q.one - ISR2) <= ERRSQ;
}
unsigned short qlg_is_minus(struct qubit q) {
  return abssq(q.zero - ISR2) <= ERRSQ && abssq(q.one + ISR2) <= ERRSQ;
}
unsigned short qlg_is_r(struct qubit q) {
  return abssq(q.zero - ISR2) <= ERRSQ &&
         abssq(q.one - ISR2 * _Complex_I) <= ERRSQ;
}
unsigned short qlg_is_l(struct qubit q) {
  return abssq(q.zero - ISR2) <= ERRSQ &&
         abssq(q.one + ISR2 * _Complex_I) <= ERRSQ;
}
unsigned short qlg_is_x_basis(struct qubit q) {
  return qlg_is_plus(q) || qlg_is_minus(q);
}
unsigned short qlg_is_y_basis(struct qubit q) {
  return qlg_is_r(q) || qlg_is_l(q);
}
unsigned short qlg_is_z_basis(struct qubit q) {
  return qlg_is_zero(q) || qlg_is_one(q);
}
