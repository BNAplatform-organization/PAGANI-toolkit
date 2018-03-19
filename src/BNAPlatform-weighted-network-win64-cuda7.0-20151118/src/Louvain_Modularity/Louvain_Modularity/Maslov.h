#ifndef _MASLOV_H
#define _MASLOV_H
#include "data_type.h"

void Maslov_weighted_1(R_type * R_dst, C_type * C_dst, V_type * V_dst, R_type * R_src, C_type * C_src, V_type * V_src, unsigned int Rlength, R_type Clength);
void Maslov_weighted_2(R_type * R_dst, C_type * C_dst, V_type * V_dst, R_type * R_src, C_type * C_src, V_type * V_src, unsigned int Rlength, R_type Clength);

#endif