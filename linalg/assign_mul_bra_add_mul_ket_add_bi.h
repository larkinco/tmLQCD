/***********************************************************************
 * Copyright (C) 2006 Thomas Chiarappa
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#ifndef _ASSIGN_MUL_BRA_ADD_MUL_KET_ADD_BI_H
#define _ASSIGN_MUL_BRA_ADD_MUL_KET_ADD_BI_H

#include "su3.h"

/* (*R) =  c2*(*R + c1*(*S)) + (*U) */
void assign_mul_bra_add_mul_ket_add_bi(bispinor * const R, bispinor * const S, bispinor * const U, const _Complex double c1, const _Complex double c2, const int N);

#endif
