/*
 * Copyright (c) The acados authors.
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#ifndef point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_COST
#define point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_COST

#ifdef __cplusplus
extern "C" {
#endif


// Cost at initial shooting node

int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_work(int *, int *, int *, int *);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_sparsity_in(int);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_sparsity_out(int);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_n_in(void);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_n_out(void);

int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_jac_ut_xt_sparsity_in(int);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_jac_ut_xt_sparsity_out(int);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_jac_ut_xt_n_in(void);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_jac_ut_xt_n_out(void);



// Cost at path shooting node

int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_work(int *, int *, int *, int *);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_sparsity_in(int);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_sparsity_out(int);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_n_in(void);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_n_out(void);

int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_jac_ut_xt_sparsity_in(int);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_jac_ut_xt_sparsity_out(int);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_jac_ut_xt_n_in(void);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_jac_ut_xt_n_out(void);



// Cost at terminal shooting node

int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_work(int *, int *, int *, int *);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_sparsity_in(int);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_sparsity_out(int);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_n_in(void);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_n_out(void);

int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_jac_ut_xt(const real_t** arg, real_t** res, int* iw, real_t* w, void *mem);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_jac_ut_xt_work(int *, int *, int *, int *);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_jac_ut_xt_sparsity_in(int);
const int *point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_jac_ut_xt_sparsity_out(int);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_jac_ut_xt_n_in(void);
int point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_jac_ut_xt_n_out(void);




#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_COST
