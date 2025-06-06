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













// standard
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
// acados
// #include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"

// example specific
#include "model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_model/model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_model.h"
#include "model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constraints/model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constraints.h"
#include "model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost/model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost.h"
#include "transition_model_8_4_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_model/transition_model_8_4_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_model.h"
#include "transition_model_8_4_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost/transition_model_8_4_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost.h"
#include "model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_model/model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_model.h"
#include "model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constraints/model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constraints.h"
#include "model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost/model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost.h"
#include "transition_model_point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_model/transition_model_point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_model.h"
#include "transition_model_point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost/transition_model_point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost.h"
#include "point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_model/point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_model.h"
#include "point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constraints/point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constraints.h"
#include "point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost/point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost.h"





#include "acados_solver_multiphase_ocp_10_1_10_1_10_20250606_034423_645279.h"


#define MULTIPHASE_OCP_10_1_10_1_10_20250606_034423_645279_N      32
#define NP_0     0
#define NP_1     0
#define NP_2     0
#define NP_3     0
#define NP_4     0


multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule * multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_capsule(void)
{
    void* capsule_mem = malloc(sizeof(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule));
    multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule *capsule = (multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule *) capsule_mem;

    return capsule;
}


int multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_free_capsule(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule *capsule)
{
    free(capsule);
    return 0;
}


int multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule)
{
    int N_shooting_intervals = MULTIPHASE_OCP_10_1_10_1_10_20250606_034423_645279_N;
    double* new_time_steps = NULL; // NULL -> don't alter the code generated time-steps
    return multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_with_discretization(capsule, N_shooting_intervals, new_time_steps);
}


/**
 * Internal function for multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create: step 1
 */
void multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_set_plan(ocp_nlp_plan_t* nlp_solver_plan, const int N)
{
    assert(N == nlp_solver_plan->N);

    /************************************************
    *  plan
    ************************************************/

    nlp_solver_plan->nlp_solver = SQP;

    nlp_solver_plan->ocp_qp_solver_plan.qp_solver = FULL_CONDENSING_HPIPM;
    nlp_solver_plan->regularization = NO_REGULARIZE;
    nlp_solver_plan->globalization = FIXED_STEP;

    nlp_solver_plan->nlp_cost[0] = NONLINEAR_LS;
    nlp_solver_plan->nlp_constraints[0] = BGH;
    for (int i = 1; i < 10; i++)
    {
        nlp_solver_plan->nlp_cost[i] = NONLINEAR_LS;
        nlp_solver_plan->nlp_constraints[i] = BGH;
    }
    for (int i = 0; i < 10; i++)
    {
        nlp_solver_plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        nlp_solver_plan->sim_solver_plan[i].sim_solver = IRK;
    }
    for (int i = 10; i < 11; i++)
    {
        nlp_solver_plan->nlp_cost[i] = NONLINEAR_LS;
        nlp_solver_plan->nlp_constraints[i] = BGH;
    }
    for (int i = 10; i < 11; i++)
    {
        nlp_solver_plan->nlp_dynamics[i] = DISCRETE_MODEL;
        // discrete dynamics does not need sim solver option, this field is ignored
        nlp_solver_plan->sim_solver_plan[i].sim_solver = INVALID_SIM_SOLVER;
    }
    for (int i = 11; i < 21; i++)
    {
        nlp_solver_plan->nlp_cost[i] = NONLINEAR_LS;
        nlp_solver_plan->nlp_constraints[i] = BGH;
    }
    for (int i = 11; i < 21; i++)
    {
        nlp_solver_plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        nlp_solver_plan->sim_solver_plan[i].sim_solver = IRK;
    }
    for (int i = 21; i < 22; i++)
    {
        nlp_solver_plan->nlp_cost[i] = NONLINEAR_LS;
        nlp_solver_plan->nlp_constraints[i] = BGH;
    }
    for (int i = 21; i < 22; i++)
    {
        nlp_solver_plan->nlp_dynamics[i] = DISCRETE_MODEL;
        // discrete dynamics does not need sim solver option, this field is ignored
        nlp_solver_plan->sim_solver_plan[i].sim_solver = INVALID_SIM_SOLVER;
    }
    for (int i = 22; i < 32; i++)
    {
        nlp_solver_plan->nlp_cost[i] = NONLINEAR_LS;
        nlp_solver_plan->nlp_constraints[i] = BGH;
    }
    for (int i = 22; i < 32; i++)
    {
        nlp_solver_plan->nlp_dynamics[i] = CONTINUOUS_MODEL;
        nlp_solver_plan->sim_solver_plan[i].sim_solver = ERK;
    }

    nlp_solver_plan->nlp_cost[N] = NONLINEAR_LS;
    nlp_solver_plan->nlp_constraints[N] = BGH;
}



/**
 * Internal function for multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create: step 2
 */
ocp_nlp_dims* multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_setup_dimensions(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule)
{
    ocp_nlp_plan_t* nlp_solver_plan = capsule->nlp_solver_plan;
    const int N = nlp_solver_plan->N;
    ocp_nlp_config* nlp_config = capsule->nlp_config;
    int i;

    /************************************************
    *  dimensions
    ************************************************/
    #define NINTNP1MEMS 18
    int* intNp1mem = (int*)malloc( (N+1)*sizeof(int)*NINTNP1MEMS );

    int* nx    = intNp1mem + (N+1)*0;
    int* nu    = intNp1mem + (N+1)*1;
    int* nbx   = intNp1mem + (N+1)*2;
    int* nbu   = intNp1mem + (N+1)*3;
    int* nsbx  = intNp1mem + (N+1)*4;
    int* nsbu  = intNp1mem + (N+1)*5;
    int* nsg   = intNp1mem + (N+1)*6;
    int* nsh   = intNp1mem + (N+1)*7;
    int* nsphi = intNp1mem + (N+1)*8;
    int* ns    = intNp1mem + (N+1)*9;
    int* ng    = intNp1mem + (N+1)*10;
    int* nh    = intNp1mem + (N+1)*11;
    int* nphi  = intNp1mem + (N+1)*12;
    int* nz    = intNp1mem + (N+1)*13;
    int* ny    = intNp1mem + (N+1)*14;
    int* nr    = intNp1mem + (N+1)*15;
    int* nbxe  = intNp1mem + (N+1)*16;
    int* np    = intNp1mem + (N+1)*17;
    for (i = 0; i < 10; i++)
    {
        // common
        nx[i] = 16;
        nu[i] = 2;
        nz[i] = 0;
        ns[i] = 0;
        np[i] = 0;
        // cost
        ny[i] = 6;
        // constraints
        nbu[i] = 0;
        nbx[i] = 0;
        ng[i] = 0;
        nh[i] = 1;
        nphi[i] = 0;
        nr[i] = 0;
        // slacks
        nsbu[i] = 0;
        nsbx[i] = 0;
        nsg[i] = 0;
        nsh[i] = 0;
        nsphi[i] = 0;
        nbxe[i] = 0;
    }
    for (i = 10; i < 11; i++)
    {
        // common
        nx[i] = 16;
        nu[i] = 0;
        nz[i] = 0;
        ns[i] = 0;
        np[i] = 0;
        // cost
        ny[i] = 16;
        // constraints
        nbu[i] = 0;
        nbx[i] = 0;
        ng[i] = 0;
        nh[i] = 0;
        nphi[i] = 0;
        nr[i] = 0;
        // slacks
        nsbu[i] = 0;
        nsbx[i] = 0;
        nsg[i] = 0;
        nsh[i] = 0;
        nsphi[i] = 0;
        nbxe[i] = 0;
    }
    for (i = 11; i < 21; i++)
    {
        // common
        nx[i] = 8;
        nu[i] = 2;
        nz[i] = 0;
        ns[i] = 0;
        np[i] = 0;
        // cost
        ny[i] = 6;
        // constraints
        nbu[i] = 0;
        nbx[i] = 0;
        ng[i] = 0;
        nh[i] = 1;
        nphi[i] = 0;
        nr[i] = 0;
        // slacks
        nsbu[i] = 0;
        nsbx[i] = 0;
        nsg[i] = 0;
        nsh[i] = 0;
        nsphi[i] = 0;
        nbxe[i] = 0;
    }
    for (i = 21; i < 22; i++)
    {
        // common
        nx[i] = 8;
        nu[i] = 0;
        nz[i] = 0;
        ns[i] = 0;
        np[i] = 0;
        // cost
        ny[i] = 8;
        // constraints
        nbu[i] = 0;
        nbx[i] = 0;
        ng[i] = 0;
        nh[i] = 0;
        nphi[i] = 0;
        nr[i] = 0;
        // slacks
        nsbu[i] = 0;
        nsbx[i] = 0;
        nsg[i] = 0;
        nsh[i] = 0;
        nsphi[i] = 0;
        nbxe[i] = 0;
    }
    for (i = 22; i < 32; i++)
    {
        // common
        nx[i] = 4;
        nu[i] = 2;
        nz[i] = 0;
        ns[i] = 0;
        np[i] = 0;
        // cost
        ny[i] = 6;
        // constraints
        nbu[i] = 0;
        nbx[i] = 0;
        ng[i] = 0;
        nh[i] = 1;
        nphi[i] = 0;
        nr[i] = 0;
        // slacks
        nsbu[i] = 0;
        nsbx[i] = 0;
        nsg[i] = 0;
        nsh[i] = 0;
        nsphi[i] = 0;
        nbxe[i] = 0;
    }

    /* initial node*/
    i = 0;
    // common
    nx[i] = 16;
    nu[i] = 2;
    nz[i] = 0;
    ns[i] = 0;
    np[i] = 0;
    // cost
    ny[i] = 6;
    // constraints
    nbu[i] = 0;
    nbx[i] = 16;
    nbxe[i] = 16;
    ng[i] = 0;
    nh[i] = 0;
    nphi[i] = 0;
    nr[i] = 0;
    // slacks
    nsbu[i] = 0;
    nsbx[i] = 0;
    nsg[i] = 0;
    nsh[i] = 0;
    nsphi[i] = 0;

    /* terminal node */
    // common
    i = N;
    nx[i] = 4;
    nu[i] = 0;
    nz[i] = 0;
    ns[i] = 0;
    np[i] = 0;
    // cost
    ny[i] = 4;
    // constraints
    nbu[i] = 0;
    nbx[i] = 0;
    ng[i] = 0;
    nh[i] = 1;
    nphi[i] = 0;
    nr[i] = 0;
    // slacks
    nsbu[i] = 0;
    nsbx[i] = 0;
    nsg[i] = 0;
    nsh[i] = 0;
    nsphi[i] = 0;
    nbxe[i] = 0;

    /* create and set ocp_nlp_dims */
    ocp_nlp_dims * nlp_dims = ocp_nlp_dims_create(nlp_config);

    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nx", nx);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nu", nu);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "nz", nz);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "ns", ns);
    ocp_nlp_dims_set_opt_vars(nlp_config, nlp_dims, "np", np);

    ocp_nlp_dims_set_global(nlp_config, nlp_dims, "np_global", 0);
    ocp_nlp_dims_set_global(nlp_config, nlp_dims, "n_global_data", 0);

    for (int i = 0; i <= N; i++)
    {
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nbx", &nbx[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nbu", &nbu[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsbx", &nsbx[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsbu", &nsbu[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "ng", &ng[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsg", &nsg[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nbxe", &nbxe[i]);
    }


    ocp_nlp_dims_set_cost(nlp_config, nlp_dims, 0, "ny", &ny[0]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, 0, "nh", &nh[0]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, 0, "nsh", &nsh[0]);


    for (int i = 1; i < 10; i++)
    {
        ocp_nlp_dims_set_cost(nlp_config, nlp_dims, i, "ny", &ny[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nh", &nh[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsh", &nsh[i]);
    }
    for (int i = 10; i < 11; i++)
    {
        ocp_nlp_dims_set_cost(nlp_config, nlp_dims, i, "ny", &ny[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nh", &nh[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsh", &nsh[i]);
    }
    for (int i = 11; i < 21; i++)
    {
        ocp_nlp_dims_set_cost(nlp_config, nlp_dims, i, "ny", &ny[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nh", &nh[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsh", &nsh[i]);
    }
    for (int i = 21; i < 22; i++)
    {
        ocp_nlp_dims_set_cost(nlp_config, nlp_dims, i, "ny", &ny[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nh", &nh[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsh", &nsh[i]);
    }
    for (int i = 22; i < 32; i++)
    {
        ocp_nlp_dims_set_cost(nlp_config, nlp_dims, i, "ny", &ny[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nh", &nh[i]);
        ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, i, "nsh", &nsh[i]);
    }



    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nh", &nh[N]);
    ocp_nlp_dims_set_constraints(nlp_config, nlp_dims, N, "nsh", &nsh[N]);
    ocp_nlp_dims_set_cost(nlp_config, nlp_dims, N, "ny", &ny[N]);

    free(intNp1mem);

    return nlp_dims;
}



/**
 * Internal function for multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create: step 3
 */
void multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_setup_functions(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule)
{
    const int N = capsule->nlp_solver_plan->N;

    /************************************************
    *  external functions
    ************************************************/

#define MAP_CASADI_FNC(__CAPSULE_FNC__, __MODEL_BASE_FNC__) do{ \
        capsule->__CAPSULE_FNC__.casadi_fun = & __MODEL_BASE_FNC__ ;\
        capsule->__CAPSULE_FNC__.casadi_n_in = & __MODEL_BASE_FNC__ ## _n_in; \
        capsule->__CAPSULE_FNC__.casadi_n_out = & __MODEL_BASE_FNC__ ## _n_out; \
        capsule->__CAPSULE_FNC__.casadi_sparsity_in = & __MODEL_BASE_FNC__ ## _sparsity_in; \
        capsule->__CAPSULE_FNC__.casadi_sparsity_out = & __MODEL_BASE_FNC__ ## _sparsity_out; \
        capsule->__CAPSULE_FNC__.casadi_work = & __MODEL_BASE_FNC__ ## _work; \
        external_function_external_param_casadi_create(&capsule->__CAPSULE_FNC__, &ext_fun_opts); \
    } while(false)

    external_function_opts ext_fun_opts;
    external_function_opts_set_to_default(&ext_fun_opts);



    ext_fun_opts.external_workspace = true;


    // nonlinear least squares function
    MAP_CASADI_FNC(cost_y_0_fun, model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun);
    MAP_CASADI_FNC(cost_y_0_fun_jac_ut_xt, model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_0_fun_jac_ut_xt);



/////////////// PATH
    int n_path, n_cost_path;
    n_path = 10;
    n_cost_path = 9;
    capsule->nl_constr_h_fun_jac_0 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++) {
        MAP_CASADI_FNC(nl_constr_h_fun_jac_0[i], model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constr_h_fun_jac_uxt_zt);
    }
    capsule->nl_constr_h_fun_0 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++) {
        MAP_CASADI_FNC(nl_constr_h_fun_0[i], model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constr_h_fun);
    }
    



    // implicit dae
    capsule->impl_dae_fun_0 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(impl_dae_fun_0[i], model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_impl_dae_fun);
    }

    capsule->impl_dae_fun_jac_x_xdot_z_0 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(impl_dae_fun_jac_x_xdot_z_0[i], model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_impl_dae_fun_jac_x_xdot_z);
    }

    capsule->impl_dae_jac_x_xdot_u_z_0 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(impl_dae_jac_x_xdot_u_z_0[i], model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_impl_dae_jac_x_xdot_u_z);
    }

    // nonlinear least squares cost
    capsule->cost_y_fun_0 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_0[i], model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun);
    }

    capsule->cost_y_fun_jac_ut_xt_0 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_jac_ut_xt_0[i], model_8_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_jac_ut_xt);
    }
    n_path = 1;
    n_cost_path = 1;


    // discrete dynamics
    capsule->discr_dyn_phi_fun_1 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++)
    {
        MAP_CASADI_FNC(discr_dyn_phi_fun_1[i], transition_model_8_4_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_dyn_disc_phi_fun);
    }

    capsule->discr_dyn_phi_fun_jac_ut_xt_1 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++)
    {
        MAP_CASADI_FNC(discr_dyn_phi_fun_jac_ut_xt_1[i], transition_model_8_4_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_dyn_disc_phi_fun_jac);
    }
    // nonlinear least squares cost
    capsule->cost_y_fun_1 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_1[i], transition_model_8_4_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun);
    }

    capsule->cost_y_fun_jac_ut_xt_1 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_jac_ut_xt_1[i], transition_model_8_4_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_jac_ut_xt);
    }
    n_path = 10;
    n_cost_path = 10;
    capsule->nl_constr_h_fun_jac_2 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++) {
        MAP_CASADI_FNC(nl_constr_h_fun_jac_2[i], model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constr_h_fun_jac_uxt_zt);
    }
    capsule->nl_constr_h_fun_2 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++) {
        MAP_CASADI_FNC(nl_constr_h_fun_2[i], model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constr_h_fun);
    }
    



    // implicit dae
    capsule->impl_dae_fun_2 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(impl_dae_fun_2[i], model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_impl_dae_fun);
    }

    capsule->impl_dae_fun_jac_x_xdot_z_2 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(impl_dae_fun_jac_x_xdot_z_2[i], model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_impl_dae_fun_jac_x_xdot_z);
    }

    capsule->impl_dae_jac_x_xdot_u_z_2 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(impl_dae_jac_x_xdot_u_z_2[i], model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_impl_dae_jac_x_xdot_u_z);
    }

    // nonlinear least squares cost
    capsule->cost_y_fun_2 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_2[i], model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun);
    }

    capsule->cost_y_fun_jac_ut_xt_2 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_jac_ut_xt_2[i], model_4_links_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_jac_ut_xt);
    }
    n_path = 1;
    n_cost_path = 1;


    // discrete dynamics
    capsule->discr_dyn_phi_fun_3 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++)
    {
        MAP_CASADI_FNC(discr_dyn_phi_fun_3[i], transition_model_point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_dyn_disc_phi_fun);
    }

    capsule->discr_dyn_phi_fun_jac_ut_xt_3 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++)
    {
        MAP_CASADI_FNC(discr_dyn_phi_fun_jac_ut_xt_3[i], transition_model_point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_dyn_disc_phi_fun_jac);
    }
    // nonlinear least squares cost
    capsule->cost_y_fun_3 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_3[i], transition_model_point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun);
    }

    capsule->cost_y_fun_jac_ut_xt_3 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_jac_ut_xt_3[i], transition_model_point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_jac_ut_xt);
    }
    n_path = 10;
    n_cost_path = 10;
    capsule->nl_constr_h_fun_jac_4 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++) {
        MAP_CASADI_FNC(nl_constr_h_fun_jac_4[i], point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constr_h_fun_jac_uxt_zt);
    }
    capsule->nl_constr_h_fun_4 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++) {
        MAP_CASADI_FNC(nl_constr_h_fun_4[i], point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constr_h_fun);
    }
    



    // explicit ode
    capsule->expl_vde_forw_4 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(expl_vde_forw_4[i], point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_expl_vde_forw);
    }

    capsule->expl_vde_adj_4 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(expl_vde_adj_4[i], point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_expl_vde_adj);
    }

    capsule->expl_ode_fun_4 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_path);
    for (int i = 0; i < n_path; i++) {
        MAP_CASADI_FNC(expl_ode_fun_4[i], point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_expl_ode_fun);
    }


    // nonlinear least squares cost
    capsule->cost_y_fun_4 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_4[i], point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun);
    }

    capsule->cost_y_fun_jac_ut_xt_4 = (external_function_external_param_casadi *) malloc(sizeof(external_function_external_param_casadi)*n_cost_path);
    for (int i = 0; i < n_cost_path; i++)
    {
        MAP_CASADI_FNC(cost_y_fun_jac_ut_xt_4[i], point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_fun_jac_ut_xt);
    }




    MAP_CASADI_FNC(nl_constr_h_e_fun_jac, point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constr_h_e_fun_jac_uxt_zt);
    MAP_CASADI_FNC(nl_constr_h_e_fun, point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_constr_h_e_fun);
    
    
    // nonlinear least square function
    MAP_CASADI_FNC(cost_y_e_fun, point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun);
    MAP_CASADI_FNC(cost_y_e_fun_jac_ut_xt, point_mass_8_4_p_10_10_10_015000000000000005_FULL_CONDENSING_HPIPM_cost_y_e_fun_jac_ut_xt);

#undef MAP_CASADI_FNC
}



/**
 * Internal function for multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create: step 4
 */
void multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_set_default_parameters(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule) {

    double* p = calloc(0, sizeof(double));
    // initialize parameters to nominal value

    for (int i = 0; i < 10; i++) {
        multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_update_params(capsule, i, p, NP_0);
    }

    for (int i = 10; i < 11; i++) {
        multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_update_params(capsule, i, p, NP_1);
    }

    for (int i = 11; i < 21; i++) {
        multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_update_params(capsule, i, p, NP_2);
    }

    for (int i = 21; i < 22; i++) {
        multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_update_params(capsule, i, p, NP_3);
    }

    for (int i = 22; i < 32; i++) {
        multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_update_params(capsule, i, p, NP_4);
    }
    free(p);
}




/**
 * Internal function for multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create: step 5
 */
void multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_setup_nlp_in(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule, int N)
{
    assert(N == capsule->nlp_solver_plan->N);
    ocp_nlp_config* nlp_config = capsule->nlp_config;
    ocp_nlp_dims* nlp_dims = capsule->nlp_dims;

    int tmp_int = 0;

    /************************************************
    *  nlp_in
    ************************************************/
    ocp_nlp_in * nlp_in = capsule->nlp_in;


    // set up time_steps

double* time_steps = malloc(N*sizeof(double));
    time_steps[0] = 0.005;
    time_steps[1] = 0.005;
    time_steps[2] = 0.004999999999999999;
    time_steps[3] = 0.005000000000000001;
    time_steps[4] = 0.005000000000000001;
    time_steps[5] = 0.005000000000000001;
    time_steps[6] = 0.005000000000000001;
    time_steps[7] = 0.0049999999999999975;
    time_steps[8] = 0.0049999999999999975;
    time_steps[9] = 0.0049999999999999975;
    time_steps[10] = 1;
    time_steps[11] = 0.004999999999999893;
    time_steps[12] = 0.004999999999999893;
    time_steps[13] = 0.004999999999999893;
    time_steps[14] = 0.004999999999999893;
    time_steps[15] = 0.004999999999999893;
    time_steps[16] = 0.004999999999999893;
    time_steps[17] = 0.004999999999999893;
    time_steps[18] = 0.004999999999999893;
    time_steps[19] = 0.004999999999999893;
    time_steps[20] = 0.004999999999999893;
    time_steps[21] = 0.9999999999999998;
    time_steps[22] = 0.004999999999999893;
    time_steps[23] = 0.004999999999999893;
    time_steps[24] = 0.004999999999999893;
    time_steps[25] = 0.004999999999999893;
    time_steps[26] = 0.004999999999999893;
    time_steps[27] = 0.004999999999999893;
    time_steps[28] = 0.004999999999999893;
    time_steps[29] = 0.004999999999999893;
    time_steps[30] = 0.004999999999999893;
    time_steps[31] = 0.004999999999999893;
    for (int i = 0; i < N; i++)
    {
        ocp_nlp_in_set(nlp_config, nlp_dims, nlp_in, i, "Ts", &time_steps[i]);
    }
    free(time_steps);
    // set cost scaling
    double* cost_scaling = malloc((N+1)*sizeof(double));
    cost_scaling[0] = 0.005;
    cost_scaling[1] = 0.005;
    cost_scaling[2] = 0.005;
    cost_scaling[3] = 0.005;
    cost_scaling[4] = 0.005;
    cost_scaling[5] = 0.005;
    cost_scaling[6] = 0.005;
    cost_scaling[7] = 0.005;
    cost_scaling[8] = 0.005;
    cost_scaling[9] = 0.005;
    cost_scaling[10] = 1;
    cost_scaling[11] = 0.005;
    cost_scaling[12] = 0.005;
    cost_scaling[13] = 0.005;
    cost_scaling[14] = 0.005;
    cost_scaling[15] = 0.005;
    cost_scaling[16] = 0.005;
    cost_scaling[17] = 0.005;
    cost_scaling[18] = 0.005;
    cost_scaling[19] = 0.005;
    cost_scaling[20] = 0.005;
    cost_scaling[21] = 1;
    cost_scaling[22] = 0.005;
    cost_scaling[23] = 0.005;
    cost_scaling[24] = 0.005;
    cost_scaling[25] = 0.005;
    cost_scaling[26] = 0.005;
    cost_scaling[27] = 0.005;
    cost_scaling[28] = 0.005;
    cost_scaling[29] = 0.005;
    cost_scaling[30] = 0.005;
    cost_scaling[31] = 0.005;
    cost_scaling[32] = 1;
    for (int i = 0; i <= N; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "scaling", &cost_scaling[i]);
    }
    free(cost_scaling);

    /* INITIAL NODE */
    double* yref_0 = calloc(6, sizeof(double));
    // change only the non-zero elements:
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "yref", yref_0);
    free(yref_0);

    double* W_0 = calloc(6*6, sizeof(double));
    // change only the non-zero elements:
    W_0[0+(6) * 0] = 4;
    W_0[1+(6) * 1] = 4;
    W_0[2+(6) * 2] = 0.05;
    W_0[3+(6) * 3] = 0.05;
    W_0[4+(6) * 4] = 0.015;
    W_0[5+(6) * 5] = 0.015;
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, 0, "W", W_0);
    free(W_0);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nls_y_fun", &capsule->cost_y_0_fun);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, 0, "nls_y_fun_jac", &capsule->cost_y_0_fun_jac_ut_xt);




    // constraints at initial node
    // x0
    int nbx_0 = 16;
    int* idxbx0 = malloc(nbx_0 * sizeof(int));
    idxbx0[0] = 0;
    idxbx0[1] = 1;
    idxbx0[2] = 2;
    idxbx0[3] = 3;
    idxbx0[4] = 4;
    idxbx0[5] = 5;
    idxbx0[6] = 6;
    idxbx0[7] = 7;
    idxbx0[8] = 8;
    idxbx0[9] = 9;
    idxbx0[10] = 10;
    idxbx0[11] = 11;
    idxbx0[12] = 12;
    idxbx0[13] = 13;
    idxbx0[14] = 14;
    idxbx0[15] = 15;

    double* lubx0 = calloc(2*nbx_0, sizeof(double));
    double* lbx0 = lubx0;
    double* ubx0 = lubx0 + nbx_0;
    // change only the non-zero elements:

    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbx", idxbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", lbx0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", ubx0);
    free(idxbx0);
    free(lubx0);
    // idxbxe_0
    int* idxbxe_0 = malloc(16 * sizeof(int));
    idxbxe_0[0] = 0;
    idxbxe_0[1] = 1;
    idxbxe_0[2] = 2;
    idxbxe_0[3] = 3;
    idxbxe_0[4] = 4;
    idxbxe_0[5] = 5;
    idxbxe_0[6] = 6;
    idxbxe_0[7] = 7;
    idxbxe_0[8] = 8;
    idxbxe_0[9] = 9;
    idxbxe_0[10] = 10;
    idxbxe_0[11] = 11;
    idxbxe_0[12] = 12;
    idxbxe_0[13] = 13;
    idxbxe_0[14] = 14;
    idxbxe_0[15] = 15;
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "idxbxe", idxbxe_0);
    free(idxbxe_0);








    /* Path related delarations */
    int i_fun;

    // cost
    double* yref;
    double* Vx;
    double* Vu;
    double* Vz;
    double* W;

    // bounds on u
    int* idxbu;
    double* lubu;
    double* lbu;
    double* ubu;

    // bounds on x
    double* lubx;
    double* lbx;
    double* ubx;
    int* idxbx;

    // general linear constraints
    double* D;
    double* C;
    double* lug;
    double* lg;
    double* ug;

    // nonlinear constraints
    double* luh;
    double* lh;
    double* uh;
    double* luphi;
    double* lphi;
    double* uphi;

    // general slack related
    double* zlumem;
    double* Zl;
    double* Zu;
    double* zl;
    double* zu;

    // specific slack types
    int* idxsbx;
    double* lusbx;
    double* lsbx;
    double* usbx;

    int* idxsbu;
    double* lusbu;
    double* lsbu;
    double* usbu;

    int* idxsg;
    double* lusg;
    double* lsg;
    double* usg;

    int* idxsh;
    double* lush;
    double* lsh;
    double* ush;

    int* idxsphi;
    double* lusphi;
    double* lsphi;
    double* usphi;

    /*********************
     *  Phase 0
     * *******************/
    /**** Dynamics ****/
    for (int i = 0; i < 10; i++)
    {
        i_fun = i - 0;
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "impl_dae_fun", &capsule->impl_dae_fun_0[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i,
                                   "impl_dae_fun_jac_x_xdot_z", &capsule->impl_dae_fun_jac_x_xdot_z_0[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i,
                                   "impl_dae_jac_x_xdot_u", &capsule->impl_dae_jac_x_xdot_u_z_0[i_fun]);
    }
    for (int i = 1; i < 10; i++)
    {
        i_fun = i - 1;

        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun", &capsule->cost_y_fun_0[i_fun]);
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun_jac", &capsule->cost_y_fun_jac_ut_xt_0[i_fun]);
    }

    /**** Cost phase 0 ****/
    yref = calloc(6, sizeof(double));
    // change only the non-zero elements:

    for (int i = 1; i < 10; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "yref", yref);
    }
    free(yref);
    W = calloc(6*6, sizeof(double));
    // change only the non-zero elements:
    W[0+(6) * 0] = 4;
    W[1+(6) * 1] = 4;
    W[2+(6) * 2] = 0.05;
    W[3+(6) * 3] = 0.05;
    W[4+(6) * 4] = 0.015;
    W[5+(6) * 5] = 0.015;

    for (int i = 1; i < 10; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "W", W);
    }
    free(W);


    /**** Constraints phase 0 ****/


    /* constraints that are the same for initial and intermediate */












    // set up nonlinear constraints for stage 1 to N-1
    luh = calloc(2*1, sizeof(double));
    lh = luh;
    uh = luh + 1;
    lh[0] = -1000000000000000;
    uh[0] = 0.2;

    for (int i = 1; i < 10; i++)
    {
        i_fun = i - 1;
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun_jac",
                                      &capsule->nl_constr_h_fun_jac_0[i_fun]);
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun",
                                      &capsule->nl_constr_h_fun_0[i_fun]);
        
        
        
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lh", lh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "uh", uh);
    }
    free(luh);



    /*********************
     *  Phase 1
     * *******************/
    /**** Dynamics ****/
    for (int i = 10; i < 11; i++)
    {
        i_fun = i - 10;
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "disc_dyn_fun", &capsule->discr_dyn_phi_fun_1[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "disc_dyn_fun_jac",
                                   &capsule->discr_dyn_phi_fun_jac_ut_xt_1[i_fun]);
    }
    for (int i = 10; i < 11; i++)
    {
        i_fun = i - 10;

        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun", &capsule->cost_y_fun_1[i_fun]);
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun_jac", &capsule->cost_y_fun_jac_ut_xt_1[i_fun]);
    }

    /**** Cost phase 1 ****/
    yref = calloc(16, sizeof(double));
    // change only the non-zero elements:

    for (int i = 10; i < 11; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "yref", yref);
    }
    free(yref);
    W = calloc(16*16, sizeof(double));
    // change only the non-zero elements:

    for (int i = 10; i < 11; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "W", W);
    }
    free(W);


    /**** Constraints phase 1 ****/


    /* constraints that are the same for initial and intermediate */















    /*********************
     *  Phase 2
     * *******************/
    /**** Dynamics ****/
    for (int i = 11; i < 21; i++)
    {
        i_fun = i - 11;
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "impl_dae_fun", &capsule->impl_dae_fun_2[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i,
                                   "impl_dae_fun_jac_x_xdot_z", &capsule->impl_dae_fun_jac_x_xdot_z_2[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i,
                                   "impl_dae_jac_x_xdot_u", &capsule->impl_dae_jac_x_xdot_u_z_2[i_fun]);
    }
    for (int i = 11; i < 21; i++)
    {
        i_fun = i - 11;

        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun", &capsule->cost_y_fun_2[i_fun]);
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun_jac", &capsule->cost_y_fun_jac_ut_xt_2[i_fun]);
    }

    /**** Cost phase 2 ****/
    yref = calloc(6, sizeof(double));
    // change only the non-zero elements:

    for (int i = 11; i < 21; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "yref", yref);
    }
    free(yref);
    W = calloc(6*6, sizeof(double));
    // change only the non-zero elements:
    W[0+(6) * 0] = 4;
    W[1+(6) * 1] = 4;
    W[2+(6) * 2] = 0.05;
    W[3+(6) * 3] = 0.05;
    W[4+(6) * 4] = 0.015;
    W[5+(6) * 5] = 0.015;

    for (int i = 11; i < 21; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "W", W);
    }
    free(W);


    /**** Constraints phase 2 ****/


    /* constraints that are the same for initial and intermediate */












    // set up nonlinear constraints for stage 1 to N-1
    luh = calloc(2*1, sizeof(double));
    lh = luh;
    uh = luh + 1;
    lh[0] = -1000000000000000;
    uh[0] = 0.2;

    for (int i = 11; i < 21; i++)
    {
        i_fun = i - 11;
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun_jac",
                                      &capsule->nl_constr_h_fun_jac_2[i_fun]);
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun",
                                      &capsule->nl_constr_h_fun_2[i_fun]);
        
        
        
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lh", lh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "uh", uh);
    }
    free(luh);



    /*********************
     *  Phase 3
     * *******************/
    /**** Dynamics ****/
    for (int i = 21; i < 22; i++)
    {
        i_fun = i - 21;
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "disc_dyn_fun", &capsule->discr_dyn_phi_fun_3[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "disc_dyn_fun_jac",
                                   &capsule->discr_dyn_phi_fun_jac_ut_xt_3[i_fun]);
    }
    for (int i = 21; i < 22; i++)
    {
        i_fun = i - 21;

        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun", &capsule->cost_y_fun_3[i_fun]);
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun_jac", &capsule->cost_y_fun_jac_ut_xt_3[i_fun]);
    }

    /**** Cost phase 3 ****/
    yref = calloc(8, sizeof(double));
    // change only the non-zero elements:

    for (int i = 21; i < 22; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "yref", yref);
    }
    free(yref);
    W = calloc(8*8, sizeof(double));
    // change only the non-zero elements:

    for (int i = 21; i < 22; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "W", W);
    }
    free(W);


    /**** Constraints phase 3 ****/


    /* constraints that are the same for initial and intermediate */















    /*********************
     *  Phase 4
     * *******************/
    /**** Dynamics ****/
    for (int i = 22; i < 32; i++)
    {
        i_fun = i - 22;
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "expl_vde_forw", &capsule->expl_vde_forw_4[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "expl_vde_adj", &capsule->expl_vde_adj_4[i_fun]);
        ocp_nlp_dynamics_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "expl_ode_fun", &capsule->expl_ode_fun_4[i_fun]);
    }
    for (int i = 22; i < 32; i++)
    {
        i_fun = i - 22;

        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun", &capsule->cost_y_fun_4[i_fun]);
        ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nls_y_fun_jac", &capsule->cost_y_fun_jac_ut_xt_4[i_fun]);
    }

    /**** Cost phase 4 ****/
    yref = calloc(6, sizeof(double));
    // change only the non-zero elements:

    for (int i = 22; i < 32; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "yref", yref);
    }
    free(yref);
    W = calloc(6*6, sizeof(double));
    // change only the non-zero elements:
    W[0+(6) * 0] = 4;
    W[1+(6) * 1] = 4;
    W[2+(6) * 2] = 0.05;
    W[3+(6) * 3] = 0.05;
    W[4+(6) * 4] = 0.00001;
    W[5+(6) * 5] = 0.00001;

    for (int i = 22; i < 32; i++)
    {
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, i, "W", W);
    }
    free(W);


    /**** Constraints phase 4 ****/


    /* constraints that are the same for initial and intermediate */












    // set up nonlinear constraints for stage 1 to N-1
    luh = calloc(2*1, sizeof(double));
    lh = luh;
    uh = luh + 1;
    lh[0] = -1000000000000000;
    uh[0] = 0.2;

    for (int i = 22; i < 32; i++)
    {
        i_fun = i - 22;
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun_jac",
                                      &capsule->nl_constr_h_fun_jac_4[i_fun]);
        ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, i, "nl_constr_h_fun",
                                      &capsule->nl_constr_h_fun_4[i_fun]);
        
        
        
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "lh", lh);
        ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, i, "uh", uh);
    }
    free(luh);





    // TERMINAL node
    double* yref_e = calloc(4, sizeof(double));
    // change only the non-zero elements:
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "yref", yref_e);
    free(yref_e);

    double* W_e = calloc(4*4, sizeof(double));
    // change only the non-zero elements:
    W_e[0+(4) * 0] = 4;
    W_e[1+(4) * 1] = 4;
    W_e[2+(4) * 2] = 0.05;
    W_e[3+(4) * 3] = 0.05;
    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, N, "W", W_e);
    free(W_e);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nls_y_fun", &capsule->cost_y_e_fun);
    ocp_nlp_cost_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nls_y_fun_jac", &capsule->cost_y_e_fun_jac_ut_xt);









    /* terminal constraints */













    // set up nonlinear constraints for last stage
    double* luh_e = calloc(2*1, sizeof(double));
    double* lh_e = luh_e;
    double* uh_e = luh_e + 1;
    lh_e[0] = -1000000000000000;
    uh_e[0] = 0.2;

    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_fun_jac", &capsule->nl_constr_h_e_fun_jac);
    ocp_nlp_constraints_model_set_external_param_fun(nlp_config, nlp_dims, nlp_in, N, "nl_constr_h_fun", &capsule->nl_constr_h_e_fun);
    
    
    
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "lh", lh_e);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, N, "uh", uh_e);
    free(luh_e);
}


/**
 * Internal function for multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create: step 6
 */
void multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_set_opts(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule)
{
    const int N = capsule->nlp_solver_plan->N;
    ocp_nlp_config* nlp_config = capsule->nlp_config;
    void *nlp_opts = capsule->nlp_opts;

    // declare
    bool tmp_bool;
    int newton_iter_val;
    double newton_tol_val;
    /************************************************
    *  opts
    ************************************************/



    int fixed_hess = 0;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "fixed_hess", &fixed_hess);
    double globalization_fixed_step_length = 1;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "globalization_fixed_step_length", &globalization_fixed_step_length);





    int with_solution_sens_wrt_params = false;
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "with_solution_sens_wrt_params", &with_solution_sens_wrt_params);

    int with_value_sens_wrt_params = false;
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "with_value_sens_wrt_params", &with_value_sens_wrt_params);

    int globalization_full_step_dual = 0;
    ocp_nlp_solver_opts_set(nlp_config, capsule->nlp_opts, "globalization_full_step_dual", &globalization_full_step_dual);

    double levenberg_marquardt = 0;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "levenberg_marquardt", &levenberg_marquardt);

    /* options QP solver */

    bool store_iterates = false;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "store_iterates", &store_iterates);

    int nlp_solver_ext_qp_res = 0;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "ext_qp_res", &nlp_solver_ext_qp_res);
    // set HPIPM mode: should be done before setting other QP solver options
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_hpipm_mode", "BALANCE");


    // set SQP specific options
    double nlp_solver_tol_stat = 0.0000001;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_stat", &nlp_solver_tol_stat);

    double nlp_solver_tol_eq = 0.0000001;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_eq", &nlp_solver_tol_eq);

    double nlp_solver_tol_ineq = 0.0000001;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_ineq", &nlp_solver_tol_ineq);

    double nlp_solver_tol_comp = 0.0000001;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "tol_comp", &nlp_solver_tol_comp);

    int nlp_solver_max_iter = 150;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "max_iter", &nlp_solver_max_iter);

    // set options for adaptive Levenberg-Marquardt Update
    bool with_adaptive_levenberg_marquardt = false;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "with_adaptive_levenberg_marquardt", &with_adaptive_levenberg_marquardt);

    double adaptive_levenberg_marquardt_lam = 5;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "adaptive_levenberg_marquardt_lam", &adaptive_levenberg_marquardt_lam);

    double adaptive_levenberg_marquardt_mu_min = 0.0000000000000001;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "adaptive_levenberg_marquardt_mu_min", &adaptive_levenberg_marquardt_mu_min);

    double adaptive_levenberg_marquardt_mu0 = 0.001;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "adaptive_levenberg_marquardt_mu0", &adaptive_levenberg_marquardt_mu0);

    bool eval_residual_at_max_iter = false;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "eval_residual_at_max_iter", &eval_residual_at_max_iter);

    int qp_solver_iter_max = 50;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "qp_iter_max", &qp_solver_iter_max);



    int print_level = 0;
    ocp_nlp_solver_opts_set(nlp_config, nlp_opts, "print_level", &print_level);



    /* Stage varying options */
    int ext_cost_num_hess;
    bool output_z_val = true;
    bool sens_algebraic_val = true;
    sim_collocation_type collocation_type;

    // set up sim_method_num_stages
    int* sim_method_num_stages = malloc(N*sizeof(int));
    sim_method_num_stages[0] = 4;
    sim_method_num_stages[1] = 4;
    sim_method_num_stages[2] = 4;
    sim_method_num_stages[3] = 4;
    sim_method_num_stages[4] = 4;
    sim_method_num_stages[5] = 4;
    sim_method_num_stages[6] = 4;
    sim_method_num_stages[7] = 4;
    sim_method_num_stages[8] = 4;
    sim_method_num_stages[9] = 4;
    sim_method_num_stages[10] = 4;
    sim_method_num_stages[11] = 4;
    sim_method_num_stages[12] = 4;
    sim_method_num_stages[13] = 4;
    sim_method_num_stages[14] = 4;
    sim_method_num_stages[15] = 4;
    sim_method_num_stages[16] = 4;
    sim_method_num_stages[17] = 4;
    sim_method_num_stages[18] = 4;
    sim_method_num_stages[19] = 4;
    sim_method_num_stages[20] = 4;
    sim_method_num_stages[21] = 4;
    sim_method_num_stages[22] = 4;
    sim_method_num_stages[23] = 4;
    sim_method_num_stages[24] = 4;
    sim_method_num_stages[25] = 4;
    sim_method_num_stages[26] = 4;
    sim_method_num_stages[27] = 4;
    sim_method_num_stages[28] = 4;
    sim_method_num_stages[29] = 4;
    sim_method_num_stages[30] = 4;
    sim_method_num_stages[31] = 4;

    // set up sim_method_num_steps
    int* sim_method_num_steps = malloc(N*sizeof(int));
    sim_method_num_steps[0] = 1;
    sim_method_num_steps[1] = 1;
    sim_method_num_steps[2] = 1;
    sim_method_num_steps[3] = 1;
    sim_method_num_steps[4] = 1;
    sim_method_num_steps[5] = 1;
    sim_method_num_steps[6] = 1;
    sim_method_num_steps[7] = 1;
    sim_method_num_steps[8] = 1;
    sim_method_num_steps[9] = 1;
    sim_method_num_steps[10] = 1;
    sim_method_num_steps[11] = 1;
    sim_method_num_steps[12] = 1;
    sim_method_num_steps[13] = 1;
    sim_method_num_steps[14] = 1;
    sim_method_num_steps[15] = 1;
    sim_method_num_steps[16] = 1;
    sim_method_num_steps[17] = 1;
    sim_method_num_steps[18] = 1;
    sim_method_num_steps[19] = 1;
    sim_method_num_steps[20] = 1;
    sim_method_num_steps[21] = 1;
    sim_method_num_steps[22] = 1;
    sim_method_num_steps[23] = 1;
    sim_method_num_steps[24] = 1;
    sim_method_num_steps[25] = 1;
    sim_method_num_steps[26] = 1;
    sim_method_num_steps[27] = 1;
    sim_method_num_steps[28] = 1;
    sim_method_num_steps[29] = 1;
    sim_method_num_steps[30] = 1;
    sim_method_num_steps[31] = 1;

    // set up sim_method_jac_reuse
    bool* sim_method_jac_reuse = malloc(N*sizeof(bool));
    sim_method_jac_reuse[0] = (bool)0;
    sim_method_jac_reuse[1] = (bool)0;
    sim_method_jac_reuse[2] = (bool)0;
    sim_method_jac_reuse[3] = (bool)0;
    sim_method_jac_reuse[4] = (bool)0;
    sim_method_jac_reuse[5] = (bool)0;
    sim_method_jac_reuse[6] = (bool)0;
    sim_method_jac_reuse[7] = (bool)0;
    sim_method_jac_reuse[8] = (bool)0;
    sim_method_jac_reuse[9] = (bool)0;
    sim_method_jac_reuse[10] = (bool)0;
    sim_method_jac_reuse[11] = (bool)0;
    sim_method_jac_reuse[12] = (bool)0;
    sim_method_jac_reuse[13] = (bool)0;
    sim_method_jac_reuse[14] = (bool)0;
    sim_method_jac_reuse[15] = (bool)0;
    sim_method_jac_reuse[16] = (bool)0;
    sim_method_jac_reuse[17] = (bool)0;
    sim_method_jac_reuse[18] = (bool)0;
    sim_method_jac_reuse[19] = (bool)0;
    sim_method_jac_reuse[20] = (bool)0;
    sim_method_jac_reuse[21] = (bool)0;
    sim_method_jac_reuse[22] = (bool)0;
    sim_method_jac_reuse[23] = (bool)0;
    sim_method_jac_reuse[24] = (bool)0;
    sim_method_jac_reuse[25] = (bool)0;
    sim_method_jac_reuse[26] = (bool)0;
    sim_method_jac_reuse[27] = (bool)0;
    sim_method_jac_reuse[28] = (bool)0;
    sim_method_jac_reuse[29] = (bool)0;
    sim_method_jac_reuse[30] = (bool)0;
    sim_method_jac_reuse[31] = (bool)0;

    // set collocation type (relevant for implicit integrators)
    collocation_type = GAUSS_LEGENDRE;
    for (int i = 0; i < 10; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_collocation_type", &collocation_type);

    // sim_method_newton_iter
    newton_iter_val = 3;
    for (int i = 0; i < 10; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_newton_iter", &newton_iter_val);

    // sim_method_newton_tol
    newton_tol_val = 0;
    for (int i = 0; i < 10; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_newton_tol", &newton_tol_val);

    // possibly varying: sim_method_num_steps, sim_method_num_stages, sim_method_jac_reuse
    for (int i = 0; i < 10; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_steps", &sim_method_num_steps[i]);

    for (int i = 0; i < 10; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_stages", &sim_method_num_stages[i]);

    for (int i = 0; i < 10; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_jac_reuse", &sim_method_jac_reuse[i]);

    // set collocation type (relevant for implicit integrators)
    collocation_type = GAUSS_LEGENDRE;
    for (int i = 11; i < 21; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_collocation_type", &collocation_type);

    // sim_method_newton_iter
    newton_iter_val = 3;
    for (int i = 11; i < 21; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_newton_iter", &newton_iter_val);

    // sim_method_newton_tol
    newton_tol_val = 0;
    for (int i = 11; i < 21; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_newton_tol", &newton_tol_val);

    // possibly varying: sim_method_num_steps, sim_method_num_stages, sim_method_jac_reuse
    for (int i = 11; i < 21; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_steps", &sim_method_num_steps[i]);

    for (int i = 11; i < 21; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_stages", &sim_method_num_stages[i]);

    for (int i = 11; i < 21; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_jac_reuse", &sim_method_jac_reuse[i]);

    // set collocation type (relevant for implicit integrators)
    collocation_type = GAUSS_LEGENDRE;
    for (int i = 22; i < 32; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_collocation_type", &collocation_type);

    // sim_method_newton_iter
    newton_iter_val = 3;
    for (int i = 22; i < 32; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_newton_iter", &newton_iter_val);

    // sim_method_newton_tol
    newton_tol_val = 0;
    for (int i = 22; i < 32; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_newton_tol", &newton_tol_val);

    // possibly varying: sim_method_num_steps, sim_method_num_stages, sim_method_jac_reuse
    for (int i = 22; i < 32; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_steps", &sim_method_num_steps[i]);

    for (int i = 22; i < 32; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_num_stages", &sim_method_num_stages[i]);

    for (int i = 22; i < 32; i++)
        ocp_nlp_solver_opts_set_at_stage(nlp_config, nlp_opts, i, "dynamics_jac_reuse", &sim_method_jac_reuse[i]);

    // free arrays
    free(sim_method_num_steps);
    free(sim_method_num_stages);
    free(sim_method_jac_reuse);
}

/**
 * Internal function for multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create: step 7
 */
void multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_set_nlp_out(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule)
{
    const int N = capsule->nlp_solver_plan->N;
    ocp_nlp_config* nlp_config = capsule->nlp_config;
    ocp_nlp_dims* nlp_dims = capsule->nlp_dims;
    ocp_nlp_out* nlp_out = capsule->nlp_out;

    int nx_max = 16;
    int nu_max = 2;

    // initialize primal solution
    double* xu0 = calloc(nx_max+nu_max, sizeof(double));
    double* x0 = xu0;

    // initialize with x0


    double* u0 = xu0 + nx_max;

    for (int i = 0; i < N; i++)
    {
        // x0
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", x0);
        // u0
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", u0);
    }
    ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, N, "x", x0);
    free(xu0);
}



/**
 * Internal function for multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create: step 9
 */
int multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_precompute(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule) {
    int status = ocp_nlp_precompute(capsule->nlp_solver, capsule->nlp_in, capsule->nlp_out);

    if (status != ACADOS_SUCCESS) {
        printf("\nocp_nlp_precompute failed!\n\n");
        exit(1);
    }

    return status;
}






int multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_with_discretization(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule, int N, double* new_time_steps)
{
    // If N does not match the number of shooting intervals used for code generation, new_time_steps must be given.
    if (new_time_steps) {
        fprintf(stderr, "multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_with_discretization: new_time_steps should be NULL " \
            "for multi-phase solver!\n", \
             N, MULTIPHASE_OCP_10_1_10_1_10_20250606_034423_645279_N);
        return 1;
    }

    // 1) create and set nlp_solver_plan; create nlp_config
    capsule->nlp_solver_plan = ocp_nlp_plan_create(N);
    multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_set_plan(capsule->nlp_solver_plan, N);
    capsule->nlp_config = ocp_nlp_config_create(*capsule->nlp_solver_plan);

    // 2) create and set dimensions
    capsule->nlp_dims = multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_setup_dimensions(capsule);

    // 3) create and set nlp_opts
    capsule->nlp_opts = ocp_nlp_solver_opts_create(capsule->nlp_config, capsule->nlp_dims);
    multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_set_opts(capsule);

    // 4) create nlp_in
    capsule->nlp_in = ocp_nlp_in_create(capsule->nlp_config, capsule->nlp_dims);

    // 5) set default parameters in functions
    multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_setup_functions(capsule);
    multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_setup_nlp_in(capsule, N);
    multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_set_default_parameters(capsule);

    // 6) create solver
    capsule->nlp_solver = ocp_nlp_solver_create(capsule->nlp_config, capsule->nlp_dims, capsule->nlp_opts, capsule->nlp_in);

    // 7) create and set nlp_out
    // 7.1) nlp_out
    capsule->nlp_out = ocp_nlp_out_create(capsule->nlp_config, capsule->nlp_dims);
    // 7.2) sens_out
    capsule->sens_out = ocp_nlp_out_create(capsule->nlp_config, capsule->nlp_dims);
    multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_set_nlp_out(capsule);

    // 8) do precomputations
    int status = multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_create_precompute(capsule);

    return status;
}




int multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_reset(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule, int reset_qp_solver_mem)
{
    // set initialization to all zeros

    const int N = capsule->nlp_solver_plan->N;
    ocp_nlp_config* nlp_config = capsule->nlp_config;
    ocp_nlp_dims* nlp_dims = capsule->nlp_dims;
    ocp_nlp_out* nlp_out = capsule->nlp_out;
    ocp_nlp_in* nlp_in = capsule->nlp_in;
    ocp_nlp_solver* nlp_solver = capsule->nlp_solver;

    double* buffer = calloc(88, sizeof(double));
    // Reset stage 0
    for (int i = 0; i < 10; i++)
    {
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "sl", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "su", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "lam", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "z", buffer);
        if (i<N)
        {
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "pi", buffer);
            ocp_nlp_set(nlp_solver, i, "xdot_guess", buffer);
            ocp_nlp_set(nlp_solver, i, "z_guess", buffer);
        }
    }
    // Reset stage 1
    for (int i = 10; i < 11; i++)
    {
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "sl", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "su", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "lam", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "z", buffer);
        if (i<N)
        {
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "pi", buffer);
        }
    }
    // Reset stage 2
    for (int i = 11; i < 21; i++)
    {
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "sl", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "su", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "lam", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "z", buffer);
        if (i<N)
        {
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "pi", buffer);
            ocp_nlp_set(nlp_solver, i, "xdot_guess", buffer);
            ocp_nlp_set(nlp_solver, i, "z_guess", buffer);
        }
    }
    // Reset stage 3
    for (int i = 21; i < 22; i++)
    {
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "sl", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "su", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "lam", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "z", buffer);
        if (i<N)
        {
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "pi", buffer);
        }
    }
    // Reset stage 4
    for (int i = 22; i < 32; i++)
    {
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "x", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "u", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "sl", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "su", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "lam", buffer);
        ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "z", buffer);
        if (i<N)
        {
            ocp_nlp_out_set(nlp_config, nlp_dims, nlp_out, i, "pi", buffer);
        }
    }

    free(buffer);
    return 0;
}




int multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_update_params(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule, int stage, double *p, int np)
{
    int solver_status = 0;
    if (stage >= 0 && stage < 10)
    {
        if (NP_0 != np)
        {
            printf("acados_update_params: trying to set %i parameters at stage %i."
                " Parameters should be of length %i. Exiting.\n", np, stage, NP_0);
            exit(1);
        }
    }
    
    if (stage >= 10 && stage < 11)
    {
        if (NP_1 != np)
        {
            printf("acados_update_params: trying to set %i parameters at stage %i."
                " Parameters should be of length %i. Exiting.\n", np, stage, NP_1);
            exit(1);
        }
    }
    
    if (stage >= 11 && stage < 21)
    {
        if (NP_2 != np)
        {
            printf("acados_update_params: trying to set %i parameters at stage %i."
                " Parameters should be of length %i. Exiting.\n", np, stage, NP_2);
            exit(1);
        }
    }
    
    if (stage >= 21 && stage < 22)
    {
        if (NP_3 != np)
        {
            printf("acados_update_params: trying to set %i parameters at stage %i."
                " Parameters should be of length %i. Exiting.\n", np, stage, NP_3);
            exit(1);
        }
    }
    
    if (stage >= 22 && stage < 32)
    {
        if (NP_4 != np)
        {
            printf("acados_update_params: trying to set %i parameters at stage %i."
                " Parameters should be of length %i. Exiting.\n", np, stage, NP_4);
            exit(1);
        }
    }
    
    ocp_nlp_in_set(capsule->nlp_config, capsule->nlp_dims, capsule->nlp_in, stage, "parameter_values", p);

    return solver_status;
}


int multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_update_params_sparse(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule * capsule, int stage, int *idx, double *p, int n_update)
{
    ocp_nlp_in_set_params_sparse(capsule->nlp_config, capsule->nlp_dims, capsule->nlp_in, stage, idx, p, n_update);

    return 0;
}


int multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_set_p_global_and_precompute_dependencies(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule, double* data, int data_len)
{

    printf("p_global is not defined, multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_set_p_global_and_precompute_dependencies does nothing.\n");
}



int multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_solve(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule)
{
    // solve NLP
    int solver_status = ocp_nlp_solve(capsule->nlp_solver, capsule->nlp_in, capsule->nlp_out);

    return solver_status;
}




void multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_print_stats(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule)
{
    int sqp_iter, stat_m, stat_n, tmp_int;
    ocp_nlp_get(capsule->nlp_solver, "sqp_iter", &sqp_iter);
    ocp_nlp_get(capsule->nlp_solver, "stat_n", &stat_n);
    ocp_nlp_get(capsule->nlp_solver, "stat_m", &stat_m);


    double stat[1800];
    ocp_nlp_get(capsule->nlp_solver, "statistics", stat);

    int nrow = sqp_iter+1 < stat_m ? sqp_iter+1 : stat_m;

    printf("iter\tres_stat\tres_eq\t\tres_ineq\tres_comp\tqp_stat\tqp_iter\talpha");
    if (stat_n > 8)
        printf("\t\tqp_res_stat\tqp_res_eq\tqp_res_ineq\tqp_res_comp");
    printf("\n");

    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < stat_n + 1; j++)
        {
            if (j == 0 || j == 5 || j == 6)
            {
                tmp_int = (int) stat[i + j * nrow];
                printf("%d\t", tmp_int);
            }
            else
            {
                printf("%e\t", stat[i + j * nrow]);
            }
        }
        printf("\n");
    }

}




int multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_free(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule)
{
    // before destroying, keep some info
    const int N = capsule->nlp_solver_plan->N;
    // free memory
    ocp_nlp_solver_opts_destroy(capsule->nlp_opts);
    ocp_nlp_in_destroy(capsule->nlp_in);
    ocp_nlp_out_destroy(capsule->nlp_out);
    ocp_nlp_out_destroy(capsule->sens_out);
    ocp_nlp_solver_destroy(capsule->nlp_solver);
    ocp_nlp_dims_destroy(capsule->nlp_dims);
    ocp_nlp_config_destroy(capsule->nlp_config);
    ocp_nlp_plan_destroy(capsule->nlp_solver_plan);

    /* free external function */
    // initial node
    external_function_external_param_casadi_free(&capsule->cost_y_0_fun);
    external_function_external_param_casadi_free(&capsule->cost_y_0_fun_jac_ut_xt);
    /* Path phase {jj} */
    // dynamics
    for (int i_fun = 0; i_fun < 10; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->impl_dae_fun_0[i_fun]);
        external_function_external_param_casadi_free(&capsule->impl_dae_fun_jac_x_xdot_z_0[i_fun]);
        external_function_external_param_casadi_free(&capsule->impl_dae_jac_x_xdot_u_z_0[i_fun]);
    }
    free(capsule->impl_dae_fun_0);
    free(capsule->impl_dae_fun_jac_x_xdot_z_0);
    free(capsule->impl_dae_jac_x_xdot_u_z_0);
    for (int i_fun = 0; i_fun < 9; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->cost_y_fun_0[i_fun]);
        external_function_external_param_casadi_free(&capsule->cost_y_fun_jac_ut_xt_0[i_fun]);
    }
    free(capsule->cost_y_fun_0);
    free(capsule->cost_y_fun_jac_ut_xt_0);

    // constraints
    for (int i_fun = 0; i_fun < 9; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->nl_constr_h_fun_jac_0[i_fun]);
        external_function_external_param_casadi_free(&capsule->nl_constr_h_fun_0[i_fun]);
    }
    free(capsule->nl_constr_h_fun_jac_0);
    free(capsule->nl_constr_h_fun_0);
    /* Path phase {jj} */
    // dynamics
    for (int i_fun = 0; i_fun < 1; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->discr_dyn_phi_fun_1[i_fun]);
        external_function_external_param_casadi_free(&capsule->discr_dyn_phi_fun_jac_ut_xt_1[i_fun]);
    }
    free(capsule->discr_dyn_phi_fun_1);
    free(capsule->discr_dyn_phi_fun_jac_ut_xt_1);
    for (int i_fun = 0; i_fun < 1; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->cost_y_fun_1[i_fun]);
        external_function_external_param_casadi_free(&capsule->cost_y_fun_jac_ut_xt_1[i_fun]);
    }
    free(capsule->cost_y_fun_1);
    free(capsule->cost_y_fun_jac_ut_xt_1);

    // constraints
    /* Path phase {jj} */
    // dynamics
    for (int i_fun = 0; i_fun < 10; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->impl_dae_fun_2[i_fun]);
        external_function_external_param_casadi_free(&capsule->impl_dae_fun_jac_x_xdot_z_2[i_fun]);
        external_function_external_param_casadi_free(&capsule->impl_dae_jac_x_xdot_u_z_2[i_fun]);
    }
    free(capsule->impl_dae_fun_2);
    free(capsule->impl_dae_fun_jac_x_xdot_z_2);
    free(capsule->impl_dae_jac_x_xdot_u_z_2);
    for (int i_fun = 0; i_fun < 10; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->cost_y_fun_2[i_fun]);
        external_function_external_param_casadi_free(&capsule->cost_y_fun_jac_ut_xt_2[i_fun]);
    }
    free(capsule->cost_y_fun_2);
    free(capsule->cost_y_fun_jac_ut_xt_2);

    // constraints
    for (int i_fun = 0; i_fun < 10; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->nl_constr_h_fun_jac_2[i_fun]);
        external_function_external_param_casadi_free(&capsule->nl_constr_h_fun_2[i_fun]);
    }
    free(capsule->nl_constr_h_fun_jac_2);
    free(capsule->nl_constr_h_fun_2);
    /* Path phase {jj} */
    // dynamics
    for (int i_fun = 0; i_fun < 1; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->discr_dyn_phi_fun_3[i_fun]);
        external_function_external_param_casadi_free(&capsule->discr_dyn_phi_fun_jac_ut_xt_3[i_fun]);
    }
    free(capsule->discr_dyn_phi_fun_3);
    free(capsule->discr_dyn_phi_fun_jac_ut_xt_3);
    for (int i_fun = 0; i_fun < 1; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->cost_y_fun_3[i_fun]);
        external_function_external_param_casadi_free(&capsule->cost_y_fun_jac_ut_xt_3[i_fun]);
    }
    free(capsule->cost_y_fun_3);
    free(capsule->cost_y_fun_jac_ut_xt_3);

    // constraints
    /* Path phase {jj} */
    // dynamics
    for (int i_fun = 0; i_fun < 10; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->expl_vde_forw_4[i_fun]);
        external_function_external_param_casadi_free(&capsule->expl_vde_adj_4[i_fun]);
        external_function_external_param_casadi_free(&capsule->expl_ode_fun_4[i_fun]);
    }
    free(capsule->expl_vde_forw_4);
    free(capsule->expl_vde_adj_4);
    free(capsule->expl_ode_fun_4);
    for (int i_fun = 0; i_fun < 10; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->cost_y_fun_4[i_fun]);
        external_function_external_param_casadi_free(&capsule->cost_y_fun_jac_ut_xt_4[i_fun]);
    }
    free(capsule->cost_y_fun_4);
    free(capsule->cost_y_fun_jac_ut_xt_4);

    // constraints
    for (int i_fun = 0; i_fun < 10; i_fun++)
    {
        external_function_external_param_casadi_free(&capsule->nl_constr_h_fun_jac_4[i_fun]);
        external_function_external_param_casadi_free(&capsule->nl_constr_h_fun_4[i_fun]);
    }
    free(capsule->nl_constr_h_fun_jac_4);
    free(capsule->nl_constr_h_fun_4);


    /* Terminal node */
    external_function_external_param_casadi_free(&capsule->nl_constr_h_e_fun_jac);
    external_function_external_param_casadi_free(&capsule->nl_constr_h_e_fun);
    external_function_external_param_casadi_free(&capsule->cost_y_e_fun);
    external_function_external_param_casadi_free(&capsule->cost_y_e_fun_jac_ut_xt);



    return 0;
}




int multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_custom_update(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule, double* data, int data_len)
{
    (void)capsule;
    (void)data;
    (void)data_len;
    printf("\ndummy function that can be called in between solver calls to update parameters or numerical data efficiently in C.\n");
    printf("nothing set yet..\n");
    return 1;

}


ocp_nlp_in *multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_get_nlp_in(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule) { return capsule->nlp_in; }
ocp_nlp_out *multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_get_nlp_out(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule) { return capsule->nlp_out; }
ocp_nlp_out *multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_get_sens_out(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule) { return capsule->sens_out; }
ocp_nlp_solver *multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_get_nlp_solver(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule) { return capsule->nlp_solver; }
ocp_nlp_config *multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_get_nlp_config(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule) { return capsule->nlp_config; }
void *multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_get_nlp_opts(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule) { return capsule->nlp_opts; }
ocp_nlp_dims *multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_get_nlp_dims(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule) { return capsule->nlp_dims; }
ocp_nlp_plan_t *multiphase_ocp_10_1_10_1_10_20250606_034423_645279_acados_get_nlp_plan(multiphase_ocp_10_1_10_1_10_20250606_034423_645279_solver_capsule* capsule) { return capsule->nlp_solver_plan; }
