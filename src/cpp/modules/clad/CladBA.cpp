// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.

#include "CladBA.h"
 
// This function must be called before any other function.
void CladBA::prepare(BAInput&& input)
{
    this->input = input;
    result = {
        std::vector<double>(2 * this->input.p),
        std::vector<double>(this->input.p),
        BASparseMat(this->input.n, this->input.m, this->input.p)
    };

    reproj_err_d = std::vector<double>(2 * (BA_NCAMPARAMS + 3 + 1));
    reproj_err_d_row = std::vector<double>(BA_NCAMPARAMS + 3 + 1);
}



BAOutput CladBA::output()
{
    return result;
}



void CladBA::calculate_objective(int times)
{
    for (int i = 0; i < times; i++)
    {
        ba_objective(
            input.n,
            input.m,
            input.p,
            input.cams.data(),
            input.X.data(),
            input.w.data(),
            input.obs.data(),
            input.feats.data(),
            result.reproj_err.data(),
            result.w_err.data()
        );
    }
}



void CladBA::calculate_jacobian(int times)
{
    for (int i = 0; i < times; i++)
    {
        result.J.clear();
        calculate_reproj_error_jacobian_part();
        calculate_weight_error_jacobian_part();
    }
}



void CladBA::calculate_reproj_error_jacobian_part()
{
    double d_feat[2];

    double* cam_gradient_part = reproj_err_d_row.data();
    double* x_gradient_part = reproj_err_d_row.data() + BA_NCAMPARAMS;
    double* weight_gradient_part = reproj_err_d_row.data() + BA_NCAMPARAMS + 3;

    for (int i = 0; i < input.p; i++)
    {
        int camIdx = input.obs[2 * i + 0];
        int ptIdx = input.obs[2 * i + 1];

        // calculate first row
        compute_reproj_error0_grad(
            &input.cams[camIdx * BA_NCAMPARAMS],
            &input.X[ptIdx * 3],
            input.w[i],
            &input.feats[i * 2],
            clad::array_ref<double>(cam_gradient_part, 11),
            clad::array_ref<double>(x_gradient_part, 3),
            weight_gradient_part,
            clad::array_ref<double>(d_feat, 2)
        );

        // fill first row elements
        for (int j = 0; j < BA_NCAMPARAMS + 3 + 1; j++)
        {
            reproj_err_d[2 * j] = reproj_err_d_row[j];
        }

        // calculate second row
        compute_reproj_error1_grad(
            &input.cams[camIdx * BA_NCAMPARAMS],
            &input.X[ptIdx * 3],
            input.w[i],
            &input.feats[i * 2],
            clad::array_ref<double>(cam_gradient_part, 11),
            clad::array_ref<double>(x_gradient_part, 3),
            weight_gradient_part,
            clad::array_ref<double>(d_feat, 2)
        );

        // fill second row elements
        for (int j = 0; j < BA_NCAMPARAMS + 3 + 1; j++)
        {
            reproj_err_d[2 * j + 1] = reproj_err_d_row[j];
        }

        result.J.insert_reproj_err_block(i, camIdx, ptIdx, reproj_err_d.data());
    }
}



void CladBA::calculate_weight_error_jacobian_part()
{
    for (int j = 0; j < input.p; j++)
    {
        double wb;              // stores calculated derivative

        double err = 0.0;       // stores fictive result
                                // (Clad doesn't calculate an original function in reverse mode)

        double errb = 1.0;      // stores dY
                                // (equals to 1.0 for derivative calculation)

        compute_zach_weight_error_grad(input.w[j], &wb);
        result.J.insert_w_err_block(j, wb);
    }
}



extern "C" DLL_PUBLIC ITest<BAInput, BAOutput>* get_ba_test()
{
    return new CladBA();
}