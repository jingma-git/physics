#ifndef SIM_linear_tri2dstvk_dq2_H
#define SIM_linear_tri2dstvk_dq2_H

#include <Eigen/Dense>
#include <EigenTypes.h>

#include <d2psi_stvk_dF2.h>

namespace sim {

template<typename HessianType, typename DefoType, typename DerivedV, typename Scalar, typename ParamType>
void linear_tri2d_stvk_dq2(Eigen::DenseBase<HessianType> &out, const Eigen::MatrixBase<DerivedV> &q, const Eigen::Ref<const Eigen::RowVectorXi> element,  
                                    const Eigen::MatrixBase<DefoType> &dphidX, const Eigen::MatrixBase<ParamType> &params, const Scalar &volume);

}

#ifndef SIM_STATIC_LIBRARY
#   include <../src/linear_tri2d_stvk_dq2.cpp>
#endif

#endif