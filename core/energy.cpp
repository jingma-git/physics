#include "energy.h"
#include "util.h"
#include "macros.h"
#include "massmatrix.h"

#ifdef DEBUG
#include <iostream>
using namespace std;
#endif

namespace egl
{

    extern "C"
    {
        void tet_linear_(double *val, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);
        void tet_linear_jac_(double *jac, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);
        void tet_linear_hes_(double *hes, const double *x, const double *Dm, const double *vol, const double *lam, const double *miu);

        void bw98_stretch_(double *val, const double *x, const double *invUV, const double *area);
        void bw98_stretch_jac_(double *jac, const double *x, const double *invUV, const double *area);
        void bw98_stretch_hes_(double *hes, const double *x, const double *invUV, const double *area);

        void bw98_shear_(double *val, const double *x, const double *invUV, const double *area);
        void bw98_shear_jac_(double *jac, const double *x, const double *invUV, const double *area);
        void bw98_shear_hes_(double *hes, const double *x, const double *invUV, const double *area);

        void calc_edge_length_(double *val, const double *x);
        void calc_edge_length_jac_(double *jac, const double *x);
        void calc_edge_length_hes_(double *hes, const double *x);

        void calc_dih_angle_(double *val, const double *x);
        void calc_dih_angle_jac_(double *jac, const double *x);
        void calc_dih_angle_hes_(double *hes, const double *x);

        void surf_bending_(double *val, const double *x, const double *d, const double *l, const double *area);
        void surf_bending_jac_(double *jac, const double *x, const double *d, const double *l, const double *area);
        void surf_bending_hes_(double *hes, const double *x, const double *d, const double *l, const double *area);
    }

    momentum_potential_imp_euler::momentum_potential_imp_euler(const matd_t &nods,
                                                               const mati_t &cell,
                                                               const double rho,
                                                               const double h,
                                                               const double w)
        : rho_(rho),
          h_(h),
          w_(w),
          dof_(nods.size())
    {
        using namespace Eigen;
        calc_mass_matrix(nods, cell, rho, nods.rows(), M_, false);
        xn_ = Map<const VectorXd>(nods.data(), dof_);
        vn_.setZero(dof_);
    }

    int momentum_potential_imp_euler::Val(const double *x, double *val) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        using namespace Eigen;
        Map<const VectorXd> X(x, dof_);
        VectorXd dv = (X - xn_) / h_ - vn_;
        double momentum_val = w_ * 0.5 * dv.dot(M_ * dv);
        *val += momentum_val;
        return 0;
    }

    int momentum_potential_imp_euler::Gra(const double *x, double *gra) const
    {
        using namespace Eigen;
        Map<const VectorXd> X(x, dof_);
        Map<VectorXd> G(gra, dof_);
        G += w_ * M_ * ((X - xn_) / h_ - vn_) / h_;
        return 0;
    }

    int momentum_potential_imp_euler::Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        using namespace Eigen;
        Map<const VectorXd> X(x, dof_);
        const double coeff = w_ / (h_ * h_);
        double momentum_hess = 0;
        for (size_t j = 0; j < M_.outerSize(); ++j)
        {
            for (SparseMatrix<double>::InnerIterator it(M_, j); it; ++it)
            {
                hes->push_back(Triplet<double>(it.row(), it.col(), coeff * it.value()));
                momentum_hess += it.value();
            }
        }
        // cout << "momentum_hess=" << momentum_hess << endl;
        return 0;
    }

    void momentum_potential_imp_euler::Init(const double *x0, const double *v0)
    {
    }

    void momentum_potential_imp_euler::Update(const double *x)
    {
        using namespace Eigen;
        Map<const VectorXd> X(x, dof_);
        vn_ = (X - xn_) / h_;
        xn_ = X;
    }

    double momentum_potential_imp_euler::QueryKineticEnergy() const
    {
    }

    elastic_potential::elastic_potential(const matd_t &nods,
                                         const mati_t &tets,
                                         Material type,
                                         const double Ym,
                                         const double Pr,
                                         const double w) : rest_(nods),
                                                           tets_(tets),
                                                           dof_(nods.size()),
                                                           type_(type),
                                                           w_(w)
    {
        using namespace Eigen;
        vol_.resize(tets_.cols());
        Dm_.resize(9, tets_.cols());

        for (size_t i = 0; i < tets.cols(); ++i)
        {
            matd_t dm(3, 3);
            ele_coord(nods, tets_.col(i), dm);
            vol_[i] = fabs(dm.determinant()) / 6.0;
            Eigen::Map<Matrix3d>(&Dm_(0, i)) = dm.inverse();
        }

        lam_ = Ym * Pr / ((1 + Pr) * (1 - 2.0 * Pr));
        miu_ = Ym / (2.0 * (1.0 + Pr));
    }

    int elastic_potential::Val(const double *x, double *val) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        using namespace Eigen;
        const matd_t X = Map<const matd_t>(x, 3, dof_ / 3);

        for (size_t i = 0; i < tets_.cols(); ++i)
        {
            matd_t vert(3, 4);
            ele_verts(X, tets_.col(i), vert);
            double value;
            switch (type_)
            {
            case LINEAR:
                tet_linear_(&value, vert.data(), &Dm_(0, i), &vol_[i], &lam_, &miu_);
                break;
            case COROTATIONAL:
                value = 0.0;
                break;
            }
            *val += w_ * value;
        }

        return 0;
    }

    int elastic_potential::Gra(const double *x, double *gra) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        using namespace Eigen;
        const matd_t X = Map<const matd_t>(x, 3, dof_ / 3);
        Map<matd_t> G(gra, 3, dof_ / 3);

        for (size_t i = 0; i < tets_.cols(); ++i)
        {
            matd_t vert(3, 4), vert0(3, 4);
            ele_verts(X, tets_.col(i), vert);
            ele_verts(rest_, tets_.col(i), vert0);

            matd_t grad(3, 4);
            grad.setZero();

            switch (type_)
            {
            case LINEAR:
                tet_linear_jac_(grad.data(), vert.data(), &Dm_(0, i), &vol_[i], &lam_, &miu_);
                break;
            case COROTATIONAL:
                Matrix3d Ds, Dm, F, Q;
                Ds.col(0) = vert.col(1) - vert.col(0);
                Ds.col(1) = vert.col(2) - vert.col(0);
                Ds.col(2) = vert.col(3) - vert.col(0);
                Dm = Map<Matrix3d>(const_cast<double *>(&Dm_(0, i)));
                F = Ds * Dm;
                polar(F, Q);
                // cout << i << "Q:\n"
                //      << Q << endl;
                // cout << "Dm:\n"
                //      << Dm << endl;
                // cout << "F:\n"
                //      << F << endl;

                matd_t H(12, 12), R(12, 12);
                R.setZero();
                R.block(0, 0, 3, 3) = Q;
                R.block(3, 3, 3, 3) = Q;
                R.block(6, 6, 3, 3) = Q;
                R.block(9, 9, 3, 3) = Q;
                tet_linear_hes_(H.data(), nullptr, &Dm_(0, i), &vol_[i], &lam_, &miu_);
                // cout << "R:\n"
                //      << R << endl;
                // cout << "H:\n"
                //      << H << endl;

                Map<VectorXd> V(vert.data(), 12), V0(vert0.data(), 12);
                Map<VectorXd>(grad.data(), 12) = R * H * (R.transpose() * V - V0); // ToDO: why?
                // cout << "grad:\n"
                //      << grad << endl;
                break;
            }
            G.col(tets_(0, i)) += w_ * grad.col(0);
            G.col(tets_(1, i)) += w_ * grad.col(1);
            G.col(tets_(2, i)) += w_ * grad.col(2);
            G.col(tets_(3, i)) += w_ * grad.col(3);
        }

        return 0;
    }

    int elastic_potential::Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const
    {
        using namespace Eigen;
        const matd_t X = Map<const matd_t>(x, 3, dof_ / 3);

        for (size_t i = 0; i < tets_.cols(); ++i)
        {
            matd_t vert(3, 4);
            ele_verts(X, tets_.col(i), vert);
            matd_t H = matd_t::Zero(12, 12);
            switch (type_)
            {
            case LINEAR:
                tet_linear_hes_(H.data(), vert.data(), &Dm_(0, i), &vol_[i], &lam_, &miu_);
                break;
            case COROTATIONAL:
                Matrix3d Ds, Dm, F, Q;
                Ds.col(0) = vert.col(1) - vert.col(0);
                Ds.col(1) = vert.col(2) - vert.col(0);
                Ds.col(2) = vert.col(3) - vert.col(0);
                Dm = Map<const Matrix3d>(&Dm_(0, i));
                F = Ds * Dm;
                polar(F, Q);

                matd_t R(12, 12);
                R.setZero();
                R.block(0, 0, 3, 3) = Q;
                R.block(3, 3, 3, 3) = Q;
                R.block(6, 6, 3, 3) = Q;
                R.block(9, 9, 3, 3) = Q;
                // cout << "R\n"
                //      << R << endl;
                tet_linear_hes_(H.data(), nullptr, Dm.data(), &vol_[i], &lam_, &miu_);

                H = (R * H * R.transpose()).eval();
                break;
            }

            for (size_t p = 0; p < 12; ++p)
            {
                for (size_t q = 0; q < 12; ++q)
                {
                    if (H(p, q) != 0.0)
                    {
                        int I = 3 * tets_(p / 3, i) + p % 3;
                        int J = 3 * tets_(q / 3, i) + q % 3;
                        hes->push_back(Triplet<double>(I, J, w_ * H(p, q)));
                    }
                }
            }
        }
        return 0;
    }

    gravitational_potential::gravitational_potential(const matd_t &nods, const mati_t &cell,
                                                     const double rho,
                                                     const double w,
                                                     const int direction)
        : dof_(nods.size()), w_(w), direction_(direction)
    {
        calc_mass_matrix(nods, cell, rho, 1, M_, true);
    }

    int gravitational_potential::Val(const double *x, double *val) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        using namespace Eigen;
        Map<const matd_t> X(x, 3, dof_ / 3);
        // double gravity = 0;
        for (size_t i = 0; i < X.cols(); ++i)
        {
            *val += w_ * 9.8 * M_.coeff(i, i) * X(direction_, i);
            // gravity += X(direction_, i);
        }
        // cout << "gravity=" << gravity << endl;
        return 0;
    }

    int gravitational_potential::Gra(const double *x, double *gra) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        using namespace Eigen;
        Map<matd_t> G(gra, 3, dof_ / 3);
        // double gravG = 0;
        // #pragma omp parallel for
        for (size_t i = 0; i < G.cols(); ++i)
        {
            G(direction_, i) += w_ * 9.8 * M_.coeff(i, i);
            // gravG += w_ * 9.8 * M_.coeff(i, i);
        }
        // cout << "gravG:" << gravG << endl;
        return 0;
    }

    positional_potential::positional_potential(const matd_t &nods, const double w)
        : dof_(nods.size()), w_(w) {}

    int positional_potential::Val(const double *x, double *val) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0 || fixed_.size() == 0);

        using namespace Eigen;
        Map<const matd_t> X(x, 3, dof_ / 3);
        for (auto &ele : fixed_)
        {
            size_t id = ele.first;
            *val += 0.5 * w_ * (X.col(id) - ele.second).squaredNorm();
        }
        return 0;
    }

    int positional_potential::Gra(const double *x, double *gra) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0 || fixed_.size() == 0);

        using namespace Eigen;
        Map<const matd_t> X(x, 3, dof_ / 3);
        Map<matd_t> G(gra, 3, dof_ / 3);
        for (auto &ele : fixed_)
        {
            size_t id = ele.first;
            G.col(id) += w_ * (X.col(id) - ele.second);
        }
        return 0;
    }

    int positional_potential::Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0 || fixed_.size() == 0);
        using namespace Eigen;

        for (auto &ele : fixed_)
        {
            size_t id = ele.first;
            hes->push_back(Triplet<double>(3 * id + 0, 3 * id + 0, w_));
            hes->push_back(Triplet<double>(3 * id + 1, 3 * id + 1, w_));
            hes->push_back(Triplet<double>(3 * id + 2, 3 * id + 2, w_));
        }
        return 0;
    }

    int positional_potential::Pin(const size_t id, const double *pos)
    {
        assert(id >= 0 && id < dof_ / 3 && "Vert does not exist");
        fixed_[id] = Eigen::Vector3d(pos);
        return 0;
    }

    int positional_potential::Release(const size_t id)
    {
        assert(id >= 0 && id < dof_ / 3 && "Vert does not exist");
        auto it = fixed_.find(id);

        if (it == fixed_.end())
        {
            cerr << "[info] vertex " << id << " is not fixed\n";
            return __LINE__;
        }
        fixed_.erase(it);
        return 0;
    }

    ext_force_energy::ext_force_energy(const matd_t &nods, const double w)
        : dof_(nods.size()), w_(w)
    {
        force_.resize(dof_);
        force_.setZero();
    }

    int ext_force_energy::Val(const double *x, double *val) const
    {
        using namespace Eigen;
        Map<const VectorXd> X(x, dof_);
        *val += -w_ * force_.dot(X);
        return 0;
    }

    int ext_force_energy::Gra(const double *x, double *gra) const
    {
        using namespace Eigen;
        Map<VectorXd> G(gra, dof_);
        G += -w_ * force_;
        return 0;
    }

    int ext_force_energy::ApplyForce(const size_t id, const double *f)
    {
        assert(id >= 0 && id < dof_ / 3 && "Vert does not exist");
        std::copy(f, f + 3, &force_[3 * id]);
        // cout << "ApplyForce:" << id << " f:" << force_.block(3 * id, 0, 3, 1).transpose() << endl;
        return 0;
    }

    int ext_force_energy::RemoveForce(const size_t id)
    {
        assert(id >= 0 && id < dof_ / 3 && "Vert does not exist");
        std::fill(&force_[3 * id + 0], &force_[3 * id + 3], 0);
        return 0;
    }

    bw98_stretch_energy::bw98_stretch_energy(const matd_t &nods, const mati_t &tris, const double w)
        : dof_(nods.size()), w_(w), tris_(tris)
    {
        Dm_.resize(4, tris.cols());
        area_.resize(tris.cols());
        matd_t O(3, tris.cols()), N(3, tris.cols()), T(3, tris.cols()), B(3, tris.cols());
        face_normal(nods, tris, N, false);
        for (size_t i = 0; i < tris.cols(); ++i)
        {
            area_[i] = 0.5 * N.col(i).norm();
            N.col(i).normalize();
            O.col(i) = (nods.col(tris(0, i)) + nods.col(tris(1, i)) + nods.col(tris(2, i))) / 3.0;
            T.col(i) = nods.col(tris(1, i)) - nods.col(tris(0, i));
            T.col(i).normalize();
            B.col(i) = cross(N.col(i), T.col(i));
            B.col(i).normalize();
        }

        for (size_t i = 0; i < tris.cols(); ++i)
        {
            matd_t uv(2, 3), base(2, 2);
            matd_t vert(3, 3);
            tri_verts(nods, tris_.col(i), vert);
            for (size_t j = 0; j < 3; ++j)
            {
                uv(0, j) = (vert.col(j) - O.col(i)).dot(T.col(i));
                uv(1, j) = (vert.col(j) - O.col(i)).dot(B.col(i));
            }
            base.col(0) = uv.col(1) - uv.col(0);
            base.col(1) = uv.col(2) - uv.col(0);
            assert(fabs(base.determinant()) > 1e-8 && "bw98_stretch: triangle degenerate!");
            base = base.inverse();
            Dm_.col(i) = Map<VectorXd>(base.data(), 4, 1);
        }
    }

    int bw98_stretch_energy::Val(const double *x, double *val) const
    {
        if (w_ == 0.0)
            return 0;
        // E = k/2 * C(i) * C(i) where C(i) = area * ( norm(Wu)-1)
        Map<const matd_t> X(x, 3, dof_ / 3);
        double value_sum = 0;
        double tmpval_sum = 0;
        for (size_t i = 0; i < tris_.cols(); ++i)
        {
            // matd_t Ds(3, 2);
            // Ds.col(0) = X.col(tris_(1, i)) - X.col(tris_(0, i));
            // Ds.col(1) = X.col(tris_(2, i)) - X.col(tris_(0, i));

            // matd_t Wu = Ds * Map<const matd_t>(&(Dm_(0, i)), 2, 2);
            // Vector2d C;
            // C(0) = (Wu.col(0).norm() - 1.0);
            // C(1) = (Wu.col(1).norm() - 1.0);
            // double value = 0.5 * area_[i] * C.dot(C);

            matd_t vert(3, 3);
            vert.col(0) = X.col(tris_(0, i));
            vert.col(1) = X.col(tris_(1, i));
            vert.col(2) = X.col(tris_(2, i));
            double tmp_val = 0;
            bw98_stretch_(&tmp_val, vert.data(), &Dm_(0, i), &area_[i]);
            value_sum += tmp_val;
            // if (value > 1e-8)
            //     cout << tris_.col(i).transpose() << " value=" << value << " tmp_val=" << tmp_val << " err=" << fabs(value - tmp_val) << endl;
        }
        *val += w_ * value_sum;
        return 0;
    }

    int bw98_stretch_energy::Gra(const double *x, double *gra) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        Map<const matd_t> X(x, 3, dof_ / 3);
        Map<matd_t> G(gra, 3, dof_ / 3);
        for (size_t i = 0; i < tris_.cols(); ++i)
        {
            matd_t vert(3, 3);
            vert.col(0) = X.col(tris_(0, i));
            vert.col(1) = X.col(tris_(1, i));
            vert.col(2) = X.col(tris_(2, i));
            matd_t Grad = matd_t::Zero(3, 3);
            bw98_stretch_jac_(Grad.data(), vert.data(), &Dm_(0, i), &area_[i]);
            G.col(tris_(0, i)) += w_ * Grad.col(0);
            G.col(tris_(1, i)) += w_ * Grad.col(1);
            G.col(tris_(2, i)) += w_ * Grad.col(2);
        }
        return 0;
    }

    int bw98_stretch_energy::Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        Map<const matd_t> X(x, 3, dof_ / 3);
        for (size_t i = 0; i < tris_.cols(); ++i)
        {
            matd_t vert(3, 3);
            vert.col(0) = X.col(tris_(0, i));
            vert.col(1) = X.col(tris_(1, i));
            vert.col(2) = X.col(tris_(2, i));
            matd_t H = matd_t::Zero(9, 9);
            bw98_stretch_hes_(H.data(), vert.data(), &Dm_(0, i), &area_[i]);

            for (int p = 0; p < 9; ++p)
            {
                for (int q = 0; q < 9; ++q)
                {
                    int I = 3 * tris_(p / 3, i) + p % 3; // nod's p%3 dimension
                    int J = 3 * tris_(q / 3, i) + q % 3;
                    hes->emplace_back(I, J, H(p, q));
                }
            }
        }
        return 0;
    }

    bw98_shear_energy::bw98_shear_energy(const matd_t &nods, const mati_t &tris, const double w)
        : dof_(nods.size()), w_(w), tris_(tris)
    {
        Dm_.resize(4, tris.cols());
        area_.resize(tris.cols());
        matd_t O(3, tris.cols()), N(3, tris.cols()), T(3, tris.cols()), B(3, tris.cols());
        face_normal(nods, tris, N, false);
        for (size_t i = 0; i < tris.cols(); ++i)
        {
            area_[i] = 0.5 * N.col(i).norm();
            N.col(i).normalize();
            O.col(i) = (nods.col(tris(0, i)) + nods.col(tris(1, i)) + nods.col(tris(2, i))) / 3.0;
            T.col(i) = nods.col(tris(1, i)) - nods.col(tris(0, i));
            T.col(i).normalize();
            B.col(i) = cross(N.col(i), T.col(i));
            B.col(i).normalize();
        }

        for (size_t i = 0; i < tris.cols(); ++i)
        {
            matd_t uv(2, 3), base(2, 2);
            matd_t vert(3, 3);
            tri_verts(nods, tris_.col(i), vert);
            for (size_t j = 0; j < 3; ++j)
            {
                uv(0, j) = (vert.col(j) - O.col(i)).dot(T.col(i));
                uv(1, j) = (vert.col(j) - O.col(i)).dot(B.col(i));
            }
            base.col(0) = uv.col(1) - uv.col(0);
            base.col(1) = uv.col(2) - uv.col(0);
            assert(fabs(base.determinant()) > 1e-8 && "bw98_stretch: triangle degenerate!");
            base = base.inverse();
            Dm_.col(i) = Map<VectorXd>(base.data(), 4, 1);
        }
    }

    int bw98_shear_energy::Val(const double *x, double *val) const
    {
        if (w_ == 0.0)
            return 0;
        //  C(i) = area*( Wu'Wv )
        Map<const matd_t> X(x, 3, dof_ / 3);
        double value_sum = 0;
        double tmpval_sum = 0;
        for (size_t i = 0; i < tris_.cols(); ++i)
        {
            matd_t vert(3, 3);
            vert.col(0) = X.col(tris_(0, i));
            vert.col(1) = X.col(tris_(1, i));
            vert.col(2) = X.col(tris_(2, i));
            double tmp_val = 0;
            bw98_shear_(&tmp_val, vert.data(), &Dm_(0, i), &area_[i]);
            value_sum += tmp_val;
        }
        *val += w_ * value_sum;
        return 0;
    }

    int bw98_shear_energy::Gra(const double *x, double *gra) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        Map<const matd_t> X(x, 3, dof_ / 3);
        Map<matd_t> G(gra, 3, dof_ / 3);
        for (size_t i = 0; i < tris_.cols(); ++i)
        {
            matd_t vert(3, 3);
            vert.col(0) = X.col(tris_(0, i));
            vert.col(1) = X.col(tris_(1, i));
            vert.col(2) = X.col(tris_(2, i));
            matd_t Grad = matd_t::Zero(3, 3);
            bw98_shear_jac_(Grad.data(), vert.data(), &Dm_(0, i), &area_[i]);
            G.col(tris_(0, i)) += w_ * Grad.col(0);
            G.col(tris_(1, i)) += w_ * Grad.col(1);
            G.col(tris_(2, i)) += w_ * Grad.col(2);
        }
        return 0;
    }

    int bw98_shear_energy::Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        Map<const matd_t> X(x, 3, dof_ / 3);
        for (size_t i = 0; i < tris_.cols(); ++i)
        {
            matd_t vert(3, 3);
            vert.col(0) = X.col(tris_(0, i));
            vert.col(1) = X.col(tris_(1, i));
            vert.col(2) = X.col(tris_(2, i));
            matd_t H = matd_t::Zero(9, 9);
            bw98_shear_hes_(H.data(), vert.data(), &Dm_(0, i), &area_[i]);

            for (int p = 0; p < 9; ++p)
            {
                for (int q = 0; q < 9; ++q)
                {
                    int I = 3 * tris_(p / 3, i) + p % 3; // nod's p%3 dimension
                    int J = 3 * tris_(q / 3, i) + q % 3;
                    hes->emplace_back(I, J, H(p, q));
                }
            }
        }
        return 0;
    }

    surf_bending_potential::surf_bending_potential(const mati_t &diams, const matd_t &nods, const double w)
        : dof_(nods.size()), diams_(diams), w_(w)
    {
        len_.resize(diams_.cols());
        angle_.resize(diams_.cols());
        area_.resize(diams_.cols());

        for (size_t i = 0; i < diams_.cols(); ++i)
        {
            matd_t vert(3, 4);
            vert.col(0) = nods.col(diams_(0, i));
            vert.col(1) = nods.col(diams_(1, i));
            vert.col(2) = nods.col(diams_(2, i));
            vert.col(3) = nods.col(diams_(3, i));
            calc_edge_length_(&len_[i], &vert(0, 1));
            calc_dih_angle_(&angle_[i], vert.data());
            area_[i] = calc_tri_area(vert.leftCols(3)) + calc_tri_area(vert.rightCols(3));
        }
    }

    int surf_bending_potential::Val(const double *x, double *val) const
    {
        using namespace Eigen;
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        Map<const matd_t> X(x, 3, dof_ / 3);
        for (size_t i = 0; i < diams_.cols(); ++i)
        {
            matd_t vert(3, 4);
            vert.col(0) = X.col(diams_(0, i));
            vert.col(1) = X.col(diams_(1, i));
            vert.col(2) = X.col(diams_(2, i));
            vert.col(3) = X.col(diams_(3, i));
            double value = 0;
            surf_bending_(&value, vert.data(), &angle_[i], &len_[i], &area_[i]);
            *val += w_ * value;
        }
        return 0;
    }

    int surf_bending_potential::Gra(const double *x, double *gra) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        Map<const matd_t> X(x, 3, dof_ / 3);
        Map<matd_t> G(gra, 3, dof_ / 3);
        for (size_t i = 0; i < diams_.cols(); ++i)
        {
            matd_t vert(3, 4);
            vert.col(0) = X.col(diams_(0, i));
            vert.col(1) = X.col(diams_(1, i));
            vert.col(2) = X.col(diams_(2, i));
            vert.col(3) = X.col(diams_(3, i));
            matd_t Grad = matd_t::Zero(3, 4);
            surf_bending_jac_(Grad.data(), vert.data(), &angle_[i], &len_[i], &area_[i]);
            G.col(diams_(0, i)) += w_ * Grad.col(0);
            G.col(diams_(1, i)) += w_ * Grad.col(1);
            G.col(diams_(2, i)) += w_ * Grad.col(2);
            G.col(diams_(3, i)) += w_ * Grad.col(3);
        }
        return 0;
    }

    int surf_bending_potential::Hes(const double *x, vector<Triplet<double>> *hes) const
    {
        RETURN_WITH_COND_TRUE(w_ == 0.0);
        Map<const matd_t> X(x, 3, dof_ / 3);
        for (size_t i = 0; i < diams_.cols(); ++i)
        {
            matd_t vert(3, 4);
            vert.col(0) = X.col(diams_(0, i));
            vert.col(1) = X.col(diams_(1, i));
            vert.col(2) = X.col(diams_(2, i));
            vert.col(3) = X.col(diams_(3, i));
            matd_t H = matd_t::Zero(12, 12);
            surf_bending_hes_(H.data(), vert.data(), &angle_[i], &len_[i], &area_[i]);
            for (int p = 0; p < 12; ++p)
            {
                for (int q = 0; q < 12; ++q)
                {
                    int I = 3 * diams_(p / 3, i) + p % 3; // nod's p%3 dimension
                    int J = 3 * diams_(q / 3, i) + q % 3;
                    hes->emplace_back(I, J, H(p, q));
                }
            }
        }
        return 0;
    }
}
