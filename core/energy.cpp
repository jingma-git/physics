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

}
