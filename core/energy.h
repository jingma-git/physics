#pragma once
#include "def.h"
#include <unordered_map>

namespace egl
{
    class momentum_potential : public Functional<double>
    {
    public:
        ~momentum_potential() {}
        virtual size_t Dof() const = 0;
        virtual int Val(const double *x, double *val) const = 0;
        virtual int Gra(const double *x, double *gra) const = 0;
        virtual int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const = 0;
        virtual void Update(const double *x) = 0;
        virtual double QueryKineticEnergy() const = 0;
        virtual const Eigen::SparseMatrix<double> &MassMatrix() const = 0;
        virtual const Eigen::VectorXd &CurrVelocity() const = 0;
        virtual double timestep() const = 0;
    };

    class momentum_potential_imp_euler : public momentum_potential
    {
    public:
        momentum_potential_imp_euler(const matd_t &nods,
                                     const mati_t &cell,
                                     const double rho,
                                     const double h,
                                     const double w = 1.0);
        size_t Dof() const { return dof_; }
        int Val(const double *x, double *val) const;
        int Gra(const double *x, double *gra) const;
        int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
        void Init(const double *x0, const double *v0);
        void Update(const double *x);
        double QueryKineticEnergy() const;
        const Eigen::SparseMatrix<double> &MassMatrix() const { return M_; }
        const Eigen::VectorXd &CurrVelocity() const { return vn_; }
        double timestep() const { return h_; }

    private:
        const double rho_, h_;
        const size_t dof_;
        double w_;
        Eigen::SparseMatrix<double> M_;
        Eigen::VectorXd xn_, vn_;
    };

    class elastic_potential : public Functional<double>
    {
    public:
        enum Material
        {
            LINEAR,
            COROTATIONAL
        };

        elastic_potential(const matd_t &nods,
                          const mati_t &tets,
                          Material type,
                          const double Ym,
                          const double Pr,
                          const double w);
        size_t Dof() const { return dof_; }
        int Val(const double *x, double *val) const;
        int Gra(const double *x, double *gra) const;
        int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;

    private:
        matd_t rest_;
        const mati_t &tets_;
        Material type_;

        const size_t dof_;
        double w_;
        double lam_, miu_;

        std::vector<double> vol_;
        matd_t Dm_;
    };

    class gravitational_potential : public Functional<double>
    {
    public:
        gravitational_potential(const matd_t &nods, const mati_t &cell,
                                const double rho, const double w, const int direction = 1);
        size_t Dof() const { return dof_; }
        int Val(const double *x, double *val) const;
        int Gra(const double *x, double *gra) const;
        int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const { return __LINE__; }

    private:
        const int direction_;
        const size_t dof_;
        double w_;
        Eigen::SparseMatrix<double> M_;
    };

    class positional_potential : public Functional<double>
    {
    public:
        positional_potential(const matd_t &nods, const double w);
        size_t Dof() const { return dof_; };
        int Val(const double *x, double *val) const;
        int Gra(const double *x, double *gra) const;
        int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const;
        int Pin(const size_t id, const double *pos);
        int Release(const size_t id);

    private:
        const size_t dof_;
        double w_;
        std::unordered_map<size_t, Eigen::Vector3d> fixed_;
    };

    class ext_force_energy : public Functional<double>
    {
    public:
        ext_force_energy(const matd_t &nods, const double w);
        size_t Dof() const { return dof_; }
        int Val(const double *x, double *val) const;
        int Gra(const double *x, double *gra) const;
        int Hes(const double *x, std::vector<Eigen::Triplet<double>> *hes) const { return __LINE__; }
        int ApplyForce(const size_t id, const double *f);
        int RemoveForce(const size_t id);

    private:
        const size_t dof_;
        double w_;
        Eigen::VectorXd force_;
    };
}
