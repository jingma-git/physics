#pragma once
#include <Eigen/Eigen>
#include <memory>
#include <exception>

#ifndef DEBUG
#define DEBUG
#endif

#ifdef DEBUG
#include <iostream>
using namespace std;
using namespace Eigen;
#endif

namespace egl
{
    typedef Eigen::Matrix<size_t, -1, -1> mati_t;
    typedef Eigen::Matrix<double, -1, -1> matd_t;
    typedef Eigen::SparseMatrix<double> spmat_t;

    template <typename T>
    class Functional
    {
    public:
        virtual ~Functional() {}
        virtual size_t Dof() const = 0;                // how many nodes/particles in the system
        virtual int Val(const T *x, T *val) const = 0; // value
        virtual int Gra(const T *x, T *gra) const = 0;
        virtual int Hes(const T *x, std::vector<Eigen::Triplet<T>> *hes) const = 0;
    };

    typedef std::shared_ptr<Functional<double>> pfunc;

    template <typename T>
    class energy_t : public Functional<T>
    {
    public:
        class null_input_exception : public std::exception
        {
        public:
            const char *what() const throw()
            {
                return "null input exception";
            }
        };

        class compatibility_exception : public std::exception
        {
        public:
            const char *what() const throw()
            {
                return "compatibility exception";
            }
        };

        energy_t(const std::vector<std::shared_ptr<Functional<T>>> &buffer)
            : buffer_(buffer), dof_(-1)
        {
            for (auto &e : buffer_)
            {
                if (e.get())
                {
                    dof_ = e->Dof();
                    break;
                }
            }

            if (dof_ == -1)
            {
                throw null_input_exception();
            }

            for (auto &e : buffer_)
            {
                if (e.get() && e->Dof() != dof_)
                {
                    throw compatibility_exception();
                }
            }
        }

        size_t Dof() const { return dof_; }

        int Val(const T *x, T *val) const
        {
            // double oval = 0;
            // int i = 0;
            for (auto &e : buffer_)
            {
                if (e.get())
                {
                    e->Val(x, val);
                    // cout << i++ << ": val=" << *val << " curVal=" << *val - oval << endl;
                    // oval = *val;
                }
            }
            return 0;
        }

        int Gra(const T *x, T *gra) const
        {
            // cout << "Grad" << endl;
            // double ograd = 0;
            // int i = 0;
            for (auto &e : buffer_)
            {
                if (e.get())
                {
                    e->Gra(x, gra);
                    // Eigen::Map<VectorXd> G(gra, dof_);
                    // double g = G.norm();
                    // cout << i++ << " g:" << g << " curG:" << g - ograd << endl;
                    // ograd = g;
                }
            }
            return 0;
        }

        int Hes(const T *x, std::vector<Eigen::Triplet<T>> *hes) const
        {
            // cout << "Hes:" << endl;
            // double ohes = 0;
            // int i = 0;
            for (auto &e : buffer_)
            {
                if (e.get())
                {
                    e->Hes(x, hes);
                    // Eigen::SparseMatrix<double> H(dof_, dof_);
                    // H.setFromTriplets(hes->begin(), hes->end());
                    // double hnorm = H.norm();
                    // cout << i++ << " hnorm=" << hnorm << " curH=" << hnorm - ohes << endl;
                    // ohes = hnorm;
                }
            }
            return 0;
        }

    protected:
        const std::vector<pfunc> &buffer_;
        size_t dof_;
    };

}