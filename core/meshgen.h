#pragma once
#include <Eigen/Eigen>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <unordered_map>
#include <random>

namespace egl
{
    template <typename matd_t, typename mati_t>
    void generate_square_mesh(const std::string strategy,
                              const size_t pts_num,
                              const size_t dim,
                              matd_t &nods,
                              mati_t &tris)
    {
        typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
        typedef CGAL::Delaunay_triangulation_2<K> Delaunay;
        typedef K::Point_2 Point;
        typedef Delaunay::Vertex_handle Vertex_handle;
        typedef Delaunay::Face_handle Face_handle;

        const size_t sq_num = static_cast<size_t>(std::sqrt(pts_num));

        std::vector<Point> pts;
        if (strategy == "random")
        {
            std::mt19937 gen;
            std::uniform_real_distribution<double> dis(0.0, 1.0);

            pts.push_back(Point(0, 0));
            pts.push_back(Point(1, 0));
            pts.push_back(Point(1, 1));
            pts.push_back(Point(0, 1));

            //-> sample others
            size_t curr_size = pts.size();
            for (size_t i = curr_size; i < pts_num; ++i)
            {

                double x = dis(gen), y = dis(gen);
                pts.emplace_back(Point(x, y));
            }
        }
        else if (strategy == "regular")
        {
            double dx = 1.0 / double(sq_num - 1);
            for (size_t i = 0; i < sq_num; ++i)
                for (size_t j = 0; j < sq_num; ++j)
                    pts.push_back(Point(i * dx, j * dx));
        }

        std::unordered_map<Vertex_handle, size_t> v2i;
        Delaunay dt;
        dt.insert(pts.begin(), pts.end());

        int count = 0;
        nods.resize(dim, pts.size());
        for (Vertex_handle v : dt.finite_vertex_handles())
        {
            Point p = v->point();
            v2i[v] = count;
            nods(0, count) = p.x();
            nods(1, count) = p.y();
            if (dim == 3)
                nods(2, count) = 0;
            ++count;
        }
        tris.resize(3, dt.number_of_faces());
        count = 0;
        for (Face_handle f : dt.finite_face_handles())
        {
            tris(0, count) = v2i[f->vertex(0)];
            tris(1, count) = v2i[f->vertex(1)];
            tris(2, count) = v2i[f->vertex(2)];
            ++count;
        }
    }
}