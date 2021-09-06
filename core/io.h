#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Eigen>

namespace egl
{
    template <typename FT, typename IT> // FT: float type, IT: index type
    inline bool tet_mesh_read_from_vtk(const char *path,
                                       std::vector<std::vector<FT>> &nods,
                                       std::vector<std::vector<IT>> &tets)
    {
        using namespace std;
        ifstream ifs(path);
        if (ifs.fail())
        {
            cerr << "Can not open file" << path << endl;
            return false;
        }

        string str;
        int point_num = 0, cell_num = 0;
        while (!ifs.eof())
        {
            ifs >> str;
            if (str == "POINTS")
            {
                ifs >> point_num >> str;
                nods.resize(point_num, std::vector<FT>(3));
                for (size_t i = 0; i < point_num; ++i)
                {
                    for (size_t j = 0; j < 3; ++j)
                        ifs >> nods[i][j];
                }
            }

            if (str == "CELLS")
            {
                ifs >> cell_num >> str;
                int point_number_of_cell = 0;
                tets.resize(cell_num, std::vector<IT>(4));
                for (size_t i = 0; i < cell_num; ++i)
                {
                    ifs >> point_number_of_cell;
                    if (point_number_of_cell != 4)
                    {
                        cerr << "Tet #vertices is not 4" << endl;
                        return false;
                    }
                    for (size_t j = 0; j < 4; ++j)
                        ifs >> tets[i][j];
                }
            }
        }
        ifs.close();
        return true;
    }

    template <typename FT, typename IT>
    bool tet_mesh_read_from_vtk(const char *path,
                                Eigen::Matrix<FT, -1, -1> &nods,
                                Eigen::Matrix<IT, -1, -1> &tets)
    {
        using namespace std;
        vector<vector<FT>> nods_list;
        vector<vector<IT>> tets_list;
        bool flag = tet_mesh_read_from_vtk(path, nods_list, tets_list);
        if (flag)
        {
            size_t num_nods = nods_list.size();
            size_t num_tets = tets_list.size();
            nods.resize(3, nods_list.size());
            tets.resize(4, tets_list.size());

            for (size_t i = 0; i < num_nods; ++i)
                std::copy(nods_list[i].begin(), nods_list[i].end(), &nods(0, i));

            for (size_t i = 0; i < num_tets; ++i)
                std::copy(tets_list[i].begin(), tets_list[i].end(), &tets(0, i));
        }
        return flag;
    }

    template <typename OS, typename FLOAT, typename INT>
    void tet2vtk(
        OS &os,
        const FLOAT *node, size_t node_num,
        const INT *tet, size_t tet_num)
    {
        os << "# vtk DataFile Version 2.0\nTET\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
        os << "POINTS " << node_num << " float\n";
        for (size_t i = 0; i < node_num; ++i)
            os << node[i * 3 + 0] << " " << node[i * 3 + 1] << " " << node[i * 3 + 2] << "\n";

        os << "CELLS " << tet_num << " " << tet_num * 5 << "\n";
        for (size_t i = 0; i < tet_num; ++i)
            os << 4 << "  "
               << tet[i * 4 + 0] << " " << tet[i * 4 + 1] << " "
               << tet[i * 4 + 2] << " " << tet[i * 4 + 3] << "\n";
        os << "CELL_TYPES " << tet_num << "\n";
        for (size_t i = 0; i < tet_num; ++i)
            os << 10 << "\n";
    }

    int read_fixed_verts(const char *filename, std::vector<size_t> &fixed);
}
