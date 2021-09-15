#pragma once

#include <Eigen/Dense>
namespace egl
{
    template <typename Scalar, typename Index>
    bool readOBJ(
        const std::string obj_file_name,
        std::vector<std::vector<Scalar>> &V,
        std::vector<std::vector<Scalar>> &TC,
        std::vector<std::vector<Scalar>> &N,
        std::vector<std::vector<Index>> &F,
        std::vector<std::vector<Index>> &FTC,
        std::vector<std::vector<Index>> &FN)
    {
        // Open file, and check for error
        FILE *obj_file = fopen(obj_file_name.c_str(), "r");
        if (NULL == obj_file)
        {
            fprintf(stderr, "IOError: %s could not be opened...\n",
                    obj_file_name.c_str());
            return false;
        }
        // File open was succesfull so clear outputs
        V.clear();
        TC.clear();
        N.clear();
        F.clear();
        FTC.clear();
        FN.clear();

        // variables an constants to assist parsing the .obj file
        // Constant strings to compare against
        std::string v("v");
        std::string vn("vn");
        std::string vt("vt");
        std::string f("f");
        std::string tic_tac_toe("#");
#ifndef IGL_LINE_MAX
#define IGL_LINE_MAX 2048
#endif

        char line[IGL_LINE_MAX];
        int line_no = 1;
        while (fgets(line, IGL_LINE_MAX, obj_file) != NULL)
        {
            char type[IGL_LINE_MAX];
            // Read first word containing type
            if (sscanf(line, "%s", type) == 1)
            {
                // Get pointer to rest of line right after type
                char *l = &line[strlen(type)];
                if (type == v)
                {
                    double x[4];
                    int count =
                        sscanf(l, "%lf %lf %lf %lf\n", &x[0], &x[1], &x[2], &x[3]);
                    if (count != 3 && count != 4)
                    {
                        fprintf(stderr,
                                "Error: readOBJ() vertex on line %d should have 3 or 4 coordinates",
                                line_no);
                        fclose(obj_file);
                        return false;
                    }
                    std::vector<Scalar> vertex(count);
                    for (int i = 0; i < count; i++)
                    {
                        vertex[i] = x[i];
                    }
                    V.push_back(vertex);
                }
                else if (type == vn)
                {
                    double x[3];
                    int count =
                        sscanf(l, "%lf %lf %lf\n", &x[0], &x[1], &x[2]);
                    if (count != 3)
                    {
                        fprintf(stderr,
                                "Error: readOBJ() normal on line %d should have 3 coordinates",
                                line_no);
                        fclose(obj_file);
                        return false;
                    }
                    std::vector<Scalar> normal(count);
                    for (int i = 0; i < count; i++)
                    {
                        normal[i] = x[i];
                    }
                    N.push_back(normal);
                }
                else if (type == vt)
                {
                    double x[3];
                    int count =
                        sscanf(l, "%lf %lf %lf\n", &x[0], &x[1], &x[2]);
                    if (count != 2 && count != 3)
                    {
                        fprintf(stderr,
                                "Error: readOBJ() texture coords on line %d should have 2 "
                                "or 3 coordinates (%d)",
                                line_no, count);
                        fclose(obj_file);
                        return false;
                    }
                    std::vector<Scalar> tex(count);
                    for (int i = 0; i < count; i++)
                    {
                        tex[i] = x[i];
                    }
                    TC.push_back(tex);
                }
                else if (type == f)
                {
                    std::vector<Index> f;
                    std::vector<Index> ftc;
                    std::vector<Index> fn;
                    // Read each "word" after type
                    char word[IGL_LINE_MAX];
                    int offset;
                    while (sscanf(l, "%s%n", word, &offset) == 1)
                    {
                        // adjust offset
                        l += offset;
                        // Process word
                        unsigned int i, it, in;
                        if (sscanf(word, "%u/%u/%u", &i, &it, &in) == 3)
                        {
                            f.push_back(i - 1);
                            ftc.push_back(it - 1);
                            fn.push_back(in - 1);
                        }
                        else if (sscanf(word, "%u/%u", &i, &it) == 2)
                        {
                            f.push_back(i - 1);
                            ftc.push_back(it - 1);
                        }
                        else if (sscanf(word, "%u//%u", &i, &in) == 2)
                        {
                            f.push_back(i - 1);
                            fn.push_back(in - 1);
                        }
                        else if (sscanf(word, "%u", &i) == 1)
                        {
                            f.push_back(i - 1);
                        }
                        else
                        {
                            fprintf(stderr,
                                    "Error: readOBJ() face on line %d has invalid element format\n",
                                    line_no);
                            fclose(obj_file);
                            return false;
                        }
                    }
                    if (
                        (f.size() > 0 && fn.size() == 0 && ftc.size() == 0) ||
                        (f.size() > 0 && fn.size() == f.size() && ftc.size() == 0) ||
                        (f.size() > 0 && fn.size() == 0 && ftc.size() == f.size()) ||
                        (f.size() > 0 && fn.size() == f.size() && ftc.size() == f.size()))
                    {
                        // No matter what add each type to lists so that lists are the
                        // correct lengths
                        F.push_back(f);
                        FTC.push_back(ftc);
                        FN.push_back(fn);
                    }
                    else
                    {
                        fprintf(stderr,
                                "Error: readOBJ() face on line %d has invalid format\n", line_no);
                        fclose(obj_file);
                        return false;
                    }
                }
                else if (strlen(type) >= 1 && (type[0] == '#' ||
                                               type[0] == 'g' ||
                                               type[0] == 's' ||
                                               strcmp("usemtl", type) == 0 ||
                                               strcmp("mtllib", type) == 0))
                {
                    //ignore comments or other shit
                }
                else
                {
                    //ignore any other lines
                    fprintf(stderr,
                            "Warning: readOBJ() ignored non-comment line %d:\n  %s",
                            line_no,
                            line);
                }
            }
            else
            {
                // ignore empty line
            }
            line_no++;
        }
        fclose(obj_file);

        assert(F.size() == FN.size());
        assert(F.size() == FTC.size());

        return true;
    }

    template <typename DerivedV, typename DerivedF>
    bool readOBJ(
        const std::string str,
        Eigen::PlainObjectBase<DerivedV> &V,
        Eigen::PlainObjectBase<DerivedF> &F)
    {
        std::vector<std::vector<double>> vV, vTC, vN;
        std::vector<std::vector<int>> vF, vFTC, vFN;
        bool success = egl::readOBJ(str, vV, vTC, vN, vF, vFTC, vFN);
        if (success)
        {
            assert(vF.size() > 0 && "readOBJ.h OBJ format is wrong!");
            size_t numV = vV.size();
            size_t dimV = vV[0].size();
            size_t numF = vF.size();
            size_t dimF = vF[0].size();
            V.resize(dimV, numV);
            F.resize(dimF, numF);
            for (size_t i = 0; i < numV; ++i)
            {
                for (size_t j = 0; j < dimV; ++j)
                {
                    V(j, i) = vV[i][j];
                }
            }

            for (size_t i = 0; i < numF; ++i)
            {
                for (size_t j = 0; j < dimF; ++j)
                {
                    F(j, i) = vF[i][j];
                }
            }
        }
        return success;
    }
}