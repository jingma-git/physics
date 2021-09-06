#include "io.h"
#include <string>
#include <iostream>
#include <fstream>

namespace egl
{
    int read_fixed_verts(const char *filename, std::vector<size_t> &fixed)
    {
        using namespace std;
        fixed.clear();
        ifstream ifs(filename);
        if (ifs.fail())
        {
            cerr << "[error] can not open " << filename << endl;
            return __LINE__;
        }
        size_t temp;
        while (ifs >> temp)
        {
            fixed.push_back(temp);
        }
        cout << "[info] fixed verts number: " << fixed.size() << endl;
        ifs.close();
        return 0;
    }
}
