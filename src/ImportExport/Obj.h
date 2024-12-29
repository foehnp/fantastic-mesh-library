#pragma once

#include "Mesh.h"
#include <string>
#include <fstream>

namespace FML {

namespace Obj {

// declarations

/**
 * @brief Imports a mesh of the 'obj' format as in the repo https://github.com/alecjacobson/common-3d-test-models
 * @param filePath input
 * @param mesh output
 * @return whether success
 */
template <typename realType>
bool importMesh(const std::string& filePath, Mesh<realType>& mesh);



// definitions

template <typename realType>
bool importMesh(const std::string& filePath, Mesh<realType>& mesh)
{
    std::ifstream ifs(filePath);

    if (!ifs.is_open())
    {
        return false;
    }

    std::string realTypeString;
    if (std::is_same_v<realType, float>)
    {
        realTypeString = 'f';
    }
    else if (std::is_same_v<realType, double>)
    {
        realTypeString = "lf";
    }
    else if (std::is_same_v<realType, int>)
    {
        realTypeString = 'i';
    }
    else if (std::is_same_v<realType, unsigned int>)
    {
        realTypeString = 'u';
    }
    else
    {
        return false;
    }
    std::string sscanfStringSTD = "v ";
    for (int i = 0; i < 3; ++i)
    {
        sscanfStringSTD += "%";
        sscanfStringSTD += realTypeString;
        if (i < 2)
        {
            sscanfStringSTD += " ";
        }
    }
    const char* sscanfString = sscanfStringSTD.c_str();

    std::vector<Eigen::Vector3<realType>> vertices;
    std::vector<std::vector<size_t>> triangles;

    std::string line;
    while (std::getline(ifs, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }
        if (line[0] == 'v')
        {
            Eigen::Vector3<realType> v;
            realType a, b ,c;
            auto debRet = sscanf(line.c_str(), sscanfString, &a, &b, &c);
            // sscanf(line.c_str(), sscanfString, &v(0), &v(1), &v(2));
            // vertices.push_back(std::move(v));
            vertices.emplace_back(a, b, c);
        }
        else if (line[0] == 'f')
        {
            unsigned int a, b, c;
            sscanf(line.c_str(), "f %u %u %u", &a, &b, &c);
            // triangles.emplace_back(a-1, b-1, c-1);
            std::vector<size_t> t(3);
            t[0] = a - 1;
            t[1] = b - 1;
            t[2] = c - 1;
            triangles.push_back(std::move(t));
        }
    }

    mesh = Mesh(vertices, triangles);

    return true;
}

} // namespace Obj

}
