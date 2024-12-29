#pragma once

#include <Eigen/Core>


namespace FML {

template <typename realType>
class Halfedge;

template <typename realType>
class Vertex : public Eigen::Vector3<realType>
{
public:
    Halfedge<realType>* halfedge = nullptr;

public:
    Vertex(const Eigen::Vector3f& vector) :
        Eigen::Vector3<realType>(vector)
    {}


    Vertex(realType x, realType y, realType z) :
        Eigen::Vector3<realType>(x, y, z)
    {
    }

};

}

