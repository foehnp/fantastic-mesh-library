#pragma once

#include <Eigen/Core>


namespace FML {


template <typename realType>
class Triangle
{
public:
    Triangle();
    Triangle(size_t i0, size_t i1, size_t i2);

    void setVertex(size_t which, size_t idx);

    size_t vertex(size_t which) const;

private:
    size_t vertices[3];

    Eigen::Vector3<realType> normal;
};

template <typename realType>
Triangle<realType>::Triangle() {}

template<typename realType>
inline Triangle<realType>::Triangle(size_t i0, size_t i1, size_t i2)
{
    vertices[0] = i0;
    vertices[1] = i1;
    vertices[2] = i2;
}

template<typename realType>
inline void Triangle<realType>::setVertex(size_t which, size_t idx)
{
    vertices[which] = idx;
}

template<typename realType>
size_t Triangle<realType>::vertex(size_t which) const
{
    return vertices[which];
}

}
