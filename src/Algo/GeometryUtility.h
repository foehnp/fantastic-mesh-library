#include <Eigen/Core>

#include <iostream>
#include <vector>

#pragma once


namespace FML {

namespace Util {

template<typename realType>
using Polygon = std::vector<Eigen::Vector3<realType>>;


template<typename realType>
bool splitPoints(const Eigen::Vector3<realType>& left, const Eigen::Vector3<realType>& right,
                 size_t dim, realType planeCoord,
                 Eigen::Vector3<realType>& res)
{
    realType diffWhole = right(dim) - left(dim);
    if (fabs(diffWhole) < 1e-5)
    {
        return false;
    }
    realType diffInter = planeCoord - left(dim);
    realType param = diffInter / diffWhole;

    res = left + param * (right - left);
    return true;
}

template<typename realType>
std::pair<Polygon<realType>, Polygon<realType>>
splitPolygonAtPlane(const Polygon<realType>& polygon, size_t dim, realType planeCoord)
{
    int numPoints = polygon.size();
    if (numPoints == 0)
    {
        return {{},{}};
    }
    bool initialSide = polygon[0](dim) > planeCoord;
    bool side = initialSide;
    int firstSplit = -1;
    int secondSplit = -1;


    for (size_t i = 1; i < polygon.size(); ++i)
    {
        bool newSide = polygon[i](dim) > planeCoord;
        if (newSide != side)
        {
            if (firstSplit == -1)
            {
                firstSplit = i;
            }
            else
            {
                secondSplit = i;
                break;
            }
            side = newSide;
        }
    }
    int secondSplitRolled;
    if (firstSplit != -1 && secondSplit == -1)
    {
        secondSplit = numPoints;
        secondSplitRolled = 0;
    }
    else
    {
        secondSplitRolled = secondSplit;
    }

    if (firstSplit == -1)
    {
        if (initialSide)
        {
            return {{}, polygon};
        }
        else
        {
            return {polygon, {}};
        }
    }

    Eigen::Vector3<realType> firstInter, secondInter;
    bool splitSuccess = splitPoints(polygon[firstSplit - 1], polygon[firstSplit], dim, planeCoord, firstInter) &&
    splitPoints(polygon[secondSplit - 1], polygon[secondSplitRolled], dim, planeCoord, secondInter);

    if (!splitSuccess)
    {
        return {{}, {}};
    }

    Polygon<realType> firstPolygon;
    firstPolygon.reserve(numPoints - secondSplit + firstSplit + 2);
    firstPolygon.push_back(secondInter);
    firstPolygon.insert(firstPolygon.end(), polygon.begin() + secondSplit, polygon.end());
    firstPolygon.insert(firstPolygon.end(), polygon.begin(), polygon.begin() + firstSplit);
    firstPolygon.push_back(firstInter);

    Polygon<realType> secondPolygon;
    secondPolygon.reserve(secondSplit - firstSplit + 2);
    secondPolygon.push_back(firstInter);
    secondPolygon.insert(secondPolygon.end(), polygon.begin() + firstSplit, polygon.begin() + secondSplit);
    secondPolygon.push_back(secondInter);

    if (initialSide)
    {
        return {std::move(secondPolygon), std::move(firstPolygon)};
    }
    else
    {
        return {std::move(firstPolygon), std::move(secondPolygon)};
    }
}

template<typename realType>
std::vector<Polygon<realType>> subPolygons(const std::vector<Eigen::Vector3<realType>>& polygon, const Eigen::Vector3<realType>& mid)
{
    std::vector<std::vector<Polygon<realType>>> lvls(4);

    lvls[0].resize(1);
    lvls[0][0] = polygon;
    lvls[1].reserve(2);
    lvls[2].reserve(4);
    lvls[3].reserve(8);

    int num = 1;
    for (int d = 0; d < 3; ++d)
    {
        for (int i = 0; i < num; ++ i)
        {
            std::pair<Polygon<realType>, Polygon<realType>> newPolygons = splitPolygonAtPlane(lvls[d][i], d, mid(d));
            lvls[d + 1].push_back(std::move(newPolygons.first));
            lvls[d + 1].push_back(std::move(newPolygons.second));
        }
        num *= 2;
    }

    return std::move(lvls[3]);
}


template<typename realType>
realType distToBoxSqrd(const std::vector<Eigen::Vector3<realType>>& box, const Eigen::Vector3<realType> &point)
{
    Eigen::Vector3<realType> effDist;
    for (int d = 0; d < 3; ++d)
    {
        if (point(d) < box[0](d))
        {
            effDist(d) = box[0](d) - point(d);
        }
        else if (point(d) <= box[1](d))
        {
            effDist(d) = 0.;
        }
        else
        {
            effDist(d) = point(d) - box[1](d);
        }
    }
    return effDist.squaredNorm();
}

template<typename realType>
realType maxDistOfBoxSqrd(const std::vector<Eigen::Vector3<realType>>& box, const Eigen::Vector3<realType> &point)
{
    Eigen::Vector3<realType> effDist;
    for (int d = 0; d < 3; ++d)
    {
        realType leftDist = fabs(point(d) - box[0](d));
        realType rightDist = fabs(point(d) - box[1](d));
        effDist(d) = (leftDist < rightDist) ? rightDist : leftDist;
    }
    return effDist.squaredNorm();
}


template<typename realType>
Eigen::Vector2<realType> perpVector2D(const Eigen::Vector2<realType>& v)
{
    Eigen::Vector2<realType> res;
    res(0) = -v(1);
    res(1) = v(0);
    return res;
    // return {-v(1), v(0)};
}

/**
 * @brief distanceToTriangleSquared2D:
 * @param triangle corner indices need to go counter-clockwise
 * @param point
 * @return distance
 */
template<typename realType>
realType distanceToTriangleSquared2D(const std::vector<Eigen::Vector2<realType>>& triangle, const Eigen::Vector2<realType>& point)
{
    std::array<Eigen::Vector2<realType>, 3> e;
    e[0] = triangle[0] - triangle[2];
    e[1] = triangle[1] - triangle[0];
    e[2] = triangle[2] - triangle[1];

    std::array<Eigen::Vector2<realType>, 3> p;
    p[0] = perpVector2D<realType>(e[0]);
    p[1] = perpVector2D<realType>(e[1]);
    p[2] = perpVector2D<realType>(e[2]);

    std::array<bool, 3> insideOfSide;
    std::array<int, 3> whereOnProjection;

    for (int s = 0; s < 3; ++s)
    {
        realType vp = p[s].dot(point - triangle[s]);
        insideOfSide[s] = (vp >= 0);
        realType proj = e[s].dot(point);
        whereOnProjection[s] = (proj > e[s].dot(triangle[s])) ? 1 :
                                (proj > e[s].dot(triangle[(s+2)%3])) ? 0 : -1;
        if (!insideOfSide[s] && whereOnProjection[s] == 0)
        {
            return vp*vp / e[s].squaredNorm();
        }
    }

    if (insideOfSide[0] && insideOfSide[1] && insideOfSide[2])
    {
        return 0.;
    }

    for (int i = 0; i < 3; ++i)
    {
        if (whereOnProjection[i] == 1 && whereOnProjection[(i+1)%3] == -1)
        {
            return (point - triangle[i]).squaredNorm();
        }
    }

    return 0.;
}

template<typename realType>
realType distanceToTriangleSquared(const Polygon<realType>& triangle, const Eigen::Vector3<realType>& point)
{
    Eigen::Vector3<realType> e1 = (triangle[1] - triangle[0]).normalized();
    Eigen::Vector3<realType> normal = ((triangle[2] - triangle[0]).cross(e1)).normalized();
    Eigen::Vector3<realType> e2 = e1.cross(normal);
    std::vector<Eigen::Vector2<realType>> projTriangle(3);
    for (int i = 0; i < 3; ++i)
    {
        projTriangle[i] = {e1.dot(triangle[i]), e2.dot(triangle[i])};
    }
    Eigen::Vector2<realType> projPoint = {e1.dot(point), e2.dot(point)};
    realType distToPlaneSqrd = pow(normal.dot(point) - normal.dot(triangle[0]), 2);
    return distanceToTriangleSquared2D(projTriangle, projPoint) + distToPlaneSqrd;
}



}
}
