#pragma once

#include "Algo/GeometryUtility.h"

#include <Eigen/Core>

#include <forward_list>
#include <vector>

namespace FML {

namespace {

template<typename realType>
using Polygon = std::vector<Eigen::Vector3<realType>>;


template<typename realType>
struct PolygonOctNode
{
    std::vector<Eigen::Vector3<realType>> bounds;
    std::vector<PolygonOctNode<realType>*> children;
    std::vector<std::pair<size_t, Polygon<realType>>> polygons;
    PolygonOctNode() :
        bounds(2)
    {}

    ~PolygonOctNode()
    {
        for (const auto child : children)
        {
            if (child)
            {
                delete child;
            }
        }
    }
};

}


/**
 * @brief The TriangleOctree class
 * Creates an octree given a list of triangles in n*log(n) time, where n is the number of triangles
 * Allowes log(n) nearest triangle search.
 */
template<typename realType>
class TriangleOctree
{
private:
    std::vector<std::vector<Eigen::Vector3<realType>>> m_polygons;

    PolygonOctNode<realType> m_rootNode;

    // settings
    realType m_tol;
    size_t m_maxPolygonsPerNode;

public:
    TriangleOctree(const std::vector<Polygon<realType>>& polygons, realType tol = 1e-2, size_t maxPolygonsPerNode = 15);


    size_t nearestPoint(const Eigen::Vector3<realType>& point, realType maxDist) const;

private:

    /**
     * @brief meiosis subdivides the box of inNode into 8 subboxes
     *  (if it contains more than maxPolygonsPerNode polygons and at least
     *  one of the sides of its cuboid is wider than tol). Cuts up all the polygons in inNode.polygons that intersect
     *  the bounding box of more than one child node and assigns the uncut and new cut polygons to the children.
     * @param inNode
     * @param tol
     * @param maxPolygonsPerNode
     * @return
     */
    bool meiosis(PolygonOctNode<realType>& inNode, realType tol, size_t maxPolygonsPerNode);

    /**
     * @brief recursiveMeiosis performs meiosis on inNode and recursively on its new children
     * @param inNode
     * @param tol
     * @param maxPolygonsPerNode
     */
    void recursiveMeiosis(PolygonOctNode<realType>& inNode, realType tol, size_t maxPolygonsPerNode);

};

template<typename realType>
TriangleOctree<realType>::TriangleOctree(const std::vector<std::vector<Eigen::Vector3<realType>>>& polygons,
                                       realType tol, size_t maxPolygonsPerNode) :
    m_tol(tol),
    m_maxPolygonsPerNode(maxPolygonsPerNode),
    m_polygons(polygons)
{
    realType infty = std::numeric_limits<realType>::max();

    std::vector<Eigen::Vector3<realType>> bb(2);
    bb[0] = {infty, infty, infty};
    bb[1] = {-infty, -infty, -infty};

    for (const auto& polygon : polygons)
    {
        for (const Eigen::Vector3<realType>& p : polygon)
        {
            for (int d = 0; d < 3; ++d)
            {
                bb[0](d) = std::min(bb[0](d), p(d));
                bb[1](d) = std::max(bb[1](d), p(d));
            }
        }
    }

    m_rootNode.bounds = bb;
    m_rootNode.polygons.resize(polygons.size());
    for (int i = 0; i < polygons.size(); ++i)
    {
        m_rootNode.polygons[i] = {i, polygons[i]};
    }
    recursiveMeiosis(m_rootNode, m_tol, m_maxPolygonsPerNode);
}



template<typename realType>
size_t TriangleOctree<realType>::nearestPoint(const Eigen::Vector3<realType> &point, realType maxDist) const
{
    std::forward_list<const PolygonOctNode<realType>*> nodes;
    nodes.push_front(&m_rootNode);
    while (true)
    {
        realType minMaxDist = std::numeric_limits<realType>::max();
        for (const PolygonOctNode<realType>* node : nodes)
        {
            minMaxDist = std::min(minMaxDist, Util::maxDistOfBoxSqrd(node->bounds, point));
        }

        std::forward_list<const PolygonOctNode<realType>*> newNodes;
        bool anyNewNodes = false;
        for (const PolygonOctNode<realType>* node : nodes)
        {
            if (Util::distToBoxSqrd(node->bounds, point) < minMaxDist)
            {
                if (node->children.empty())
                {
                    newNodes.push_front(node);
                }
                else
                {
                    anyNewNodes = true;
                    for (auto child : node->children)
                    {
                        if (child != nullptr)
                        {
                            newNodes.push_front(child);
                        }
                    }
                }
            }
        }

        if (!anyNewNodes)
        {
            //deb
            std::vector<const PolygonOctNode<realType>*> debugNodes;
            for (const PolygonOctNode<realType>* node : nodes)
            {
                debugNodes.push_back(node);
            }
            realType debClosest = 1e10;
            size_t debClosestIdx = 9999999999;
            for (const auto debNode : debugNodes)
            {
                for (const std::pair<size_t, Polygon<realType>>& currPolygonPair : debNode->polygons)
                {
                    for (const auto& point : currPolygonPair.second)
                    {
                        if (point[0] < debClosest)
                        {
                            debClosest = point[0];
                            debClosestIdx = currPolygonPair.first;
                        }
                    }
                }
            }
            //end deb

            size_t closestPolygon;
            realType closestDist = std::numeric_limits<realType>::max();
            for (const PolygonOctNode<realType>* node : nodes)
            {
                for (const std::pair<size_t, Polygon<realType>>& currPolygonPair : node->polygons)
                {
                    size_t currPolygonIdx = currPolygonPair.first;
                    realType dist = Util::distanceToTriangleSquared(m_polygons[currPolygonIdx], point);
                    if (dist < closestDist)
                    {
                        closestPolygon = currPolygonIdx;
                        closestDist = dist;
                    }
                }
            }
            return closestPolygon;
        }

        nodes = std::move(newNodes);
    }
    return SIZE_MAX;
}

template<typename realType>
bool TriangleOctree<realType>::meiosis(PolygonOctNode<realType> &inNode, realType tol, size_t maxPolygonsPerNode)
{
    if (inNode.polygons.size() <= maxPolygonsPerNode)
    {
        // No need to further split node
        return false;
    }

    // Check if node is too small to split
    unsigned char tooSmallDims = 0;
    for (int d = 0; d < 3; ++d)
    {
        if (fabs(inNode.bounds[1](0) - inNode.bounds[0](0)) < tol)
        {
            ++tooSmallDims;
        }
    }
    if (tooSmallDims == 3)
    {
        return false;
    }

    Eigen::Vector3<realType> seps{0.5 * (inNode.bounds[0](0) + inNode.bounds[1](0)),
                                  0.5 * (inNode.bounds[0](1) + inNode.bounds[1](1)),
                                  0.5 * (inNode.bounds[0](2) + inNode.bounds[1](2)),};

    std::vector<std::vector<std::pair<size_t, Polygon<realType>>>> polygonsPerChild(8);
    for (const std::pair<size_t, Polygon<realType>>& polygonIdxPair : inNode.polygons)
    {
        std::vector<Polygon<realType>> childPolygons = Util::subPolygons(polygonIdxPair.second, seps);
        for (int cat = 0; cat < 8; ++cat)
        {
            if (!childPolygons[cat].empty())
            {
                // auto pair = std::make_pair(polygonIdxPair.first, (childPolygons[cat]));
                polygonsPerChild[cat].push_back(std::make_pair(polygonIdxPair.first, (childPolygons[cat])));
            }
        }
    }

    bool changed = false;
    inNode.children.resize(8);
    for (size_t cat = 0; cat < 8; ++cat)
    {
        if (polygonsPerChild[cat].empty())
        {
            inNode.children[cat] = nullptr;
        }
        else
        {
            inNode.children[cat] = new PolygonOctNode<realType>;
            size_t iCopy = cat;
            for (int d = 2; d >= 0; --d)
            {
                if (iCopy % 2 == 0)
                {
                    inNode.children[cat]->bounds[0](d) = inNode.bounds[0](d);
                    inNode.children[cat]->bounds[1](d) = seps(d);
                }
                else
                {
                    inNode.children[cat]->bounds[0](d) = seps(d);
                    inNode.children[cat]->bounds[1](d) = inNode.bounds[1](d);
                }
                iCopy /= 2;
            }
            inNode.children[cat]->polygons = polygonsPerChild[cat];
        }
    }

    return true;
}

template<typename realType>
void TriangleOctree<realType>::recursiveMeiosis(PolygonOctNode<realType> &inNode, realType tol, size_t maxPolygonsPerNode)
{
    if (!meiosis(inNode, tol, maxPolygonsPerNode))
    {
        return;
    }
    for (PolygonOctNode<realType>* newChild : inNode.children)
    {
        if (newChild)
        {
            recursiveMeiosis(*newChild, tol, maxPolygonsPerNode);
        }
    }
}

}
