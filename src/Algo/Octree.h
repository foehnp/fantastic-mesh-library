#ifndef OCTREE_H
#define OCTREE_H

#include "Algo/GeometryUtility.h"

#include <Eigen/Core>

#include <queue>
#include <vector>


namespace FML {

namespace {

template<typename realType>
struct OctNode
{
    bool isPoint = false; // if this is true, the box is just a point, with coords bound[0]
    std::vector<Eigen::Vector3<realType>> bounds;
    std::vector<OctNode<realType>*> children;
    std::vector<size_t> points;
    OctNode() :
        bounds(2)
    {}

    void destroyChildren()
    {
        for (const auto child : children)
        {
            if (child)
            {
                child->destroyChildren();
                delete child;
            }
        }
    }
};

}


/**
 * @brief The Octree class
 * Creates an octree given a list of points in n*log(n) time, where n is the number of points
 * Allowes log(n) nearest neighbor search.
 */
template<typename realType>
class Octree
{
private:
    std::vector<Eigen::Vector3<realType>> m_points;

    OctNode<realType> m_rootNode;

    // settings
    realType m_tol;

public:
    Octree(const std::vector<Eigen::Vector3<realType>>& points, realType tol);

    ~Octree();


    size_t nearestPoint(const Eigen::Vector3<realType>& point, realType maxDist) const;

private:

    /**
     * @brief meiosis subdivides the box of inNode into 8 subboxes
     *  (if at least one of the sides of its cuboid is wider than tol).
     *  Assigns each point to the box it is contained in. A child node is created
     *  for each non-empty sub-box.
     * @param inNode
     * @param tol
     * @return
     */
    bool meiosis(OctNode<realType>& inNode, realType tol);

    /**
     * @brief recursiveMeiosis performs meiosis on inNode and recursively on its new children
     * @param inNode
     * @param tol
     */
    void recursiveMeiosis(OctNode<realType>& inNode, realType tol);

};

template<typename realType>
Octree<realType>::Octree(const std::vector<Eigen::Vector3<realType>> &points, realType tol) :
    m_tol(tol),
    m_points(points)
{
    realType infty = std::numeric_limits<realType>::max();

    std::vector<Eigen::Vector3<realType>> bb(2);
    bb[0] = {infty, infty, infty};
    bb[1] = {-infty, -infty, -infty};

    for (const Eigen::Vector3<realType>& p : points)
    {
        for (int d = 0; d < 3; ++d)
        {
            bb[0](d) = std::min(bb[0](d), p(d));
            bb[1](d) = std::max(bb[1](d), p(d));
        }
    }

    m_rootNode.bounds = bb;
    m_rootNode.points.resize(m_points.size());
    for (int i = 0; i < m_points.size(); ++i)
    {
        m_rootNode.points[i] = i;
    }
    recursiveMeiosis(m_rootNode, m_tol);
}

template<typename realType>
inline Octree<realType>::~Octree()
{
    m_rootNode.destroyChildren();
}


namespace {

template<typename realType>
realType distToNodeSqrd(const OctNode<realType>& node, const Eigen::Vector3<realType>& point)
{
    if (node.isPoint)
    {
        return (node.bounds[0]-point).squaredNorm();
    }
    else
    {
        return Util::distToBoxSqrd(node.bounds, point);
    }
}

template<typename realType>
realType maxDistOfNodeSqrd(const OctNode<realType>& node, const Eigen::Vector3<realType>& point)
{
    if (node.isPoint)
    {
        return (node.bounds[0]-point).squaredNorm();
    }
    else
    {
        return Util::maxDistOfBoxSqrd(node.bounds, point);
    }
}

}


template<typename realType>
size_t Octree<realType>::nearestPoint(const Eigen::Vector3<realType> &point, realType maxDist) const
{
    auto compare = [](const std::pair<realType, const OctNode<realType>*>& left, const std::pair<realType, const OctNode<realType>*>& right)
    {return left.first > right.first;};
    std::priority_queue<std::pair<realType, const OctNode<realType>*>,
                        std::vector<std::pair<realType, const OctNode<realType>*>>,
                        decltype(compare)>
        prQueue(compare);
    prQueue.emplace(distToNodeSqrd(m_rootNode, point), &m_rootNode);
    while (!prQueue.top().second->isPoint)
    {
        const OctNode<realType>* node = prQueue.top().second;
        prQueue.pop();
        for (const OctNode<realType>* child : node->children)
        {
            if (child)
            {
                prQueue.emplace(distToNodeSqrd(*child, point), child);
            }
        }
    }
    return prQueue.top().second->points[0];
}

template<typename realType>
bool Octree<realType>::meiosis(OctNode<realType> &inNode, realType tol)
{
    size_t inNodeNumPoints = inNode.points.size();
    if (inNode.points.size() == 0)
    {
        inNode.isPoint = true;
        inNode.bounds[0] = m_points[inNode.points[0]];
        return false;
    }
    if (inNode.points.size() < 9)
    {
        // We make a point child node for each point instead of splitting the box into subboxes
        inNode.children.resize(8);
        for (size_t i = 0; i < inNodeNumPoints; ++i)
        {
            size_t point = inNode.points[i];
            const Eigen::Vector3<realType>& pointCoord = m_points[point];
            int matchedChildIdx = -1;
            for (size_t j = 0; j < i; ++j)
            {
                if (pointCoord == m_points[inNode.points[j]])
                {
                    matchedChildIdx = j;
                    break;
                }
            }
            if (matchedChildIdx == -1)
            {
                inNode.children[i] = new OctNode<realType>();
                inNode.children[i]->isPoint = true;
                inNode.children[i]->bounds[0] = pointCoord;
                inNode.children[i]->points.push_back(point);
            }
            else
            {
                inNode.children[i] = nullptr;
                inNode.children[matchedChildIdx]->points.push_back(point);
            }
        }
        // No need to further split node
        return false;
    }

    // Check if node is too small to split
    unsigned char tooSmallDims = 0;
    for (int d = 0; d < 3; ++d)
    {
        if (inNode.bounds[1](0) - inNode.bounds[0](0) < tol)
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

    std::vector<std::vector<size_t>> pointsPerChild(8);
    for (size_t pointIdx : inNode.points)
    {
        const Eigen::Vector3<realType>& point = m_points[pointIdx];
        size_t cat = (point(0) > seps(0) ? 1 : 0) + (point(1) > seps(1) ? 2 : 0) + (point(2) > seps(2) ? 4 : 0);
        pointsPerChild[cat].push_back(pointIdx);
    }

    bool changed = false;
    inNode.children.resize(8);
    for (size_t cat = 0; cat < 8; ++cat)
    {
        if (pointsPerChild[cat].empty())
        {
            inNode.children[cat] = nullptr;
        }
        else
        {
            inNode.children[cat] = new OctNode<realType>;
            size_t iCopy = cat;
            for (int d = 0; d < 3; ++d)
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
            inNode.children[cat]->points = pointsPerChild[cat];
        }
    }

    return true;
}

template<typename realType>
void Octree<realType>::recursiveMeiosis(OctNode<realType> &inNode, realType tol)
{
    if (!meiosis(inNode, tol))
    {
        return;
    }
    for (OctNode<realType>* newChild : inNode.children)
    {
        if (newChild)
        {
            recursiveMeiosis(*newChild, tol);
        }
    }
}

}

#endif // OCTREE_H
