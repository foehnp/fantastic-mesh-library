#ifndef OCTREE_H
#define OCTREE_H

#include <Eigen/Core>

#include <forward_list>
#include <vector>

namespace FML {

namespace {

template<typename realType>
struct OctNode
{
    std::vector<Eigen::Vector3<realType>> bounds;
    std::vector<OctNode<realType>*> children;
    std::vector<size_t> points;
    OctNode() :
        bounds(2)
    {}

    ~OctNode()
    {
        for (const auto child : children)
        {
            delete child;
        }
    }
};

}



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


    size_t nearestPoint(const Eigen::Vector3<realType>& point, realType maxDist) const;

private:

    /* subdivides inNode into 8 subnodes if it contains more than one points and at least
     *  one of the sides of its cuboid is wider than tol
     *
     * returns  false if no subdivision had to be done
     **/
    bool meiosis(OctNode<realType>& inNode, realType tol);

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

namespace {

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

}

template<typename realType>
size_t Octree<realType>::nearestPoint(const Eigen::Vector3<realType> &point, realType maxDist) const
{
    std::forward_list<const OctNode<realType>*> nodes;
    nodes.push_front(&m_rootNode);
    while (true)
    {
        realType minMaxDist = std::numeric_limits<realType>::max();
        for (const OctNode<realType>* node : nodes)
        {
            minMaxDist = std::min(minMaxDist, maxDistOfBoxSqrd(node->bounds, point));
        }

        std::forward_list<const OctNode<realType>*> newNodes;
        bool anyNewNodes = false;
        for (const OctNode<realType>* node : nodes)
        {
            if (distToBoxSqrd(node->bounds, point) < minMaxDist)
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
            size_t closestPoint;
            realType closestDist = std::numeric_limits<realType>::max();
            for (const OctNode<realType>* node : nodes)
            {
                size_t currPointIdx = node->points[0];
                realType dist = (m_points[currPointIdx] - point).squaredNorm();
                if (dist < closestDist)
                {
                    closestPoint = currPointIdx;
                    closestDist = dist;
                }
            }
            return closestPoint;
        }

        nodes = std::move(newNodes);
    }
    return SIZE_MAX;
}

template<typename realType>
bool Octree<realType>::meiosis(OctNode<realType> &inNode, realType tol)
{
    if (inNode.points.size() < 2)
    {
        // No need to further split node
        return false;
    }

    // Check if node is too small to split
    unsigned char tooSmallDims = 0;
    for (int d = 0; d < 3; ++d)
    {
        if (inNode.bounds[1](0) - inNode.bounds[0](0) < tol)
        {
            ++d;
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