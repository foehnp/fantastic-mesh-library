#ifndef MESH_H
#define MESH_H

#include "Halfedge.h"
#include "Triangle.h"
#include "Vertex.h"
#include "LineCol.h"

#include <unordered_set>
#include <limits>


namespace FML {


template <typename realType>
class Mesh
{
private:
    std::vector<Vertex<realType>> vertices;
    std::vector<Triangle<realType>> triangles;
    std::vector<Halfedge<realType>> halfedges;
    bool halfedgesComputed = false;
    bool normalsComputed = false;


public:
    Mesh(){};

    Mesh(const std::vector<Eigen::Vector3<realType>>& vertexVec, const std::vector<std::vector<size_t>>& triangleVec);

    size_t numTriangles() const;

    Triangle<realType>* getTriangle(size_t idx);

    const Triangle<realType>* getTriangle(size_t idx) const;

    size_t numVertices() const;

    Vertex<realType>* getVertex(size_t idx);

    const Vertex<realType>* getVertex(size_t idx) const;

    std::vector<Eigen::Vector3<realType>> getVerticesAsVectors() const;

    void computeNormals();

    void computeHalfEdges();

    LineCol<realType> computeBoundary() const;

    std::vector<Eigen::Vector3<realType>> computeBoundingBox() const;
};


template<typename realType>
Mesh<realType>::Mesh(const std::vector<Eigen::Vector3<realType> > &vertexVec, const std::vector<std::vector<size_t> > &triangleVec)
{
    vertices.reserve(vertexVec.size());
    for (const auto& vec : vertexVec)
    {
        vertices.emplace_back(vec);
    }
    triangles.reserve(triangleVec.size());
    for (size_t i = 0; i < triangleVec.size(); ++i)
    {
        const std::vector<size_t>& t = triangleVec[i];
        triangles.emplace_back(t[0], t[1], t[2]);
    }
}

template<typename realType>
inline size_t Mesh<realType>::numTriangles() const
{
    return triangles.size();
}

template<typename realType>
inline Triangle<realType> *Mesh<realType>::getTriangle(size_t idx)
{
    return &triangles[idx];
}

template<typename realType>
inline const Triangle<realType> *Mesh<realType>::getTriangle(size_t idx) const
{
    return &triangles[idx];
}

template<typename realType>
inline size_t Mesh<realType>::numVertices() const
{
    return triangles.size();
}

template<typename realType>
inline Vertex<realType> *Mesh<realType>::getVertex(size_t idx)
{
    return &vertices[idx];
}

template<typename realType>
inline const Vertex<realType> *Mesh<realType>::getVertex(size_t idx) const
{
    return &vertices[idx];
}

template<typename realType>
inline std::vector<Eigen::Vector3<realType> > Mesh<realType>::getVerticesAsVectors() const
{
    std::vector<Eigen::Vector3<realType> > res(vertices.size());
    for (int i = 0; i < vertices.size(); ++i)
    {
        res[i] = vertices[i];
    }
    return res;
}

template<typename realType>
void Mesh<realType>::computeNormals()
{
    if (normalsComputed)
    {
        return;
    }

    for (const Triangle<realType>& t : triangles)
    {
        t.computeNormal();
    }

    normalsComputed = true;
}

template<typename realType>
void Mesh<realType>::computeHalfEdges()
{
    if (halfedgesComputed)
    {
        return;
    }

    halfedges.clear();
    halfedges.reserve(vertices.size()*6);
    // (startVertex, (endVertex, halfEdge)):
    std::vector<std::unordered_map<size_t, Halfedge<realType>*>> perVertexSecondVertexToHE(vertices.size());
    for (Triangle<realType>& t : triangles)
    {
        Vertex<realType>* tVertices[3];
        Halfedge<realType>* tHalfedges[3];
        Halfedge<realType>* tNeighborHalfedges[3];
        for (size_t i = 0; i < 3; ++i)
        {
            size_t thisVertIdx = t.vertex(i);
            Vertex<realType>* thisVert = &vertices[thisVertIdx];
            size_t nextVertIdx = t.vertex((i+1)%3);
            Vertex<realType>* nextVert = &vertices[nextVertIdx];

            tVertices[i] = thisVert;
            tHalfedges[i] = &(halfedges.emplace_back(nextVert));

            auto finder = perVertexSecondVertexToHE[nextVertIdx].find(thisVertIdx);
            if (finder != perVertexSecondVertexToHE[nextVertIdx].end())
            {
                tNeighborHalfedges[i] = finder->second;
                finder->second->neighbor = tHalfedges[i];
            }
            else
            {
                tNeighborHalfedges[i] = nullptr;
            }
            (perVertexSecondVertexToHE[thisVertIdx])[nextVertIdx] = tHalfedges[i];
        }
        for (size_t i = 0; i < 3; ++i)
        {
            tHalfedges[i]->next = tHalfedges[(i+1)%3];
            tHalfedges[i]->prev = tHalfedges[(i+2)%3];
            tHalfedges[i]->neighbor = tNeighborHalfedges[i];
        }
    }
    halfedgesComputed = true;
}

template<typename realType>
LineCol<realType> Mesh<realType>::computeBoundary() const
{
    LineCol<realType> lineCol;

    std::unordered_set<const Halfedge<realType>*> halfedgeUsed;
    for (int i = 0; i < halfedges.size(); ++i)
    {
        if ((halfedges[i].neighbor != nullptr)  || halfedgeUsed.count(&halfedges[i]) > 0)
        {
            continue;
        }
        std::vector<Eigen::Vector3<realType>> line;
        const Halfedge<realType>* startHe = &halfedges[i];
        const Halfedge<realType>* currHe = startHe;
        halfedgeUsed.insert(currHe);
        do
        {
            line.push_back(*currHe->destination);

            const Halfedge<realType>* inpointing = currHe->next;
            const Halfedge<realType>* neighbor = inpointing->neighbor;
            while (neighbor)
            {
                inpointing = neighbor->next;
                neighbor = inpointing->neighbor;
            }
            halfedgeUsed.insert(currHe);
            currHe = inpointing;
        } while (currHe != startHe);

        lineCol.addLineMove(line, true);
    }
    return lineCol;
}

template<typename realType>
std::vector<Eigen::Vector3<realType>> Mesh<realType>::computeBoundingBox() const
{
    realType infty = std::numeric_limits<realType>::max();

    std::vector<Eigen::Vector3<realType>> res(2);
    res[0] = {infty, infty, infty};
    res[1] = {-infty, -infty, -infty};

    for (const Vertex<realType>& v : vertices)
    {
        for (int d = 0; d < 3; ++d)
        {
            res[0](d) = std::min(res[0](d), v(d));
            res[1](d) = std::max(res[1](d), v(d));
        }
    }
    return res;
}

}

#endif // MESH_H
