#pragma once

namespace FML {

template <typename realType>
class Vertex;

template <typename realType>
class Halfedge
{
public:
    Vertex<realType>* destination = nullptr;
    Halfedge<realType>* next = nullptr;
    Halfedge<realType>* prev = nullptr;
    Halfedge<realType>* neighbor = nullptr;

public:
    Halfedge()
    {}

    Halfedge(Vertex<realType>* destination) :
        destination(destination)
    {}

};

}
