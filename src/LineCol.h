#ifndef LINECOL_H
#define LINECOL_H

#include <Eigen/Dense>
#include <memory>
#include <vector>


namespace FML {

template <typename realType>
class LineCol
{
private:
    std::vector<std::pair<std::vector<Eigen::Vector3<realType>>, bool>> lines;
    bool closed = false;

public:
    LineCol(){}

    void addLine(const std::vector<Eigen::Vector3<realType>>& points, bool closed = false);

    void addLineMove(std::vector<Eigen::Vector3<realType>>& points, bool closed = false);

    size_t numLines() const;

    std::pair<std::vector<Eigen::Vector3<realType>>, bool> getLine(size_t idx) const;

};

template<typename realType>
inline void LineCol<realType>::addLine(const std::vector<Eigen::Vector3<realType>>& points, bool closed)
{
    lines.emplace_back(points, closed);
}

template<typename realType>
inline void LineCol<realType>::addLineMove(std::vector<Eigen::Vector3<realType>>& points, bool closed)
{
    lines.emplace_back(std::move(points), closed);
}

template<typename realType>
inline size_t LineCol<realType>::numLines() const
{
    return lines.size();
}

template<typename realType>
inline std::pair<std::vector<Eigen::Vector3<realType>>, bool> LineCol<realType>::getLine(size_t idx) const
{
    return lines.at(idx);
}

}

#endif // LINECOL_H
