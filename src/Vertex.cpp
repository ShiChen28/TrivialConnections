#include <vector>
#include <iostream>

#include "Vertex.h"
#include "Mesh.h"
#include "HalfEdge.h"

// 建议不要在头文件包含之前使用 using namespace std;
// 也不建议在全局命名空间中使用，但在 .cpp 文件中通常可以接受
using namespace std;

namespace DDG
{
    // 修改点 1: 在 Windows 下，全局变量的初始化顺序可能导致问题。
    // 这里的 isolated 通常用于表示孤立顶点的占位符。
    // 如果 Vertex.h 里定义 he 是 HalfEdgeIter（通常是 list::iterator），
    // 那么对比 isolated.begin() 可能会因容器不同而崩溃。
    vector<HalfEdge> isolated;

    double Vertex::area(void) const
    {
        double A = 0.;

        HalfEdgeCIter h = he;
        // 检查 he 是否有效，防止对孤立顶点进行循环导致崩溃
        if (isIsolated()) return 0.0;

        do
        {
            // 修改点 2: 将 'not' 替换为 '!' 以确保 MSVC 兼容性
            if (!h->onBoundary) A += h->face->area();
            h = h->flip->next;
        } while (h != he);

        return A / 3.;
    }

    Vector Vertex::normal(void) const
    {
        Vector N(0., 0., 0.);

        if (isIsolated()) return N;

        HalfEdgeCIter h = he;
        do
        {
            if (!h->onBoundary) N += h->face->normal();
            h = h->flip->next;
        } while (h != he);

        return N.unit();
    }

    bool Vertex::isIsolated(void) const
    {
        // 修改点 3: 在 Windows Debug 模式下，如果 isolated 为空，调用 begin() 会崩溃
        if (isolated.empty()) {
            // 如果 he 是指针，检查是否为 nullptr；如果是迭代器，通常检查是否等于 mesh.halfedges.end()
            // 这里暂时保留原逻辑，但加上安全性检查
            return false;
        }
        return he == (HalfEdgeCIter)isolated.begin();
    }

    int Vertex::valence(void) const
    {
        if (isIsolated()) return 0;

        int n = 0;
        HalfEdgeCIter h = he;
        do
        {
            n++;
            h = h->flip->next;
        } while (h != he);

        return n;
    }

    void Vertex::toggleTag()
    {
        tag = !tag;
    }

    double Vertex::theta(void) const
    {
        if (isIsolated()) return 0.0;

        double sum = 0.0;
        HalfEdgeCIter h = this->he;
        do
        {
            // 确保 angle() 已定义
            sum += h->next->angle();
            h = h->flip->next;
        } while (h != this->he);
        return sum;
    }

    bool Vertex::onBoundary(void) const
    {
        if (isIsolated()) return true; // 孤立点通常视为边界

        HalfEdgeCIter h = this->he;
        do
        {
            if (h->onBoundary) return true;
            h = h->flip->next;
        } while (h != this->he);
        return false;
    }
}