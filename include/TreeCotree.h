#pragma once // 确保 Windows 下头文件只包含一次

#include <queue>
#include <vector>
#include <algorithm> // 用于 std::swap
#include <cassert>   // 用于 assert

// 如果在很旧的 VS 版本上编译，取消下面这行的注释
// #include <iso646.h> 

namespace DDG
{
    class TreeCotree
    {
    public:
        // 假设这些类型在 Mesh 类中已定义
        // 为了兼容，确保 Mesh 结构体在外部正确实现
        void build(Mesh& mesh)
        {
            mesh.generators.clear();
            buildPrimalSpanningTree(mesh);
            buildDualSpanningCoTree(mesh);
            buildCycles(mesh);
        }

    protected:
        void buildCycles(Mesh& mesh)
        {
            for (EdgeIter e = mesh.edges.begin();
                e != mesh.edges.end();
                e++)
            {
                // 将 'or' 改为 '||'，'and' 改为 '&&' 以增强 MSVC 兼容性
                if (e->he->onBoundary || e->he->flip->onBoundary) continue;
                if (inPrimalSpanningTree(e->he)) continue;
                if (inDualSpanningTree(e->he)) continue;

                // path to root from right face
                Mesh::Generator c1;
                FaceIter f = e->he->flip->face;
                while (f != f->parent)
                {
                    c1.push_back(sharedHalfEdge(f, f->parent));
                    f = f->parent;
                }

                // path to root from left face
                Mesh::Generator c2;
                f = e->he->face;
                while (f != f->parent)
                {
                    c2.push_back(sharedHalfEdge(f, f->parent));
                    f = f->parent;
                }

                // build cycle
                Mesh::Generator g;
                g.push_back(e->he);

                // 使用 int 确保减法安全
                int m = (int)c1.size() - 1;
                int n = (int)c2.size() - 1;

                while (m >= 0 && n >= 0 && c1[m] == c2[n]) { m--; n--; }

                for (int i = 0; i <= m; i++) g.push_back(c1[i]);
                for (int i = n; i >= 0; i--) g.push_back(c2[i]->flip);

                // make sure that boundary loops wind around the boundary
                if (isDualBoundaryLoop(g))
                {
                    if (g[0]->next->vertex->onBoundary())
                    {
                        // reverse cycle
                        size_t num_elements = g.size();
                        for (size_t i = 0; i < num_elements; i++)
                        {
                            g[i] = g[i]->flip;
                        }
                        // 使用标准 std::swap
                        for (size_t i = 0; i < num_elements / 2; i++)
                        {
                            std::swap(g[i], g[num_elements - 1 - i]);
                        }
                    }
                }

                mesh.generators.push_back(g);
            }
        }

        void buildPrimalSpanningTree(Mesh& mesh)
        {
            for (VertexIter v = mesh.vertices.begin();
                v != mesh.vertices.end();
                v++)
            {
                v->parent = v;
            }

            VertexIter root = mesh.vertices.begin();
            // 寻找第一个非边界顶点作为根
            while (root != mesh.vertices.end() && root->onBoundary()) root++;
            if (root == mesh.vertices.end()) return;

            std::queue<VertexIter> Q;
            Q.push(root);
            while (!Q.empty())
            {
                VertexIter v = Q.front();
                Q.pop();

                HalfEdgeIter he = v->he;
                do
                {
                    VertexIter w = he->flip->vertex;
                    if (w->parent == w && w != root && !w->onBoundary())
                    {
                        w->parent = v;
                        Q.push(w);
                    }
                    he = he->flip->next;
                } while (he != v->he);
            }
        }

        void buildDualSpanningCoTree(Mesh& mesh)
        {
            for (FaceIter f = mesh.faces.begin();
                f != mesh.faces.end();
                f++)
            {
                f->parent = f;
            }

            FaceIter root = mesh.faces.begin();
            while (root != mesh.faces.end() && root->isBoundary()) root++;
            if (root == mesh.faces.end()) return;

            std::queue<FaceIter> Q;
            Q.push(root);
            while (!Q.empty())
            {
                FaceIter f = Q.front();
                Q.pop();

                HalfEdgeIter he = f->he;
                do
                {
                    FaceIter g = he->flip->face;
                    if (g->parent == g && g != root && !g->isBoundary() &&
                        !inPrimalSpanningTree(he))
                    {
                        g->parent = f;
                        Q.push(g);
                    }
                    he = he->next;
                } while (he != f->he);
            }
        }

        bool inPrimalSpanningTree(HalfEdgeIter he) const
        {
            VertexIter v = he->vertex;
            VertexIter w = he->flip->vertex;
            return v->parent == w || w->parent == v;
        }

        bool inDualSpanningTree(HalfEdgeIter he) const
        {
            FaceIter f = he->face;
            FaceIter g = he->flip->face;
            return f->parent == g || g->parent == f;
        }

        HalfEdgeIter sharedHalfEdge(FaceIter f, FaceIter g) const
        {
            HalfEdgeIter he = f->he;
            do
            {
                if (he->flip->face == g) return he;
                he = he->next;
            } while (he != f->he);
            assert(false);
            return f->he; // 避免编译器警告
        }

        bool isDualBoundaryLoop(const Mesh::Generator& cycle) const
        {
            if (cycle.empty()) return false;
            return (cycle[0]->vertex->onBoundary() ||
                cycle[0]->flip->vertex->onBoundary());
        }
    };
}