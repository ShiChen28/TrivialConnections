#ifndef DDG_DIRECTION_H
#define DDG_DIRECTION_H

// 1. 确保在 Windows 上定义 M_PI
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <iostream>
#include <stack>      // 引入栈用于迭代处理，防止栈溢出
#include <algorithm>  // std::abs

#include "Mesh.h"
#include "Quaternion.h"

// 兼容某些旧版 MSVC 编译器对 and/or 的支持
#include <iso646.h>

namespace DDG
{
    class DirectionField
    {
    public:
        void generate(Mesh& mesh, double angle, bool debugMode = false)
        {
            setFaceTag(mesh, false);

            if (mesh.faces.empty()) return;

            FaceIter f = mesh.faces.begin();
            init(f, angle);

            // 2. 将递归 propagate 改为迭代，防止 Windows 栈溢出
            std::stack<std::pair<HalfEdgeIter, Vector>> s;

            HalfEdgeIter he = f->he;
            do
            {
                s.push(std::make_pair(he->flip, f->vector));
                he = he->next;
            } while (he != f->he);

            while (!s.empty())
            {
                std::pair<HalfEdgeIter, Vector> curr = s.top();
                s.pop();
                propagateIterative(mesh, curr.first, curr.second, s);
            }

            if (debugMode) debug(mesh);
        }

    protected:
        void setFaceTag(Mesh& mesh, bool flag)
        {
            for (FaceIter f = mesh.faces.begin();
                f != mesh.faces.end();
                f++)
            {
                f->tag = flag;
            }
        }

        void init(FaceIter f, double angle)
        {
            VertexIter v0 = f->he->vertex;
            VertexIter v1 = f->he->next->vertex;
            Vector v = (v1->position - v0->position).unit();

            f->vector = rotate(v, f->normal(), angle);
            f->tag = true;
        }

        // 3. 迭代版本的 propagate
        void propagateIterative(const Mesh& mesh, HalfEdgeIter h, const Vector& v,
            std::stack<std::pair<HalfEdgeIter, Vector>>& s)
        {
            if (h->onBoundary) return;
            if (h->face->tag) return;

            FaceIter f = h->face;
            f->vector = transport(mesh, v, h);
            f->tag = true;

            HalfEdgeIter he = f->he;
            do
            {
                // 将后续任务压入栈，而不是递归调用
                if (!he->flip->face->tag) {
                    s.push(std::make_pair(he->flip, f->vector));
                }
                he = he->next;
            } while (he != f->he);
        }

        Vector transport(const Mesh& mesh, const Vector& v, HalfEdgeIter h) const
        {
            // 4. 将 or 改为 ||
            if (h->onBoundary || h->flip->onBoundary) return v;

            Vector p1 = h->flip->vertex->position;
            Vector p0 = h->vertex->position;
            Vector e = (p1 - p0).unit();

            Vector nR = h->flip->face->normal();
            Vector nL = h->face->normal();
            Vector nLxnR = cross(nL, nR);

            // 5. 使用 std::atan2
            double dihedral = std::atan2(nLxnR.norm(), dot(nL, nR));
            if (dot(e, nLxnR) < 0.0) dihedral = -dihedral;

            Vector u = rotate(v, nR, -mesh.connectionOneForm(h));
            return rotate(u, e, dihedral);
        }

        Vector rotate(const Vector& v, const Vector& axis, double angle) const
        {
            angle *= 0.5;
            Quaternion q(std::cos(angle), std::sin(angle) * axis);
            return (q.conj() * Quaternion(v) * q).im();
        }

        void debug(Mesh& mesh) const
        {
            std::cout << "DEBUG-BEGIN" << std::endl;

            for (VertexIter v = mesh.vertices.begin();
                v != mesh.vertices.end();
                v++)
            {
                if (v->onBoundary()) continue;

                Vector n = v->he->flip->face->normal();
                Vector u0 = v->he->flip->face->vector;

                Vector u = u0;
                HalfEdgeIter h = v->he;
                do
                {
                    u = transport(mesh, u, h);
                    h = h->next->next->flip;
                } while (h != v->he);

                double offset = std::atan2(cross(u, u0).norm(), dot(u, u0));
                if (dot(n, cross(u, u0)) > 0.0) offset = -offset;

                // 6. 确保 M_PI 可用
                while (offset < 0.0) offset += 2. * M_PI;
                while (offset >= 2.0 * M_PI) offset -= 2. * M_PI;
                if (std::abs(offset - 2.0 * M_PI) < 1.0e-6) offset = 0.0;

                double sing = mesh.vertexHolonomy(v);

                while (sing < 0.0) sing += 2. * M_PI;
                while (sing >= 2.0 * M_PI) sing -= 2. * M_PI;
                if (std::abs(sing - 2.0 * M_PI) < 1.0e-6) sing = 0.0;

                if (std::abs(offset - sing) > 1.0e-8)
                {
                    std::cout << "Vertex" << v->index << ": "
                        << "offset = " << (180. / M_PI) * offset << " ; "
                        << "sing = " << (180. / M_PI) * sing
                        << std::endl;
                }
            }

            // 7. 使用 size_t 避免警告
            for (size_t i = 0; i < mesh.generators.size(); ++i)
            {
                const Mesh::Generator& cycle = mesh.generators[i];
                if (cycle.empty()) continue;

                Vector n = cycle[0]->face->normal();
                Vector u0 = cycle[0]->face->vector;

                Vector u = u0;
                // 8. 注意 size() 是 unsigned，减法要小心溢出
                for (int k = (int)cycle.size() - 1; k >= 0; k--)
                {
                    u = transport(mesh, u, cycle[k]);
                }

                double offset = std::atan2(cross(u, u0).norm(), dot(u, u0));
                if (dot(n, cross(u, u0)) > 0.0) offset = -offset;

                while (offset < 0.0) offset += 2. * M_PI;
                while (offset >= 2.0 * M_PI) offset -= 2. * M_PI;
                if (std::abs(offset - 2.0 * M_PI) < 1.0e-6) offset = 0.0;

                double sing = mesh.generatorHolonomy(cycle);

                while (sing < 0.0) sing += 2. * M_PI;
                while (sing >= 2.0 * M_PI) sing -= 2. * M_PI;
                if (std::abs(sing - 2.0 * M_PI) < 1.0e-6) sing = 0.0;

                if (std::abs(offset - sing) > 1.0e-8)
                {
                    std::cout << "Generator" << (unsigned)i << ": "
                        << "offset = " << (180. / M_PI) * offset << " ; "
                        << "sing = " << (180. / M_PI) * sing
                        << std::endl;
                }
            }

            std::cout << "DEBUG-END" << std::endl;
        }
    };
}

#endif