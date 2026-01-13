#ifndef DDG_HARMONIC_BASES_H
#define DDG_HARMONIC_BASES_H

#pragma once // Windows 推荐使用的头文件保护

#include "Mesh.h"
#include "Real.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "DiscreteExteriorCalculus.h"

// 包含此头文件以确保 'and', 'or', 'not' 在旧版编译器中可用，
// 或者直接将它们替换为 &&, ||, !
#include <iso646.h> 

namespace DDG
{
    class HarmonicBases
    {
    public:
        void compute(Mesh& mesh)
        {
            cleanHalfEdges(mesh);
            if (mesh.generators.empty()) return;

            SparseMatrix<Real> star1, d0, div;
            HodgeStar1Form<Real>::build(mesh, star1);
            ExteriorDerivative0Form<Real>::build(mesh, d0);

            // 矩阵转置和乘法：确保 SparseMatrix 类在 Windows 下正确实现了这些重载
            div = d0.transpose() * star1;

            bool skipBoundaryLoop = true;
            // 将 unsigned 改为 size_t 以匹配 .size() 的返回类型，防止警告
            for (size_t i = 0; i < mesh.generators.size(); ++i)
            {
                const Mesh::Generator& cycle = mesh.generators[i];

                // 将 'and' 改为 '&&'
                if (skipBoundaryLoop && mesh.isBoundaryGenerator(cycle))
                {
                    skipBoundaryLoop = false;
                    continue;
                }

                DenseMatrix<Real> w;
                buildClosedPrimalOneForm(mesh, cycle, w);

                DenseMatrix<Real> u;
                DenseMatrix<Real> divw = div * w;

                // 注意：backsolvePositiveDefinite 通常依赖于线性代数库（如 Cholmod 或 Eigen）
                // 在 Windows 上，请确保对应的库（如 SuiteSparse）已正确链接
                backsolvePositiveDefinite(mesh.L, u, divw);

                // 计算谐波形式: h = w - du
                // 这里使用星号算子进行对偶转换
                DenseMatrix<Real> h = star1 * (w - (d0 * u));
                storeHarmonicForm(h, mesh);
            }
        }

    protected:
        void buildClosedPrimalOneForm(const Mesh& mesh,
            const Mesh::Generator& cycle,
            DenseMatrix<Real>& oneform) const
        {
            // 确保 DenseMatrix 的构造函数在 Windows 下工作正常
            oneform = DenseMatrix<Real>((int)mesh.edges.size());
            for (size_t i = 0; i < cycle.size(); ++i)
            {
                Real value = (Real)1.0;
                HalfEdgeIter he = cycle[i];

                // 判断半边方向与边方向是否一致
                if (he->edge->he != he) value = -value;

                // oneform(index) 必须在 DenseMatrix 中有定义
                oneform((int)he->edge->index) = value;
            }
        }

        void storeHarmonicForm(const DenseMatrix<Real>& oneform,
            Mesh& mesh)
        {
            for (EdgeIter e = mesh.edges.begin();
                e != mesh.edges.end();
                e++)
            {
                Real value = oneform((int)e->index);
                // 确保 harmonicBases 是 std::vector<Real>
                e->he->harmonicBases.push_back(value);
                e->he->flip->harmonicBases.push_back(-value);
            }
        }

        void cleanHalfEdges(Mesh& mesh)
        {
            for (EdgeIter e = mesh.edges.begin();
                e != mesh.edges.end();
                e++)
            {
                e->he->harmonicBases.clear();
                e->he->flip->harmonicBases.clear();
            }
        }
    };
}

#endif