#ifndef DDG_APPLICATION_H
#define DDG_APPLICATION_H

// 1. 确保定义 M_PI
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <cmath>
#include <ctime>
#include <iostream>
#include <algorithm> // 必须包含以使用 std::min

#include "Mesh.h"
#include "Real.h"
#include "Quaternion.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "DiscreteExteriorCalculus.h"
#include "Utility.h"

// 2. 增加逻辑运算符兼容性
#include <iso646.h> 

namespace DDG
{
    class Application
    {
    public:
        bool solveForConnection(Mesh& mesh)
        {
            std::cout << "      ------------------" << std::endl;
            std::cout << "      computing connection" << std::endl;

            bool ok = checkGaussBonnet(mesh);
            if (!ok) // 将 not 改为 !
            {
                std::cout << "Gauss-Bonnet thm does not hold" << std::endl;
                return false;
            }

            std::cout << "      ------------------" << std::endl;
            std::cout << "      computing trivial holonomy" << std::endl;
            clock_t t0 = clock();
            solveForTrivialHolonomy(mesh);
            clock_t t1 = clock();
            // 修正时间计算，防止 seconds 函数未定义或不兼容
            std::cout << "[trivial] time: " << (double)(t1 - t0) / CLOCKS_PER_SEC << "s" << "\n";

            std::cout << "      ------------------" << std::endl;
            std::cout << "      computing nontrivial holonomy" << std::endl;
            t0 = clock();
            solveForNonTrivialHolonomy(mesh);
            t1 = clock();
            std::cout << "[nontrivial] time: " << (double)(t1 - t0) / CLOCKS_PER_SEC << "s" << "\n";

            return true;
        }

        double solveForGeodesic(double dt, Mesh& mesh)
        {
            std::cout << "      ------------------" << std::endl;
            std::cout << "      initial condition" << std::endl;
            DenseMatrix<Real> u0;
            int nb = builImpulseSignal(mesh, u0);
            if (nb == 0) return 1.0;

            std::cout << "      ------------------" << std::endl;
            std::cout << "      Build Hodge * 0Form" << std::endl;
            SparseMatrix<Real> star0;
            HodgeStar0Form<Real>::build(mesh, star0);

            std::cout << "      ------------------" << std::endl;
            std::cout << "      Build Hodge * 1Form" << std::endl;
            SparseMatrix<Real> star1;
            HodgeStar1Form<Real>::build(mesh, star1);

            std::cout << "      ------------------" << std::endl;
            std::cout << "      Build Exterior Derivative 0Form" << std::endl;
            SparseMatrix<Real> d0;
            ExteriorDerivative0Form<Real>::build(mesh, d0);

            SparseMatrix<Real> L = d0.transpose() * star1 * d0;
            L += Real(1.0e-8) * star0;

            // 修正：sqr 可能在某些环境下未定义，可以改为 x*x
            double meanLen = mesh.meanEdgeLength();
            dt *= (meanLen * meanLen);
            SparseMatrix<Real> A = star0 + Real(dt) * L;

            DenseMatrix<Real> u;
            solvePositiveDefinite(A, u, u0);

            computeVectorField(u, mesh);

            DenseMatrix<Real> div;
            computeDivergence(mesh, div);

            DenseMatrix<Real> phi;
            solvePositiveDefinite(L, phi, div);

            setMinToZero(phi);
            assignDistance(phi, mesh);
            return (double)phi.norm();
        }

    protected:

        bool checkGaussBonnet(const Mesh& mesh) const
        {
            int k = 0;
            for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
            {
                k += v->singularity;
            }
            k += (int)mesh.firstGeneratorIndex;

            return (mesh.getEulerCharacteristicNumber() == k);
        }

        void solveForTrivialHolonomy(Mesh& mesh)
        {
            DenseMatrix<Real> b((int)mesh.vertices.size());
            for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
            {
                double value = 0.0;
                if (!v->onBoundary()) // 将 not 改为 !
                {
                    value -= (2. * M_PI - v->theta());
                    value += 2. * M_PI * v->singularity;
                }
                b((int)v->index) = value;
            }

            DenseMatrix<Real> u((int)mesh.vertices.size());
            if (b.norm() > 1.0e-8) backsolvePositiveDefinite(mesh.L, u, b);

            for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
            {
                v->potential = (double)u((int)v->index);
            }
        }

        void solveForNonTrivialHolonomy(Mesh& mesh)
        {
            unsigned nb = mesh.numberHarmonicBases();
            if (nb == 0) return;
            mesh.harmonicCoefs = std::vector<double>(nb, 0.0);

            DenseMatrix<Real>  b((int)nb);
            SparseMatrix<Real> H((int)nb, (int)nb);

            int row = 0;
            bool skipBoundaryLoop = true;
            for (size_t i = 0; i < mesh.generators.size(); ++i)
            {
                const Mesh::Generator& cycle = mesh.generators[i];
                if (skipBoundaryLoop && mesh.isBoundaryGenerator(cycle)) // 将 and 改为 &&
                {
                    skipBoundaryLoop = false;
                    continue;
                }

                for (size_t j = 0; j < cycle.size(); ++j)
                {
                    HalfEdgeIter he = cycle[j];
                    for (unsigned col = 0; col < nb; ++col)
                    {
                        H(row, (int)col) += he->harmonicBases[col];
                    }
                }

                double value = -mesh.generatorHolonomy(cycle);
                if (row == 0)
                {
                    value += 2.0 * M_PI * mesh.firstGeneratorIndex;
                }
                b(row) = value;
                row++;
            }

            DenseMatrix<Real> x((int)nb);
            if (b.norm() > 1.0e-8) solve(H, x, b);

            for (unsigned i = 0; i < nb; ++i)
                mesh.harmonicCoefs[i] = (double)x((int)i);
        }

        int builImpulseSignal(const Mesh& mesh, DenseMatrix<Real>& x) const
        {
            int count = 0;
            x = DenseMatrix<Real>((int)mesh.vertices.size());
            for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
            {
                x((int)v->index) = 0.0;
                if (v->tag)
                {
                    x((int)v->index) = 1.0;
                    count++;
                }
            }
            return count;
        }

        void computeVectorField(const DenseMatrix<Real>& u, Mesh& mesh)
        {
            for (FaceIter f = mesh.faces.begin(); f != mesh.faces.end(); f++)
            {
                if (f->isBoundary()) continue;

                HalfEdgeIter hij = f->he;
                HalfEdgeIter hjk = hij->next;
                HalfEdgeIter hki = hjk->next;

                double ui = u((int)hij->vertex->index);
                double uj = u((int)hjk->vertex->index);
                double uk = u((int)hki->vertex->index);

                Vector eij90 = hij->rotatedEdge();
                Vector ejk90 = hjk->rotatedEdge();
                Vector eki90 = hki->rotatedEdge();

                Vector X = 0.5 * (ui * ejk90 + uj * eki90 + uk * eij90) / f->area();
                f->vector = -X.unit();
            }
        }

        void computeDivergence(const Mesh& mesh, DenseMatrix<Real>& div) const
        {
            div = DenseMatrix<Real>((int)mesh.vertices.size());
            for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
            {
                double sum = 0.0;
                HalfEdgeIter he = v->he;
                do
                {
                    if (!he->onBoundary) // 将 not 改为 !
                    {
                        Vector n = he->next->rotatedEdge();
                        Vector vec = he->face->vector;
                        sum += dot(n, vec);
                    }
                    he = he->flip->next;
                } while (he != v->he);
                div((int)v->index) = sum;
            }
        }

        void assignDistance(const DenseMatrix<Real>& phi, Mesh& mesh)
        {
            for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++)
            {
                v->distance = (double)phi((int)v->index);
            }
        }

        void setMinToZero(DenseMatrix<Real>& phi) const
        {
            double minValue = 1.0e100;
            int rows = phi.nRows();
            int cols = phi.nColumns();
            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    minValue = (std::min)(minValue, (double)phi(i, j)); // 使用 (std::min)

            for (int i = 0; i < rows; ++i)
                for (int j = 0; j < cols; ++j)
                    phi(i, j) -= minValue;
        }
    };
}

#endif