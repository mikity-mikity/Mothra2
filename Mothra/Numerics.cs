using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ShoNS.Array;
using Rhino.Geometry;
namespace mikity.ghComponents
{
    public partial class Mothra2 : Grasshopper.Kernel.GH_Component
    {
        private DoubleArray[] baseFunction;
        private DoubleArray[] coeff;
        Func<double, double, double>[] Function;
        Action<double, double, double[]>[] dFunction;
        Action<double, double, double[,]>[] ddFunction;
        private static double epsilon = 0.2;
        Func<double, double, double, double, double> F = (x1, x2, y1, y2) => { return Math.Sqrt(1 + epsilon * (((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)))); };
        Action<double, double, double, double, double[]> dF = (x1, x2, y1, y2, grad) =>
        {
            grad[0] = epsilon * (x1 - x2) / Math.Sqrt(1 + epsilon * (((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2))));
            grad[1] = epsilon * (y1 - y2) / Math.Sqrt(1 + epsilon * (((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2))));
        };
        Action<double, double, double, double, double[,]> ddF = (x1, x2, y1, y2, grad) =>
        {
            grad[0, 0] = epsilon / Math.Sqrt(1 + epsilon * (((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)))) - epsilon * epsilon * (x1 - x2) * (x1 - x2) * Math.Pow(1 + epsilon * (((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2))), -3 / 2d);
            grad[0, 1] = -epsilon * epsilon * (x1 - x2) * (y1 - y2) * Math.Pow(1 + epsilon * (((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2))), -3 / 2d);
            grad[1, 0] = -epsilon * epsilon * (x1 - x2) * (y1 - y2) * Math.Pow(1 + epsilon * (((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2))), -3 / 2d);
            grad[1, 1] = epsilon / Math.Sqrt(1 + epsilon * (((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2)))) - epsilon * epsilon * (y1 - y2) * (y1 - y2) * Math.Pow(1 + epsilon * (((x1 - x2) * (x1 - x2)) + ((y1 - y2) * (y1 - y2))), -3 / 2d);
        };
        public void computeBaseFunction(int lastComputed)
        {
            if (lastComputed >= nOutterSegments)
            {
                computeBaseFunction2(lastComputed);
            }
            else
            {
                computeBaseFunction1(lastComputed);
            }
        }

        private void computeBaseFunctionCommon(int lastComputed,DoubleArray origX)
        {
            var M = (shiftArray.T.Multiply(Laplacian) as SparseDoubleArray) * shiftArray as SparseDoubleArray;
            int T1 = n - fixedPoints.Count;
            int T2 = n;
            var slice1 = new SparseDoubleArray(T1, T2);
            var slice2 = new SparseDoubleArray(T2, T2 - T1);
            for (int i = 0; i < T1; i++)
            {
                slice1[i, i] = 1;
            }
            for (int i = 0; i < T2 - T1; i++)
            {
                slice2[i + T1, i] = 1;
            }
            var DIB = (slice1.Multiply(M) as SparseDoubleArray).Multiply(slice2) as SparseDoubleArray;
            var DII = (slice1.Multiply(M) as SparseDoubleArray).Multiply(slice1.T) as SparseDoubleArray;
            var solver = new SparseLU(DII);
            origX = shiftArray.T * origX;
            var fixX = origX.GetSlice(T1, T2 - 1, 0, 0);
            var B = -DIB * fixX;
            var dx = solver.Solve(B);
            var ret = DoubleArray.Zeros(n, 1);
            for (int i = 0; i < T1; i++)
            {
                ret[i, 0] = dx[i, 0];
            }
            for (int i = T1; i < T2; i++)
            {
                ret[i, 0] = fixX[i - T1, 0];
            }
            if (lastComputed < baseFunction.Length)
            {
                baseFunction[lastComputed] = shiftArray * ret;
            }
            //vertices...x,y
            //baseFunction[]...z
            var MM = DoubleArray.Zeros(n, n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    MM[i, j] = F(vertices[i].X, vertices[j].X, vertices[i].Y, vertices[j].Y);
                }
            }
            var V = DoubleArray.Zeros(n, 1);
            for (int i = 0; i < n; i++)
            {
                V[i, 0] = baseFunction[lastComputed][i, 0];
            }
            //var solver2 = new LU(MM);
            var solver2 = new SVD(MM);
            coeff[lastComputed] = solver2.Solve(V);
            Function[lastComputed] = (x, y) =>
            {
                double val = 0;
                for (int j = 0; j < n; j++)
                {
                    val += coeff[lastComputed][j, 0] * F(x, vertices[j].X, y, vertices[j].Y);
                }
                return val;
            };
            dFunction[lastComputed] = (x, y, val) =>
            {
                val[0] = 0;
                val[1] = 0;
                double[] r = new double[2] { 0, 0 };
                for (int j = 0; j < n; j++)
                {
                    dF(x, vertices[j].X, y, vertices[j].Y, r);
                    val[0] += coeff[lastComputed][j, 0] * r[0];
                    val[1] += coeff[lastComputed][j, 0] * r[1];
                }

            };
            ddFunction[lastComputed] = (x, y, val) =>
            {
                val[0, 0] = 0;
                val[0, 1] = 0;
                val[1, 0] = 0;
                val[1, 1] = 0;
                double[,] r = new double[2, 2] { {0, 0}, {0, 0} };
                for (int j = 0; j < n; j++)
                {
                    ddF(x, vertices[j].X, y, vertices[j].Y, r);
                    val[0, 0] += coeff[lastComputed][j, 0] * r[0, 0];
                    val[0, 1] += coeff[lastComputed][j, 0] * r[0, 1];
                    val[1, 0] += coeff[lastComputed][j, 0] * r[1, 0];
                    val[1, 1] += coeff[lastComputed][j, 0] * r[1, 1];
                }

            };
            if (lastComputed == 0) resultToPreview(0);
        }
        private void computeBaseFunction2(int lastComputed)
        {
            DoubleArray origX = new DoubleArray(n, 1);
            for (int i = 0; i < n; i++)
            {
                origX[i, 0] = 0;
            }
            var b = bbIn[lastComputed-nOutterSegments];
            for (int i = 0; i < b.Count; i++)
            {
                origX[b[i].P0, 0] = 5d;
            }
            computeBaseFunctionCommon(lastComputed, origX);

        }
        private void computeBaseFunction1(int lastComputed)
        {
            DoubleArray origX = new DoubleArray(n, 1);
            for (int i = 0; i < n; i++)
            {
                origX[i, 0] = 0;
            }
            var b = bbOut[lastComputed];
            for (int i = 0; i < b.Count; i++)
            {
                double y = Math.Pow((i - ((double)b.Count / 2d)) / ((double)b.Count / 2d), 2) - 1d;
                origX[b[i].P0, 0] = -3d * y;
                if (i == b.Count - 1)
                {
                    y = Math.Pow(((i + 1d) - ((double)b.Count / 2d)) / ((double)b.Count / 2d), 2) - 1d;
                    origX[b[i].P1, 0] = -3d * y;
                }
            }
            computeBaseFunctionCommon(lastComputed, origX);
        }
        public void resultToPreview(int num)
        {
            result = new List<Rhino.Geometry.Line>();
            if (baseFunction[num] == null) { result = null; return; }
            var xx = baseFunction[num];
            foreach (var edge in edges)
            {
                Rhino.Geometry.Point3d P = new Rhino.Geometry.Point3d(vertices[edge.P0].X, vertices[edge.P0].Y, xx[edge.P0]);
                Rhino.Geometry.Point3d Q = new Rhino.Geometry.Point3d(vertices[edge.P1].X, vertices[edge.P1].Y, xx[edge.P1]);
                result.Add(new Rhino.Geometry.Line(P, Q));
            }
            for (int i = 0; i < a.Count; i++)
            {
                Point3d P = new Point3d(a[i].X, a[i].Y, Function[num](a2[i].X, a2[i].Y));
                a.RemoveAt(i);
                a.Insert(i, P);
            }
        }
        public SparseDoubleArray computeLaplacian(int[,] lines, int nP)
        {
            if (lines.GetLength(1) != 2) return null;
            SparseDoubleArray C = new SparseDoubleArray(lines.GetLength(0), nP);
            int m = lines.GetLength(0);
            for (int i = 0; i < m; i++)
            {
                int P = lines[i, 0];
                int Q = lines[i, 1];
                C[i, P] = 1;
                C[i, Q] = -1;
            }
            SparseDoubleArray D = (C.T * C) as SparseDoubleArray;
            return D;
        }
        public SparseDoubleArray computeShiftArray(List<int> fixedPoints,int n)
        {
            int[] shift = new int[n];

            int NN = fixedPoints.Count();
            int DD = n - NN; //number of free nodes
            for (int i = 0; i < n; i++)
            {
                shift[i] = -100;
            }
            for (int i = 0; i < NN; i++)
            {
                shift[fixedPoints[i]] = DD + i;
            }
            int S = 0;
            for (int i = 0; i < n; i++)
            {
                if (shift[i] == -100)
                {
                    shift[i] = S;
                    S++;
                }
            }
            if (S != DD) AddRuntimeMessage(Grasshopper.Kernel.GH_RuntimeMessageLevel.Error, "contradiction!");
            var shiftArray = new SparseDoubleArray(n, n);
            for (int i = 0; i < n; i++)
            {
                shiftArray[i, shift[i]] = 1;
            }
            return shiftArray;
        }
    }
}
