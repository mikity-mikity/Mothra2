using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ShoNS.Array;
namespace mikity.ghComponents
{
    public partial class Mothra2 : Grasshopper.Kernel.GH_Component
    {
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
                C[i, Q] = -11;
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
