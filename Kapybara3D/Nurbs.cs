using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Minilla3D.Elements
{
    public class nurbsElement:element
    {
        public class tuple
        {
            public int index;  //element index
            public double area;
            int N;

            public double ou, ov, u, v;   //scaled coordinate, coordinate on Rhino, local coordinate on an element
            public double[] lo;
            public double[] f;
            public double[][] df;
            public double[][,] ddf;
            public double[] nf;
            public double[][] ndf;
            public double[][,] nddf;
            public double[,] kernel;
            public double x, y, z;
            public double[][] gi;
            public double[][] Gi;
            public double[,] gij;
            public double[,] Gij;
            public double[][] gi2;
            public double[][] Gi2;
            public double[,] gij2;
            public double[,] Gij2;
            public double[, ,] Gammaijk;
            public double[,][] second;
            public double[, ,] Gammaijk2;
            public double[,][] second2;
            public double[,] shape;
            public double[, ,] C;
            public double[, , ,] B;
            public double[, , ,] D;
            public double dv, refDv;
            public void CtoB(int elemDim,int nDV)
            {
                for (int i = 0; i < elemDim; i++)
                {
                    for (int j = 0; j < elemDim; j++)
                    {
                        for (int u = 0; u < nDV; u++)
                        {
                            for (int v = 0; v < nDV; v++)
                            {
                                B[i, j, u, v] = 0;
                                for (int r = 0; r < __DIM; r++)
                                {
                                    //対称化
                                    B[i, j, u, v] += C[i, r, u] * C[j, r, v];
                                    B[i, j, u, v] += C[j, r, u] * C[i, r, v];
                                }
                            }
                        }
                    }
                }
            }

            public tuple(int _N, double _ou, double _ov, double _u, double _v, int _index, double _loU, double _loV, double _area)
            {

                N = _N;
                ou = _ou;
                ov = _ov;
                u = _u;
                v = _v;
                lo=new double[2]{_loU,_loV};

                index = _index;
                area = _area;
                f = new double[N];
                df = new double[N][];
                ddf = new double[N][,];
                nf = new double[N];
                ndf = new double[N][];
                nddf = new double[N][,];
                kernel = new double[N, N];
                for (int i = 0; i < N; i++)
                {
                    df[i] = new double[2];
                    ndf[i] = new double[2];
                    ddf[i] = new double[2, 2];
                    nddf[i] = new double[2, 2];
                }
                x = 0;
                y = 0;
                z = 0;
                gi = new double[2][] { new double[3], new double[3] };
                Gi = new double[2][] { new double[3], new double[3] };
                gij = new double[2, 2];
                Gij = new double[2, 2];
                gi2 = new double[2][] { new double[3], new double[3] };
                Gi2 = new double[2][] { new double[3], new double[3] };
                gij2 = new double[2, 2];
                Gij2 = new double[2, 2];
                Gammaijk = new double[2, 2, 2];
                second = new double[2, 2][] { { new double[3], new double[3] }, { new double[3], new double[3] } };
                Gammaijk2 = new double[2, 2, 2];
                second2 = new double[2, 2][] { { new double[3], new double[3] }, { new double[3], new double[3] } };
            }
        }
        public void precompute(tuple tup)
        {
            //Assume M[0] and M[1] are precomputed,

            //double[][,] M=new double[2][,];
		    //M[0]=fM(uNum,_uDim,_uDim-1,uKnot);
		    //M[1]=fM(vNum,_vDim,_vDim-1,vKnot);

            tup.shape = new double[__DIM, nDV];                        //Global coordinate *coefficient*
            tup.C = new double[elemDim, __DIM, nDV];                //Base vectors *coefficient*
            tup.B = new double[elemDim, elemDim, nDV, nDV];          //Metric *coefficient*
            tup.D = new double[elemDim, elemDim, __DIM, nDV];   //Hessian coefficient

			//Shape functions  [N] (for global coordinate)
			for(int j=0;j<elemDim;j++)
			{
				double t=tup.lo[j];
				for(int k=0;k<dim[j];k++)
				{
					hh[j][k]=Math.Pow(t,(dim[j]-k-1));
				}
				for(int k=0;k<dim[j];k++)
				{
					double val=0;
					for(int l=0;l<dim[j];l++)
					{
						val+=hh[j][l]*M[j][l,k];
					}
					tt[j][k]=val;
				}
			}
            for (int j = 0; j < __DIM; j++)
            {
                for (int k = 0; k < nDV; k++)
                {
                    tup.shape[j, k] = 0;
                }
            }
            for (int k = 0; k < nNode; k++)
            {
                //Shape functinos
                double shape = 1.0;
                for (int j = 0; j < elemDim; j++)
                {
                    shape *= tt[j][dd[k, j]];
                }
                for (int j = 0; j < __DIM; j++)
                {
                    tup.shape[j, k * __DIM + j] = shape;
                }
            }
            //Create [C]  (for base vectors)
            for (int m = 0; m < elemDim; m++)
            {
                for (int j = 0; j < elemDim; j++)
                {
                    double t = tup.lo[j];
                    if (j != m)
                    {
                        for (int k = 0; k < dim[j]; k++)
                        {
                            hh[j][k] = Math.Pow(t, (dim[j] - k - 1));
                        }
                    }
                    else
                    {
                        for (int k = 0; k < dim[j] - 1; k++)
                        {
                            hh[j][k] = (dim[j] - k - 1) * Math.Pow(t, (dim[j] - k - 2));
                        }
                        hh[j][dim[j] - 1] = 0;
                    }
                    for (int k = 0; k < dim[j]; k++)
                    {
                        double val = 0;
                        for (int l = 0; l < dim[j]; l++)
                        {
                            val += hh[j][l] * M[j][l, k];
                        }
                        tt[j][k] = val;
                    }
                }
                for (int jj = 0; jj < __DIM; jj++)
                {
                    for (int j = 0; j < nDV; j++)
                    {
                        tup.C[m, jj, j] = 0;
                    }
                }
                for (int k = 0; k < nNode; k++)
                {
                    //[C]
                    double C = 1.0;
                    for (int j = 0; j < elemDim; j++)
                    {
                        C *= tt[j][dd[k, j]];
                    }
                    for (int j = 0; j < __DIM; j++)
                    {
                        tup.C[m, j, k * __DIM + j] = C;
                    }
                }
            }
            //Create [B]  (for metric)
            tup.CtoB(elemDim, nDV);
            
            //Create [D] (for second derivative)
            for (int m = 0; m < elemDim; m++)
            {
                for (int n = 0; n < elemDim; n++)
                {
                    for (int j = 0; j < elemDim; j++)
                    {
                        double t = tup.lo[j];
                        if (j != m && j != n)
                        {
                            for (int k = 0; k < dim[j]; k++)
                            {
                                hh[j][k] = Math.Pow(t, (dim[j] - k - 1));
                            }
                        }
                        if ((j != m && j == n) || (j == m && j != n))
                        {
                            for (int k = 0; k < dim[j] - 1; k++)
                            {
                                hh[j][k] = (dim[j] - k - 1) * Math.Pow(t, (dim[j] - k - 2));
                            }
                            hh[j][dim[j] - 1] = 0;
                        }
                        if (j == m && j == n)
                        {
                            for (int k = 0; k < dim[j] - 1; k++)
                            {
                                hh[j][k] = (dim[j] - k - 1) * (dim[j] - k - 2) * Math.Pow(t, (dim[j] - k - 3));
                            }
                            hh[j][dim[j] - 1] = 0;
                            hh[j][dim[j] - 2] = 0;
                        }

                        for (int k = 0; k < dim[j]; k++)
                        {
                            double val = 0;
                            for (int l = 0; l < dim[j]; l++)
                            {
                                val += hh[j][l] * M[j][l, k];
                            }
                            tt[j][k] = val;
                        }
                    }
                    for (int jj = 0; jj < __DIM; jj++)
                    {
                        for (int j = 0; j < nDV; j++)
                        {
                            tup.D[m, n, jj, j] = 0;
                        }
                    }
                    for (int k = 0; k < nNode; k++)
                    {
                        //[D]
                        double D = 1.0;
                        for (int j = 0; j < elemDim; j++)
                        {
                            D *= tt[j][dd[k, j]];
                        }
                        for (int j = 0; j < __DIM; j++)
                        {
                            tup.D[m, n, j, k * __DIM + j] = D;
                        }
                    }
                }
            }
            //covariant base vectors
            for (int n = 0; n < elemDim; n++)
            {
                double fx = 0, fy = 0;
                for (int i = 0; i < nDV; i++)
                {
                    fx += tup.C[n, 0, i] * node[i];
                    fy += tup.C[n, 1, i] * node[i];
                }
                tup.gi[n][0] = fx;
                tup.gi[n][1] = fy;
                tup.gi[n][2] = 0;
            }
            for (int n = 0; n < elemDim; n++)
            {
                for (int m = 0; m < elemDim; m++)
                {
                    tup.gij[n, m] = tup.gi[n][0] * tup.gi[m][0] + tup.gi[n][1] * tup.gi[m][1]+ tup.gi[n][2] * tup.gi[m][2];
                }
            }
            if (elemDim == 1)
            {
                _inv1(tup.gij, tup.Gij);
            }
            else if (elemDim == 2)
            {
                _inv2(tup.gij, tup.Gij);
            }
            else if (elemDim == 3)
            {
                _inv3(tup.gij, tup.Gij);
            }
            if (elemDim == 1)
            {
                tup.dv = Math.Sqrt(_det1(tup.gij));
            }
            else if (elemDim == 2)
            {
                tup.dv = Math.Sqrt(_det2(tup.gij));
            }
            else if (elemDim == 3)
            {
                tup.dv = Math.Sqrt(_det3(tup.gij));
            }
            tup.refDv = tup.dv;

            //contravatiant base vectors
            for (int n = 0; n < elemDim; n++)
            {
                double Fx = 0, Fy = 0;
                for (int m = 0; m < elemDim; m++)
                {
                    Fx += tup.gi[m][0] * tup.Gij[m, n];
                    Fy += tup.gi[m][1] * tup.Gij[m, n];
                }
                tup.Gi[n][0] = Fx;
                tup.Gi[n][1] = Fy;
            }
            //Connection coefficients
            for (int n = 0; n < elemDim; n++)
            {
                for (int m = 0; m < elemDim; m++)
                {
                    double gx = 0, gy = 0;
                    for (int i = 0; i < nDV; i++)
                    {
                        gx += tup.D[n, m, 0, i] * node[i];
                        gy += tup.D[n, m, 1, i] * node[i];
                    }
                    tup.second[n, m][0] = gx;
                    tup.second[n, m][1] = gy;
                    tup.second[n, m][2] = 0;
                    for (int k = 0; k < elemDim; k++)
                    {
                        tup.Gammaijk[n, m, k] = gx * tup.Gi[k][0] + gy * tup.Gi[k][1];
                    }
                }
            }
            /*
            //Create gradient of hessian with computed connection coefficients
            for (int m = 0; m < elemDim; m++)
            {
                for (int n = 0; n < elemDim; n++)
                {
                    for (int k = 0; k < nNode; k++)
                    {
                        gradient[m, n][k] = D[m, n, 2, k * __DIM + 2];
                        for (int i = 0; i < elemDim; i++)
                        {
                            gradient[m, n][k] -= Gamma[m, n, i] * C[i, 2, k * 3 + 2];
                        }
                    }
                }
            }
            //Raising index
            double[,] tmp = new double[elemDim, elemDim];
            for (int k = 0; k < nNode; k++)
            {
                for (int m = 0; m < elemDim; m++)
                {
                    for (int n = 0; n < elemDim; n++)
                    {
                        double val = 0;
                        for (int s = 0; s < elemDim; s++)
                        {
                            val += gradient[m, s][k] * invMetric[s, n];
                        }
                        gradient2[m, n][k] = val;
                    }
                }
            }*/
        }
/*
 * 
 * 
                        }
 * */
        double[][,] M;
        int[] dim;
		double[][] hh;
		double[][] tt;
        int[,] dd;
            
        double[][] _cu;
		double[][] _pu;
        public int uDim,vDim;
        private void _inv1(double[,] from, double[,] to)
        {
            if (from[0, 0] <= 0)
            {
                to[0, 0] = 0;
                from[0, 0] = 0;
            }
            else
            {
                to[0, 0] = 1 / from[0, 0];
            }
        }
        private void _inv2(double[,] from, double[,] to)
        {
            double det = _det2(from);
            if (det <= 0)
            {
                to[0, 0] = 0;
                to[1, 0] = 0;
                to[0, 1] = 0;
                to[1, 1] = 0;
                from[0, 0] = 0;
                from[1, 0] = 0;
                from[0, 1] = 0;
                from[1, 1] = 0;
            }
            else
            {
                to[0, 0] = from[1, 1] / det;
                to[1, 1] = from[0, 0] / det;
                to[1, 0] = -from[1, 0] / det;
                to[0, 1] = -from[0, 1] / det;
            }
        }
        private void _inv3(double[,] from, double[,] to)
        {
            double det = _det3(from);
            if (det <= 0)
            {
                to[0, 0] = 0;
                to[0, 1] = 0;
                to[0, 2] = 0;
                to[1, 0] = 0;
                to[1, 1] = 0;
                to[1, 2] = 0;
                to[1, 3] = 0;
                to[2, 0] = 0;
                to[2, 1] = 0;
                to[2, 2] = 0;
                from[0, 0] = 0;
                from[0, 1] = 0;
                from[0, 2] = 0;
                from[1, 0] = 0;
                from[1, 1] = 0;
                from[1, 2] = 0;
                from[1, 3] = 0;
                from[2, 0] = 0;
                from[2, 1] = 0;
                from[2, 2] = 0;
            }
            else
            {
                to[0, 0] = (from[1, 1] * from[2, 2] - from[1, 2] * from[2, 1]) / det;
                to[1, 0] = (from[2, 1] * from[0, 2] - from[2, 2] * from[0, 1]) / det;
                to[2, 0] = (from[0, 1] * from[1, 2] - from[0, 2] * from[1, 1]) / det;
                to[0, 1] = (from[1, 2] * from[2, 0] - from[1, 0] * from[2, 2]) / det;
                to[1, 1] = (from[2, 2] * from[0, 0] - from[2, 0] * from[0, 2]) / det;
                to[2, 1] = (from[0, 2] * from[1, 0] - from[0, 0] * from[1, 2]) / det;
                to[0, 2] = (from[1, 0] * from[2, 1] - from[1, 1] * from[2, 0]) / det;
                to[1, 2] = (from[2, 0] * from[0, 1] - from[2, 1] * from[0, 0]) / det;
                to[2, 2] = (from[0, 0] * from[1, 1] - from[0, 1] * from[1, 0]) / det;
            }
        }
        private double _det1(double[,] m)
        {
            return m[0, 0];
        }
        private double _det2(double[,] m)
        {
            return m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0];
        }
        private double _det3(double[,] m)
        {
            return m[0, 0] * m[1, 1] * m[2, 2]
            + m[1, 0] * m[2, 1] * m[0, 2]
            + m[2, 0] * m[0, 1] * m[1, 2]
            - m[2, 0] * m[1, 1] * m[0, 2]
            - m[1, 0] * m[0, 1] * m[2, 2]
            - m[0, 0] * m[2, 1] * m[1, 2];
        }

        double[] fN(int _i, int _k, int _dim, int dim, double[] knot)
		{
		    if (_dim==1)
			{
		        double[] F=new double[dim]; 
				for (int i=0;i<dim;i++)
				{
					F[i]=0;
				}
		        if (_k==_i)
				{
					F[dim-1]=1;
				}
		        return F;
			}
		    double[] S1=fN(_i,_k,_dim-1,dim,knot);
			double[] S2=fN(_i,_k+1,_dim-1,dim,knot);
			double E1=knot[_k+_dim-2]-knot[_k-1];
			double E2=knot[_k+_dim-1]-knot[_k];
			double[] D1=new double[2]{0,0};
			double[] D2=new double[2]{0,0};
		    if (E1>0)
			{
				D1[0]=1d/E1;
				D1[1]=-knot[_k-1]/E1;
			}
			if (E2>0)
			{
				D2[0]=-1d/E2;
				D2[1]=knot[_k+_dim-1]/E2;
			}
		    double[] F2=new double[dim]; 
			for (int i=0;i<dim;i++)
			{
				F2[i]=0;
			}
			for(int i=1;i<dim;i++)
			{
				F2[i-1]=F2[i-1]+S1[i]*D1[0];
				F2[i]=F2[i]+S1[i]*D1[1];
				F2[i-1]=F2[i-1]+S2[i]*D2[0];
				F2[i]=F2[i]+S2[i]*D2[1];
			}
			return F2;
		}
        double[,] fM(int shift, int dim, int ddim, double[] knot){
			double[,] M=new double[dim,dim];
            for (int i = 0; i < dim; i++)
            {
                for (int j = 0; j < dim; j++)
                {
                    M[i,j] = 0;
                }
            }
		    for(int k=shift;k<dim+shift;k++)
			{
		        double[] D=fN(shift+ddim,k,dim,dim,knot);
				for (int n =0;n<dim;n++)
				{
					M[n,k-shift]=D[n];
				}
			}

			double[,] S=new double[dim,dim];
            for(int i=0;i<dim;i++)
			{
                for(int j=0;j<dim;j++)
			    {
				    S[i,j]=0;
			    }
            }
			for(int  n =1;n<dim+1;n++)
			{
				for (int k=1+n;k<dim+2;k++)
				{
					if (n==dim)
					{
						for (int t =0;t<n-1;t++)
						{
						   S[t,n-1]=0;
						}
						S[(n-1),n-1]=1;
					}else
					{
						S[(k-2),n-1]=binominal(dim-n,dim+1-k)*Math.Pow(shift-1,k-1-n);
					}
				}
			}
			double[,] G=new double[dim,dim];
			for (int j=0;j<dim;j++)
			{
				for(int k=0;k<dim;k++)
				{
					double v=0;
					for(int l=0;l<dim;l++)
					{
						v+=S[j,l]*M[l,k];
					}
					G[j,k]=v;
				}
			}
			return G;
		}
        private static int Factorial(int x)
        {
            if (x == 0) return 1;
            if (x == 1) return 1;
            if (x == 2) return 2;
            if (x == 3) return 6;
            if (x == 4) return 24;
            if (x == 5) return 120;
            int val = 1;
            for (int i = 2; i <= x; i++)
            {
                val *= i;
            }
            return val;
        }

        public static double binominal(int N, int k)
        {
            return Factorial(N) / Factorial(N - k) / Factorial(k);
        }
        bool exceptional()
        {
            if (this.typeOfBorder == (border.Left | border.Top))
            {
                return true;
            }
            if (this.typeOfBorder == (border.Left | border.Bottom))
            {
                return true;
            }
            if (this.typeOfBorder == (border.Right | border.Top))
            {
                return true;
            }
            if (this.typeOfBorder == (border.Right | border.Bottom))
            {
                return true;
            }
            return false;
        }
        public int numberOfConstraintConditions()
        {
            if (exceptional()) return 0;
            return nBIntPoint * 2;
        }
        public int numberOfConstraintConditions2(bool T)
        {
            if (T&&exceptional()) return 0;
            return nIntPoint * 3;
        }
        public int mergeResidual(ShoNS.Array.DoubleArray residual, int i)
        {
            if (exceptional()) return i;
            for (int j = 0; j < nBIntPoint; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    double resid = bIntP[j].getResidualOfBoundaryCondition(node,k);
                    residual[i + j * 2 + k] = resid;
                }
            }
            return i + nBIntPoint * 2;
        }
        public int mergeJacobian(ShoNS.Array.SparseDoubleArray jacobian, int i)
        {
            if (exceptional()) return i;
            for (int j = 0; j < nBIntPoint; j++)
            {
                for (int k = 0; k < 2; k++)
                {
                    var grad = bIntP[j].getGradientOfBoundaryCondition(k);
                    for(int f=0;f<nNode;f++)
                    {
                        jacobian[i + j*2 + k, index[f]] = grad[f];
                    }
                }
            }
            return i + nBIntPoint * 2;
        }
        public int mergeJacobianOfCurvature(ShoNS.Array.SparseDoubleArray jacobH, int num,bool T)
        {
            if (T&&exceptional()) return num;
            for (int i = 0; i < nIntPoint; i++)
            {
                var grad = intP[i].getGradientOfH(0, 0);
                for (int f = 0; f < nNode; f++)
                {
                    jacobH[num + i * 3 + 0, index[f]] = grad[f];
                }
                grad = intP[i].getGradientOfH(1, 1);
                for (int f = 0; f < nNode; f++)
                {
                    jacobH[num + i * 3 + 1, index[f]] = grad[f];
                }
                grad = intP[i].getGradientOfH(0, 1);
                for (int f = 0; f < nNode; f++)
                {
                    jacobH[num + i * 3 + 2, index[f]] = grad[f];
                }
            }
            return num + nIntPoint*3;
        }

        double __cu(int n, int dim)
        {
            switch (dim)
            {
                case 6:
                    switch (n)
                    {
                        case 0:
                            return (-0.9491079123427585245261897) * 0.5 + 0.5;
                        case 1:
                            return (-0.7415311855993944398638648) * 0.5 + 0.5;
                        case 2:
                            return (-0.4058451513773971669066064) * 0.5 + 0.5;
                        case 3:
                            return 0.5;
                        case 4:
                            return (0.4058451513773971669066064) * 0.5 + 0.5;
                        case 5:
                            return (0.7415311855993944398638648) * 0.5 + 0.5;
                        case 6:
                            return (0.9491079123427585245261897) * 0.5 + 0.5;
                        default:
                            return 0;
                    }
                case 5:
                    switch (n)
                    {
                        case 0:
                            return (-0.9324695142031521) * 0.5 + 0.5;
                        case 1:
                            return (-0.6612093864662645) * 0.5 + 0.5;
                        case 2:
                            return (-0.2386191860831969) * 0.5 + 0.5;
                        case 3:
                            return (0.2386191860831969) * 0.5 + 0.5;
                        case 4:
                            return (0.6612093864662645) * 0.5 + 0.5;
                        case 5:
                            return (0.9324695142031521) * 0.5 + 0.5;
                        default:
                            return 0;
                    }
                case 4:
                    switch (n)
                    {
                        case 0:
                            return (-0.9061798459386640)*0.5+0.5;
                        case 1:
                            return (-0.5384693101056831)*0.5+0.5;
                        case 2:
                            return (0.0000000000000000)*0.5+0.5;
                        case 3:
                            return (0.5384693101056831)*0.5+0.5;
                        case 4:
                            return (0.9061798459386640)*0.5+0.5;
				        default:
                            return 0;
                    }
                case 3:
                    switch (n)
                    {
                        case 0:
                            return (-0.8611363115940526)*0.5+0.5;
                        case 1:
                            return (-0.3399810435848563)*0.5+0.5;
                        case 2:
                            return (0.3399810435848563)*0.5+0.5;
                        case 3:
                            return (0.8611363115940526)*0.5+0.5;
                        default:
                            return 0;
                    }
                case 2:
                    switch (n)
                    {
                        case 0:
                            return (-0.7745966692414834) * 0.5 + 0.5;
                        case 1:
                            return (0.0000000000000000) * 0.5 + 0.5;
                        case 2:
                            return (0.7745966692414834) * 0.5 + 0.5;
                        default:
                            return 0;
                    }
                default:
                    return 0;

            }
        }
        double __pu(int n, int dim)
        {
            switch (dim)
            {
                case 6:
                    switch (n)
                    {
                        case 0:
                            return 0.1294849661688696932706114 * 0.5;
                        case 1:
                            return 0.2797053914892766679014678 * 0.5;;
                        case 2:
                            return 0.3818300505051189449503698 * 0.5;;
                        case 3:
                            return 0.4179591836734693877551020 * 0.5;;
                        case 4:
                            return 0.3818300505051189449503698 * 0.5;;
                        case 5:
                            return 0.2797053914892766679014678 * 0.5;;
                        case 6:
                            return 0.1294849661688696932706114 * 0.5;;
                        default:
                            return 0;
                    }
                case 5:
                    switch (n)
                    {
                        case 0:
                            return 0.1713244923791704 * 0.5;
                        case 1:
                            return 0.3607615730481386 * 0.5;
                        case 2:
                            return 0.4679139345726910 * 0.5;
                        case 3:
                            return 0.4679139345726910 * 0.5;
                        case 4:
                            return 0.3607615730481386 * 0.5;
                        case 5:
                            return 0.1713244923791704 * 0.5;
                        default:
                            return 0;
                    }
                case 4:
                    switch (n)
                    {
                        case 0:
                            return 0.2369268850561891 * 0.5;
                        case 1:
                            return 0.4786286704993665 * 0.5;
                        case 2:
                            return 0.5688888888888889 * 0.5;
                        case 3:
                            return 0.4786286704993665 * 0.5;
                        case 4:
                            return 0.2369268850561891 * 0.5;
                        default:
                            return 0;
                    }
                case 3:
                    switch (n)
                    {
                        case 0:
                            return 0.3478548451374538*0.5;
                        case 1:
                            return 0.6521451548625461*0.5;
                        case 2:
                            return 0.6521451548625461*0.5;
                        case 3:
                            return 0.3478548451374538 * 0.5;
                        default:
                            return 0;
                    }
                case 2:
                    switch (n)
                    {
                        case 0:
                            return 0.5555555555555556*0.5;
                        case 1:
                            return 0.8888888888888888*0.5;
                        case 2:
                            return 0.5555555555555556 * 0.5;
                        default:
                            return 0;
                    }
                default:
                    return 0;
            }
        }
        public void stitch(params nurbsCurve[] elems)
        {
            int count = 0;
            foreach (var e in elems)
            {
                for (int i = 0; i < e.nIntPoint; i++, count++)
                {
                    bIntP[count].refIntP = e.intP[i];
                    e.intP[i].refIntP = bIntP[count];
                }
            }
        }
        public void setPlane(double a, double b, double c, double d)
        {
            foreach (var p in intP)
            {
                p.refIntP.setPlane(a, b, c, d);
            }
        }
        public void computeAngle()
        {
            foreach (var p in bIntP)
            {
                p.computeAngle(this.node);
                p.transferAngleToTension();
            }
        }
        public nurbsElement(int _uDim, int _vDim, int[] _index, int uNum, int vNum, double[] uKnot, double[] vKnot, border _border = border.None)
            : base(_index, _uDim * _vDim, 2, 4)
        {

            this.typeOfBorder = _border;
            int bVdim = _vDim +1;
            int bUdim = _uDim +1;
            switch (_border)
            {
                case border.Left | border.Top|border.Right:
                    this.bIntP = new integratingPoint[bVdim * 2 + bUdim];
                    nBIntPoint = bVdim * 2 + bUdim;
                    break;
                case border.Left | border.Bottom | border.Right:
                    this.bIntP = new integratingPoint[bVdim * 2 + bUdim];
                    nBIntPoint = bVdim * 2 + bUdim;
                    break;
                case border.Left | border.Top | border.Bottom:
                    this.bIntP = new integratingPoint[bUdim * 2 + bVdim];
                    nBIntPoint = bUdim * 2 + bVdim;
                    break;
                case border.Right | border.Top | border.Bottom:
                    this.bIntP = new integratingPoint[bUdim * 2 + bVdim];
                    nBIntPoint = bUdim * 2 + bVdim;
                    break;
                case border.Left | border.Right:
                    this.bIntP = new integratingPoint[bVdim + bVdim];
                    nBIntPoint = bVdim + bVdim;
                    break;
                case border.Top | border.Bottom:
                    this.bIntP = new integratingPoint[bUdim + bUdim];
                    nBIntPoint = bUdim + bUdim;
                    break;
                case border.Left | border.Top:
                    this.bIntP = new integratingPoint[bUdim + bVdim];
                    nBIntPoint = bUdim + bVdim;
                    break;
                case border.Left | border.Bottom:
                    this.bIntP = new integratingPoint[bUdim + bVdim];
                    nBIntPoint = bUdim + bVdim;
                    break;
                case border.Right | border.Top:
                    this.bIntP = new integratingPoint[bUdim + bVdim];
                    nBIntPoint = bUdim + bVdim;
                    break;
                case border.Right | border.Bottom:
                    this.bIntP = new integratingPoint[bUdim + bVdim];
                    nBIntPoint = bUdim + bVdim;
                    break;
                case border.Left:
                    this.bIntP = new integratingPoint[bVdim];
                    nBIntPoint = bVdim;
                    break;
                case border.Right:
                    this.bIntP = new integratingPoint[bVdim];
                    nBIntPoint = bVdim;
                    break;
                case border.Top:
                    this.bIntP = new integratingPoint[bUdim];
                    nBIntPoint = bUdim;
                    break;
                case border.Bottom:
                    this.bIntP = new integratingPoint[bUdim];
                    nBIntPoint = bUdim;
                    break;
                default:
                    this.bIntP = new integratingPoint[0];
                    nBIntPoint = 0;
                    break;
            }
            for (int i = 0; i < nBIntPoint; i++)
            {
                bIntP[i] = new integratingPoint(nNode, elemDim);
            }
            uDim = _uDim;
            vDim=_vDim;
            _cu = new double[elemDim][];
            _cu[0] = new double[2/*uDim + 1*/];
            _cu[1] = new double[2/*vDim + 1*/];
            _pu = new double[elemDim][];
            _pu[0] = new double[2/*uDim + 1*/];
            _pu[1] = new double[2/*vDim + 1*/];
            int[,] ss = new int[nIntPoint, elemDim];		//Indeces for integrating points
			dd=new int[nNode,elemDim];		    //Indeces for nodes
            dim=new int[2]{uDim,vDim};
			hh=new double[elemDim][];
		    tt=new double[elemDim][];

            //For polynominal
			for(int i=0;i<elemDim;i++)
			{
				hh[i]=new double[dim[i]];
				tt[i]=new double[dim[i]];
			}
					
			//Weight coefficient distribution
			//Coordinates distribution
            /*for (int i = 0; i < uDim+1; i++)
            {
                _cu[0][i] = __cu(i, uDim);
                _pu[0][i] = __pu(i, uDim);
            }*/
            /*for (int i = 0; i < vDim+1; i++)
            {
                _cu[1][i] = __cu(i, vDim);
                _pu[1][i] = __pu(i, vDim);
            }*/
            _cu[0][0] = 0.3;
            _cu[0][1] = 0.7;
            _cu[1][0] = 0.3;
            _cu[1][1] = 0.7;
            _pu[0][0] = 0.5;
            _pu[0][1] = 0.5;
            _pu[1][0] = 0.5;
            _pu[1][1] = 0.5;
            //Indeces for integrating points
			for(int i=0;i<elemDim;i++)
			{
				ss[0,i]=0;
			}
			for(int i=1;i<nIntPoint;i++)
			{
				for(int j=0;j<elemDim;j++)
				{
					ss[i,j]=ss[i-1,j];
				}
				for(int j=0;j<elemDim;j++)
				{

					if(ss[i,j]<1/*dim[j]*/)
					{
						ss[i,j]++;
						for(int k=0;k<j;k++)
						{
							ss[i,k]=0;
						}
						break;
					}
				}
			}

			//Indices for nodes
            for(int i=0;i<elemDim;i++)
			{
				dd[0,i]=0;
			}
			for(int i=1;i<nNode;i++)
			{
				for(int j=0;j<elemDim;j++)
				{
					dd[i,j]=dd[i-1,j];
				}
				for(int j=0;j<elemDim;j++)
				{
					if(dd[i,j]<dim[j]-1)
					{
						dd[i,j]++;
						for(int k=0;k<j;k++)
						{
							dd[i,k]=0;
						}
						break;
					}
				}
			}
			//weight coefficients and loacl coordinates for integrating points
			for(int i=0;i<nIntPoint;i++)
			{
				intP[i].weight=1.0;
				for(int j=0;j<elemDim;j++)
				{
					intP[i].localCoord[j]=_cu[j][ss[i,j]];
					intP[i].weight*=_pu[j][ss[i,j]];
				}
			}
            //local coordinates for integrating points on border
            int count = 0;
            switch (_border)
            {
                case border.Left | border.Top | border.Right:
                    count = 0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bUdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    for (int i = 0; i < bVdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim - 1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    break;
                case border.Left | border.Bottom | border.Right:
                    count = 0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i,bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    break;
                case border.Left | border.Top | border.Bottom:
                    count = 0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    for (int i = 0; i <bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                case border.Right | border.Top | border.Bottom:
                    count = 0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim - 1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim - 1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                case border.Left | border.Right:
                    count = 0;
                    for (int i = 0; i < bVdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bVdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    break;
                case border.Top | border.Bottom:
                    count = 0;
                    for (int i = 0; i < bUdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    for (int i = 0; i < bUdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                case border.Left | border.Top:
                    count = 0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    break;
                case border.Left | border.Bottom:
                    count=0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                case border.Right | border.Top:
                    count=0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    break;
                case border.Right | border.Bottom:
                    count=0;
                    for (int i = 0; i < bVdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    for (int i = 0; i < bUdim; i++,count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                case border.Left:
                    count=0;
                    for (int i = 0; i < bVdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = 0;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, -1 };
                    }
                    break;
                case border.Right:
                    count=0;
                    for (int i = 0; i < bVdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = 1;
                        bIntP[count].localCoord[1] = __cu(i, bVdim-1);
                        bIntP[count].weight = __pu(i, bVdim - 1);
                        bIntP[count].edge = new double[2] { 0, 1 };
                    }
                    break;
                case border.Top:
                    count=0;
                    for (int i = 0; i < bUdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 0;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { 1, 0 };
                    }
                    break;
                case border.Bottom:
                    count=0;
                    for (int i = 0; i < bUdim; i++, count++)
                    {
                        bIntP[count].localCoord[0] = __cu(i, bUdim-1);
                        bIntP[count].localCoord[1] = 1;
                        bIntP[count].weight = __pu(i, bUdim - 1);
                        bIntP[count].edge = new double[2] { -1, 0 };
                    }
                    break;
                default:
                    break;
            }
            List<Minilla3D.Elements.integratingPoint> allIntP = new List<integratingPoint>();
            allIntP.AddRange(intP);
            allIntP.AddRange(bIntP);
            int nAllIntP = nIntPoint + nBIntPoint;
            M=new double[2][,];
		    M[0]=fM(uNum,_uDim,_uDim-1,uKnot);
		    M[1]=fM(vNum,_vDim,_vDim-1,vKnot);
			//Shape functions  [N] (for global coordinate)
			for(int i=0;i<nAllIntP;i++)
			{
				for(int j=0;j<elemDim;j++)
				{
					double t=allIntP[i].localCoord[j];
					for(int k=0;k<dim[j];k++)
					{
						hh[j][k]=Math.Pow(t,(dim[j]-k-1));
					}
					for(int k=0;k<dim[j];k++)
					{
						double val=0;
						for(int l=0;l<dim[j];l++)
						{
							val+=hh[j][l]*M[j][l,k];
						}
						tt[j][k]=val;
					}
				}
				for(int j=0;j<__DIM;j++)
				{
                    for(int k=0;k<nDV;k++)
                    {
					    allIntP[i].N[j,k]=0;
                    }
				}
				for(int k=0;k<nNode;k++)
				{
					//Shape functinos
					double N=1.0;
					for(int j=0;j<elemDim;j++)
					{
						N*=tt[j][dd[k,j]];
					}
					for(int j=0;j<__DIM;j++)
					{
						allIntP[i].N[j,k*__DIM+j]=N;
					}
				}

                //Create [C]  (for base vectors)
				for (int m=0;m<elemDim;m++)
				{
					for(int j=0;j<elemDim;j++)
					{
						double t=allIntP[i].localCoord[j];
						if(j!=m)
						{
							for(int k=0;k<dim[j];k++)
							{
								hh[j][k]=Math.Pow(t,(dim[j]-k-1));
							}
						}else
						{
							for(int k=0;k<dim[j]-1;k++)
							{
								hh[j][k]=(dim[j]-k-1)*Math.Pow(t,(dim[j]-k-2));
							}
							hh[j][dim[j]-1]=0;
						}
						for(int k=0;k<dim[j];k++)
						{
							double val=0;
							for(int l=0;l<dim[j];l++)
							{
								val+=hh[j][l]*M[j][l,k];
							}
							tt[j][k]=val;
						}
					}
                    for(int jj=0;jj<__DIM;jj++)
                    {
					    for(int j=0;j<nDV;j++)
					    {
						    allIntP[i].C[m,jj,j]=0;
					    }
                    }
					for(int k=0;k<nNode;k++)
					{
						//[C]
						double C=1.0;
						for(int j=0;j<elemDim;j++)
						{
							C*=tt[j][dd[k,j]];
						}
						for(int j=0;j<__DIM;j++)
						{
							allIntP[i].C[m,j,k*__DIM+j]=C;
						}
					}
                }
				//Create [B]  (for metric)
                allIntP[i].CtoB();

                //Create [D] (for second derivative)
                for (int m = 0; m < elemDim; m++)
                {
                    for (int n = 0; n < elemDim; n++)
                    {
                        for (int j = 0; j < elemDim; j++)
                        {
                            double t = allIntP[i].localCoord[j];
                            if (j != m&&j!=n)
                            {
                                for (int k = 0; k < dim[j]; k++)
                                {
                                    hh[j][k] = Math.Pow(t, (dim[j] - k - 1));
                                }
                            }
                            if((j!=m&&j==n)||(j==m&&j!=n))
                            {
                                for (int k = 0; k < dim[j] - 1; k++)
                                {
                                    hh[j][k] = (dim[j] - k - 1) * Math.Pow(t, (dim[j] - k - 2));
                                }
                                hh[j][dim[j] - 1] = 0;
                            }
                            if (j == m && j == n)
                            {
                                for (int k = 0; k < dim[j] - 1; k++)
                                {
                                    hh[j][k] = (dim[j] - k - 1) *(dim[j] - k - 2) * Math.Pow(t, (dim[j] - k - 3));
                                }
                                hh[j][dim[j] - 1] = 0;
                                hh[j][dim[j] - 2] = 0;
                            }

                            for (int k = 0; k < dim[j]; k++)
                            {
                                double val = 0;
                                for (int l = 0; l < dim[j]; l++)
                                {
                                    val += hh[j][l] * M[j][l, k];
                                }
                                tt[j][k] = val;
                            }
                        }
                        for (int jj = 0; jj < __DIM; jj++)
                        {
                            for (int j = 0; j < nDV; j++)
                            {
                                allIntP[i].D[m, n,jj, j] = 0;
                            }
                        }
                        for (int k = 0; k < nNode; k++)
                        {
                            //[D]
                            double D = 1.0;
                            for (int j = 0; j < elemDim; j++)
                            {
                                D *= tt[j][dd[k, j]];
                            }
                            for (int j = 0; j < __DIM; j++)
                            {
                                allIntP[i].D[m,n, j, k * __DIM + j] = D;
                            }
                        }
                    }
                }
            }
        }

    }
}
