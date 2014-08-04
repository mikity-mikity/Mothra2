using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Rhino.Geometry;

using Mesher;
using Mesher.Geometry;
using Mesher.Tools;

using Mothra.UI;

using ShoNS.Array;

namespace mikity.ghComponents
{
    public partial class Mothra2 : Grasshopper.Kernel.GH_Component
    {
        public struct tuple
        {
            public double u, v;
            public double area;
            int N;
            public double[] f;
            public double[][] df;
            public double[][,] ddf;
            public double[] nf;
            public double[][] ndf;
            public double[][,] nddf;
            public double[,] kernel;
            public tuple(int _N,double _u,double _v,double _area)
            {
                N = _N;
                u = _u;
                v = _v;
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
            }
        }
        ControlBox myControlBox = new ControlBox();

        List<Point3d> a;
        List<Point3d> a2;
        List<Point3d> b;
        List<Curve> c;
        List<Point3d> d;
        List<Point3d> d2;
        List<Line> f;
        List<Point3d> g;
        List<Line>[] boundaries;
        List<Line>[] holes;
        List<Line> result;
        int nOutterSegments = 0;
        int nInnerLoops = 0;
        Rhino.Geometry.Mesh gmesh = new Rhino.Geometry.Mesh();
        int lastComputed = -1;

        SparseDoubleArray Laplacian;
        SparseDoubleArray shiftArray;
        List<Mesher.Geometry.Edge>[] bbOut;
        List<Mesher.Geometry.Edge>[] bbIn;
        int n, m, r;  //Number of vertices, edges and triangles.

        Mesher.Data.Vertex[] vertices;
        Mesher.Geometry.Edge[] edges;
        Mesher.Data.Triangle[] triangles;
        tuple[] tuples;
        List<int> fixedPoints;
            
        private void init()
        {
            a = new List<Point3d>();
            a2 = new List<Point3d>();
            b = new List<Point3d>();
            c = new List<Curve>();
            d = new List<Point3d>();
            d2 = new List<Point3d>();
            f = new List<Line>();
            g = new List<Point3d>();
            lastComputed = -1;
        }
        public Mothra2()
            : base("Mothra2", "Mothra2", "Mothra2", "Kapybara3D", "Computation")
        {
        }
        public override Guid ComponentGuid
        {
            get { return new Guid("191d64c1-2ed5-4bfc-96e4-800fe372ad0a"); }
        }
        protected override void RegisterInputParams(Grasshopper.Kernel.GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddBrepParameter("TrimmedSurface", "trmSrf", "a trimmed surface, outer trim and holes are supported", Grasshopper.Kernel.GH_ParamAccess.item);
        }
        protected override void RegisterOutputParams(Grasshopper.Kernel.GH_Component.GH_OutputParamManager pManager)
        {
        }
        public override void AddedToDocument(Grasshopper.Kernel.GH_Document document)
        {
            base.AddedToDocument(document);
            myControlBox.Show();
            myControlBox.setFunctionToCompute2(() => { computeF(); });
        }
        void computeF()
        {
            if (lastComputed == nInnerLoops + nOutterSegments - 1) {
                int N=nInnerLoops+nOutterSegments;
                tuples = new tuple[r];
                for(int i=0;i<r;i++)
                {
                    var tri = triangles[i];
                    var A = tri.GetVertex(0);
                    var B = tri.GetVertex(1);
                    var C = tri.GetVertex(2);
                    double centerU = (A.X + B.X + C.X) / 3d;
                    double centerV = (A.Y + B.Y + C.Y) / 3d;
                    tuples[i] = new tuple(nInnerLoops + nOutterSegments, centerU, centerV,tri.Area);
                    for(int s=0;s<N;s++)
                    {
                        tuples[i].f[s] = Function[s](centerU, centerV);
                        dFunction[s](centerU, centerV, tuples[i].df[s]);
                        ddFunction[s](centerU, centerV, tuples[i].ddf[s]);
                    }
                    //computes kernel...
                    double squaredLength = 0;
                    for (int s = 0; s < N; s++)
                    {
                        squaredLength+=tuples[i].f[s] * tuples[i].f[s];
                    }
                    for (int s = 0; s < N; s++)
                    {
                        for (int t = 0; t < N; t++)
                        {
                            tuples[i].kernel[s, t] = -tuples[i].f[s] * tuples[i].f[t] / squaredLength;
                            if (s == t) tuples[i].kernel[s, t] += 1d;
                        }
                    }
                    double norm=Math.Sqrt(squaredLength);
                    //compute normalized base function
                    for(int s=0;s<N;s++)
                    {
                        tuples[i].nf[s]=tuples[i].f[s]/norm;
                    }
                    // compute normalized first derivative
                    for(int v=0;v<2;v++)
                    {
                        for (int s = 0; s < N; s++)
                        {
                            double val = 0;
                            for (int t = 0; t < N; t++)
                            {
                                val += tuples[i].df[s][v] * tuples[i].kernel[s, t] / norm;
                            }
                            tuples[i].ndf[s][v] = val;
                        }
                    }
                    //compute normalized second derivative
                    for (int v = 0; v < 2; v++)
                    {
                        for (int w = 0; w < 2; w++)
                        {
                            for (int s = 0; s < N; s++)
                            {
                                double val = 0;
                                for (int t = 0; t < N; t++)
                                {
                                    val += tuples[i].ddf[s][v,w] * tuples[i].kernel[s, t] / norm;
                                    val -= 3 * tuples[i].df[s][v] * tuples[i].df[s][w] * tuples[i].kernel[s, t] / norm / squaredLength;
                                }
                                tuples[i].nddf[s][v,w] = val;
                            }
                        }
                    }
                }
            } else { System.Windows.Forms.MessageBox.Show("Not Ready."); }
        }
        protected override void SolveInstance(Grasshopper.Kernel.IGH_DataAccess DA)
        {
            Brep brep = null;
            init();
            if (!DA.GetData(0, ref brep)) { return; }
            var face=brep.Faces[0];
            var domU=face.Domain(0);
            var domV=face.Domain(1);
            int nPt = 50;
            int nPt2 = 30;

            for (int i = 0; i <= nPt; i++)
            {
                double u = domU[0] + (domU[1] - domU[0]) / ((double)nPt) * i;
                var C = face.TrimAwareIsoCurve(1, u);
                foreach (var curve in C)
                {
                    c.Add(curve);
                }
            }
            for (int i = 0; i <= nPt; i++)
            {
                double v = domV[0] + (domV[1] - domV[0]) / ((double)nPt) * i;
                var D = face.TrimAwareIsoCurve(0, v);
                foreach (var curve in D)
                {
                    c.Add(curve);
                }
            }
            for (int i = 0; i <= nPt; i++)
            {
                for (int j = 0; j <= nPt; j++)
                {
                    double u = domU[0] + (domU[1] - domU[0]) / ((double)nPt) * i;
                    double v = domV[0] + (domV[1] - domV[0]) / ((double)nPt) * j;
                    Point3d P;
                    Vector3d[] tmp;
                    var flag=face.Evaluate(u,v,0,out P,out tmp);
                    var C = face.TrimAwareIsoCurve(0, v);
                    var D = face.TrimAwareIsoCurve(1, u);
                    bool flagU=false,flagV=false;
                    if (C.Length > 0)
                    {
                        foreach (var s in C)
                        {
                            var domC = s.Domain;
                            double f;
                            s.ClosestPoint(P, out f);
                            var P2=s.PointAt(f);
                            if((P-P2).Length<0.00000001) flagU = true;
                        }
                    }
                    if (D.Length > 0)
                    {
                        foreach (var s in D)
                        {
                            var domD = s.Domain;
                            double f;
                            s.ClosestPoint(P, out f);
                            var P2 = s.PointAt(f);
                            if ((P - P2).Length < 0.00000001) flagV = true;
                        }
                    }
                    if (flagU && flagV)
                    {
                        a.Add(P);
                        a2.Add(new Point3d(u, v, 0));
                    }
                    else
                    {
                        b.Add(P);
                    }
                }
            }
            InputGeometry input = new InputGeometry();
            int tmpN = 0;
            int N = 0;
            int ss = 100;
            foreach (var loop in face.Loops)
            {
                var _edges3D = loop.To3dCurve();
                var _edges2D = loop.To2dCurve();
                if (_edges3D is PolyCurve)
                {
                    var edges3D = _edges3D as PolyCurve;
                    var edges2D = _edges2D as PolyCurve;
                    nOutterSegments = edges3D.SegmentCount;
                    for (int s = 0; s < edges3D.SegmentCount; s++)
                    {
                        var edge3D = edges3D.SegmentCurve(s);
                        var edge2D = edges2D.SegmentCurve(s);
                        var dom2D = edge2D.Domain;
                        var dom3D=edge3D.Domain;
                        for (int _t = 0; _t <= nPt2; _t++)
                        {
                            double t = dom2D[0] + (dom2D[1] - dom2D[0]) / ((double)nPt2) * _t;
                            var P2D=edge2D.PointAt(t);
                            var P3D=face.PointAt(P2D.X,P2D.Y);
                            d.Add(P3D);
                            d2.Add(P2D);
                            if (_t == nPt2-1 && s == edges3D.SegmentCount - 1)
                            {
                                input.AddPoint(P2D.X, P2D.Y);
                                N++;
                                input.AddSegment(N - 1, tmpN,s+1);
                            }
                            else if(_t<nPt2-1)
                            {
                                if (_t == 0)
                                {
                                    input.AddPoint(P2D.X, P2D.Y);
                                }
                                else
                                {
                                    input.AddPoint(P2D.X, P2D.Y);
                                }
                                N++;
                                input.AddSegment(N - 1, N,s+1);
                            }
                        }
                    }
                    tmpN = N;
                }
                else if(_edges3D is NurbsCurve)
                {
                    var edges3D = _edges3D as NurbsCurve;
                    var edges2D = _edges2D as NurbsCurve;
                    var dom2D=edges2D.Domain;
                    var dom3D=edges3D.Domain;
                    var center = new Point3d(0, 0, 0);
                    for (int _t = 0; _t <= nPt2; _t++)
                    {
                        double t = dom2D[0] + (dom2D[1] - dom2D[0]) / ((double)nPt2) * _t;
                        var P2D = edges2D.PointAt(t);
                        var P3D = face.PointAt(P2D.X, P2D.Y);
                        d.Add(P3D);
                        d2.Add(P2D);
                        if (_t < nPt2-1)
                        {
                            input.AddPoint(P2D.X, P2D.Y);
                            center += P2D;
                            N++;
                            input.AddSegment(N - 1, N,ss);
                        }
                        else if(_t==nPt2-1)
                        {
                            input.AddPoint(P2D.X, P2D.Y);
                            center += P2D;
                            N++;
                            input.AddSegment(N - 1, tmpN,ss);
                        }
                    }
                    ss++;
                    center /= nPt2;
                    input.AddHole(center.X, center.Y);
                    tmpN = N;
                }
            }
            nInnerLoops = ss-100;
            foreach (var l in input.Holes)
            {
                g.Add(new Point3d(l.X, l.Y, 0));
            }
            foreach (var l in input.Segments)
            {
                var P = input.Points.ElementAt(l.P0);
                var Q = input.Points.ElementAt(l.P1);
                f.Add(new Line(new Point3d(P.X, P.Y, 0), new Point3d(Q.X, Q.Y, 0)));
            }


            myControlBox.setNumF(nInnerLoops + nOutterSegments);
            baseFunction = new DoubleArray[nInnerLoops + nOutterSegments];
            coeff = new DoubleArray[nInnerLoops + nOutterSegments];
            Function = new Func<double, double, double>[nInnerLoops + nOutterSegments];
            dFunction = new Action<double, double, double[]>[nInnerLoops + nOutterSegments];
            ddFunction = new Action<double, double, double[,]>[nInnerLoops + nOutterSegments];
            myControlBox.setFunctionToCompute(() =>
            {
                if (lastComputed == nInnerLoops + nOutterSegments - 1) return;
                lastComputed++;
                computeBaseFunction(lastComputed);
                this.ExpirePreview(true);
                myControlBox.EnableRadio(lastComputed, (i) => { resultToPreview(i); this.ExpirePreview(true); });
            }
                );
            Mesher.Mesh mesh = new Mesher.Mesh();
            
            mesh.Behavior.UseBoundaryMarkers = true;
            mesh.Behavior.MaxArea = Math.Pow(Math.Min(input.Bounds.Width, input.Bounds.Height) / 10d, 2);
            mesh.Behavior.Convex = false;
            mesh.Behavior.Algorithm = TriangulationAlgorithm.SweepLine;
            mesh.Behavior.ConformingDelaunay = true;
            mesh.Triangulate(input);
/*            Mesher.Tools.Statistic statistic = new Statistic();
            statistic.Update(mesh, 0);
            mesh.Behavior.MaxArea = statistic.SmallestArea * 1.2;            
            */
            mesh.Behavior.Quality = true;
            mesh.Behavior.MinAngle = 30;
            mesh.Behavior.MaxAngle = 100;
            
            mesh.Refine();
            
            mesh.Smooth();
            mesh.Smooth();
            foreach(var P in mesh.Vertices)
            {
                if (P.Attributes == null)
                {
                    gmesh.Vertices.Add(new Point3d(P.X, P.Y, 0));
                }
                else
                {
                    gmesh.Vertices.Add(new Point3d(P.X, P.Y, P.Attributes[0]));
                }
            }
            foreach(var tri in mesh.Triangles)
            {
                gmesh.Faces.AddFace(tri.P0, tri.P1, tri.P2);
            }

            boundaries = new List<Line>[nOutterSegments];
            holes = new List<Line>[nInnerLoops];
            bbOut = new List<Edge>[nOutterSegments];
            bbIn = new List<Edge>[nInnerLoops];
            foreach (var edge in mesh.Edges)
            {
                var f = edge.Boundary;
                if (f == 0) continue;
                if (f < 100)
                {
                    if (boundaries[f-1] == null)
                    {
                        boundaries[f-1] = new List<Line>();
                    }
                    if (bbOut[f - 1] == null)
                    {
                        bbOut[f - 1] = new List<Edge>();
                    }
                    bbOut[f - 1].Add(new Edge(edge.P0, edge.P1));
                    var P = mesh.Vertices.ElementAt(edge.P0);
                    var Q = mesh.Vertices.ElementAt(edge.P1);
                    boundaries[f-1].Add(new Line(new Point3d(P.X, P.Y, 0), new Point3d(Q.X, Q.Y, 0)));
                }
                else
                {
                    if (holes[f-100] == null)
                    {
                        holes[f-100] = new List<Line>();
                    }
                    if (bbIn[f - 100] == null)
                    {
                        bbIn[f - 100] = new List<Edge>();
                    }
                    bbIn[f - 100].Add(new Edge(edge.P0, edge.P1));
                    var P = mesh.Vertices.ElementAt(edge.P0);
                    var Q = mesh.Vertices.ElementAt(edge.P1);
                    holes[f-100].Add(new Line(new Point3d(P.X, P.Y, 0), new Point3d(Q.X, Q.Y, 0)));
                }
            }
            n = mesh.Vertices.Count();
            m = mesh.Edges.Count();
            r = mesh.Triangles.Count();
            int[,] lines = new int[m, 2];
            int _i = 0;
            foreach (var edge in mesh.Edges)
            {
                lines[_i,0]=edge.P0;
                lines[_i,1]=edge.P1;
                _i++;
            }
            fixedPoints=new List<int>();
            vertices=mesh.Vertices.ToArray();
            edges=mesh.Edges.ToArray();
            triangles = mesh.Triangles.ToArray();
            foreach(var V in vertices)
            {
                if(V.Boundary>0)fixedPoints.Add(Array.IndexOf(vertices,V));
            }
            Laplacian = computeLaplacian(lines, n);
            shiftArray = computeShiftArray(fixedPoints, n);

            //Arranging bb
            foreach (var bbb in bbOut)
            {
                //Look for the first end
                Mesher.Geometry.Edge first=null;
                bool reverse = false;
                foreach (var bbbb in bbb)
                {
                    int P = bbbb.P0;
                    int count = 0;
                    foreach (var bbbbb in bbb)
                    {
                        if (bbbbb.P0 == P || bbbbb.P1 == P) count++;
                    }
                    if (count == 1) { first = bbbb; break; }
                    P = bbbb.P1;
                    count = 0;
                    foreach (var bbbbb in bbb)
                    {
                        if (bbbbb.P0 == P || bbbbb.P1 == P) count++;
                    }
                    if (count == 1) { first = bbbb; reverse = true; break; }                    
                }
                bbb.Remove(first);
                if (reverse)
                {
                    first = new Edge(first.P1, first.P0);
                }
                bbb.Insert(0, first);
                for (int i = 0; i < bbb.Count-1; i++)
                {
                    Mesher.Geometry.Edge next = null;
                    reverse = false;
                    int P = bbb[i].P1;
                    for (int j = i + 1; j < bbb.Count; j++)
                    {
                        if (bbb[j].P0 == P) { next = bbb[j]; break; }
                        if (bbb[j].P1 == P) { next = bbb[j]; reverse = true; break; }
                    }
                    bbb.Remove(next);
                    if (reverse) next = new Edge(next.P1, next.P0);
                    bbb.Insert(i + 1, next);
                }
            }
            foreach (var bbb in bbIn)
            {
                //For inner loops, any edge can be the first edge
                Mesher.Geometry.Edge first = bbb[0];
                bool reverse = false;
                bbb.Remove(first);
                if (reverse)
                {
                    first = new Edge(first.P1, first.P0);
                }
                bbb.Insert(0, first);
                for (int i = 0; i < bbb.Count - 1; i++)
                {
                    Mesher.Geometry.Edge next = null;
                    reverse = false;
                    int P = bbb[i].P1;
                    for (int j = i + 1; j < bbb.Count; j++)
                    {
                        if (bbb[j].P0 == P) { next = bbb[j]; break; }
                        if (bbb[j].P1 == P) { next = bbb[j]; reverse = true; break; }
                    }
                    bbb.Remove(next);
                    if (reverse) next = new Edge(next.P1, next.P0);
                    bbb.Insert(i + 1, next);
                }
            }
        }
    }
}
