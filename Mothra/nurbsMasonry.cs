using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using Rhino.Geometry;

using Mesher;
using Mesher.Geometry;
using Mesher.Tools;

namespace mikity.ghComponents
{
    public partial class Mothra2 : Grasshopper.Kernel.GH_Component
    {
        List<Point3d> a;
        List<Point3d> a2;
        List<Point3d> b;
        List<Curve> c;
        List<Point3d> d;
        List<Point3d> d2;
        Rhino.Geometry.Mesh gmesh = new Rhino.Geometry.Mesh();
        private void init()
        {
            a = new List<Point3d>();
            a2 = new List<Point3d>();
            b = new List<Point3d>();
            c = new List<Curve>();
            d = new List<Point3d>();
            d2 = new List<Point3d>();
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
        protected override void SolveInstance(Grasshopper.Kernel.IGH_DataAccess DA)
        {
            Brep brep = null;
            init();
            if (!DA.GetData(0, ref brep)) { return; }
            var face=brep.Faces[0];
            var domU=face.Domain(0);
            var domV=face.Domain(1);

            for (int i = 0; i <= 50; i++)
            {
                double u = domU[0] + (domU[1] - domU[0]) / 50d * i;
                var C = face.TrimAwareIsoCurve(1, u);
                foreach (var curve in C)
                {
                    c.Add(curve);
                }
            }
            for (int i = 0; i <= 50; i++)
            {
                double v = domV[0] + (domV[1] - domV[0]) / 50d * i;
                var D = face.TrimAwareIsoCurve(0, v);
                foreach (var curve in D)
                {
                    c.Add(curve);
                }
            }
            for (int i = 0; i <= 50; i++)
            {
                for (int j = 0; j <= 50; j++)
                {
                    double u = domU[0] + (domU[1] - domU[0]) / 50d * i;
                    double v = domV[0] + (domV[1] - domV[0]) / 50d * j;
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
                            if (f > domC[0] && f < domC[1]) flagU = true;
                        }
                    }
                    if (D.Length > 0)
                    {
                        foreach (var s in D)
                        {
                            var domD = s.Domain;
                            double f;
                            s.ClosestPoint(P, out f);
                            if (f > domD[0] && f < domD[1]) flagV = true;
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
            InputGeometry input = new InputGeometry(12);
            foreach (var loop in face.Loops)
            {
                var _edges3D = loop.To3dCurve();
                var _edges2D = loop.To2dCurve();
                int tmpN = 0;
                int N = 0;
                if (_edges3D is PolyCurve)
                {
                    var edges3D = _edges3D as PolyCurve;
                    var edges2D = _edges2D as PolyCurve;
                    for (int s = 0; s < edges3D.SegmentCount; s++)
                    {
                        var edge3D = edges3D.SegmentCurve(s);
                        var edge2D = edges2D.SegmentCurve(s);
                        var dom2D = edge2D.Domain;
                        var dom3D=edge3D.Domain;
                        for (int _t = 0; _t <= 50; _t++)
                        {
                            double t = dom2D[0] + (dom2D[1] - dom2D[0]) / 50d * _t;
                            var P2D=edge2D.PointAt(t);
                            var P3D=face.PointAt(P2D.X,P2D.Y);
                            d.Add(P3D);
                            d2.Add(P2D);
                            if (_t == 49 && s == edges3D.SegmentCount - 1)
                            {
                                input.AddPoint(P2D.X, P2D.Y);
                                N++;
                                input.AddSegment(N - 1, tmpN);
                            }
                            else if(_t<49)
                            {
                                input.AddPoint(P2D.X, P2D.Y);
                                N++;
                                input.AddSegment(N - 1, N);
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
                    for (int _t = 0; _t <= 50; _t++)
                    {
                        double t = dom2D[0] + (dom2D[1] - dom2D[0]) / 50d * _t;
                        var P2D = edges2D.PointAt(t);
                        var P3D = face.PointAt(P2D.X, P2D.Y);
                        d.Add(P3D);
                        d2.Add(P2D);
                        if (_t < 49)
                        {
                            input.AddPoint(P2D.X, P2D.Y);
                            center += P2D;
                            N++;
                            input.AddSegment(N - 1, N);
                        }
                        else if(_t==49)
                        {
                            input.AddPoint(P2D.X, P2D.Y);
                            center += P2D;
                            N++;
                            input.AddSegment(N - 1, tmpN);
                        }
                    }
                    center /= 50;
                    //input.AddHole(center.X, center.Y);
                    tmpN = N;
                }
            }
            //input.AddHole(-23.9180922915,-17.2936274735);
            //input.AddHole(-17.6846493856, -18.4197444922);




            Mesher.Mesh mesh = new Mesher.Mesh();
            //mesh.Behavior.MaxArea = Math.Pow(Math.Min(input.Bounds.Width, input.Bounds.Height) / 10d, 2);
            //mesh.Behavior.Quality = true;
            //mesh.Behavior.MinAngle = 25;
            //mesh.Behavior.MaxAngle = 110;
            //mesh.Behavior.Convex = false;
            mesh.Behavior.Algorithm = TriangulationAlgorithm.Dwyer;
            mesh.Triangulate(input);
            /*mesh.Behavior.MaxArea = mesh.Behavior.MaxArea*0.5;
            mesh.Refine();
            mesh.Behavior.MaxArea = mesh.Behavior.MaxArea * 0.5;
            mesh.Refine();
            mesh.Behavior.MaxArea = mesh.Behavior.MaxArea * 0.5;
            mesh.Refine();
            Mesher.Tools.Statistic statistic = new Statistic();
            statistic.Update(mesh, 0);
            mesh.Behavior.MaxArea = statistic.SmallestArea * 5;
            mesh.Refine();
            mesh.Smooth();
            mesh.Smooth();*/
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                gmesh.Vertices.Add(new Point3d(mesh.Vertices.ElementAt(i).X, mesh.Vertices.ElementAt(i).Y, 0));
            }
            for (int i = 0; i < mesh.Triangles.Count; i++)
            {
                gmesh.Faces.AddFace(mesh.Triangles.ElementAt(i).P0, mesh.Triangles.ElementAt(i).P1, mesh.Triangles.ElementAt(i).P2);
            }
        }
    }
}
