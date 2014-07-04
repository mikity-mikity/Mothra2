using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mikity.ghComponents
{
    public partial class Mothra2 : Grasshopper.Kernel.GH_Component
    {

        public override void DrawViewportWires(Grasshopper.Kernel.IGH_PreviewArgs args)
        {
            if (Hidden)
            {
                return;
            }
            if (boundaries != null)
            {
                for (int i = 0; i < boundaries.Count(); i++)
                {
                    var boundary = boundaries[i];
                    System.Drawing.Color color = System.Drawing.Color.FromArgb(i * 30, 255 - i * 30, i * 60);
                    foreach (var bb in boundary)
                    {
                        args.Display.DrawLine(bb, color);
                    }
                }
            }

            if (f != null)
            {
                foreach (var l in f)
                {
                    args.Display.DrawLine(l, System.Drawing.Color.Bisque);
                }
            }
            if (g != null)
            {
                foreach (var l in g)
                {
                    args.Display.DrawPoint(l, System.Drawing.Color.Bisque);
                }
            }
            if (gmesh != null)
            {
                args.Display.DrawMeshWires(gmesh, System.Drawing.Color.Red);
            }
            if (a != null)
            {
                foreach (var point in a)
                {
                    args.Display.DrawPoint(point, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Blue);
                }
            }
            if (a2 != null)
            {
                foreach (var point in a2)
                {
                    args.Display.DrawPoint(point, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Blue);
                }
            }
            if (b != null)
            {
                foreach (var point in b)
                {
                    args.Display.DrawPoint(point, Rhino.Display.PointStyle.ControlPoint, 2, System.Drawing.Color.Orange);
                }
            }
            if (d != null)
            {
                foreach (var point in d)
                {
                    args.Display.DrawPoint(point, Rhino.Display.PointStyle.Simple, 2, System.Drawing.Color.Brown);
                }
            }
            if (d2 != null)
            {
                foreach (var point in d2)
                {
                    args.Display.DrawPoint(point, Rhino.Display.PointStyle.Simple, 2, System.Drawing.Color.Brown);
                }
            }
            if (c != null)
            {
                foreach(var curve in c)
                {
                    args.Display.DrawCurve(curve, System.Drawing.Color.Red);
                }
            }
        }
    }
}
