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
            if (result != null)
            {
                args.Display.DrawLines(result, System.Drawing.Color.Red);
            }
            if (boundaries != null)
            {
                for (int i = 0; i < boundaries.Count(); i++)
                {
                    var boundary = boundaries[i];
                    if (boundary != null)
                    {
                        System.Drawing.Color color = System.Drawing.Color.Red;
                        if (i == 0) color = System.Drawing.Color.Red;
                        if (i == 1) color = System.Drawing.Color.Blue;
                        if (i == 2) color = System.Drawing.Color.Green;
                        if (i == 3) color = System.Drawing.Color.Yellow;
                        if (i == 4) color = System.Drawing.Color.Orange;
                        if (i == 5) color = System.Drawing.Color.Purple;
                        if (i == 6) color = System.Drawing.Color.Brown;
                        foreach (var bb in boundary)
                        {
                            args.Display.DrawLine(bb, color, 5);

                        }
                    }
                }
            }
            if (holes != null)
            {
                for (int i = 0; i < holes.Count(); i++)
                {
                    var hole = holes[i];
                    if (hole != null)
                    {
                        System.Drawing.Color color = System.Drawing.Color.Red;
                        if (i == 0) color = System.Drawing.Color.Pink;
                        if (i == 1) color = System.Drawing.Color.LightBlue;
                        if (i == 2) color = System.Drawing.Color.LightGreen;
                        if (i == 3) color = System.Drawing.Color.LightYellow;
                        if (i == 4) color = System.Drawing.Color.LightCoral;
                        if (i == 5) color = System.Drawing.Color.LightSeaGreen;
                        if (i == 6) color = System.Drawing.Color.LightCyan;
                        foreach (var bb in hole)
                        {
                            args.Display.DrawLine(bb, color, 5);

                        }
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
                args.Display.DrawPoints(g, Rhino.Display.PointStyle.ActivePoint, 2, System.Drawing.Color.Bisque);
            }
            if (gmesh != null)
            {
                args.Display.DrawMeshWires(gmesh, System.Drawing.Color.Gray);
            }
            if (a != null)
            {
                args.Display.DrawPoints(a, Rhino.Display.PointStyle.X, 2, System.Drawing.Color.Blue);

            }
            if (a2 != null)
            {
                args.Display.DrawPoints(a2, Rhino.Display.PointStyle.X, 10, System.Drawing.Color.Blue);
            }
            if (b != null)
            {
                args.Display.DrawPoints(b, Rhino.Display.PointStyle.ControlPoint, 2, System.Drawing.Color.Orange);
            }
            if (b2 != null)
            {
                args.Display.DrawPoints(b2, Rhino.Display.PointStyle.ControlPoint, 2, System.Drawing.Color.Orange);
            }
            if (d != null)
            {
                args.Display.DrawPoints(d, Rhino.Display.PointStyle.Simple, 2, System.Drawing.Color.Brown);
            }
            //Integrating Points
            if (dd2 != null)
            {
                args.Display.DrawPoints(dd2, Rhino.Display.PointStyle.X, 8, System.Drawing.Color.Violet);
            }
            if (dd != null)
            {
                args.Display.DrawPoints(dd, Rhino.Display.PointStyle.ControlPoint, 4, System.Drawing.Color.LimeGreen);
            }
            //base vectors
            if (basis2 != null)
            {
                args.Display.DrawLines(basis2, System.Drawing.Color.Red, 4);
            }
            if (basis != null)
            {
                args.Display.DrawLines(basis, System.Drawing.Color.Blue, 2);
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
