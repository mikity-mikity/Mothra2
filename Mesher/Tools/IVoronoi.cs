﻿// -----------------------------------------------------------------------
// <copyright file="IVoronoi.cs" company="">
// Triangle.NET code by Christian Woltering, http://triangle.codeplex.com/
// </copyright>
// -----------------------------------------------------------------------

namespace Mesher.Tools
{
    using System.Collections.Generic;
    using Mesher.Geometry;

    /// <summary>
    /// TODO: Update summary.
    /// </summary>
    public interface IVoronoi
    {
        /// <summary>
        /// Gets the list of Voronoi vertices.
        /// </summary>
        Point[] Points { get; }

        /// <summary>
        /// Gets the list of Voronoi regions.
        /// </summary>
        List<VoronoiRegion> Regions { get; }
    }
}
