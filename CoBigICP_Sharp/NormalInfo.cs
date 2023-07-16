using StarMathLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using TVGL;
using TVGL.PointCloud;

namespace CoBigICP_Sharp
{
    public static partial class CoBigICP
    {

        static NormalInfo[] CalculateNormalInfo(KDTree<Vector3> data)
        {
            var numPoints = data.Count;
            var normals = new NormalInfo[numPoints];
            var numNeighbors = Math.Min(20, numPoints - 1);
            var S = new Matrix3x3(0.001, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
            //var mean =Vector3.Zero;
            //for (int i = 0; i < numPoints; i++)
            //    mean += data.Points[i];
            //mean /= numPoints;
            //var centeredData = new List<Vector3>(numPoints);
            //for (int i = 0; i < numPoints; i++)
            //{
            //    centeredData.Add(data.Points[i] - mean);
            //}
            for (int i = 0; i < numPoints; i++)
            {
                var point = data.OriginalPoints[i];
                var neighbors = data.FindNearest(point, numNeighbors).ToList();
                Matrix3x3 covariance = MakeCovarianceMatrix(neighbors);
                var eigenValues = covariance.GetEigenValuesAndVectors(out var eigenVectors);
                int[] orderOfEigens = OrderEigenValues(eigenValues[0], eigenValues[1], eigenValues[2]);

                var u = new Matrix3x3(eigenVectors[orderOfEigens[0]].X, eigenVectors[orderOfEigens[1]].X,eigenVectors[orderOfEigens[2]].X,
                    eigenVectors[orderOfEigens[0]].Y, eigenVectors[orderOfEigens[1]].Y, eigenVectors[orderOfEigens[2]].Y,
                    eigenVectors[orderOfEigens[0]].Z, eigenVectors[orderOfEigens[1]].Z, eigenVectors[orderOfEigens[2]].Z);
                covariance = u * S * (u.Transpose());
                Matrix3x3.Invert(covariance, out var omega);

                normals[i] = new NormalInfo { U = u, Normal = eigenVectors[orderOfEigens[0]], Omega = omega };
            }
            return normals;
        }

        private static int[] OrderEigenValues(ComplexNumber e1, ComplexNumber e2, ComplexNumber e3)
        {
            var e1Mag = e1.LengthSquared();
            var e2Mag = e2.LengthSquared();
            var e3Mag = e3.LengthSquared();
            if (e1Mag < e2Mag)
            {
                if (e1Mag < e3Mag)
                {
                    if (e2Mag < e3Mag) return new int[] { 0, 1, 2 };
                    else return new int[] { 0, 2, 1 };
                }
                else return new int[] { 2, 0, 1 };
            }
            else if (e2Mag < e3Mag)
            {
                if (e1Mag < e3Mag) return new int[] { 1, 0, 2 };
                else return new int[] { 1, 2, 0 };
            }
            else return new int[] { 2, 1, 0 };
        }
    

    private static Matrix3x3 MakeCovarianceMatrix(List<Vector3> vectors)
    {
        var num = vectors.Count;
        var m11 = 0.0;
        var m12 = 0.0;
        var m13 = 0.0;
        var m22 = 0.0;
        var m23 = 0.0;
        var m33 = 0.0;
        var average = Vector3.Zero;
        foreach (var neigbor in vectors)
            average += neigbor;
        average /= num;
        var centeredPoints = vectors.Select(n => n - average).ToList();
        foreach (var point in centeredPoints)
        {
            m11 += point.X * point.X;
            m12 += point.X * point.Y;
            m13 += point.X * point.Z;
            m22 += point.Y * point.Y;
            m23 += point.Y * point.Z;
            m33 += point.Z * point.Z;
        }
        var covariance = new Matrix3x3(m11, m12, m13, m12, m22, m23, m13, m23, m33);
        covariance *= 1.0 / (num - 1);
        return covariance;
    }

    internal readonly struct NormalInfo
    {
        public Vector3 Normal { get; init; }
        public Matrix3x3 U { get; init; }
        public Matrix3x3 Omega { get; init; }
    }
}
}
