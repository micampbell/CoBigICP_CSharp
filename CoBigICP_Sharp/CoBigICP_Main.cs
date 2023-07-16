using HelixToolkit.SharpDX.Core.Core2D;
using Newtonsoft.Json.Linq;
using System.Data;
using System.Diagnostics.Metrics;
using System.Runtime.CompilerServices;
using System.Transactions;
using System.Windows.Controls;
using System.Windows.Documents;
using System.Windows.Media.Imaging;
using System.Windows.Media.Media3D;
using TVGL;
using TVGL.PointCloud;

namespace CoBigICP_Sharp
{
    public static partial class CoBigICP
    {
        internal const double cosd30 = 0.8660254037844387;
        //        function[R, T, Flag] = CoBigICP_fun(Md, Mo, MovData, RefData, MovInfo, RefInfo, Tf0, DistThr, AngThr, sigma_times)
        public static Matrix4x4 Run(IList<Vector3> MovData, IList<Vector3> RefData,
            double DistThr = double.PositiveInfinity, double AngThr = cosd30)
        {
            return Run(KDTree.Create(MovData, Enumerable.Range(0, MovData.Count).ToList()),
                KDTree.Create(RefData, Enumerable.Range(0, MovData.Count).ToList()), DistThr, AngThr);
        }

        private static Matrix4x4 Run(KDTree<Vector3, int> MovCloud, KDTree<Vector3, int> RefCloud,
                      double DistThr = double.PositiveInfinity, double AngThr = cosd30)
        {
            return Run(MovCloud, CalculateNormalInfo(MovCloud),
                RefCloud, CalculateNormalInfo(RefCloud), DistThr, AngThr);
        }

        private static Matrix4x4 Run(KDTree<Vector3, int> movCloud, NormalInfo[] movingNormalInfo, KDTree<Vector3, int> refCloud,
            NormalInfo[] refNormalInfo, double distThr, double angThr)
        {
            var T = MakeInitialTranslationMatrix(movCloud, refCloud);
            var R = Matrix4x4.Identity;
            var dR = Matrix4x4.Identity;
            var dT = Matrix4x4.Identity;
            var success = false;
            var sigma_itr = 31.62;
            var DecayPram = 0.97;  // 1.03;
            var U0 = ConcatenateAlongColumns(refNormalInfo.Select(n => n.U), refNormalInfo.Length);
            var U1 = ConcatenateAlongColumns(movingNormalInfo.Select(n => n.U), movingNormalInfo.Length);
            var Norm_Ref = ConcatenateAlongColumns(refNormalInfo.Select(n => n.Normal), refNormalInfo.Length);//cat(2, RefInfo(:).normal);
            var Norm_Mov = ConcatenateAlongColumns(movingNormalInfo.Select(n => n.Normal), movingNormalInfo.Length);//cat(2, RefInfo(:).normal);
            var epsilon = 1e-3;
            var MaxIter = 50;
            var JArray = new List<double>();
            for (int nIter = 0; nIter < MaxIter; nIter++)
            {
                var T_curr = R*T; // [R T; 0 0 0 1];
                var RTranspose = R.Transpose();
                Matrix4x4.Invert(T_curr, out var T_curr_inv);
                //            R_inv = T_curr_inv(1:3, 1:3); T_inv = T_curr_inv(1:3, end);
                var AftData = movCloud.OriginalPoints.Select(p => p.Transform(RTranspose*T)).ToList();  // apply transformation to move points.

                var refAftData = refCloud.OriginalPoints.Select(p => p.Transform(T_curr_inv)).ToList();  // apply transformation to ref points.
                                                                                                         //                AftData = Loc2Glo(MovData, R', T );   % apply transformation to move points.

                //  refAftData = Loc2Glo(RefData, R_inv', T_inv );   % apply transformation to ref points.

                var NNidx = AftData.SelectMany(p => refCloud.FindNearest(p, 1)).ToList();
                //  [NNIdx, DD] = knnsearch(refCloud, AftData' ); % establish correspondence.
                //  % bidirectional correspondence
                var NNIdx_inverse = refAftData.SelectMany(p => movCloud.FindNearest(p, 1)).ToList();
                //  [NNIdx_inverse, ~] = knnsearch(movCloud, refAftData' );
                var bi_eff = Enumerable.Range(0, AftData.Count)
                    .Where(i => AftData[i].Distance(AftData[NNIdx_inverse[i].Item2]) < 0.01).ToList();
                //var Angle = sum(Norm_Ref(:, NNIdx).*(R * Norm_Mov));
                //                EffIdx_sim = find(abs(Angle) > AngThr & DD' < DistThr );

                //                EffIdx = intersect(EffIdx_sim, bi_eff);
                //                MovIdx = EffIdx;
                //                RefIdx = NNIdx(EffIdx);

                //                tmpAft = AftData(:, MovIdx);
                //                tmpRef = RefData(:, RefIdx);


                //    %%%%%%%%%% after transformation, we need align AftData to RefData.
                //    %%%%%%%%%%% obtain rotation and translation via solving quadratic programming.

                //                    [H, b, J] = CalHbCobig_Gabor(tmpAft, tmpAft - tmpRef, R, RefInfo(RefIdx), MovInfo(MovIdx), sigma_itr);
                //                dx = -pinv(H) * b;
                //                dR = expm(SkewFun(dx(1:3)));
                //                dT = dx(4:6);
                //                R = R * dR;
                //                T = T + dT;
                //                bTest = 1;
                //                
                //                Err = max(norm(dR - eye(3)), norm(dT));
                //                JArray(end + 1) = J / length(MovIdx);
                //                if length(JArray) >= 2
                //                    if abs((JArray(end) - JArray(end - 1)) / JArray(end - 1)) <= 1e-4
                //            break;
                //                end
                //            end
                //    str = sprintf('nIter = %02d/%02d, J = %.6f', nIter, MaxIter, JArray(end));
                //                disp(str);
                //                bTest = 1;
                //                sigma_itr = sigma_itr * DecayPram;
            }

            //                Tf_est = [R T; 0 0 0 1];
            //% Tf_gt
            //bTest = 1;
            //                end
            return R * T;
        }

        private static IEnumerable<Vector3> Loc2Glo(IEnumerable<Vector3> points, Matrix4x4 transform)
        {
            return points.Select(p => p.Transform(transform));
            // the input pcData must be Dim* N format.
            //function[Pts, varargout] = Loc2Glo(pcData, varargin)
            //    [Dim, DataLen] = size(pcData);
            //                if Dim ~= 2 && Dim ~= 3
            //        error('the input pcData must be 2 * N or 3 * N format!');
            //                end
            //    if (nargin - 1) == 2
            //        R = varargin{ 1};
            //                T = varargin{ 2};
            //                end
            //    if (nargin - 1) == 1
            //        Pose = varargin{ 1};
            //                [R T] = ExtractRT(Pose);
            //                end
            //                Pts = bsxfun(@plus, inv(R) * pcData, T);
            //                if (nargout - 1) == 1
            //        varargout{ 1} = [R T];
            //                end
            //    if (nargout - 1) == 2
            //        varargout{ 1} = R;
            //                varargout{ 2} = T;
            //                end
            //            end


        }

        private static double[,] ConcatenateAlongColumns(IEnumerable<Vector3> vectors, int length)
        {
            var result = new double[3, length];
            var i = 0;
            foreach (var v in vectors)
            {
                result[0, i] = v.X;
                result[1, i] = v.Y;
                result[2, i] = v.Z;
                i++;
            }
            return result;
        }

        private static double[,] ConcatenateAlongColumns(IEnumerable<Matrix3x3> matrices, int length)
        {
            var result = new double[3, 3 * length];
            var i = 0;
            foreach (var matrix in matrices)
            {
                result[0, i] = matrix.M11;
                result[1, 1] = matrix.M21;
                result[2, i] = matrix.M31;
                i++;
                result[0, i] = matrix.M12;
                result[1, i] = matrix.M22;
                result[2, i] = matrix.M32;
                i++;
                result[0, i] = matrix.M13;
                result[1, i] = matrix.M23;
                result[2, i] = matrix.M33;
                i++;
            }
            return result;
        }

        private static Matrix4x4 MakeInitialTranslationMatrix(KDTree<Vector3> movCloud, KDTree<Vector3> refCloud)
        {
            var difference = Vector3.Zero;
            var n = Math.Min(movCloud.Count, refCloud.Count);
            for (int i = 0; i < n; i++)
                difference += refCloud.OriginalPoints[i] - movCloud.OriginalPoints[i];
            difference /= n;
            return new Matrix4x4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, difference.X, difference.Y, difference.Z, 1);
        }
    }
}
