﻿using StarMathLib;
using System.Data;
using TVGL;
using TVGL.PointCloud;

namespace CoBigICP_Sharp
{
    public static partial class CoBigICP
    {
        internal const double cosd30 = 0.8660254037844387;
        //        function[R, T, Flag] = CoBigICP_fun(Md, Mo, MovData, RefData, MovInfo, RefInfo, Tf0, DistThr, AngThr, sigma_times)
        public static Matrix4x4 Run(IList<Vector3> MovData, IList<Vector3> RefData,
            double DistThr = double.PositiveInfinity, double AngThr = 0.0)
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
            var forwardTransform = GetTranslationMatrix(movCloud.OriginalPoints, refCloud.OriginalPoints);
            var R = Matrix4x4.Identity;
            //var dR = Matrix4x4.Identity;
            //var dT = Matrix4x4.Identity;
            var success = false;
            var sigma_itr = 31.62;
            var DecayPram = 0.97;  // 1.03;
            var U0 = ConcatenateAlongColumns(refNormalInfo.Select(n => n.U), refNormalInfo.Length);
            var U1 = ConcatenateAlongColumns(movingNormalInfo.Select(n => n.U), movingNormalInfo.Length);
            var Norm_Ref = refNormalInfo.Select(n => n.Normal).ToList();
            //var Norm_Ref = ConcatenateAlongColumns(refNormalInfo.Select(n => n.Normal), refNormalInfo.Length);//cat(2, RefInfo(:).normal);
            var Norm_Mov = movingNormalInfo.Select(n => n.Normal).ToList();
            //var Norm_Mov = ConcatenateAlongColumns(movingNormalInfo.Select(n => n.Normal), movingNormalInfo.Length);//cat(2, RefInfo(:).normal);
            var epsilon = 1e-3;
            var MaxIter = 50;
            var JArray = new List<double>();
            for (int nIter = 0; nIter < MaxIter; nIter++)
            {
                //var forwardTransform = AddTtoR(R, TVector);
                //Matrix4x4.Invert(R, out var R_inv);
                Matrix4x4.Invert(forwardTransform, out var backTransform);
                //Matrix4x4 backTransform = AddTtoR(R_inv, -TVector);
                var refAftData = refCloud.OriginalPoints.Select(p => p.Transform(forwardTransform)).ToList();  // apply transformation to ref points.
                var AftData = movCloud.OriginalPoints.Select(p => p.Transform(backTransform)).ToList();  // apply transformation to move points.

                var NNidx = AftData.SelectMany(p => refCloud.FindNearest(p, 1)).ToList();
                //  [NNIdx, DD] = knnsearch(refCloud, AftData' ); % establish correspondence.
                //  % bidirectional correspondence
                var NNIdx_inverse = refAftData.SelectMany(p => movCloud.FindNearest(p, 1)).ToList();
                //  [NNIdx_inverse, ~] = knnsearch(movCloud, refAftData' );
                var bi_eff = Enumerable.Range(0, AftData.Count)
                    .Where(i => AftData[i].Distance(AftData[NNIdx_inverse[NNidx[i].Item2].Item2]) < 0.1).ToList();
                var Angle = GetAngle(Norm_Ref, Norm_Mov, NNidx, R.Transpose());
                //var Angle = sum(Norm_Ref(:, NNIdx).*(R * Norm_Mov));
                var EffIdx_sim = Enumerable.Range(0, Angle.Count)
                    .Where(i => Math.Abs(Angle[i]) > angThr && NNidx[i].Item1.Distance(AftData[i]) < distThr).ToList();
                //                EffIdx_sim = find(abs(Angle) > AngThr & DD' < DistThr );
                var EffIdx = EffIdx_sim.Intersect(bi_eff).ToList();
                //                EffIdx = intersect(EffIdx_sim, bi_eff);
                var MovIdx = EffIdx;
                //                MovIdx = EffIdx;
                var RefIdx = EffIdx.Select(index => NNidx[index].Item2);
                //                RefIdx = NNIdx(EffIdx);
                var tmpAft = MovIdx.Select(index => AftData[index]).ToList();
                //                tmpAft = AftData(:, MovIdx);
                var tmpRef = RefIdx.Select(index => refCloud.OriginalPoints[index]);
                //                tmpRef = RefData(:, RefIdx);
                var diffTemp = tmpAft.Zip(tmpRef, (a, b) => a - b).ToList();
                var relevantRefNormalInfo = RefIdx.Select(index => refNormalInfo[index]).ToList();
                var relevantMovNormalInfo = MovIdx.Select(index => movingNormalInfo[index]).ToList();
                //    %%%%%%%%%% after transformation, we need align AftData to RefData.
                //    %%%%%%%%%%% obtain rotation and translation via solving quadratic programming.
                CalHbCobig_Gabor(tmpAft, diffTemp, R, relevantRefNormalInfo, relevantMovNormalInfo, sigma_itr, out var H, out var b, out var J);
                //                    [H, b, J] = CalHbCobig_Gabor(tmpAft, tmpAft - tmpRef, R, RefInfo(RefIdx), MovInfo(MovIdx), sigma_itr);
                StarMathLib.StarMath.solve(H, b, out var dx);
                dx = dx.Select(x => -x).ToArray();
                //                dx = -pinv(H) * b;
                double[,] dRMatrix = StarMath.ExpMatrix(SkewFun(dx.Take(3).ToArray()));
                //                dR = expm(SkewFun(dx(1:3)));
                var dR = new Matrix4x4(dRMatrix[0, 0], dRMatrix[0, 1], dRMatrix[0, 2], 0,
                    dRMatrix[1, 0], dRMatrix[1, 1], dRMatrix[1, 2], 0,
                    dRMatrix[2, 0], dRMatrix[2, 1], dRMatrix[2, 2], 0,
                    0, 0, 0, 1);
                R = R * dR;
                //var dTVector = new Vector3(dx[3], dx[4], dx[5]);
                //TVector -= dTVector;
                var Transl = Matrix4x4.CreateTranslation(forwardTransform.M41- dx[3], forwardTransform.M42 - dx[4], forwardTransform.M43 - dx[5]);
                forwardTransform = R * Transl;
                var t2 =GetTranslationVector(movCloud.OriginalPoints,
                    refCloud.OriginalPoints.Select(p => p.Transform(R)));  
                //Console.WriteLine(t2.ToString());
                //dT = Matrix4x4.CreateTranslation(dTVector);
                //                dT = dx(4:6);
                //T = T * dT;
                //                T = T + dT;
                var bTest = 1;
                //                bTest = 1;
                //var Err = Math.Max(norm(dR - Matrix4x4.Identity), dTVector.Length());
                //                Err = max(norm(dR - eye(3)), norm(dT));
                JArray.Add(J / MovIdx.Count);
                //                JArray(end + 1) = J / length(MovIdx);
                if (JArray.Count >= 2 && Math.Abs((JArray[^1] - JArray[^2]) / JArray[^2]) <= 1e-4)
                    break;
                //                if length(JArray) >= 2
                //                    if abs((JArray(end) - JArray(end - 1)) / JArray(end - 1)) <= 1e-4
                //            break;
                //                end
                //            end
                Console.WriteLine("nIter = " + nIter + "/" + MaxIter + ", J = " + JArray[^1], 1);
                //    str = sprintf('nIter = %02d/%02d, J = %.6f', nIter, MaxIter, JArray(end));
                //                disp(str);
                bTest = 1;
                //                bTest = 1;
                sigma_itr = sigma_itr * DecayPram;
                //                sigma_itr = sigma_itr * DecayPram;
            }

            //                Tf_est = [R T; 0 0 0 1];
            //% Tf_gt
            //bTest = 1;
            //                end
            return forwardTransform;  // R * T;
            //return AddTtoR(R, TVector);  // R * T;
        }

        private static Matrix4x4 AddTtoR(Matrix4x4 m, Vector3 t)
        {
            return new Matrix4x4(m.M11, m.M12, m.M13, m.M14, m.M21, m.M22, m.M23, m.M24,
                m.M31, m.M32, m.M33, m.M34, t.X, t.Y, t.Z, 1);
        }

        private static double norm(Matrix4x4 a4x4)
        {
            var aM = new double[,] { { a4x4.M11, a4x4.M12, a4x4.M13, a4x4.M14 },
                {a4x4.M21, a4x4.M22, a4x4.M23, a4x4.M24  }, {a4x4.M31, a4x4.M32, a4x4.M33, a4x4.M34  }, {a4x4.M41, a4x4.M42, a4x4.M43, a4x4.M44  } };
            var aMt = aM.Transpose();
            var maxEigenValue = 0.0;
            try
            {
                var eigenValues = StarMathLib.StarMath.GetEigenValues(aMt.multiply(aM));
                foreach (var eigenValue in eigenValues)
                {
                    var magSqd = eigenValue.Real * eigenValue.Real + eigenValue.Imaginary * eigenValue.Imaginary;
                    if (magSqd > maxEigenValue)
                        maxEigenValue = magSqd;
                }
            }
            catch { }
            return Math.Sqrt(Math.Sqrt(maxEigenValue));
        }

        private static void CalHbCobig_Gabor(IList<Vector3> Mov, List<Vector3> PtsDiff, Matrix4x4 R,
            List<NormalInfo> RefNormInfo, List<NormalInfo> MovNormInfo, double sigma, out double[,] H, out double[] b, out double J)
        {
            var p1 = Mov.Select(p => p.X).ToList();
            var p2 = Mov.Select(p => p.Y).ToList();
            var p3 = Mov.Select(p => p.Z).ToList();

            var v1 = PtsDiff.Select(p => p.X).ToList();
            var v2 = PtsDiff.Select(p => p.Y).ToList();
            var v3 = PtsDiff.Select(p => p.Z).ToList();

            var infoRMat = ConcatenateAlongColumns(RefNormInfo.Select(rni => rni.Omega), RefNormInfo.Count);
            var infoR1 = RefNormInfo.Select(rni => rni.Omega.M11).ToList();
            //infoR1 = infoRMat(1, 1:3:end);
            var infoR2 = RefNormInfo.Select(rni => rni.Omega.M12).ToList();
            //infoR2 = infoRMat(1, 2:3:end);
            var infoR3 = RefNormInfo.Select(rni => rni.Omega.M13).ToList();
            //infoR3 = infoRMat(1, 3:3:end);
            var infoR4 = RefNormInfo.Select(rni => rni.Omega.M22).ToList();
            //infoR4 = infoRMat(2, 2:3:end);
            var infoR5 = RefNormInfo.Select(rni => rni.Omega.M23).ToList();
            //infoR5 = infoRMat(2, 3:3:end);
            var infoR6 = RefNormInfo.Select(rni => rni.Omega.M33).ToList();
            //infoR6 = infoRMat(3, 3:3:end);
            //infoMMat = cat(2, MovNormInfo(:).omega);
            var infoM1 = MovNormInfo.Select(mni => mni.Omega.M11).ToList();
            //infoM1 = infoMMat(1, 1:3:end);
            var infoM2 = MovNormInfo.Select(mni => mni.Omega.M12).ToList();
            //infoM2 = infoMMat(1, 2:3:end);
            var infoM3 = MovNormInfo.Select(mni => mni.Omega.M13).ToList();
            //infoM3 = infoMMat(1, 3:3:end);
            var infoM4 = MovNormInfo.Select(mni => mni.Omega.M22).ToList();
            //infoM4 = infoMMat(2, 2:3:end);
            var infoM5 = MovNormInfo.Select(mni => mni.Omega.M23).ToList();
            //infoM5 = infoMMat(2, 3:3:end);
            var infoM6 = MovNormInfo.Select(mni => mni.Omega.M33).ToList();
            //infoM6 = infoMMat(3, 3:3:end);
            var R1 = R.M11;
            var R2 = R.M12;
            var R3 = R.M13;
            var R4 = R.M21;
            var R5 = R.M22;
            var R6 = R.M23;
            var R7 = R.M31;
            var R8 = R.M32;
            var R9 = R.M33;

            var infobi = new List<double[]>();
            var m1 = ZipCoefficientAdd((1, infoR1), (R1 * R1, infoM1), (R2 * R2, infoM4), (R3 * R3, infoM6), (R1 * R2 * 2, infoM2), (R1 * R3 * 2, infoM3), (R2 * R3 * 2, infoM5));
            infobi.Add(m1);
            var m2 = ZipCoefficientAdd((1, infoR2), (R4 * R1, infoM1), (R4 * R2, infoM2), (R4 * R3, infoM3), (R5 * R1, infoM2), (R5 * R2, infoM4), (R5 * R3, infoM5),
                (R6 * R1, infoM3), (R6 * R2, infoM5), (R6 * R3, infoM6));
            //infobi(end + 1, :) = infoR2 + R4.* (R1.* infoM1 + R2.* infoM2 + R3.* infoM3) + R5.* (R1.* infoM2 + R2.* infoM4 + R3.* infoM5) + R6.* (R1.* infoM3 + R2.* infoM5 + R3.* infoM6);
            var m3 = ZipCoefficientAdd((1, infoR3), (R7 * R1, infoM1), (R7 * R2, infoM2), (R7 * R3, infoM3), (R8 * R1, infoM2), (R8 * R2, infoM4), (R8 * R3, infoM5), (R9 * R1, infoM3), (R9 * R2, infoM5),
                (R9 * R3, infoM6));
            //infobi(end + 1, :) = infoR3 + R7.* (R1.* infoM1 + R2.* infoM2 + R3.* infoM3) + R8.* (R1.* infoM2 + R2.* infoM4 + R3.* infoM5) + R9.* (R1.* infoM3 + R2.* infoM5 + R3.* infoM6);
            var m4 = ZipCoefficientAdd((1, infoR4), (R4 * R4, infoM1), (R5 * R5, infoM4), (R6 * R6, infoM6), (R4 * R5 * 2, infoM2), (R4 * R6 * 2, infoM3), (R5 * R6 * 2, infoM5));
            //infobi(end + 1, :) = infoR4 + (R4.* R4).* infoM1 + (R5.* R5).* infoM4 + (R6.* R6).* infoM6 + R4.* R5.* infoM2 * 2.0 + R4.* R6.* infoM3 * 2.0 + R5.* R6.* infoM5 * 2.0;
            var m5 = ZipCoefficientAdd((1, infoR5), (R7 * R4, infoM1), (R7 * R5, infoM2), (R7 * R6, infoM3), (R8 * R4, infoM2), (R8 * R5, infoM4), (R8 * R6, infoM5), (R9 * R4, infoM3),
                (R9 * R5, infoM5), (R9 * R6, infoM6));
            //infobi(end + 1, :) = infoR5 + R7.* (R4.* infoM1 + R5.* infoM2 + R6.* infoM3) + R8.* (R4.* infoM2 + R5.* infoM4 + R6.* infoM5) + R9.* (R4.* infoM3 + R5.* infoM5 + R6.* infoM6);
            var m6 = ZipCoefficientAdd((1, infoR6), (R7 * R7, infoM1), (R8 * R8, infoM4), (R9 * R9, infoM6), (R7 * R8 * 2, infoM2), (R7 * R9 * 2, infoM3), (R8 * R9 * 2, infoM5));
            //infobi(end + 1, :) = infoR6 + (R7.* R7).* infoM1 + (R8.* R8).* infoM4 + (R9.* R9).* infoM6 + R7.* R8.* infoM2 * 2.0 + R7.* R9.* infoM3 * 2.0 + R8.* R9.* infoM5 * 2.0;

            // MArray = sum(infobi, 2);
            /*
                        m1 = infobi(1,:);
                        m2 = infobi(2,:);
                        m3 = infobi(3,:);
                        m4 = infobi(4,:);
                        m5 = infobi(5,:);
                        m6 = infobi(6,:);
            */
            var w = new double[v1.Count];
            for (int i = 0; i < w.Length; i++)
            {
                w[i] = -(v1[i] * (m1[i] * v1[i] + m2[i] * v2[i] + m3[i] * v3[i])
                    + v2[i] * (m2[i] * v1[i] + m4[i] * v2[i] + m5[i] * v3[i])
                    + v3[i] * (m3[i] * v1[i] + m5[i] * v2[i] + m6[i] * v3[i]));
                // w = -exp(-(v1.* (m1.* v1 + m2.* v2 + m3.* v3) + v2.* (m2.* v1 + m4.* v2 + m5.* v3) + v3.* (m3.* v1 + m5.* v2 + m6.* v3)) / (2.* sigma ^ 2)) / sigma ^ 2;
                w[i] = -Math.Exp(w[i] / (2 * sigma * sigma));
                w[i] /= sigma * sigma;
            }
            //w = w.Select(wi => -Math.Exp(wi / (2 * sigma * sigma)) / sigma * sigma).ToArray();
            //var sez = ZipMultiplyAdd((v1.Select(x => -x).ToList(), ZipMultiplyAdd((m1, v1), (m2, v2), (m3, v3))), (v2, ZipMultiplyAdd((m2, v1), (m4, v2), (m5, v3))),
            //    (v3, ZipMultiplyAdd((m3, v1), (m5, v2), (m6, v3))));
            //var w = ZipMultiplyAdd((v1.Select(x => -x).ToList(), ZipMultiplyAdd((m1, v1), (m2, v2), (m3, v3))), (v2, ZipMultiplyAdd((m2, v1), (m4, v2), (m5, v3))),
            //    (v3, ZipMultiplyAdd((m3, v1), (m5, v2), (m6, v3)))).Select(y => -Math.Exp(y / (2 * sigma * sigma)) / sigma * sigma).ToList();
            // w = -exp(-(v1.* (m1.* v1 + m2.* v2 + m3.* v3) + v2.* (m2.* v1 + m4.* v2 + m5.* v3) + v3.* (m3.* v1 + m5.* v2 + m6.* v3)) / (2.* sigma ^ 2)) / sigma ^ 2;
            H = new double[6, 6];
            for (int i = 0; i < w.Length; i++)
            {
                H[0, 0] += w[i] * (m4[i] * (p3[i] * p3[i]) + m6[i] * (p2[i] * p2[i]) - m5[i] * p2[i] * p3[i] * 2.0);
                H[1, 0] += -p3[i] * (m2[i] * p3[i] * w[i] - m3[i] * p2[i] * w[i]) + p1[i] * (m5[i] * p3[i] * w[i] - m6[i] * p2[i] * w[i]);
                H[2, 0] += p2[i] * (m2[i] * p3[i] * w[i] - m3[i] * p2[i] * w[i]) - p1[i] * (m4[i] * p3[i] * w[i] - m5[i] * p2[i] * w[i]);
                H[3, 0] += -w[i] * (m2[i] * p3[i] - m3[i] * p2[i]);
                H[4, 0] += -w[i] * (m4[i] * p3[i] - m5[i] * p2[i]);
                H[5, 0] += -w[i] * (m5[i] * p3[i] - m6[i] * p2[i]);
                H[0, 1] += -p3[i] * (m2[i] * p3[i] * w[i] - m5[i] * p1[i] * w[i]) + p2[i] * (m3[i] * p3[i] * w[i] - m6[i] * p1[i] * w[i]);
                H[1, 1] += w[i] * (m1[i] * (p3[i] * p3[i]) + m6[i] * (p1[i] * p1[i]) - m3[i] * p1[i] * p3[i] * 2.0);
                H[2, 1] += -p2[i] * (m1[i] * p3[i] * w[i] - m3[i] * p1[i] * w[i]) + p1[i] * (m2[i] * p3[i] * w[i] - m5[i] * p1[i] * w[i]);
                H[3, 1] += w[i] * (m1[i] * p3[i] - m3[i] * p1[i]);
                H[4, 1] += w[i] * (m2[i] * p3[i] - m5[i] * p1[i]);
                H[5, 1] += w[i] * (m3[i] * p3[i] - m6[i] * p1[i]);
                H[0, 2] += p3[i] * (m2[i] * p2[i] * w[i] - m4[i] * p1[i] * w[i]) - p2[i] * (m3[i] * p2[i] * w[i] - m5[i] * p1[i] * w[i]);
                H[1, 2] += -p3[i] * (m1[i] * p2[i] * w[i] - m2[i] * p1[i] * w[i]) + p1[i] * (m3[i] * p2[i] * w[i] - m5[i] * p1[i] * w[i]);
                H[2, 2] += w[i] * (m1[i] * (p2[i] * p2[i]) + m4[i] * (p1[i] * p1[i]) - m2[i] * p1[i] * p2[i] * 2.0);
                H[3, 2] += -w[i] * (m1[i] * p2[i] - m2[i] * p1[i]);
                H[4, 2] += -w[i] * (m2[i] * p2[i] - m4[i] * p1[i]);
                H[5, 2] += -w[i] * (m3[i] * p2[i] - m5[i] * p1[i]);
                H[0, 3] += -w[i] * (m2[i] * p3[i] - m3[i] * p2[i]);
                H[1, 3] += w[i] * (m1[i] * p3[i] - m3[i] * p1[i]);
                H[2, 3] += -w[i] * (m1[i] * p2[i] - m2[i] * p1[i]);
                H[3, 3] += m1[i] * w[i];
                H[4, 3] += m2[i] * w[i];
                H[5, 3] += m3[i] * w[i];
                H[0, 4] += -w[i] * (m4[i] * p3[i] - m5[i] * p2[i]);
                H[1, 4] += w[i] * (m2[i] * p3[i] - m5[i] * p1[i]);
                H[2, 4] += -w[i] * (m2[i] * p2[i] - m4[i] * p1[i]);
                H[3, 4] += m2[i] * w[i];
                H[4, 4] += m4[i] * w[i];
                H[5, 4] += m5[i] * w[i];
                H[0, 5] += -w[i] * (m5[i] * p3[i] - m6[i] * p2[i]);
                H[1, 5] += w[i] * (m3[i] * p3[i] - m6[i] * p1[i]);
                H[2, 5] += -w[i] * (m3[i] * p2[i] - m5[i] * p1[i]);
                H[3, 5] += m3[i] * w[i];
                H[4, 5] += m5[i] * w[i];
                H[5, 5] += m6[i] * w[i];
                //H = reshape(sum(H, 2), 6, 6);
            }
            b = new double[6];
            for (int i = 0; i < w.Length; i++)
            {
                b[0] += -v1[i] * (m2[i] * p3[i] * w[i] - m3[i] * p2[i] * w[i]) - v2[i] * (m4[i] * p3[i] * w[i] - m5[i] * p2[i] * w[i]) - v3[i] * (m5[i] * p3[i] * w[i] - m6[i] * p2[i] * w[i]);
                b[1] += v1[i] * (m1[i] * p3[i] * w[i] - m3[i] * p1[i] * w[i]) + v2[i] * (m2[i] * p3[i] * w[i] - m5[i] * p1[i] * w[i]) + v3[i] * (m3[i] * p3[i] * w[i] - m6[i] * p1[i] * w[i]);
                b[2] += -v1[i] * (m1[i] * p2[i] * w[i] - m2[i] * p1[i] * w[i]) - v2[i] * (m2[i] * p2[i] * w[i] - m4[i] * p1[i] * w[i]) - v3[i] * (m3[i] * p2[i] * w[i] - m5[i] * p1[i] * w[i]);
                b[3] += w[i] * (m1[i] * v1[i] + m2[i] * v2[i] + m3[i] * v3[i]);
                b[4] += w[i] * (m2[i] * v1[i] + m4[i] * v2[i] + m5[i] * v3[i]);
                b[5] += w[i] * (m3[i] * v1[i] + m5[i] * v2[i] + m6[i] * v3[i]);
                //b += sum(b, 2);
            }
            J = 0.0;
            for (int i = 0; i < w.Length; i++)
            {
                J += m1[i] * (v1[i] * v1[i]) + m4[i] * (v2[i] * v2[i]) + m6[i] * (v3[i] * v3[i]) + m2[i] * v1[i] * v2[i] * 2.0
                    + m3[i] * v1[i] * v3[i] * 2.0 + m5[i] * v2[i] * v3[i] * 2.0;
            }
            var maxValue = b.Max(x => Math.Abs(x));
            for (int i = 0; i < 6; i++)
                for (int j = 0; j < 6; j++)
                    if (maxValue < Math.Abs(H[i, j]))
                        maxValue = Math.Abs(H[i, j]);
            if (maxValue < 0.1e-6)
            {
                var k = 1 / maxValue;
                b = b.multiply(k);
                H = H.multiply(k);
            }
        }
        private static double[] ZipCoefficientAdd(params (double, IList<double>)[] terms)
        {
            var length = terms[0].Item2.Count;
            var result = new double[length];

            for (int i = 0; i < result.Length; i++)
            {
                var sum = 0.0;
                foreach (var term in terms)
                    sum += term.Item1 * term.Item2[i];
                result[i] = sum;
            }
            return result;
        }
        private static double[] ZipMultiplyAdd(params (IList<double>, IList<double>)[] terms)
        {
            var length = terms[0].Item2.Count;
            var result = new double[length];

            for (int i = 0; i < result.Length; i++)
            {
                var sum = 0.0;
                foreach (var term in terms)
                    sum += term.Item1[i] * term.Item2[i];
                result[i] = sum;
            }
            return result;
        }
        private static List<double> GetAngle(IList<Vector3> norm_Ref, IList<Vector3> norm_Mov,
            List<(Vector3, int)> nNidx, Matrix4x4 r)
        {
            var angles = new List<double>(norm_Ref.Count);
            for (int i = 0; i < norm_Ref.Count; i++)
                angles.Add(norm_Ref[nNidx[i].Item2].Dot(norm_Mov[i].Transform(r)));
            return angles;
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

        private static double[,] SkewFun(double[] a)
        {
            if (a.Length == 3)
                return new double[,]  {
                { 0, -a[2], a[1]},
                { a[2],  0, -a[0]},
                               { -a[1], a[0], 0}
            };
            if (a.Length == 2)
                return new double[,]  {
                { a[1]},
                { -a[0]}
            };
            else throw new ArgumentException("a must be 2 or 3 dimensional");
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

        private static Matrix4x4 GetTranslationMatrix(IEnumerable<Vector3> targetPoints, IEnumerable<Vector3> startingPoints)
        {
            return Matrix4x4.CreateTranslation(GetTranslationVector(targetPoints, startingPoints));
        }

        private static Vector3 GetTranslationVector(IEnumerable<Vector3> targetPoints, IEnumerable<Vector3> startingPoints)
        {
            var targetCom = GetCenterOfMassPoints(targetPoints);
            var startingCom = GetCenterOfMassPoints(startingPoints);
            return targetCom - startingCom;
        }

        private static Vector3 GetCenterOfMassPoints(IEnumerable<Vector3> points)
        {
            var numPoints = 0;
            var center = Vector3.Zero;
            foreach (var p in points)
            {
                center += p;
                numPoints++;
            }
            center /= numPoints;
            return center;
        }
    }
}
