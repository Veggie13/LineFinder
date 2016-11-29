using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;

namespace LineFinder
{
    public static class Class1
    {
        public static List<PointF> GetPeaks(Bitmap bm)
        {
            List<PointF> results = new List<PointF>();

            double[,] img = GetBrightness(bm);

            double[] sample;
            for (int x = 1; x < bm.Width - 1; x++)
            {
                for (int y = 1; y < bm.Height - 1; y++)
                {
                    //sample = SampleAt(img, x, y).ToArray();
                }
            }
            return null;
        }

        static double[,] transform = MMult(0.25, Matrix(new double[,] {
            {0, 0, 0, 0, 4, 0, 0, 0, 0},
            {0, 0, 0, -2, 0, 2, 0, 0, 0},
            {0, 0, 0, 2, -4, 2, 0, 0, 0},
            {0, -2, 0, 0, 0, 0, 0, 2, 0},
            {1, 0, -1, 0, 0, 0, -1, 0, 1},
            {-1, 2, -1, 0, 0, 0, 1, -2, 1},
            {0, 2, 0, 0, -4, 0, 0, 2, 0},
            {-1, 0, 1, 2, 0, -2, -1, 0, 1},
            {1, -2, 1, -2, 4, -2, 1, -2, 1}
        }));
        static IMatrix trans = T(Matrix(transform));

        static double[,] transform2 = new double[,] {
                {1, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 1, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 1, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0},
                {0, 0, 0, 0, 0, 0, 0, 0, 0}
            };
        static IMatrix trans2 = T(Matrix(transform2));

        public class Peak
        {
            public PointF Loc;
            public double Val;
        }

        public static Peak FindPeak(double[] sample)
        {
            var samp = Matrix(sample);
            var coeffs = MMult(trans, samp);

            if (coeffs[0, 1] != 0)
            {
                double Z = coeffs[0, 3] / coeffs[0, 1];
                trans2[3, 1] = trans2[4, 2] = trans2[5, 3] = Z;
                trans2[6, 2] = trans2[7, 3] = trans2[8, 4] = Z * Z;

                var coeffs2 = MMult(Matrix(MMult(trans2, trans)), samp);

                var roots = GetRealCubicRoots(4 * coeffs2[0, 4], 3 * coeffs2[0, 3], 2 * coeffs2[0, 2], coeffs2[0, 1]);
                var peaks = roots.Select(r => MMult(Matrix(new double[,] { { 1, r, r * r, r * r * r, r * r * r * r } }).T(), Matrix(coeffs2))[0, 0]).ToList();

                double t = roots[peaks.IndexOf(peaks.Max())];
                return new Peak()
                { 
                    Loc = new PointF((float)t, (float)(Z * t)), 
                    Val = peaks.Max()
                };
            }
            else if (coeffs[0, 6] < 0)
            {
                // z = a0 + a3*y + a6*y^2
                // z' = a3 + 2*a6*y = 0
                // y = -a3/(2*a6)
                // z = a0 - a3^2/(4*a6)

                return new Peak()
                {
                    Loc = new PointF(0f, (float)(-coeffs[0, 3] / (2 * coeffs[0, 6]))),
                    Val = coeffs[0, 0] - coeffs[0, 3] * coeffs[0, 3] / (4 * coeffs[0, 6])
                };
            }
            else if (coeffs[0, 6] == 0)
            {
                return new Peak()
                {
                    Loc = new PointF(0f, Math.Sign(coeffs[0, 3])),
                    Val = coeffs[0, 0] + Math.Abs(coeffs[0, 3])
                };
            }
            else
            {
                return new Peak()
                {
                    Loc = new PointF(0f, -Math.Sign(coeffs[0, 3])),
                    Val = coeffs[0, 0] - Math.Abs(coeffs[0, 3]) + coeffs[0, 6]
                };
            }
        }

        public static IEnumerable<double> SampleAt(IMatrix img, int x, int y)
        {
            yield return img[x - 1, y - 1];
            yield return img[x, y - 1];
            yield return img[x + 1, y - 1];
            yield return img[x - 1, y];
            yield return img[x, y];
            yield return img[x + 1, y];
            yield return img[x - 1, y + 1];
            yield return img[x, y + 1];
            yield return img[x + 1, y + 1];
        }

        static double[,] GetBrightness(Bitmap bm)
        {
            double[,] result = new double[bm.Width, bm.Height];
            for (int x = 0; x < bm.Width; x++)
            {
                for (int y = 0; y < bm.Height; y++)
                {
                    result[x, y] = bm.GetPixel(x, y).GetBrightness();
                }
            }

            return result;
        }

        public interface IMatrix
        {
            double this[int x, int y] { get; set; }
            int Width { get; }
            int Height { get; }
        }

        class ArrayMatrix : IMatrix
        {
            private double[,] _arr;
            public ArrayMatrix(double[,] arr)
            {
                _arr = arr;
            }

            public double this[int x, int y]
            {
                get
                {
                    return _arr[x, y];
                }
                set
                {
                    _arr[x, y] = value;
                }
            }

            public int Width
            {
                get { return _arr.GetLength(0); }
            }

            public int Height
            {
                get { return _arr.GetLength(1); }
            }
        }

        class VectorMatrix : IMatrix
        {
            private double[] _vec;
            public VectorMatrix(double[] vec)
            {
                _vec = vec;
            }

            public double this[int x, int y]
            {
                get
                {
                    if (x != 0)
                        throw new ArgumentOutOfRangeException("x");
                    return _vec[y];
                }
                set
                {
                    if (x != 0)
                        throw new ArgumentOutOfRangeException("x");
                    _vec[y] = value;
                }
            }

            public int Width
            {
                get { return 1; }
            }

            public int Height
            {
                get { return _vec.Length; }
            }
        }

        class Transpose : IMatrix
        {
            private IMatrix _mat;
            public Transpose(IMatrix mat)
            {
                _mat = mat;
            }

            public double this[int x, int y]
            {
                get
                {
                    return _mat[y, x];
                }
                set
                {
                    _mat[y, x] = value;
                }
            }

            public int Width
            {
                get { return _mat.Height; }
            }

            public int Height
            {
                get { return _mat.Width; }
            }
        }

        public static IMatrix Matrix(this double[,] mat)
        {
            return new ArrayMatrix(mat);
        }

        public static IMatrix Matrix(this double[] vec)
        {
            return new VectorMatrix(vec);
        }

        public static IMatrix T(this IMatrix mat)
        {
            return new Transpose(mat);
        }

        static double[,] MMult(IMatrix a, IMatrix b)
        {
            if (a.Width != b.Height)
            {
                throw new ArgumentException("Array sizes won't multiply.");
            }

            double[,] result = new double[b.Width, a.Height];

            for (int x = 0; x < result.GetLength(0); x++)
                for (int y = 0; y < result.GetLength(1); y++)
                {
                    for (int i = 0; i < a.Width; i++)
                    {
                        result[x, y] += a[i, y] * b[x, i];
                    }
                }

            return result;
        }

        static double[,] MMult(double v, IMatrix a)
        {
            double[,] result = new double[a.Width, a.Height];

            for (int x = 0; x < result.GetLength(0); x++)
                for (int y = 0; y < result.GetLength(1); y++)
                {
                    result[x, y] = v * a[x, y];
                }

            return result;
        }

        static double[] GetRealCubicRoots(double a, double b, double c, double d)
        {
            double q = (3 * a * c - b * b) / (9 * a * a);
            double r = (9 * a * b * c - 27 * a * a * d - 2 * b * b * b) / (54 * a * a * a);

            double disc = q * q * q + r * r;
            double s, t;

            if (disc > 0)
            {
                s = Math.Pow(r + Math.Sqrt(disc), 1.0 / 3) + Math.Pow(r - Math.Sqrt(disc), 1.0 / 3);
                t = 0;
            }
            else
            {
                double rho = Math.Sqrt(-q * q * q);
                double theta = Math.Acos(r / rho) / 3;
                s = Math.Sqrt(-q) * Math.Cos(theta);
                t = Math.Sqrt(-q) * Math.Sin(theta);
            }

            double r1, r2, r3;
            r1 = 2 * s - b / (3 * a);
            if (t == 0)
                return new double[] { r1 };

            r2 = -s - b / (3 * a) + Math.Sqrt(3) * t;
            r3 = -s - b / (3 * a) - Math.Sqrt(3) * t;
            return new double[] { r1, r2, r3 };
        }
    }
}
