using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using LineFinder;

namespace Test
{
    class Program
    {
        static void Main(string[] args)
        {
            double[,] data = new double[,] {
                {0, 0, 1},
                {0, 0, 1},
                {0, 1, 0}
            };

            var sample = Class1.SampleAt(data.Matrix().T(), 1, 1);
            var peak = Class1.FindPeak(sample.ToArray());
        }
    }
}
