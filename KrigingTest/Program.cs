using KrigingAlgo;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace KrigingTest
{
    class Program
    {
        static void Main(string[] args)
        {
            string file = @"D:\Works\Kriging\Kriging\Data\ZoneA.csv";
            List<double> x = new List<double>();
            List<double> y = new List<double>();
            List<double> z = new List<double>();
            var lines =  File.ReadAllLines(file);
            for (int i = 1; i < lines.Length; i++)
            {
                string[] items = lines[i].Split(',');
                if (items.Length > 3)
                {
                    x.Add(Convert.ToDouble(items[0]));
                    y.Add(Convert.ToDouble(items[1]));
                    z.Add(Convert.ToDouble(items[3]));
                }
            }

            RasterContext rasterContext = new RasterContext(0, 0, 100, 100, 200, 160, false);
            Kriging kriging = new Kriging(x.ToArray(), y.ToArray(), z.ToArray(), rasterContext);

            kriging.Initialize();
            float[] buffer = new float[200*160];
            var timer = System.Diagnostics.Stopwatch.StartNew();
            kriging.OrdinaryKrige(Kriging.Model.Spherical, buffer, 0.0, kriging.GetEstimatedSill(), 4000, 16, 16, 0);
            Console.WriteLine(timer.ElapsedMilliseconds);
            timer.Restart();
            Kriging2 kriging2 = new Kriging2(x.ToArray(), y.ToArray(), z.ToArray(), rasterContext);
            kriging2.SimpleKrige(Kriging2.Model.Spherical, buffer, 0.0, kriging.GetEstimatedSill(), 4000);
            Console.WriteLine(timer.ElapsedMilliseconds);
            Console.ReadKey();
            Console.WriteLine("");
        }
    }
}
