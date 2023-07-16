using System;
using System.Collections.Generic;
using System.Linq;
using TVGL;
using System.Text;
using System.Threading.Tasks;
using SharpDX.Mathematics.Interop;

var RefData = readCloudCsv("Hokuyo_0.csv", 0.1);
var MovData = readCloudCsv("Hokuyo_1.csv", 0.1);



CoBigICP_Sharp.CoBigICP.Run(MovData, RefData);

//function[rawPoints] = readCloudCsv(filename, subSampleRatio)
//% READCLOUDCSV read csv files of ETH data set
//% INPUTS:
//% -filename: CSV file name of the point cloud to be ploted. It can be
//%                   local point cloud (ex.: Hokuyo_0.csv) or global map (ex.:
//% PointCloud0.csv).
//% -subSampleRatio: Ratio of point to be removed (used for performance
//%                               reasons)
static List<Vector3> readCloudCsv(string fileName, double subSampleRatio)
{
    //% Load data
    //cloud = importdata(filename);
    using var reader = new StreamReader(FindInParentDirectory(fileName).FullName);
    var rawPoints = new List<Vector3>();
    var line = reader.ReadLine();
    while (line != null)
    {
        var values = line.Split(',');
        //% Extract coordinates
        //x = cloud.data(:, strcmp('x', cloud.colheaders));
        //y = cloud.data(:, strcmp('y', cloud.colheaders));
        //z = cloud.data(:, strcmp('z', cloud.colheaders));
        if (double.TryParse(values[1], out var x) &&
            double.TryParse(values[2], out var y) &&
            double.TryParse(values[3], out var z))
        {
            rawPoints.Add(new Vector3(x, y, z));
        }
        line = reader.ReadLine();
    }
    var numPoints = rawPoints.Count;
    var numToKeep = 1000;
    rawPoints.RemoveRange(numToKeep, numPoints - numToKeep);
    return rawPoints;


    var numToRemove = (int)(numPoints * (1 - subSampleRatio));
    var random = new Random();
    //% Subsample
    //f = rand(1, length(x));
    //x = x(f > subSampleRatio);
    //y = y(f > subSampleRatio);
    //z = z(f > subSampleRatio);
    for (int i = 0; i < numToRemove; i++)
    {
        var randIndex = random.Next(rawPoints.Count);
        rawPoints.RemoveAt(randIndex);
    }
    //rawPoints = [x y z];
    return rawPoints;
}


static FileInfo FindInParentDirectory(string fileName)
{
    var dir = new DirectoryInfo(".");
    while (!File.Exists(Path.Combine(dir.FullName, fileName)))
    {
        if (dir == null) throw new FileNotFoundException("Target folder not found", fileName);
        dir = dir.Parent;
    }
    return new FileInfo(Path.Combine(dir.FullName, fileName));
}
