// See https://aka.ms/new-console-template for more information

using om38to13;
using System.Globalization;

Thread.CurrentThread.CurrentCulture = new CultureInfo("en-US");

// Root directory of the project - inside are data in a directory named "data".
var root = Path.GetFullPath(Path.Combine(AppContext.BaseDirectory, @"..\..\..\..\"));

Tool.TestPossition(root, 1, 143310164, 143361930, DataType.AllData);
