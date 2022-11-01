// <copyright file="Tool.cs" company="Marek Behalek, marek.behalek@vsb.cz ,VSB-TU Ostrava">
// Copyright (c) Marek Behalek, marek.behalek@vsb.cz ,VSB-TU Ostrava. All rights reserved.
// </copyright>

namespace Experiments
{
    using System.Text;

    /// <summary>
    /// Tool filter.
    /// </summary>
    public class Tool
    {
        private const double Epsilon = 10000;

        private static string[] messages = new string[]
        {
            "(Reversed)", // just note
            "no mapping from HG38 to CHM13", // Empty
            "new possition in CHM13", // Simple
            "alternative mapping in CHM13", // Alternative
            "in CHM13 have no mapping to HG38", // Empty backward
            "in CHM13 have multiple sources in HG38",
        };

        private enum IntervalType
        {
            Simple,
            Alternative,
            Empty,
        }
        /// <summary>
        /// Start the computing.
        /// </summary>
        public static void DoWork()
        {
            // var smap = LoadSmap(@"d:\OpticData\tool\test1.smap");
            // var smap = LoadSmap(@"d:\OpticData\tool\test2.smap");
            var smap = LoadSmap(@"d:\OpticData\tool\exp_refineFinal1_merged_filter_inversions.smap");
            var dataFromExperimentsFiles = new (string FromHGtoCHM, string FromCHMToHG)[] {
                (@"d:\OpticData\tool\fromHG38toCHM13", @$"d:\OpticData\tool\fromCHM13toHG38"),
                (@"d:\OpticData\tool\fromHG38toCHM13-alignments", @$"d:\OpticData\tool\fromCHM13toHG38-alignments"),
                (@"d:\OpticData\tool\fromHG38toCHM13-assemblies", @$"d:\OpticData\tool\fromCHM13toHG38-assemblies"),
            };
            var dataFromExperimens = new List<(List<(IntervalType Type, Position HG, Position? CHM)> All, List<Position> Empty, List<(Position CHM, Position HG1, Position HG2)> Alternatives)>();

            // selects from CHM -> HG mapping empty and alternative posiotions. No needed in final tool. use just created files.
            foreach (var (_, fileName) in dataFromExperimentsFiles)
            {
                if (!File.Exists(fileName + "-empty"))
                {
                    using var outputForEmpty = new StreamWriter(fileName + "-empty");
                    var chmToHg = LoadExperimentData(fileName);
                    var empty = chmToHg.Where(x => x.Type == IntervalType.Empty).Select(x => x.From);
                    foreach (var item in empty)
                    {
                        outputForEmpty.WriteLine(item);
                    }
                }
                if (!File.Exists(fileName + "-alternatives") || true)
                {
                    using var outputForEmpty = new StreamWriter(fileName + "-alternatives");
                    var chmToHg = LoadExperimentData(fileName);
                    var alternatives = chmToHg.Where(x => x.Type == IntervalType.Alternative);
                    foreach (var alternativeItem in alternatives)
                    {
                        var simple = chmToHg.Where(x => x.From.OverlapWith(alternativeItem.From) && x.Type == IntervalType.Simple).ToList();
                        if (simple.Count == 0) throw new Exception("Expecting just one alternative.");
                        foreach (var mainItem in simple)
                        {
                            var part = alternativeItem.From.IntersectionWith(mainItem.From);
                            var target1 = alternativeItem.To?.GetMapping(alternativeItem.From, part);
                            var target2 = mainItem.To?.GetMapping(mainItem.From, part);
                            outputForEmpty.WriteLine($"{part}\t{target1}\t{target2}");
                        }
                    }
                }
            }

            foreach (var item in dataFromExperimentsFiles)
            {
                var all = LoadExperimentData(item.FromHGtoCHM);
                var empty = LoadEmptyIntervals(item.FromCHMToHG + "-empty");
                var alternatives = LoadAlternativeIntervals(item.FromCHMToHG + "-alternatives");
                dataFromExperimens.Add((all, empty, alternatives));
            }

            var dataFromArticle = LoadArticleData(@"d:\OpticData\tool\prediction_38.bed");

            void Process(
                StreamWriter output,
                string indentation,
                Position position,
                List<(string Name, Position HG, Position CHM)> dataFromArticle,
                List<(IntervalType Type, Position HG, Position? CHM)> all,
                List<Position> empty,
                List<(Position CHM, Position HG1, Position HG2)> alternatives)
            {
                output.Write(ProcessArticleData(indentation + "\t", position, dataFromArticle));
                output.Write(ProcessExperimentData(indentation + "\t", position, all, empty, alternatives));
                output.WriteLine();
            }

            using var output = new StreamWriter(@"d:\OpticData\tool\result");
            foreach (var item in smap.Common)
            {
                output.WriteLine($"{item.Id} {item.Position} {smap.Original[item.Id]}");
                Process(output, string.Empty, item.Position, dataFromArticle, dataFromExperimens[0].All, dataFromExperimens[0].Empty, dataFromExperimens[0].Alternatives);
            }
            foreach (var item in smap.Translocations)
            {
                output.WriteLine($"{item.Id} TRANSLOCATION A={item.PositionA.Chromosome}:{item.PositionA.End - Epsilon}\tB={item.PositionB.Chromosome}:{item.PositionB.End - Epsilon}");
                output.WriteLine($"\tA: {item.PositionA}");
                Process(output, "\t", item.PositionA, dataFromArticle, dataFromExperimens[0].All, dataFromExperimens[0].Empty, dataFromExperimens[0].Alternatives);
                output.WriteLine($"\tB: {item.PositionB}");
                Process(output, "\t", item.PositionB, dataFromArticle, dataFromExperimens[0].All, dataFromExperimens[0].Empty, dataFromExperimens[0].Alternatives);
                output.WriteLine();
            }
            foreach (var item in smap.Inversions)
            {
                output.WriteLine($"{item.IdA} {item.IdB} INVERSION {item.Position}");
                Process(output, string.Empty, item.Position, dataFromArticle, dataFromExperimens[0].All, dataFromExperimens[0].Empty, dataFromExperimens[0].Alternatives);
            }
        }

        private static List<(string Name, Position HG, Position CHM)> LoadArticleData(string fileName)
        {
            var result = new List<(string Name, Position HG, Position CHM)>();
            var input = new StreamReader(fileName);
            var line = input.ReadLine();
            while (line != null)
            {
                var splitedLine = line.Split("\t;".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                if (splitedLine.Length != 3) throw new Exception("Unexpected data format.");
                result.Add((splitedLine[0], Position.ReadFrom(splitedLine[1]), Position.ReadFrom(splitedLine[2])));
                line = input.ReadLine();
            }
            return result;
        }

        private static List<(IntervalType Type, Position From, Position? To)> LoadExperimentData(string fileName)
        {
            var result = new List<(IntervalType Type, Position From, Position? To)>();
            var input = new StreamReader(fileName);
            var line = input.ReadLine();
            while (line != null)
            {
                var splitedLine = line.Split("\t;".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                if (splitedLine.Length == 2 && splitedLine[0].StartsWith('E'))
                {
                    result.Add((IntervalType.Empty, Position.ReadFrom(splitedLine[1]), null));

                }
                else if (splitedLine.Length == 3)
                {
                    var type = splitedLine[0].StartsWith('S') ? IntervalType.Simple : IntervalType.Alternative;
                    result.Add((type, Position.ReadFrom(splitedLine[1]), Position.ReadFrom(splitedLine[2])));
                }
                else
                {
                    throw new Exception("Unexpected data format.");

                }
                line = input.ReadLine();
            }
            return result;
        }

        private static List<Position> LoadEmptyIntervals(string fileName)
        {
            var result = new List<Position>();
            var input = new StreamReader(fileName);
            var line = input.ReadLine();
            while (line != null)
            {
                result.Add(Position.ReadFrom(line));
                line = input.ReadLine();
            }
            return result;
        }

        private static List<(Position CHM, Position HG1, Position HG2)> LoadAlternativeIntervals(string fileName)
        {
            var result = new List<(Position CHM, Position HG1, Position HG2)>();
            var input = new StreamReader(fileName);
            var line = input.ReadLine();
            while (line != null)
            {
                var splitedLine = line.Split("\t;".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                result.Add((Position.ReadFrom(splitedLine[0]), Position.ReadFrom(splitedLine[1]), Position.ReadFrom(splitedLine[2])));
                line = input.ReadLine();
            }
            return result;
        }

        private static (
            Dictionary<int, string> Original,
            List<(int Id, Position Position)> Common,
            List<(int Id, Position PositionA, Position PositionB)> Translocations,
            List<(int IdA, int IdB, Position Position)> Inversions) LoadSmap(string fileName)
        {
            var original = new Dictionary<int, string>();
            var common = new List<(int Id, Position Position)>();
            var translocations = new List<(int Id, Position PositionA, Position PositionB)>();
            var inversionsPartial = new List<(int Id, int LinkId, int Chromosome, double X, double Y)>();
            var inversions = new List<(int IdA, int IdB, Position Position)>();

            var input = new StreamReader(fileName);
            var line = input.ReadLine();
            while (line != null)
            {
                if (string.IsNullOrEmpty(line)) break;
                if (!line.StartsWith("#"))
                {
                    var splitedLine = line.Split("\t".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                    original.Add(int.Parse(splitedLine[0]), line);
                    if (int.Parse(splitedLine[12]) != -1)
                    {
                        // int a = int.Parse(splitedLine[0]);
                        // int b = int.Parse(splitedLine[12]);
                        inversionsPartial.Add((int.Parse(splitedLine[0]), int.Parse(splitedLine[12]), int.Parse(splitedLine[2]), double.Parse(splitedLine[6]), double.Parse(splitedLine[7])));

                    } else if (line.Contains("trans"))
                    {
                        var positionA = double.Parse(splitedLine[6]);
                        var positionB = double.Parse(splitedLine[7]);
                        var positionAStart = (positionA - Epsilon < 0) ? 0 : positionA - Epsilon;
                        var positionBStart = (positionB - Epsilon < 0) ? 0 : positionB - Epsilon;
                        var aroundA = new Position(int.Parse(splitedLine[2]), positionAStart, positionA + Epsilon);
                        var aroundB = new Position(int.Parse(splitedLine[3]), positionBStart, positionB + Epsilon);
                        translocations.Add((int.Parse(splitedLine[0]), aroundA, aroundB));
                    }
                    else
                    {
                        common.Add((int.Parse(splitedLine[0]), new Position(int.Parse(splitedLine[2]), double.Parse(splitedLine[6]), double.Parse(splitedLine[7]))));
                    }
                }

                line = input.ReadLine();
            }

            foreach (var itemA in inversionsPartial)
            {
                if (itemA.Y != -1 && itemA.X > itemA.Y) throw new Exception("We assume that x is always smaller than y.");
                if (itemA.X < 0) throw new Exception("Only the second one can be -1");

                if (!inversions.Any(x => x.IdA == itemA.Id || x.IdB == itemA.Id))
                {
                    var itemB = inversionsPartial.First(x => x.Id == itemA.LinkId);
                    if (itemA.Chromosome != itemB.Chromosome) throw new Exception("Inversion is on the same chromosome.");
                    if (itemA.Id != itemB.LinkId) throw new Exception("Wrong link id in data.");

                    if (itemA.Y == -1)
                    {
                        if (itemA.X < itemB.X) inversions.Add((itemA.Id, itemB.Id, new Position(itemA.Chromosome, itemA.X, itemB.X)));
                        else inversions.Add((itemA.Id, itemB.Id, new Position(itemA.Chromosome, itemB.Y, itemA.X)));
                    }
                    else if (itemB.Y == -1)
                    {
                        if (itemB.X < itemA.X) inversions.Add((itemA.Id, itemB.Id, new Position(itemA.Chromosome, itemB.X, itemA.X)));
                        else inversions.Add((itemA.Id, itemB.Id, new Position(itemA.Chromosome, itemA.Y, itemB.X)));
                    }
                    else
                    {
                        if (itemA.X < itemB.X) inversions.Add((itemA.Id, itemB.Id, new Position(itemA.Chromosome, itemA.Y, itemB.X)));
                        else inversions.Add((itemA.Id, itemB.Id, new Position(itemA.Chromosome, itemB.Y, itemA.X)));
                    }
                }
            }


            return (original, common, translocations, inversions);
        }

        private static string ProcessArticleData(string indentation, Position position, List<(string Name, Position HG, Position CHM)> dataFromArticle)
        {
            var result = new StringBuilder();
            foreach (var item in dataFromArticle)
            {
                if (item.HG.OverlapWith(position))
                {
                    result.AppendLine($"{indentation}{item.Name} {item.HG} {item.CHM}");
                }
            }
            return result.ToString();
        }

        private static string ProcessExperimentData(string indentation, Position position, List<(IntervalType Type, Position HG, Position? CHM)> all, List<Position> empty, List<(Position CHM, Position HG1, Position HG2)> alternatives)
        {
            var result = new StringBuilder();
            bool atLeastOne = false;
            foreach (var item in all.Where(x => x.HG.OverlapWith(position)))
            {
                // result.AppendLine($"{indentation}{item.Type} {item.HG} {item.CHM}");
                var intersection = position.IntersectionWith(item.HG);
                if (item.Type == IntervalType.Simple || item.Type == IntervalType.Alternative)
                {
                    var target = item.CHM?.GetMapping(item.HG, intersection);
                    var message = (item.Type == IntervalType.Simple) ? 2 : 3;
                    result.AppendLine($"{indentation}{intersection} {messages[message]} {target}{((item.CHM is not null && item.CHM.IsReversed) ? messages[0] : string.Empty)}");
                    if (target is not null)
                    {
                        if (target.Start > target.End) target = new Position(target.Chromosome, target.End, target.Start);
                        foreach (var emptyInterval in empty.Where(x => x.OverlapWith(target)))
                        {
                            result.AppendLine($"{indentation + "\t"}{emptyInterval.IntersectionWith(target)} {messages[4]}");
                        }
                        foreach (var alternativeInterval in alternatives.Where(x => x.CHM.OverlapWith(target)))
                        {
                            var part = alternativeInterval.CHM.IntersectionWith(target);
                            result.AppendLine($"{indentation + "\t"}{part} {messages[5]} {alternativeInterval.HG1.GetMapping(alternativeInterval.CHM, part)} {alternativeInterval.HG2.GetMapping(alternativeInterval.CHM, part)}");
                        }
                    }
                }
                else if (item.Type == IntervalType.Empty)
                {
                    result.AppendLine($"{indentation}{intersection} {messages[1]}");
                }

                atLeastOne = true;
            }

            if (!atLeastOne)
            {
                Console.WriteLine(position);
            }
            /*
            foreach (var item in complicated)
            {
                if (item.HG.OverlapWith(position))
                {
                    result.AppendLine($"{indentation}{item.HG} {item.CHM}");
                }
            }*/
            return result.ToString();
        }

        private class Position
        {
            public Position(int chromosome, double start, double end)
            {
                if (start < 0 || end < 0) throw new Exception("Unfiltered -1 from SMAP");
                this.Chromosome = chromosome;
                this.Start = start;
                this.End = end;
            }

            public double End { get; private set; }
            public int Chromosome { get; private set; }
            public bool IsReversed => this.Start > this.End;

            public double Start { get; private set; }
            public static Position ReadFrom(string input)
            {
                var elements = input.Split(":-".ToCharArray(), StringSplitOptions.RemoveEmptyEntries);
                if (elements.Length != 3) throw new Exception("Unexpected position format.");
                return new Position(int.Parse(elements[0]), double.Parse(elements[1]), double.Parse(elements[2]));
            }
            public Position GetMapping(Position whole, Position part)
            {
                double Get(double percent)
                {
                    return Math.Round(this.Start + ((this.End - this.Start) * percent));
                }
                double GetReversed(double percent)
                {
                    return Math.Round(this.Start - ((this.Start - this.End) * percent));
                }

                if (whole.Start > whole.End || part.Start > part.End) throw new Exception("Sorce interval should be from left to right.");
                // just a point
                if (this.Start == this.End || whole.Start == whole.End)
                {
                    return this;
                }

                var startPercent = (part.Start - whole.Start) / (whole.End - whole.Start);
                var endPercent = (part.End - whole.Start) / (whole.End - whole.Start);
                if (this.Start < this.End)
                {
                    return new Position(this.Chromosome, Get(startPercent), Get(endPercent));
                }
                else
                {
                    return new Position(this.Chromosome, GetReversed(startPercent), GetReversed(endPercent));
                }
            }

            public Position IntersectionWith(Position position)
            {
                if (this.Start > this.End || position.Start > position.End) throw new Exception("Reverse order is special case");
                if (this.Chromosome != position.Chromosome) throw new Exception("Intersection only valid for the same chromosome.");
                if (!this.OverlapWith(position)) throw new Exception("No ontersection.");
                return new Position(this.Chromosome, Math.Max(this.Start, position.Start), Math.Min(this.End, position.End));
            }

            public bool OverlapWith(Position position)
            {
                if (this.Start > this.End || position.Start > position.End) throw new Exception("Reverse order is special case");
                if (this.Chromosome != position.Chromosome) return false;
                return this.Start <= position.End && position.Start <= this.End;
            }
            public override string ToString()
            {
                return string.Format("{0}:{1}-{2}", this.Chromosome, this.Start, this.End);
            }
        }
    }
}
