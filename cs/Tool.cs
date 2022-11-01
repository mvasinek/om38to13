// <copyright file="Tool.cs" company="Marek Behalek, marek.behalek@vsb.cz ,VSB-TU Ostrava">
// Copyright (c) Marek Behalek, marek.behalek@vsb.cz ,VSB-TU Ostrava. All rights reserved.
// </copyright>

namespace om38to13
{
    using System.Collections.Generic;

    /// <summary>
    /// Tool om38to13 transforming position from GRCH38 to T2T-CHM13.
    /// </summary>
    public class Tool
    {
        private static string[] messages = new string[]
        {
            "(Reversed)", // just note
            "no mapping from HG38 to CHM13", // Empty
            "position in CHM13", // Simple
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

        /// <summary>Start the computing.</summary>
        /// <param name="dataDirectory">Path to directory data from the source github project.</param>
        /// <param name="chromosome">Chormosome number, use 23 for chromosome X, and 24 for chromosome y.</param>
        /// <param name="start">Start position.</param>
        /// <param name="end">End position.</param>
        /// <param name="dataType"></param>
        public static void TestPossition(string dataDirectory, int chromosome, double start, double end, DataType dataType)
        {
            var dataFromExperimentsFiles = new (string FromHGtoCHM, string FromCHMToHG)[] {
                (dataDirectory + @"data\fromHG38toCHM13", dataDirectory + @"data\fromCHM13toHG38"),
                (dataDirectory + @"data\fromHG38toCHM13-alignments", dataDirectory + @"data\fromCHM13toHG38-alignments"),
                (dataDirectory + @"data\fromHG38toCHM13-assemblies", dataDirectory + @"data\fromCHM13toHG38-assemblies"),
            };

            var all = LoadExperimentData(dataFromExperimentsFiles[(int)dataType].FromHGtoCHM);
            var empty = LoadEmptyIntervals(dataFromExperimentsFiles[(int)dataType].FromCHMToHG + "-empty");
            var alternatives = LoadAlternativeIntervals(dataFromExperimentsFiles[(int)dataType].FromCHMToHG + "-alternatives");

            var dataFromArticle = LoadArticleData(@"d:\OpticData\tool\prediction_38.bed");

            var query = new Position(chromosome, start, end);
            Console.WriteLine("Structural variants induced by transition from HG38 to CHM13-T2T");
            var resultArticle = ProcessArticleData(query, dataFromArticle);
            foreach (var item in resultArticle)
            {
                Console.WriteLine($"\t{item.Name} {item.HG} {item.CHM}");
            }

            Console.WriteLine("Ambigous and other mapping events");
            var resultExperiments = ProcessExperimentData(query, all, empty, alternatives);
            foreach (var item in resultExperiments)
            {
                Console.WriteLine($"\t{(item.Detail ? "\t" : string.Empty)} {item.Source} {item.Message} {(item.Target1 is null ? string.Empty : item.Target1.ToString())} {(item.Target2 is null ? string.Empty : item.Target2.ToString())}");
            }
            Console.WriteLine();
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

        private static List<(string Name, Position HG, Position CHM)> ProcessArticleData(Position position, List<(string Name, Position HG, Position CHM)> dataFromArticle)
        {
            var result = new List<(string, Position, Position)>();
            foreach (var item in dataFromArticle)
            {
                if (item.HG.OverlapWith(position))
                {
                    result.Add(item);
                }
            }
            return result;
        }

        private static List<(Position Source, string Message, Position? Target1, Position? Target2, bool Detail)> ProcessExperimentData(Position position, List<(IntervalType Type, Position HG, Position? CHM)> all, List<Position> empty, List<(Position CHM, Position HG1, Position HG2)> alternatives)
        {
            var result = new List<(Position HG, string Message, Position? CHM1, Position? CHM2, bool Detail)>();
            bool atLeastOne = false;
            foreach (var item in all.Where(x => x.HG.OverlapWith(position)))
            {
                var intersection = position.IntersectionWith(item.HG);
                if (item.Type == IntervalType.Simple || item.Type == IntervalType.Alternative)
                {
                    var target = item.CHM?.GetMapping(item.HG, intersection);
                    var message = (item.Type == IntervalType.Simple) ? 2 : 3;
                    result.Add((intersection, messages[message], target, null, false));
                    if (target is not null)
                    {
                        if (target.Start > target.End) target = new Position(target.Chromosome, target.End, target.Start);
                        foreach (var emptyInterval in empty.Where(x => x.OverlapWith(target)))
                        {
                            result.Add((emptyInterval.IntersectionWith(target), messages[4], null, null, true));
                        }
                        foreach (var alternativeInterval in alternatives.Where(x => x.CHM.OverlapWith(target)))
                        {
                            var part = alternativeInterval.CHM.IntersectionWith(target);
                            result.Add((part, messages[5], alternativeInterval.HG1.GetMapping(alternativeInterval.CHM, part), alternativeInterval.HG2.GetMapping(alternativeInterval.CHM, part), true));
                        }
                    }
                }
                else if (item.Type == IntervalType.Empty)
                {
                    result.Add((intersection, messages[1], null, null, false));
                }

                atLeastOne = true;
            }

            if (!atLeastOne)
            {
                throw new Exception("There should be at least one result (check the query, if it is valid).");
            }

            return result;
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