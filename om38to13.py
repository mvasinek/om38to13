import argparse
import enum
import math
import os
import re

def ExitWithPrint(msg: str):
    print(msg)
    exit(1)    

messages = ["(Reversed)", 
            "no mapping from HG38 to CHM13", 
            "site positions in CHM13", 
            "alternative site mapping in CHM13", 
            "in CHM13 have no mapping to HG38", 
            "in CHM13 have multiple sources in HG38"]

class Position:
    def __init__(self, chromosome: int, start: float, end: float):
        self.__chromosome = chromosome
        self.__start = start
        self.__end = end

    def chromosome(self) -> int:
        return self.__chromosome

    def start(self) -> float:
        return self.__start

    def end(self) -> float:
        return self.__end

    def diff(self) -> float:
        return abs(self.__end - self.__start + 1)

    def isReversed(self) -> bool:
        return self.__start > self.__end

    def getMapping(self, whole: 'Position', part: 'Position') -> 'Position':
        def get(percent: float):
            return round(self.__start + (self.__end - self.__start) * percent)

        def getReversed(percent: float):
            return round(self.__start - (self.__start - self.__end) * percent)

        if whole.start() > whole.end() or part.start() > part.end():
            ExitWithPrint("Unexpected coordinate start and end positions, start should be always less than or equal to end position.")                  

        if self.__start == self.__end or whole.start() == whole.end():
            return self

        start_percent = (part.start() - whole.start()) / (whole.end() - whole.start())
        end_percent = (part.end() - whole.start()) / (whole.end() - whole.start())

        if(self.__start < self.__end):
            return Position(self.__chromosome, get(start_percent), get(end_percent))
        else:
            return Position(self.__chromosome, getReversed(start_percent), getReversed(end_percent))

    def intersectionWith(self, position: 'Position') -> 'Position':
        return Position(self.__chromosome, max(self.__start, position.start()), min(self.__end, position.end()))

    def overlapWith(self, position: 'Position') -> bool:
        if(self.__chromosome != position.chromosome()):
            return False
        return self.__start <= position.end() and position.start() <= self.__end


    def tString(self, epsilon):
        return "%s:%d" % (Position.ChrFromBionano(self.__chromosome), self.__end - epsilon)

    def __str__(self) -> str:
        if self.__start == self.__end:
            return "%s:%d" % (Position.ChrFromBionano(self.__chromosome), self.__start)
        else:
            return "%s:%d-%d" % (Position.ChrFromBionano(self.__chromosome), self.__start, self.__end)                

    @staticmethod
    def FromString(s: str) -> 'Position':
        chr_pos_reobj = re.match(r'(.+):(((\d+)-(\d+))|((\d+)))', s)
        if chr_pos_reobj == None:
            ExitWithPrint("Incorrectly specified position: " + s + ". It should match either chr:pos or chr:start-end.")

        pos12_obj = re.match(r'(\d+)-(\d+)', chr_pos_reobj.group(2))
        if pos12_obj == None:
            #single position
            s_pos = int(chr_pos_reobj.group(2))
            return Position(Position.ChrToBionano(chr_pos_reobj.group(1)), s_pos, s_pos)
        else:
            #start and end positions differs
            start_pos = int(pos12_obj.group(1))
            end_pos = int(pos12_obj.group(2))
            return Position(Position.ChrToBionano(chr_pos_reobj.group(1)), start_pos, end_pos)

    @staticmethod
    def ChrToBionano(chromosome: str) -> int:
        chr_low = chromosome.lower()
        if chr_low == "chrx":
            return 23
        elif chr_low == "chry":
            return 24
        else:
            return int(chr_low.replace("chr",""))

    @staticmethod
    def ChrFromBionano(chromosome: int) -> str:
        if chromosome == 23:
            return "chrX"
        elif chromosome == 24:
            return "chrY"
        else:
            return "chr%d" % chromosome

class SMAP:
    def __init__(self, h_lines, original, common, translocations, inversions):
        self.__h_lines = h_lines
        self.__original = original
        self.__common = common
        self.__translocations = translocations
        self.__inversions = inversions

    def serialize(self, outpath, smap_ids):
        f = open(outpath, "w")
        f.write("\n".join(self.__h_lines) + "\n")
        f.write("\n".join(list(map(lambda x: self.__original[x], smap_ids))))
        f.close()

    def common(self):
        return self.__common

    def translocations(self):
        return self.__translocations

    def inversions(self):
        return self.__inversions

    def original(self):
        return self.__original

    @staticmethod
    def FromFile(smap_path, epsilon):
        if not os.path.exists(smap_path):
            ExitWithPrint("The specified smap input file does not exists.")            

        h_lines = []
        common = []
        translocations = []
        inversions = []
        inversions_partial = []
        original = {}

        with open(smap_path) as f:            
            for line in f.readlines():
                if line.startswith("#"):
                    h_lines.append(line.rstrip())
                    continue

                l = line.strip()
                if len(l) == 0:
                    continue

                values = l.split("\t")
                smap_id = int(values[0])
                original[smap_id] = l

                if int(values[12]) != -1:
                    inversions_partial.append((int(values[0]), int(values[12]), int(values[2]), float(values[6]), float(values[7])))
                elif "trans" in l:
                    pos_a = float(values[6])
                    pos_b = float(values[7])
                    pos_a_start = 0 if pos_a - epsilon < 0 else pos_a - epsilon
                    pos_b_start = 0 if pos_b - epsilon < 0 else pos_b - epsilon
                    aroundA = Position(int(values[2]), pos_a_start, pos_a + epsilon)
                    aroundB = Position(int(values[3]), pos_b_start, pos_b + epsilon)
                    translocations.append((smap_id, aroundA, aroundB))
                else:
                    common.append((int(values[0]), Position(int(values[2]), float(values[6]), float(values[7]))))

        for id, link_id, chromosome, x, y in inversions_partial:
            if y != -1 and x > y: ExitWithPrint("Inversion: x must be smaller than y.") 
            if x < 0: ExitWithPrint("Inversion: the first coordinate cannot be negative")

            def existsInversions(id: int) -> bool:
                for item_a, item_b, _ in inversions:
                    if item_a == id or item_b == id: return True
                return False

            def locatePartial(linkID: int) -> bool:
                for inversion in inversions_partial:
                    if inversion[0] == linkID: return inversion
                return None

            if not existsInversions(id):
                itemB = locatePartial(link_id)

                if y == -1:
                    if x < itemB[3]: inversions.append((id, itemB[0], Position(chromosome, x, itemB[3])))
                    else: inversions.append((id, itemB[0], Position(chromosome, itemB[4], x)))
                elif itemB[4] == -1:
                    if itemB[3] < x: inversions.append((id, itemB[0], Position(chromosome, itemB[3], x)))
                    else: inversions.append((id, itemB[0], Position(chromosome, y, itemB[3])))
                else:
                    if x < itemB[3]: inversions.append((id, itemB[0], Position(chromosome, y, itemB[3])))
                    else: inversions.append((id, itemB[0], Position(chromosome, itemB[4], x)))

        return SMAP(h_lines, original, common, translocations, inversions)

class IntervalType(enum.Enum):
        Simple = 1
        Alternative = 2
        Empty = 3

class SMAPFilterType(enum.Enum):
        Mapping = 1
        Induced = 2
        Both = 3

class OMGenomeTools:    
    @staticmethod
    def Filter(smap_path, options, sample_path_g12, sample_path_g21, predictions_path, result_path, epsilon=10000):
        """
        Three options - mapping, induced, both
        """
        smap = SMAP.FromFile(smap_path, epsilon)
        smap_ids = []
        data_from_genome = OMGenomeTools.LoadArticleData(predictions_path)

        for item in smap.common():
            if OMGenomeTools.ProcessArticleData("", item[1], data_from_genome) == "":
                smap_ids.append(item[0])

        for item in smap.translocations():            
            if OMGenomeTools.ProcessArticleData("", item[1], data_from_genome) == "" and OMGenomeTools.ProcessArticleData("", item[2], data_from_genome) == "":
                smap_ids.append(item[0])

        for item in smap.inversions():
            if OMGenomeTools.ProcessArticleData("", item[2], data_from_genome) == "":
                smap_ids.append(item[0])
        
        smap_ids.sort()
        smap.serialize(result_path, smap_ids)

    @staticmethod
    def View(interval, sample_path_g12, sample_path_g21, predictions_path):
        position = Position.FromString(interval)
        all = OMGenomeTools.LoadExperimentData(sample_path_g12)
        empty = OMGenomeTools.LoadEmptyIntervals(sample_path_g21 + "-empty")
        alternatives = OMGenomeTools.LoadAlternativeIntervals(sample_path_g21 + "-alternatives")
        data_from_genome = OMGenomeTools.LoadArticleData(predictions_path)

        print("Position: ", str(position))
        print("\tStructural variants induced by transition from HG38 to CHM13-T2T")
        print("\t\tType|HG38 coordinate (size)|CHM13 coordinate (size)")
        sv_data = OMGenomeTools.ProcessArticleData("\t\t", position, data_from_genome)
        if sv_data == "":
            print("\t\tNone")
        else:
            print(sv_data)

        print("\tAmbigous and other mapping events")
        print(OMGenomeTools.ProcessExperimentData("\t\t", position, all, empty, alternatives))
        print("\n")

    @staticmethod
    def Annotate(smap_path, sample_path_g12, sample_path_g21, predictions_path, result_path, epsilon=10000):
        smap = SMAP.FromFile(smap_path, epsilon)

        data_from_experiments = []
        all = OMGenomeTools.LoadExperimentData(sample_path_g12)
        empty = OMGenomeTools.LoadEmptyIntervals(sample_path_g21 + "-empty")
        alternatives = OMGenomeTools.LoadAlternativeIntervals(sample_path_g21 + "-alternatives")

        data_from_experiments.append((all, empty, alternatives))
        data_from_genome = OMGenomeTools.LoadArticleData(predictions_path)

        def Process(fstream, indentation, position, d_from_article, d_all, d_empty, d_alternatives):
            fstream.write("\tStructural variants induced by transition from HG38 to CHM13-T2T\n")
            sv_data = OMGenomeTools.ProcessArticleData(indentation + "\t\t", position, d_from_article)
            if sv_data == "":
                fstream.write("\t\tNone\n")
            else:
                fstream.write(sv_data)

            fstream.write("\tAmbigous and other mapping events\n")
            fstream.write(OMGenomeTools.ProcessExperimentData(indentation + "\t\t", position, d_all, d_empty, d_alternatives))
            fstream.write("\n")

        ostream = open(result_path, "w")
        orig = smap.original()
        for item in smap.common():
            ostream.write(str(item[0]) + " " + str(item[1]) + " " + orig[item[0]] + "\n")
            Process(ostream, "", item[1], data_from_genome, data_from_experiments[0][0], data_from_experiments[0][1], data_from_experiments[0][2])

        for item in smap.translocations():
            ostream.write(str(item[0]) + " TRANSLOCATION A=" + item[1].tString(epsilon) + "\t" + item[2].tString(epsilon) + "\n")
            ostream.write("\tA: " + str(item[1]))
            Process(ostream, "\t", item[1], data_from_genome, data_from_experiments[0][0], data_from_experiments[0][1], data_from_experiments[0][2])
            ostream.write("\tB: " + str(item[2]))
            Process(ostream, "\t", item[2], data_from_genome, data_from_experiments[0][0], data_from_experiments[0][1], data_from_experiments[0][2])
            ostream.write("\n")

        for item in smap.inversions():
            ostream.write(str(item[0]) + " " + str(item[1]) + " INVERSION " + str(item[2]) + "\n")
            Process(ostream, "", item[2], data_from_genome, data_from_experiments[0][0], data_from_experiments[0][1], data_from_experiments[0][2])

        ostream.close()

    @staticmethod
    def LoadArticleData(fileName: str):
        result = []
        with open(fileName) as f:
            for line in f.readlines():
                ls = line.strip()
                if len(ls) == 0:
                    continue

                values = ls.replace(";","\t").split("\t")
                if(len(values) != 3): ExitWithPrint("Unexpected data format of: %s" % fileName)

                result.append((values[0], Position.FromString(values[1]), Position.FromString(values[2])))

        return result

    @staticmethod
    def LoadExperimentData(fileName: str):
        result = []

        with open(fileName) as f:
            for line in f.readlines():
                ls = line.strip()
                if len(ls) == 0:
                    continue

                values = ls.replace(";","\t").split("\t")
                if len(values) == 2 and values[0].startswith("E"):
                    result.append((IntervalType.Empty, Position.FromString(values[1]), None))
                elif len(values) == 3:
                    itype = IntervalType.Simple if values[0].startswith("S") else IntervalType.Alternative
                    result.append((itype, Position.FromString(values[1]), Position.FromString(values[2])))
                else:
                    ExitWithPrint("Unexpected input data for %s" % fileName)

        return result

    @staticmethod
    def LoadEmptyIntervals(fileName: str):
        result = []

        with open(fileName) as f:
            for line in f.readlines():
                ls = line.strip()
                if len(ls) == 0:
                    continue

                result.append(Position.FromString(ls))
        
        return result

    @staticmethod
    def LoadAlternativeIntervals(fileName: str):
        result = []

        with open(fileName) as f:
            for line in f.readlines():
                ls = line.strip()
                if len(ls) == 0:
                    continue
                
                values = ls.split("\t")
                result.append((Position.FromString(values[0]), Position.FromString(values[1]), Position.FromString(values[2])))

        return result

    @staticmethod
    def ProcessArticleData(indentation: str, position: 'Position', dataFromArticle) -> str:
        result = ""

        for item in dataFromArticle:
            if item[1].overlapWith(position):
                result += indentation + item[0] + " " + str(item[1]) + " (" + str(item[1].diff()) + ") " + str(item[2]) + " (" + str(item[2].diff()) + ")" + "\n"

        return result

    @staticmethod
    def ProcessExperimentData(indentation: str, position: 'Position', all, empty, alternatives) -> str:
        result = ""
        atLeastOne = False

        all_overlaped = list(filter(lambda x: x[1].overlapWith(position), all))

        for itype, ihg, ichm in all_overlaped:
            intersection = position.intersectionWith(ihg)
            if itype == IntervalType.Simple or itype == IntervalType.Alternative:
                target = ichm.getMapping(ihg, intersection)      
                message = 2 if itype == IntervalType.Simple else 3
                
                msg_end = messages[0] if ichm != None and ichm.isReversed() else ""
                result += indentation + str(intersection) + " " + messages[message] + " " + str(target) + msg_end + "\n"

                if target != None:
                    if target.start() > target.end():
                        target = Position(target.chromosome(), target.end(), target.start())
                    empty_overlaped = list(filter(lambda x: x.overlapWith(target), empty))
                    for emptyival in empty_overlaped:
                        result += indentation + "\t" + str(emptyival.intersectionWith(target)) + " " + messages[4] + "\n"
                    alternative_overlaped = list(filter(lambda x: x[0].overlapWith(target), alternatives))
                    for alt_ival in alternative_overlaped:                                             
                        part = alt_ival[0].intersectionWith(target)
                        result += indentation + "\t" + str(part) + " " + messages[5] + " " + str(alt_ival[1].getMapping(alt_ival[0], part)) + " " + str(alt_ival[2].getMapping(alt_ival[0], part)) + "\n"
            elif itype == IntervalType.Empty:
                result += indentation + str(intersection) + " " + messages[1] + "\n"

            atLeastOne = True

        return result

def run(args):    
    op = args.op
    input = args.input
    
    #verify the input and op correctness
    if op in ["annotate" ,"filter"]:
        if not input.lower().endswith(".smap"):
            ExitWithPrint("Annotate and filter commands requires smap file as input.")            
        if op == "annotate":
            if args.output == None:
                args.output = args.input.replace(".smap",".annotated.txt")
            OMGenomeTools.Annotate(args.input, "data/fromHG38toCHM13-alignments", "data/fromCHM13toHG38-alignments", "data/prediction_38.bed", args.output)
        elif op == "filter":
            if args.output == None:
                args.output = args.input.replace(".smap",".filtered.smap")
            OMGenomeTools.Filter(args.input, None, "data/fromHG38toCHM13-alignments", "data/fromCHM13toHG38-alignments", "data/prediction_38.bed", args.output)
    elif op == "view":
        OMGenomeTools.View(args.input,"data/fromHG38toCHM13-alignments", "data/fromCHM13toHG38-alignments", "data/prediction_38.bed")
        

    #OMGenomeTools.DoWork(args.input, "data/fromHG38toCHM13-alignments", "data/fromCHM13toHG38-alignments", "data/prediction_hg38.bed", args.output)
    #OMGenomeTools.Annotate("test_data/exp_refineFinal1_merged_filter_inversions.smap", "data/fromHG38toCHM13-alignments", "data/fromCHM13toHG38-alignments", "data/prediction_38.bed", "results/out.txt")
    #OMGenomeTools.Filter("test_data/exp_refineFinal1_merged_filter_inversions.smap", None, "data/fromHG38toCHM13-alignments", "data/fromCHM13toHG38-alignments", "data/prediction_38.bed", "results/out.smap")
    #OMGenomeTools.View("chr1:1000000-2000000","data/fromHG38toCHM13-alignments", "data/fromCHM13toHG38-alignments", "data/prediction_38.bed")

if __name__ == "__main__":

    
    parser = argparse.ArgumentParser(prog="om38to13", description="The program for providing additional information to structural variants found by Bionano software tools. For more details please refer to README.txt file.")
    
    parser.add_argument("op", choices=["annotate","filter","view"], help="Select one of the supported operations annotate, filter or view")
    parser.add_argument("input", help="Specify path to the smap file for annotate and filter commands or a region in a standard genomic location format for view command: e.g. chr1:1000, chr2:1000-2000")
    parser.add_argument("-d", "--distance", type=int, default=10000, help="Specify the distance for acceptation of translocation. It is a distance from the translocation breakpoint.")
    parser.add_argument("-o", "--output", help="Specify the output file path, if not provided, .smap will be replaced by .filter.smap or .annotated.smap")    
    parser.add_argument("-w", "--data_dir", help="Specify directory with default data")
    
    args = parser.parse_args()    

    run(args)
    exit(0)