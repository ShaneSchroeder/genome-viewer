import sys  # for getting command line arguments
from collections import OrderedDict  # for the ordered dictionary in Dinucleotide counter
import re
import matplotlib.pyplot as plt
import editdistance

FASTA_LINE_NUM = 50
NON_FASTA_LINE_NUM = 100
FASTQ_FILE = False


def readFASTQ(fileName: str) -> list:
    """
    Reads a FASTQ file and retrieves the name, description, sequence, and quality scores for the sequence(s).

    :param fileName: name of the FASTQ file
    :return: list containting 4 other lists. [ [names], [descriptions], [sequences], [quality scores] ]
    """
    file_handler = open(fileName, 'r')
    my_file = file_handler.read()
    numberOfSequences = my_file.count("\n") // 4

    seq_name_list = []
    description_list = []
    sequence_list = []
    quality_list = []
    for x in range(0, numberOfSequences * 4, 4):
        curFASTQData = my_file.split('\n')[x:x + 4]
        defLine = curFASTQData[0][1:].split(" ", maxsplit=1)

        name = defLine[0]
        description = defLine[1]
        sequence = curFASTQData[1]
        quality = curFASTQData[3]

        seq_name_list.append(name)
        description_list.append(description)
        sequence_list.append(sequence)
        quality_list.append(quality)

    file_handler.close()
    return [seq_name_list, description_list, sequence_list, quality_list]


def readFASTA(fileName: str) -> list:
    """
    Takes a FASTA file name and parses the file to extract the name, description, and sequence.

    :param fileName: name of the file to process
    :return: list: in the format [name, description, sequence]
    """
    file_handler = open(fileName, 'r')
    my_file = file_handler.read()

    seq_name_list = []
    description_list = []
    sequence_list = []

    numberOfSequencesInFile = my_file.count(">")
    indices = [i for i, x in enumerate(my_file.split('\n')) if x.find(">") != -1]

    for x in range(0, numberOfSequencesInFile):
        if x == numberOfSequencesInFile - 1:
            curFASTAData = my_file.split('\n')[indices[x]:]
        else:
            curFASTAData = my_file.split('\n')[indices[x]:indices[x + 1]]

        defLine = curFASTAData[0][1:].split(" ", maxsplit=1)

        name = defLine[0]
        description = defLine[1]
        sequence = ''.join(curFASTAData[1:])

        seq_name_list.append(name)
        description_list.append(description)
        sequence_list.append(sequence)

    file_handler.close()
    return [seq_name_list, description_list, sequence_list]


def printInFASTA(SeqName: str = "", SeqDescription: str = "", Sequence: str = "",
                 Quality: str = ""):
    """
    Displays file in FASTA format, along with quality values if the file is FASTQ.

    :param SeqName: name of the sequence in the FASTA file
    :param SeqDescription: description of sequence in the FASTA file
    :param Sequence: sequence of the FASTA file
    :param Quality: If file is FASTQ, will print the quality values below the sequence.
    """
    print(">" + SeqName + " " + SeqDescription)
    # Print with a set amount of characters per line
    for charPosition in range(len(Sequence.strip())):
        if charPosition % FASTA_LINE_NUM == 0 and charPosition != 0:
            print("\n", end="")
        print(Sequence[charPosition], end="")
    print()
    print()
    if Quality != "":
        for charPosition in range(len(Quality.strip())):
            if charPosition % FASTA_LINE_NUM == 0 and charPosition != 0:
                print("\n", end="")
            print(Quality[charPosition], end="")
        print()
    print()


# Really really really ugly code, but it works, time complexity is O(n) though so thats decent.
def printWithRuler(SeqName="", SeqDescription="", Sequence="", Spacer="", Quality: str = ""):
    """
    Prints out the FASTA or FASTQ file with a ruler, and a spacer if desired.

    :param SeqName: Name of the sequence
    :param SeqDescription: Description of the sequence
    :param Sequence: The sequence itself
    :param Spacer: An extra space between blocks of the sequence if desired.
    :param Quality: quality values of the sequence, only given if file is FASTQ.
    """
    if SeqName != "" and SeqDescription != "":
        print(">" + SeqName + " " + SeqDescription + "\n")

    # prints top ruler that numbers each group of 10
    if Spacer == "":
        print("  ", end="")
    print("   " + Spacer, end="")
    for topruler in range(10):
        if topruler == 9:
            if Spacer == " ":
                print("          0", end="")
            else:
                print("         0", end="")
        else:
            print("         " + Spacer + str(topruler + 1), end="")
        if topruler == 9:
            print()

    print("Line ", end="")
    horizontalCount = 0
    for lineRulerHorizontal in range(100):
        if lineRulerHorizontal % 10 == 0 and lineRulerHorizontal != 0:
            print("0" + Spacer, end="")
            horizontalCount = 0
        elif horizontalCount != 0:
            print(horizontalCount, end="")
        horizontalCount += 1
    print("0")

    lineRulerVertical = 1
    for charPosition in range(len(Sequence.strip())):
        if charPosition == 0:
            if Spacer == "":
                print("   " + str(lineRulerVertical) + " ", end="")
            else:
                print("   " + str(lineRulerVertical), end="")
        if charPosition % NON_FASTA_LINE_NUM == 0 and charPosition != 0:
            if lineRulerVertical >= 999:
                if Spacer == " ":
                    print("\n" + str(lineRulerVertical + 1), end="")
                else:
                    print("\n" + str(lineRulerVertical + 1) + " ", end="")
            else:
                if lineRulerVertical >= 99:
                    if Spacer == " ":
                        print("\n" + " " + str(lineRulerVertical + 1), end="")
                    else:
                        print("\n" + " " + str(lineRulerVertical + 1) + " ", end="")
                else:
                    if lineRulerVertical >= 9:
                        if Spacer == " ":
                            print("\n" + "  " + str(lineRulerVertical + 1), end="")
                        else:
                            print("\n" + "  " + str(lineRulerVertical + 1) + " ", end="")
                    else:
                        if Spacer == " ":
                            print("\n" + "   " + str(lineRulerVertical + 1), end="")
                        else:
                            print("\n" + "   " + str(lineRulerVertical + 1) + " ", end="")
            lineRulerVertical += 1
        if charPosition % 10 == 0:
            print(Spacer, end="")
        print(Sequence[charPosition], end="")

    # just for formatting
    print()
    print()

    if Quality != "":

        # prints top ruler that numbers each group of 10
        if Spacer == "":
            print("  ", end="")
        print("   " + Spacer, end="")
        for topruler in range(10):
            if topruler == 9:
                if Spacer == " ":
                    print("          0", end="")
                else:
                    print("         0", end="")
            else:
                print("         " + Spacer + str(topruler + 1), end="")
            if topruler == 9:
                print()

        print("Line ", end="")
        horizontalCount = 0
        for lineRulerHorizontal in range(100):
            if lineRulerHorizontal % 10 == 0 and lineRulerHorizontal != 0:
                print("0" + Spacer, end="")
                horizontalCount = 0
            elif horizontalCount != 0:
                print(horizontalCount, end="")
            horizontalCount += 1
        print("0")

        lineRulerVertical = 1
        for charPosition in range(len(Quality.strip())):
            if charPosition == 0:
                if Spacer == "":
                    print("   " + str(lineRulerVertical) + " ", end="")
                else:
                    print("   " + str(lineRulerVertical), end="")
            if charPosition % NON_FASTA_LINE_NUM == 0 and charPosition != 0:
                if lineRulerVertical >= 999:
                    if Spacer == " ":
                        print("\n" + str(lineRulerVertical + 1), end="")
                    else:
                        print("\n" + str(lineRulerVertical + 1) + " ", end="")
                else:
                    if lineRulerVertical >= 99:
                        if Spacer == " ":
                            print("\n" + " " + str(lineRulerVertical + 1), end="")
                        else:
                            print("\n" + " " + str(lineRulerVertical + 1) + " ", end="")
                    else:
                        if lineRulerVertical >= 9:
                            if Spacer == " ":
                                print("\n" + "  " + str(lineRulerVertical + 1), end="")
                            else:
                                print("\n" + "  " + str(lineRulerVertical + 1) + " ", end="")
                        else:
                            if Spacer == " ":
                                print("\n" + "   " + str(lineRulerVertical + 1), end="")
                            else:
                                print("\n" + "   " + str(lineRulerVertical + 1) + " ", end="")
                lineRulerVertical += 1
            if charPosition % 10 == 0:
                print(Spacer, end="")
            print(Quality[charPosition], end="")

    print()
    print()
    print()


def displayMode(SpcSeq: int, cleanedSeq: str):
    """
    Asks the user a series of questions to figure out in what format the file should be displayed in.
    """
    # Read the data
    if not FASTQ_FILE:
        extractedData = readFASTA(str(sys.argv[1]))
        seq_name_list = extractedData[0]
        seq_description_list = extractedData[1]
        sequence_list = extractedData[2]
    else:
        extractedData = readFASTQ(str(sys.argv[1]))
        seq_name_list = extractedData[0]
        seq_description_list = extractedData[1]
        sequence_list = extractedData[2]
        quality_list = extractedData[3]

    print("\n[Part I]: Display Mode\n")
    print("Question 1: Do you want to view the sequence in FASTA format or not (Y|N)?\n")
    viewInFasta = input("Enter Y or N: ")
    if viewInFasta == "Y" or viewInFasta == "y":
        if not FASTQ_FILE:
            printInFASTA(seq_name_list[SpcSeq], seq_description_list[SpcSeq], cleanedSeq)
        else:
            for x in range(len(seq_name_list)):
                printInFASTA(seq_name_list[x], seq_description_list[x], sequence_list[x], quality_list[x])
    else:
        print("Question 2: Do you need a spacer for viewing nucleotide positions (Y|N)?\n")
        viewWithSpacer = input("Enter Y or N: ")
        if viewWithSpacer == "Y" or viewWithSpacer == "y":
            if not FASTQ_FILE:
                printWithRuler(seq_name_list[SpcSeq], seq_description_list[SpcSeq], cleanedSeq, " ")
            else:
                for x in range(len(seq_name_list)):
                    printWithRuler(seq_name_list[x], seq_description_list[x], sequence_list[x], " ", quality_list[x])
        else:
            if not FASTQ_FILE:
                printWithRuler(seq_name_list[SpcSeq], seq_description_list[SpcSeq], cleanedSeq, "")
            else:
                for x in range(len(seq_name_list)):
                    printWithRuler(seq_name_list[x], seq_description_list[x], sequence_list[x], "", quality_list[x])
    main()


def nucleotideCounter(Sequence: str) -> tuple:
    """
    Counts each nucleotide A, T, G, C, N, or other that is present in the sequence.

    :param Sequence: The sequence to be analyzed.
    :return: tuple: (ACount, TCount, GCount, CCount, NCount, OtherCount)
    """
    ACount = Sequence.count("A")
    TCount = Sequence.count("T")
    GCount = Sequence.count("G")
    CCount = Sequence.count("C")
    NCount = Sequence.count("N")
    OtherCount = len(Sequence) - (ACount + TCount + GCount + CCount + NCount)
    return ACount, TCount, GCount, CCount, NCount, OtherCount


def gcContent(Sequence: str) -> float:
    """
    Calculates the percentage of the nucleotides G and C in the sequence.

    :param Sequence: The sequence to be analyzed.
    :return: float: percentage of sequence is made up of C and G.
    """
    nucleotideCounts = nucleotideCounter(Sequence)
    return round((nucleotideCounts[2] + nucleotideCounts[3]) / len(Sequence) * 100, 2)


def CpGIsland(Sequence: str) -> dict:
    """
    Finds sequences of CpG islands. It will store the starting and ending indices of the island as well as the length.

    A CpG Island is defined as the minimum continuous CpG di-nucleotides of ≥6 CpG dinucleotides.

    :param Sequence: The sequence to be analyzed for CpG islands.
    :return: dict: The dictionary containting CpG islands, if none were found will be an empty dictionary.
            { island# : "startIndex-endIndex_length" }
    """
    CpGDict = {}
    finds = 0
    lengthOfCpGIsland = 12
    startIndex = Sequence.find("CGCGCGCGCGCG")
    while startIndex != -1:
        finds += 1
        startOfAdditionalIsland = Sequence[startIndex + 12:]
        for curIndex in range(0, len(startOfAdditionalIsland), 2):
            nextTwoChars = startOfAdditionalIsland[curIndex:curIndex + 2]
            if nextTwoChars == "CG":
                lengthOfCpGIsland += 2
            else:
                break

        CpGDict[finds] = str(startIndex) + "-" + str(startIndex + lengthOfCpGIsland - 1) + "_" + str(lengthOfCpGIsland)
        # new startIndex is at the sliced strings next match of "CGCGCGCGCGCG" after the just found island + old startIndex.
        if Sequence[startIndex + lengthOfCpGIsland:].find("CGCGCGCGCGCG") != -1:
            startIndex = (startIndex + lengthOfCpGIsland) + Sequence[startIndex + lengthOfCpGIsland:].find(
                "CGCGCGCGCGCG")
            lengthOfCpGIsland = 12
        else:
            break
    return CpGDict


def printDiNucleotide(od: OrderedDict):
    """
    Prints the ordered dictionary in the specified format.

    :param od: An ordered dictionary containting the dinucleotide counts for the sequence.
    """
    nucleotides = ['T', 'C', 'A', 'G']
    print("\n\t\t   2nd")
    print("       ---------------------------")
    print("1st      T      C      A      G", end="")
    counter = 0
    index = 0
    for key, value in od.items():
        if counter % 4 == 0:
            print('\n' + nucleotides[index] + '      ', end="")
            index += 1
            if counter != 0:
                counter = 0
        print("{0}={1: >3} ".format(key, value), end="")
        counter += 1
    print()


def diNucleotideCounter(Sequence: str) -> OrderedDict:
    """
    Generates an ordered map containing all dinucleotides combinations found in the sequence.

    :param Sequence: The selected DNA sequence.
    :return: an ordered map of the counts found in the sequence.
    """
    od = OrderedDict()
    od['TT'] = Sequence.count('TT')
    od['TC'] = Sequence.count('TC')
    od['TA'] = Sequence.count('TA')
    od['TG'] = Sequence.count('TG')
    od['CT'] = Sequence.count('CT')
    od['CC'] = Sequence.count('CC')
    od['CA'] = Sequence.count('CA')
    od['CG'] = Sequence.count('CG')
    od['AT'] = Sequence.count('AT')
    od['AC'] = Sequence.count('AC')
    od['AA'] = Sequence.count('AA')
    od['AG'] = Sequence.count('AG')
    od['GT'] = Sequence.count('GT')
    od['GC'] = Sequence.count('GC')
    od['GA'] = Sequence.count('GA')
    od['GG'] = Sequence.count('GG')
    return od


def codonProfile(Sequence: str) -> OrderedDict:
    """
    Will construct an ordered dictionary containting the counts for all codons in the sequence.

    :param Sequence: The DNA sequence
    :return: OrdredDict conataining the codon counts for the sequence in order.
    """
    nucleotides = ['T', 'C', 'A', 'G']
    od = OrderedDict()
    for counter1 in range(0, 4):
        for counter2 in range(0, 4):
            for counter3 in range(0, 4):
                codon = nucleotides[counter1] + nucleotides[counter3] + nucleotides[counter2]
                od[codon] = Sequence.count(codon)
    return od


def printCodonProfile(od: OrderedDict):
    """
    Will print the codon profile in the format specified by the assignment.

    :param od: The codon profile as an ordered dictionary
    """
    nucleotides = ['T', 'C', 'A', 'G']
    print("                   2nd")
    print("        ------------------------------- ")
    print("1st      T      C       A      G     3rd", end="")
    counter = 0
    index1 = 0
    index3 = 0
    for key, value in od.items():
        if counter % 4 == 0:
            if counter % 16 == 0:
                print('  ' + nucleotides[index3] + '\n\n' + nucleotides[index1] + '     ', end="")
                index1 += 1
                index3 = 0
            else:
                print('  ' + nucleotides[index3] + '\n' + '      ', end="")
                index3 += 1
            if index3 == 4:
                index3 = 0
        print("{0}={1: >3} ".format(key, value), end="")
        counter += 1
    print("  G")
    print()


def reverseComplement(cleanSeq: str) -> str:
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(nn[n] for n in reversed(cleanSeq))


def translation(cleanSeq: str) -> str:
    trans_dic = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    protein = []
    for i in range(0, len(cleanSeq)-2, 3):
        protein.append(trans_dic[cleanSeq[i:i+3]])
    return "".join(protein)


def detectHomopolymer(Seq: str) -> dict:
    homopolymerDict = {"PolyA": [], "PolyT": [], "PolyG": [], "PolyC": []}
    baseList = ['A', 'T', 'G', 'C']
    for base in baseList:
        counter = 1
        pattern = "(["+base+"]{6,})+"
        matches = re.finditer(pattern, Seq)
        for match in matches:
            lengthOfPoly = match.end() - match.start()
            resultString = str(counter) + "=" + str(match.start()) + "-" + str(match.end()) + "_" + str(lengthOfPoly) + " "
            curKey = "Poly" + base
            homopolymerDict[curKey].append(resultString)
            counter += 1
    return homopolymerDict


def detectORFs(seq: str) -> dict:
    ORFDict = {"+1": [], "+2": [], "+3": [], "-1": [], "-2": [], "-3": []}
    # first three reading frames
    for x in range(3):
        translated = translation(seq[x:])
        pattern = "(M[^\*]{8,}?\*)"
        matches = re.finditer(pattern, translated)
        endIndexOfLastMatch = 0
        curKey = "+" + str(x+1)
        for match in matches:
            ORFDict[curKey].append(translated[match.start():match.end()])
            endIndexOfLastMatch = match.end()

        # now check if there is an M that does not reach * before end of amino acid sequence.
        indexOfNonEndedStart = translated.find("M", endIndexOfLastMatch, len(translated))
        if indexOfNonEndedStart != -1 and translated[indexOfNonEndedStart:].find("*") == -1:
            ORFDict[curKey].append(translated[indexOfNonEndedStart:])

    # last three reading frames on reverse compliment string
    for x in range(3):
        reverseComplementSeq = reverseComplement(seq)
        translated = translation(reverseComplementSeq[x:])
        pattern = "(M[^\*]{8,}?\*)"
        matches = re.finditer(pattern, translated)
        endIndexOfLastMatch = 0
        curKey = "-" + str(x+1)
        for match in matches:
            ORFDict[curKey].append(translated[match.start():match.end()])
            endIndexOfLastMatch = match.end()

        # now check if there is an M that does not reach * before end of amino acid sequence.
        indexOfNonEndedStart = translated.find("M", endIndexOfLastMatch, len(translated))
        if indexOfNonEndedStart != -1 and translated[indexOfNonEndedStart:].find("*") == -1:
            ORFDict[curKey].append(translated[indexOfNonEndedStart:])

    return ORFDict


def analysisMode(seq):
    """
    Displays length, GC content, nucleotide counts, and CpG islands of the full DNA sequence.
    """
    if not FASTQ_FILE:
        if str.isnumeric(str(seq)):
            sequenceString = readFASTA(str(sys.argv[1]))[2][seq]
        else:
            sequenceString = seq
    else:
        sequenceString_list = readFASTQ(str(sys.argv[1]))[2]

    print("\n[Part II]: Analysis Mode\n")
    if FASTQ_FILE:
        for x in range(len(sequenceString_list)):
            print("Sequence Length=" + str(len(sequenceString_list[x])))
            print("GC Content=" + str(gcContent(sequenceString_list[x])) + "%")
            nCount = nucleotideCounter(sequenceString_list[x])
            print("Nucleotide Counts: A=[{0}] T=[{1}] G=[{2}] C=[{3}] N=[{4}] Other=[{5}]".format(*nCount))
            CpGIslandDict = CpGIsland(sequenceString_list[x])
            print("CpG Islands: ", end="")
            for entryNum in range(len(CpGIslandDict)):
                print(str(entryNum + 1) + "=" + str(CpGIslandDict[entryNum + 1]) + " ", end="")
            print()
            print()
    else:
        print("(2.1) Sequence Length=" + str(len(sequenceString)))

        print("(2.2) GC Content=" + str(gcContent(sequenceString)) + "%")
        nCount = nucleotideCounter(sequenceString)

        print("(2.3) Nucleotide Counts: A=[{0}] T=[{1}] G=[{2}] C=[{3}] N=[{4}] Other=[{5}]".format(*nCount))
        CpGIslandDict = CpGIsland(sequenceString)

        print("(2.4) CpG Islands: ", end="")
        for entryNum in range(len(CpGIslandDict)):
            print(str(entryNum + 1) + "=" + str(CpGIslandDict[entryNum + 1]) + " ", end="")
        print()

        print("(2.5) Homopolymers:")
        homopolymerDict = detectHomopolymer(sequenceString)
        for key, value in homopolymerDict.items():
            print("      {0} {1}".format(key, "".join(value)))

        print("(2.6) Di-nucleotide Counts:")
        diNucleotideDict = diNucleotideCounter(sequenceString)
        printDiNucleotide(diNucleotideDict)

        print("\n(2.7) Codon Profile:\n")
        codonProfileDict = codonProfile(sequenceString)
        printCodonProfile(codonProfileDict)

        print("(2.8) detected ORFs: ")
        ORFDict = detectORFs(sequenceString)
        for key, value in ORFDict.items():
            if len(value) == 0:
                print("ORF ({0}):".format(key))
            else:
                for i in range(len(value)):
                    if i == 0:
                        print("ORF ({0}): {1}".format(key, value[0]))
                    else:
                        print("          {0}".format(value[i]))

    main()


def motifFinder(Sequence: str) -> dict:
    """
    Will find motif(s) in the given sequence and store the data in a dictionary in the following format:
        {
            motif : instanceNumber=startPosition-endPosition_length ...
            ...
        }

    :param Sequence: The DNA sequence to process
    :return: A dictionary containting all apropriate data about the user entered motif(s)
    """
    # GCGA, GAAAA, GAAA
    print("Please enter the motif that you want to find (e.g., TATATTATA) ")
    motifs = input("Enter your motif sequences here: ")
    motifDict = {}
    if motifs == "":
        return motifDict
    elif motifs == "homopolymer":
        homopolymerDict = detectHomopolymer(Sequence)
        return homopolymerDict
    else:
        motifList = motifs.split(", ")
        for motif in motifList:
            finds = re.finditer(motif, Sequence)
            instance = 1
            completeResult = ""
            for result in finds:
                start = result.start()
                end = result.end()
                length = end - start
                completeResult += str(instance) + "=" + str(start) + "-" + str(end) + "_" + str(length) + " "
                instance += 1
            motifDict[motif] = completeResult
        return motifDict


def diNucleotideCompare(firstDinuc: OrderedDict, secondDinuc: OrderedDict, firstName: str, secondName: str):
    """
    Will compare di-nucleotide counts of two specified DNA sequences. It will display a helpful chart to
    see how the sequences differ in count.

    :param firstDinuc: Di-nucleotide ordered dictionary corresponding to the first sequence.
    :param secondDinuc: Di-nucleotide ordered dictionary corresponding to the second sequence.
    :param firstName: The FASTA name of the first sequence.
    :param secondName: The FASTA name of the second sequence.
    """
    if firstName == secondName:
        print("{0} and {1} are the same in di-nucleotides and codon profiles:".format(firstName, secondName))
    else:
        print("{0} and {1} are different in di-nucleotides and codon profiles:".format(firstName, secondName))

    nucleotides = ['T', 'C', 'A', 'G']
    print("                                2nd")
    print("       --------------------------------------------------")
    print("1st     T            C            A            G")
    counter = 0
    index = 0
    for (k, v), (k2, v2) in zip(firstDinuc.items(), secondDinuc.items()):
        if counter % 4 == 0:
            print('\n' + nucleotides[index] + '      ', end="")
            index += 1
            if counter != 0:
                counter = 0
        print("{0}[{1: >3}:{2: >3}]  ".format(k, v, v2), end="")
        counter += 1
    print()


def codonProfileCompare(firstCodons: OrderedDict, secondCodons: OrderedDict, firstName: str, secondname: str):
    """
    Will compare the codon profiles of two DNA sequences. Will display a helpful chart to understand how the two
    codon profiles differ.

    :param firstCodons: Codon ordered dictionary corresponding to the first sequence.
    :param secondCodons: Codon ordered dictionary corresponding to the second sequence.
    :param firstName: FASTA name of the first sequence.
    :param secondname: FASTA name of the second sequence.
    """
    nucleotides = ['T', 'C', 'A', 'G']
    print("\n                   2nd")
    print("      ----------------------------------- ")
    print("1st    T       C       A      G       3rd", end="")
    counter = 0
    index1 = 0
    index3 = 0
    for (k, v), (k2, v2) in zip(firstCodons.items(), secondCodons.items()):
        if counter % 4 == 0:
            if counter % 16 == 0:
                print('  ' + nucleotides[index3] + '\n\n' + nucleotides[index1] + '     ', end="")
                index1 += 1
                index3 = 0
            else:
                print('  ' + nucleotides[index3] + '\n' + '      ', end="")
                index3 += 1
            if index3 == 4:
                index3 = 0
        symbolToPrint = ""
        if v > v2:
            if v - v2 >= 10:
                symbolToPrint = ">"
            if v - v2 >= 20:
                symbolToPrint = ">>"
            if v - v2 >= 30:
                symbolToPrint = ">>>"
        elif v < v2:
            if v - v2 <= -10:
                symbolToPrint = "<"
            if v - v2 <= -20:
                symbolToPrint = "<<"
            if v - v2 <= -30:
                symbolToPrint = "<<<"
        print("{0} {1: >3} ".format(k, symbolToPrint), end="")
        counter += 1
    print("  G")
    print()
    createImage(firstCodons, secondCodons)


def fragment_show(SequenceChunk: str, startIndex: int, endIndex: int):
    SequenceChunkLength = len(SequenceChunk)

    ORFDict = {"+1": "", "+2": "", "+3": "", "-1": "", "-2": "", "-3": ""}


    reverseSequenceChunk = reverseComplement(SequenceChunk)[::-1]
    for x in range(3):
        translated = translation(reverseSequenceChunk[x:])
        curKey = "-" + str(x+1)
        ORFDict[curKey] += translated

    # TOP
    print("   ", end="")
    for char in ORFDict["-3"]:
        print("{0}  ".format(char), end="")
    print()

    # MIDDLE
    print(" ", end="")
    for char in ORFDict["-2"]:
        print("{0}  ".format(char), end="")
    print()

    # BOTTOM
    print("  ", end="")
    for char in ORFDict["-1"]:
        print("{0}  ".format(char), end="")
    print()

    print(reverseSequenceChunk)

    for pipe in range(SequenceChunkLength):
            print("|", end="")
    print()
    indexPlacesPrinted = len(str(startIndex) + str(endIndex))
    for nucleotideNumber in range(SequenceChunkLength - indexPlacesPrinted):
        if nucleotideNumber == 0:
            print("<" + str(startIndex), end="")
        elif nucleotideNumber == SequenceChunkLength - indexPlacesPrinted - 1:
            print(str(endIndex) + ">", end="")
        else:
            print("-", end="")
    print()
    for pipe in range(SequenceChunkLength):
        print("|", end="")
    print()

    print(SequenceChunk)

    for x in range(3):
        translated = translation(SequenceChunk[x:])
        curKey = "+" + str(x+1)
        ORFDict[curKey] += translated

    # TOP
    print(" ", end="")
    for char in ORFDict["+1"]:
        print("{0}  ".format(char), end="")
    print()

    # MIDDLE
    print("  ", end="")
    for char in ORFDict["+2"]:
        print("{0}  ".format(char), end="")
    print()

    # BOTTOM
    print("   ", end="")
    for char in ORFDict["+3"]:
        print("{0}  ".format(char), end="")
    print()

    print()
    print()


def createImage(firstCodonDict: OrderedDict, secondCodonDict: OrderedDict):
    x1 = firstCodonDict.keys()
    y1 = firstCodonDict.values()

    plt.subplot(2, 1, 1)
    plt.xlabel("64 codons")
    plt.ylabel("frequency")
    plt.bar(x1, y1)
    plt.legend(["seq1"])

    x2 = secondCodonDict.keys()
    y2 = secondCodonDict.values()
    plt.subplot(2, 1, 2)
    plt.xlabel("64 codons")
    plt.ylabel("frequency")
    plt.bar(x2, y2, color="orange")
    plt.legend(["seq2"])

    plt.show()


def inquiry(SeqNumber: int, seq: str):
    """
    Takes user input to analyze part of the DNA Sequence in sizes of 10-100 nucleotides.
    Will print data about the sequence chunk such as, length, GC Content, nucelotide counts,
    and CpG Islands out to the user.

    :param seq: the cleaned sequence.
    :param SeqNumber: The DNA Sequence number.
    """
    Sequence = seq

    print("(3.1) Extract a DNA fragment")
    print("(3.2) Find a motif")
    print("(3.3) Compare two sequences\n")

    inquirySelect = input("Your selection is: ")
    if inquirySelect == "3.2":
        motifDict = motifFinder(Sequence)
        newSequence = Sequence
        if "PolyA" in motifDict:
            for key, value in motifDict.items():
                print("      {0} {1}".format(key, "".join(value)))
                # nasty stuff to find the span of homopolymers
                for match in value:
                    rightOfEquals = match.split("=")
                    start = int(rightOfEquals[1].split("-")[0])
                    end = int(rightOfEquals[1].split("-")[1].split("_")[0])
                    newSequence = newSequence[:start] + newSequence[start:end].lower() + newSequence[end:]
            print()
            printWithRuler("", "", newSequence, " ")
        else:
            print("Motif\t\t\t Frequency\t\t\t Detail")
            for key in motifDict:
                frequency = motifDict[key].count("=")
                print("{0}\t\t\t {1}\t\t\t {2}".format(key, frequency, motifDict[key]))
                # To print with the motifs lowercase AND allow overlapping motifs to be fully lower case.
                finds = re.finditer(key, Sequence)  # redundant I know but easiest way to get this data.
                for find in finds:
                    start = find.start()
                    end = find.end()
                    newSequence = newSequence[:start] + newSequence[start:end].lower() + newSequence[end:]
            print()
            printWithRuler("", "", newSequence, " ")
    elif inquirySelect == "3.3":
        sequenceNames = readFASTA(str(sys.argv[1]))[0]
        currentSeqName = sequenceNames[SeqNumber]
        print("Which Sequence do you want to examime? [", end="")
        for x in range(1, numberOfSequencesInInputFile + 1):
            if x == numberOfSequencesInInputFile:
                print(str(x), end="")
            else:
                print(str(x) + "|", end="")
        print("]\n")
        sequenceSelect = int(input("Your selection is: ")) - 1
        print("(" + sequenceNames[sequenceSelect] + ")")
        secondSequenceName = sequenceNames[sequenceSelect]
        secondSequence = readFASTA(str(sys.argv[1]))[2][sequenceSelect]

        # I ASSUME THESE SHOULD BE TRIMMED OF ADAPTERS EVEN FOR COMPARING? NOT LISTED IN ASSIGNMENT.
        adapterFileName = input("Please enter the adapter file name: ")
        trimmingData, secondSequence = removeAdapter(sequenceSelect, adapterFileName)

        firstDinuc = diNucleotideCounter(Sequence)
        secondDinuc = diNucleotideCounter(secondSequence)
        diNucleotideCompare(firstDinuc, secondDinuc, currentSeqName, secondSequenceName)

        firstCodons = codonProfile(Sequence)
        secondCodons = codonProfile(secondSequence)
        codonProfileCompare(firstCodons, secondCodons, currentSeqName, secondSequenceName)
        print()

    else:
        rangeInput = input("Please enter the start and end positions (e.g., 19-48): ")

        # at least one dash in the input, will basically ignore extra dashes,
        # including any information entered after an additional dash.
        while rangeInput.find("-") == -1:
            print("Invalid format, please enter a range in the form x-y.")
            rangeInput = input("Please enter the start and end positions separated by a dash: ")

        # Check if the starting index input is:
        # 1. a number
        # 2. greater or equal to 0
        # 3. not out of bounds.
        startIndex = rangeInput.split("-")[0]
        while not startIndex.isnumeric():
            startIndex = input("Invalid Positions or format! Please enter a new starting position: ")
            while int(startIndex) < 0:
                startIndex = input("Invalid Positions or format! Please enter a new starting position: ")
                while int(startIndex) >= len(Sequence):
                    startIndex = input("Invalid Positions or format! Please enter a new starting position: ")
        startIndex = int(startIndex)

        # Check if the ending index input is:
        # 1. a number
        # 2. greater than starting index
        # 3. not out of bounds
        # 4. greater than start.
        endIndex = rangeInput.split("-")[1]
        while not endIndex.isnumeric():
            endIndex = input("Invalid Positions or format! Please enter a new ending position: ")
            while int(endIndex) < 0:
                endIndex = input("Invalid Positions or format! Please enter a new ending position: ")
                while int(endIndex) >= len(Sequence):
                    endIndex = input("Invalid Positions or format! Please enter a new ending position: ")
                    while int(endIndex) <= startIndex:
                        endIndex = input("Invalid Positions or format! Please enter a new ending position: ")
        endIndex = int(endIndex)  # At this point endIndex is guaranteed to be a proper index.

        if endIndex - startIndex + 1 <= 0 or endIndex - startIndex + 1 < 10 or endIndex - startIndex > 100 \
                or endIndex - startIndex + 1 > len(Sequence):
            print("Invalid Fragment Length. It should be between 10 and 100!")
            inquiry(SeqNumber, Sequence)

        SequenceChunk = Sequence[startIndex:endIndex]
        SequenceChunkLength = len(SequenceChunk)
        print("\nThe fragment you selected has a length of " + str(SequenceChunkLength) + " nucleotides:\n")

        fragment_show(SequenceChunk, startIndex, endIndex)

        print("Fragment Length=" + str(endIndex - startIndex))
        print("Fragment GC content=" + str(gcContent(SequenceChunk)) + "%")
        nucleotideCounts = nucleotideCounter(SequenceChunk)
        print("Nucleotide Counts: A=[{0}] T=[{1}] G=[{2}] C=[{3}] N=[{4}] Other=[{5}]".format(*nucleotideCounts))
        CpGIslandDict = CpGIsland(SequenceChunk)
        print("CpG Islands: ", end="")
        for entryNum in range(len(CpGIslandDict)):
            print(str(entryNum + 1) + "=" + str(CpGIslandDict[entryNum + 1]) + " ", end="")
        print()

        print("Homopolymers:")
        homopolymerDict = detectHomopolymer(SequenceChunk)
        for key, value in homopolymerDict.items():
            print("{0} {1}".format(key, "".join(value)))
        print("Di-nucleotide Counts:")
        diNucleotideDict = diNucleotideCounter(SequenceChunk)
        printDiNucleotide(diNucleotideDict)
        print("\nCodon Profile:\n")
        codonProfileDict = codonProfile(SequenceChunk)
        printCodonProfile(codonProfileDict)

        print("detected ORFs: ")
        ORFDict = detectORFs(SequenceChunk)
        for key, value in ORFDict.items():
            if len(value) == 0:
                print("ORF ({0}):".format(key))
            else:
                for i in range(len(value)):
                    if i == 0:
                        print("ORF ({0}): {1}".format(key, value[0]))
                    else:
                        print("          {0}".format(value[i]))

    main()


def getNumberOfSequences(fileName: str) -> int:
    """
    My own helper method, fetches the total amount of sequences in the supplied input file.

    :param fileName: the name of the input file.
    :return: int, the number of sequences in the file.
    """
    file_handler = open(fileName, 'r')
    my_file = file_handler.read()
    counts = my_file.count(">")
    file_handler.close()
    return counts


def removeAdapter(seqNumber, adapterFileName):
    Sequence = readFASTA(str(sys.argv[1]))[2][seqNumber]
    cleanedSequence = Sequence[:]
    if adapterFileName == '':
        adapterFileName = "Adapter.fasta"
    file_handler = open(adapterFileName, 'r')
    my_file = file_handler.read()
    adapterList = []
    eachLine = my_file.split("\n")
    for x in range(1, len(eachLine), 2):
        adapterList.append(eachLine[x])

    trimmingData = {}
    # 5' end
    for adapter in adapterList:
        fivePrimeEnd = Sequence[0:41]
        distance = editdistance.distance(fivePrimeEnd, adapter)
        if distance <= len(fivePrimeEnd) - len(adapter):
            startIndex = cleanedSequence.find(adapter, 0, 41)
            startEndIndices = [startIndex, startIndex+len(adapter)]
            fivePrimeDict = {"5'" : startEndIndices}
            cleanedSequence = cleanedSequence.replace(adapter, '', 1)
            if adapter in trimmingData:
                trimmingData[adapter] += fivePrimeDict
            else:
                trimmingData[adapter] = fivePrimeDict
            break

    # 3' end
    for adapter in adapterList:
        threePrimeEnd = Sequence[-40:]
        distance = editdistance.distance(threePrimeEnd, adapter)
        if distance <= len(threePrimeEnd) - len(adapter):
            startIndex = cleanedSequence.find(adapter, (len(cleanedSequence)-40), len(cleanedSequence))
            startEndIndices = [startIndex, startIndex+len(adapter)]
            cleanedSequence = cleanedSequence.replace(adapter, '', 1)  # assuming the adapter occurs once in the entire sequence
            threePrimeDict = {"3'" : startEndIndices}
            if adapter in trimmingData:
                # if already exists, need to merge dictionaries
                trimmingData[adapter] = {**trimmingData[adapter], **threePrimeDict}
            else:
                trimmingData[adapter] = threePrimeDict
            break

    return trimmingData, cleanedSequence


def main():
    """
    Main function of the program. Displays the different modes that the user can enter and calls the
    corresponding function based on the input.
    """
    print("Which Sequence do you want to examime? [", end="")
    for x in range(1, numberOfSequencesInInputFile + 1):
        if x == numberOfSequencesInInputFile:
            print(str(x), end="")
        else:
            print(str(x) + "|", end="")
    print("]\n")
    sequenceNames = readFASTA(str(sys.argv[1]))[0]
    sequenceSelect = int(input("Your selection is: ")) - 1
    print("The sequence name that you selected is: (" + sequenceNames[sequenceSelect] + ")")
    adapterFileName = input("Please enter the adapter file name: ")
    print("Now clearing the sequence ...")
    trimmingData, cleanedSequence = removeAdapter(sequenceSelect, adapterFileName)
    print("Sequence cleaning is finished!\n")

    if trimmingData == {}:
        print("Case A: No adapter was detected at all")
    else:
        for adapter, dictValue in trimmingData.items():
            if "5'" in dictValue:
                print("Case B: Adapter sequence ({0}) was found in 5’ end, [start:end]=[{1}:{2}]".format(adapter, dictValue["5'"][0], dictValue["5'"][1]))
            if "3'" in dictValue:
                print("Case C: Adapter sequence ({0}) was found in 3’ end, [start:end]=[{1}:{2}]".format(adapter, dictValue["3'"][0], dictValue["3'"][1]))

    print("\nThere are four options for you:\n")
    print("[1] Display Mode")
    print("[2] Analysis Mode")
    print("[3] Inquiry Mode")
    print("[4] Quit\n")
    modeselect = (input("Your selection is: "))

    while modeselect != "1" and modeselect != "2" and modeselect != "3" and modeselect != "4":
        print("Invalid input, please enter a number 1-4")
        modeselect = (input("Your selection is: "))

    if modeselect == "1":
        displayMode(sequenceSelect, cleanedSequence)

    if modeselect == "2":
        analysisMode(cleanedSequence)

    if modeselect == "3":
        print("[Part III]: Inquiry Mode\n")
        inquiry(sequenceSelect, cleanedSequence)

    if modeselect == "4":
        quit()


print("Welcome Sequence Viewer! Programmer: schroe51\n")

if len(sys.argv) == 1:
    print("No supplied file to process. Bye.")
    quit()

if sys.argv[1].find(".fasta") == -1 and sys.argv[1].find(".fastq") == -1:
    print("Invalid input file, please include a .fasta or .fastq file as an argument. Goodbye.")
    quit()

numberOfSequencesInInputFile = getNumberOfSequences(sys.argv[1])

print("There is " + str(numberOfSequencesInInputFile) + " sequence(s) detected in the file: " + str(sys.argv[1]) + "\n")
if sys.argv[1].split(".")[1] == "fastq":
    FASTQ_FILE = True
    print("!!! Input file is in FASTQ format !!!")
if __name__ == "__main__":
    main()
