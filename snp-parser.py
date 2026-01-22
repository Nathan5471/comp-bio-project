import csv
import os
import json


def parseSnpFile():
    with open("gwas-association.tsv") as file:
        reader = csv.DictReader(file, delimiter="\t")
        index = 0
        SNPs = {}
        Genes = {}
        for row in reader:
            snpData = row["SNPS"].split(", ")
            for snp in snpData:
                SNPs[snp] = SNPs[snp] + 1 if snp in SNPs else 1
            geneData = row["REPORTED GENE(S)"].split(", ")
            for gene in geneData:
                Genes[gene] = Genes[gene] + 1 if gene in Genes else 1
            index += 1
    filteredSNPs = {k: v for k, v in SNPs.items() if v > 1}
    filteredGenes = {k: v for k, v in Genes.items() if v > 1}
    with open("snp.txt", "w") as snpFile:
        for snp in filteredSNPs:
            snpFile.write(snp + "\n")
    with open("gene.txt", "w") as geneFile:
        for gene in filteredGenes:
            geneFile.write(gene + "\n")
    print(
        f"Processed {index} entries resulting in {len(filteredSNPs)} SNPs and {len(filteredGenes)} genes."
    )


def compareGenesToTfs():
    genes = set()
    with open("gene.txt") as geneFile:
        for line in geneFile:
            genes.add(line.strip())
    tfs = {}
    with open("all-tfs.tsv") as tfFile:
        reader = csv.DictReader(tfFile, delimiter="\t")
        index = 0
        for row in reader:
            tf = row["Name.TF"]
            if tf in genes:
                tfs[tf] = tfs.get(tf, []) + [row["Name.Target"]]
            index += 1
            if (index % 100000) == 0:
                print(f"Processed {index} TF entries.")
    commonGenes = genes.intersection(tfs)
    with open("tfs.txt", "w") as outputFile:
        json.dump(tfs, outputFile)
    print(f"Found {len(commonGenes)} TFs in the gene list out of {index} total TFs.")


def analyzeGeoData(fileName: str):
    geoData = []
    with open(f"geo2r/{fileName}") as geoFile:
        reader = csv.DictReader(geoFile, delimiter="\t")
        pValueTag = None
        adjPValueTag = None
        if "p.value" in reader.fieldnames:
            pValueTag = "p.value"
            adjPValueTag = "adj.p.value"
        elif "P.value" in reader.fieldnames:
            pValueTag = "P.value"
            adjPValueTag = "adj.P.value"
        elif "P.Value" in reader.fieldnames:
            pValueTag = "P.Value"
            if "adj.P.Value" in reader.fieldnames:
                adjPValueTag = "adj.P.Value"
            elif "adj.P.Val" in reader.fieldnames:
                adjPValueTag = "adj.P.Val"
        elif "pvalue" in reader.fieldnames:
            pValueTag = "pvalue"
            adjPValueTag = "padj"
        else:
            print("No recognizable p-value fields found:", reader.fieldnames)
            return
        index = 0
        for row in reader:
            if (row[pValueTag] != "NA" and float(row[pValueTag]) < 0.05) and (
                row[adjPValueTag] != "NA" and float(row[adjPValueTag]) < 0.05
            ):
                geneSymbol = row["Symbol"]
                if geneSymbol and not geneSymbol in geoData:
                    geoData.append(geneSymbol)
            index += 1
            if (index % 10000) == 0:
                print(f"Processed {index} GEO entries.")
    print(f"Total significant GEO entries: {len(geoData)}")
    tfData = None
    with open("tfs.txt") as tfFile:
        tfData: dict = json.load(tfFile)
    if tfData is None:
        print("TF data not found.")
        return
    impactedTfs = {}
    for gene in geoData:
        for tf in tfData:
            if gene in tfData[tf]:
                currentGeneList = impactedTfs.get(tf, [])
                if gene not in currentGeneList:
                    impactedTfs[tf] = impactedTfs.get(tf, []) + [gene]
    with open(f"output/impacted-genes-{fileName}.txt", "w") as outputFile:
        json.dump(impactedTfs, outputFile)


def findGeneData(fileName: str):  # Only make a list, don't fetch summary
    genes = set()
    with open(f"output/{fileName}") as geneFile:
        data = json.load(geneFile)
        for tf in data:
            for gene in data[tf]:
                if gene not in genes:
                    genes.add(gene)
    print(f"Total unique genes found: {len(genes)}")

    with open(f"output/list-{fileName}", "w") as outputFile:
        json.dump(",".join(genes), outputFile)


def identifyBrainGenes(fileName: str):
    brainGenes = set()
    with open(f"output/{fileName}") as geneFile:
        data = json.load(geneFile)
        brainGenes = data.replace('"', "").split(",")
        if len(brainGenes) == 0 or (len(brainGenes) == 1 and brainGenes[0] == ""):
            print("No genes found in the list.")
            return
    brainCategories = [
        "GOBP_AMINO_ACID_NEUROTRANSMITTER_REUPTAKE",
        "GOBP_ANTEROGRADE_DENRITIC_TRANSPORT_OF_NEUROTRANSMITTER_RECEPTOR_COMPLEX",
        "GOBP_CALCIUM_ION_REGULATED_EXOCYTOSIS_OF_NEUROTRANSMITTER",
        "GOBP_NEGATIVE_REGULATION_OF_NEUROTRANSMITTER_SECRETION",
        "GOBP_NEGATIVE_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT",
        "GOBP_NEUROTRANSMITTER_GATED_ION_CHANNEL_CLUSTERING",
        "GOBP_NEUROTRANSMITTER_LOADING_INTO_SYNAPTIC_VESICLE",
        "GOBP_NEUROTRANSMITTER_RECEPTOR_INTERNALIZATION",
        "GOBP_NEUROTRANSMITTER_RECEPTOR_LOCALIZATION_TO_POSTSYNAPTIC_SPECIALIZATION_MEMBRANE",
        "GOBP_NEUROTRANSMITTER_RECEPTOR_TRANSPORT",
        "GOBP_NEUROTRANSMITTER_RECEPTOR_TRANSPORT_ENDOSOME_TO_POSTSYNAPTIC_MEMBRANE",
        "GOBP_NEUROTRANSMITTER_RECEPTOR_TRANSPORT_TO_PLASMA_MEMBRANE",
        "GOBP_NEUROTRANSMITTER_REUPTAKE",
        "GOBP_NEUROTRANSMITTER_SECRETION",
        "GOBP_NEUROTRANSMITTER_TRANSPORT",
        "GOBP_NEUROTRANSMITTER_UPTAKE",
        "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_SECRETION",
        "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT",
        "GOBP_POSITIVE_REGULATION_OF_NEUROTRANSMITTER_UPTAKE",
        "GOBP_REGULATION_OF_CALCIUM_ION_DEPENDENT_EXOCYTOSIS_OF_NEUROTRANSMITTER",
        "GOBP_REGULATION_OF_NEUROTRANSMITTER_RECEPTOR_ACTIVITY",
        "GOBP_REGULATION_OF_NEUROTRANSMITTER_RECEPTOR_LOCALIZATION_TO_POSTSYNAPTIC_SPECIALIZATION_MEMBRANE",
        "GOBP_REGULATION_OF_NEUROTRANSMITTER_SECRETION",
        "GOBP_REGULATION_OF_NEUROTRANSMITTER_TRANSPORT",
        "GOBP_REGULATION_OF_NEUROTRANSMITTER_UPTAKE",
        "GOBP_REGULATION_OF_POSTSYNAPTIC_MEMBRANE_NEUROTRANSMITTER_RECEPTOR_LEVELS",
        "GOBP_REGULTION_OF_POSTSYNAPTIC_NEUROTRANSMITTER_RECEPTOR_INTERNALIZATION",
        "GOBP_SPONTANEOUS_NEUROTRANSMITTER_SECRETION",
        "GOBP_CENTRAL_NERVOUS_SYSTEM_MYELIN_FORMATION",
        "GOBP_CENTRAL_NERVOUS_SYSTEM_MYELIN_MAINTENANCE",
        "GOBP_MYELIN_ASSEMBLY",
        "GOBP_MYELIN_MAINTENANCE",
        "GOBP_NEGATIVE_REGULATION_OF_MYELINATION",
        "GOBP_PERIPHERAL_NERVOUS_SYSTEM_MYELIN_MAINTENANCE",
        "GOBP_POSITIVE_REGULATION_OF_MYELINATION",
        "GOBP_REGULATION_MYELINATION",
        "GOBP_SPHINGOMYELIN_BIOSYNTHETIC_PROCESS",
        "GOBP_SPHINGOMYELIN_CATABOLIC_PROCESS",
        "GOBP_SPHINGOMYELIN_METABOLIC_PROCESS",
    ]


if not os.path.exists("snp.txt") and not os.path.exists(
    "gene.txt"
):  # Skip if already done to save time
    parseSnpFile()
if not os.path.exists("tfs.txt"):
    compareGenesToTfs()
for fileName in os.listdir("geo2r"):
    if fileName.endswith(".tsv") and not os.path.exists(
        f"output/impacted-genes-{fileName}.txt"
    ):
        print(f"Analyzing GEO data file: {fileName}")
        analyzeGeoData(fileName)
for fileName in os.listdir("output"):
    if (
        fileName.startswith("impacted-genes-")
        and fileName.endswith(".txt")
        and not os.path.exists(f"output/list-{fileName}")
    ):
        print(f"Finding gene data in file: {fileName}")
        findGeneData(fileName)
for fileName in os.listdir("output"):
    if (
        fileName.startswith("list-impacted-genes-")
        and fileName.endswith(".txt")
        and not os.path.exists(f"output/brain-genes-{fileName}")
    ):
        print(f"Identifying brain genes in file: {fileName}")
        identifyBrainGenes(fileName)
