import csv
import os
import json
import time
from Bio import Entrez
from dotenv import load_dotenv


load_dotenv()


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


def findGeneData(fileName: str):
    genes = set()
    with open(f"output/{fileName}") as geneFile:
        data = json.load(geneFile)
        for tf in data:
            for gene in data[tf]:
                if gene not in genes:
                    genes.add(gene)
    print(f"Total unique genes found: {len(genes)}")
    Entrez.email = os.getenv("email", "")
    Entrez.api_key = os.getenv("api", "")
    Entrez.sleep_between_tries = 0.2  # 5/second rate limit

    geneInfo = {}
    for gene in genes:
        handle = Entrez.esearch(
            db="gene",
            term=f"{gene}[Gene Name] AND Homo sapiens[Organism]",
            retmode="xml",
        )
        record = Entrez.read(handle)
        handle.close()
        print(f"Searching for gene: {gene}, found {record['Count']} entries.")
        if not record["IdList"]:
            continue
        geneId = record["IdList"][0]
        handle = Entrez.esummary(db="gene", id=geneId, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        summary = records["DocumentSummarySet"]["DocumentSummary"][0]["Summary"]
        geneInfo[gene] = summary
    with open(f"output/gene-data-{fileName}", "w") as outputFile:
        json.dump(geneInfo, outputFile)


def checkBrainRelation(fileName: str, save=False):
    approvedData = {}
    data = None
    savedIndex = 0
    with open(f"output/{fileName}") as inputFile:
        data = json.load(inputFile)
    if save:
        with open(f"saves/{fileName}.txt") as saveFile:
            savedData = json.load(saveFile)
            approvedData = savedData["data"]
            savedIndex = savedData["index"]
    if not data:
        return
    index = 0
    for gene in data:
        if save and index < savedIndex:
            index += 1
            continue
        print("Gene Name:", gene)
        print(data[gene])
        while True:
            response = input("Is this gene related to brain function? (y/n): ")
            if response.lower() == "y":
                approvedData[gene] = data[gene]
                break
            elif response.lower() == "n":
                break
            elif response.lower() == "save":
                with open(f"saves/{fileName}.txt", "w") as saveFile:
                    json.dump({"data": approvedData, "index": index}, saveFile)
                print(f"Progress saved at index {index}.")
                quit()
            else:
                print("Please input 'y' or 'n'.")
        print("")
        print("")
        print("--------------------------------")
        print("")
        print("")
        index += 1
    with open(f"diffGenes-{fileName}", "w") as outputFile:
        json.dump(approvedData, outputFile)


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
        and not os.path.exists(f"output/gene-data-{fileName}")
    ):
        print(f"Finding gene data in file: {fileName}")
        findGeneData(fileName)
for fileName in os.listdir("output"):
    if fileName.startswith("gene-data-") and fileName.endswith(".txt"):
        print(f"Checking brain relation for file: {fileName}")
        if os.path.exists(f"saves/{fileName}.txt"):
            checkBrainRelation(fileName, save=True)
        else:
            checkBrainRelation(fileName)
