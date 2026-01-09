import csv


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
    tfs = set()
    with open("all-tfs.tsv") as tfFile:
        reader = csv.DictReader(tfFile, delimiter="\t")
        index = 0
        for row in reader:
            tfs.add(row["Name.TF"])
            index += 1
            if (index % 100000) == 0:
                print(f"Processed {index} TF entries.")
    commonGenes = genes.intersection(tfs)
    with open("tfs.txt", "w") as outputFile:
        for gene in commonGenes:
            outputFile.write(gene + "\n")
    print(f"Found {len(commonGenes)} TFs in the gene list out of {index} total TFs.")


def compareDiffGeneToTfs():  # Compare the differently expressed genes to the TFs that were impacted by the SNPs
    diffGenes = set()
    with open("gene.txt") as geneFile:
        for line in geneFile:
            diffGenes.add(line.strip())
    tfs = set()
    with open("tfs.txt") as tfFile:
        for line in tfFile:
            tfs.add(line.strip())
    tfDict = {}
    with open("all-tfs.tsv") as tfFile:
        reader = csv.DictReader(tfFile, delimiter="\t")
        index = 0
        for row in reader:
            if row["Name.TF"] in tfs:
                tfDict[row["Name.TF"]] = tfDict.get(row["Name.TF"], []) + [
                    row["Name.Target"]
                ]
            index += 1
            if (index % 100000) == 0:
                print(f"Processed {index} TF entries.")
    impactedGenes = set()
    for gene in diffGenes:
        if any(gene in targets for targets in tfDict.values()):
            impactedGenes.add(gene)
    print(f"Identified {len(impactedGenes)} impacted genes from the TFs.")
    with open("impacted-genes.txt", "w") as outputFile:
        for gene in impactedGenes:
            outputFile.write(gene + "\n")


parseSnpFile()
compareGenesToTfs()
compareDiffGeneToTfs()
