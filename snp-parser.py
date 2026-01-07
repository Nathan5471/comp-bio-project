import csv


def parseSnpFile():
    with open("gwas-association.tsv") as file:
        reader = csv.DictReader(file, delimiter="\t")
        index = 0
        SNPs = {}
        Genes = {}
        for row in reader:
            if index == 0 or index == 1 or index == 2:
                print(row)
            index += 1
        print(len(reader))


parseSnpFile()
