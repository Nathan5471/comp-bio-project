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
        "GOBP_ANTEROGRADE_NEURONAL_DENSE_CORE_VESICLE_TRANSPORT",
        "GOBP_BRANCHIOMOTOR_NEURON_AXON_GUIDANCE",
        "GOBP_CELL_MORPHOGENESIS_INVOLVED_IN_NEURON_DIFFERENTIATION",
        "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_AXONOGENESIS",
        "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DEVELOPMENT",
        "GOBP_CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION",
        "GOBP_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS",
        "GOBP_CEREBRAL_CORTEX_GABAERGIC_INTERNEURON_DEVELOPMENT",
        "GOBP_CEREBRAL_CORTEX_GABAERGIC_INTERNEURON_DIFFERENTIATION",
        "GOBP_CEREBRAL_CORTEX_GABAERGIC_INTERNEURON_MIGRATION",
        "GOBP_CEREBRAL_CORTEX_NEURON_DIFFERENTIATION",
        "GOBP_COMMISSURAL_NEURON_AXON_GUIDANCE",
        "GOBP_COMMITMENT_OF_NEURONAL_CELL_TO_SPECIFIC_NEURON_TYPE_IN_FOREBRAIN",
        "GOBP_DOPAMINERGIC_NEURON_AXON_GUIDANCE",
        "GOBP_DOPAMINERGIC_NEURON_DIFFERENTIATION",
        "GOBP_ENSHEATHMENT_OF_NEURONS",
        "GOBP_FOREBRAIN_GENERATION_OF_NEURONS",
        "GOBP_FOREBRAIN_NEURON_DEVELOPMENT",
        "GOBP_FOREBRAIN_NEURON_FATE_COMMITMENT",
        "GOBP_GABAERGIC_NEURON_DIFFERENTIATION",
        "GOBP_GENERATION_OF_NEURONS",
        "GOBP_GLUTAMATERGIC_NEURON_DIFFERENTIATION",
        "GOBP_GONADOTROPHIN_RELEASING_HORMONE_NEURONAL_MIGRATION_TO_THE_HYPOTHALAMUS",
        "GOBP_HIPPOCAMPAL_NEURON_APOPTOTIC_PROCESS",
        "GOBP_HYPOTHALAMUS_GONADOTROPHIN_RELEASING_HORMONE_NEURON_DIFFERENTIATION",
        "GOBP_INTERNEURON_MIGRATION",
        "GOBP_INTERNEURON_MIGRATION_FROM_THE_SUBPALLIUM_TO_THE_CORTEX",
        "GOBP_MIDBRAIN_DOPAMINERGIC_NEURON_DIFFERENTIATION",
        "GOBP_MOTOR_NEURON_APOPTOTIC_PROCESS",
        "GOBP_MOTOR_NEURON_AXON_GUIDANCE",
        "GOBP_MOTOR_NEURON_MIGRATION",
        "GOBP_NEGATIVE_REGULATION_OF_HIPPOCAMPAL_NEURON_APOPTOTIC_PROCESS",
        "GOBP_NEGATIVE_REGULATION_OF_MOTOR_NEURON_APOPTOTIC_PROCESS",
        "GOBP_NEGATIVE_REGULATION_OF_NEURON_APOPTOTIC_PROCESS",
        "GOBP_NEGATIVE_REGULATION_OF_NEURON_DIFFERENTIATION",
        "GOBP_NEGATIVE_REGULATION_OF_NEURON_MIGRATION",
        "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT",
        "GOBP_NEGATIVE_REGULATION_OF_NEURON_PROJECTION_REGENERATION",
        "GOBP_NEGATIVE_REGULATION_OF_OXIDATIVE_STRESS_INDUCED_NEURON_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
        "GOBP_NEURONAL_ACTION_POTENTIAL",
        "GOBP_NEURONAL_DENSE_CORE_VESICLE_EXOCYTOSIS",
        "GOBP_NEURONAL_ION_CHANNEL_CLUSTERING",
        "GOBP_NEURONAL_SIGNAL_TRANSDUCTION",
        "GOBP_NEURONAL_STEM_CELL_POPULATION_MAINTENANCE",
        "GOBP_NEURON_APOPTOTIC_PROCESS",
        "GOBP_NEURON_CELLULAR_HOMEOSTASIS",
        "GOBP_NEURON_CELL_CELL_ADHESION",
        "GOBP_NEURON_DEVELOPMENT",
        "GOBP_NEURON_FATE_COMMITMENT",
        "GOBP_NEURON_FATE_DETERMINATION",
        "GOBP_NEURON_FATE_SPECIFICATION",
        "GOBP_NEURON_GLIAL_CELL_SIGNALING",
        "GOBP_NEURON_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
        "GOBP_NEURON_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_OXIDATIVE_STRESS",
        "GOBP_NEURON_MATURATION",
        "GOBP_NEURON_MIGRATION",
        "GOBP_NEURON_NEURON_SYNAPTIC_TRANSMISSION",
        "GOBP_NEURON_PROJECTION_ARBORIZATION",
        "GOBP_NEURON_PROJECTION_EXTENSION",
        "GOBP_NEURON_PROJECTION_EXTENSION_INVOLVED_IN_NEURON_PROJECTION_GUIDANCE",
        "GOBP_NEURON_PROJECTION_GUIDANCE",
        "GOBP_NEURON_PROJECTION_MAINTENANCE",
        "GOBP_NEURON_PROJECTION_ORGANIZATION",
        "GOBP_NEURON_PROJECTION_REGENERATION",
        "GOBP_NEURON_RECOGNITION",
        "GOBP_NEURON_REMODELING",
        "GOBP_NORADRENERGIC_NEURON_DEVELOPMENT",
        "GOBP_NORADRENERGIC_NEURON_DIFFERENTIATION",
        "GOBP_OLFACTORY_BULB_INTERNEURON_DEVELOPMENT",
        "GOBP_OLFACTORY_BULB_INTERNEURON_DIFFERENTIATION",
        "GOBP_PERIPHERAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION",
        "GOBP_POSITIVE_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY",
        "GOBP_POSITIVE_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY",
        "GOBP_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION",
        "GOBP_POSITIVE_REGULATION_OF_NEURON_MIGRATION",
        "GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_ARBORIZATION",
        "GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT",
        "GOBP_POSITIVE_REGULATION_OF_NEURON_PROJECTION_REGENERATION",
        "GOBP_PYRAMIDAL_NEURON_DEVELOPMENT",
        "GOBP_PYRAMIDAL_NEURON_DIFFERENTIATION",
        "GOBP_RADIAL_GLIA_GUIDED_PYRAMIDAL_NEURON_MIGRATION",
        "GOBP_REGULATION_OF_DOPAMINERGIC_NEURON_DIFFERENTIATION",
        "GOBP_REGULATION_OF_HIPPOCAMPAL_NEURON_APOPTOTIC_PROCESS",
        "GOBP_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY",
        "GOBP_REGULATION_OF_MOTOR_NEURON_APOPTOTIC_PROCESS",
        "GOBP_REGULATION_OF_NEURONAL_ACTION_POTENTIAL",
        "GOBP_REGULATION_OF_NEURONAL_SYNAPTIC_PLASTICITY",
        "GOBP_REGULATION_OF_NEURON_APOPTOTIC_PROCESS",
        "GOBP_REGULATION_OF_NEURON_DIFFERENTIATION",
        "GOBP_REGULATION_OF_NEURON_MATURATION",
        "GOBP_REGULATION_OF_NEURON_MIGRATION",
        "GOBP_REGULATION_OF_NEURON_PROJECTION_ARBORIZATION",
        "GOBP_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT",
        "GOBP_REGULATION_OF_NEURON_PROJECTION_REGENERATION",
        "GOBP_REGULATION_OF_SHORT_TERM_NEURONAL_SYNAPTIC_PLASTICITY",
        "GOBP_RETINAL_BIPOLAR_NEURON_DIFFERENTIATION",
        "GOBP_RETROGRADE_NEURONAL_DENSE_CORE_VESICLE_TRANSPORT",
        "GOBP_SEROTONERGIC_NEURON_AXON_GUIDANCE",
        "GOBP_SOMATIC_MOTOR_NEURON_DIFFERENTIATION",
        "GOBP_SPINAL_CORD_ASSOCIATION_NEURON_DIFFERENTIATION",
        "GOBP_SPINAL_CORD_MOTOR_NEURON_CELL_FATE_SPECIFICATION",
        "GOBP_SPINAL_CORD_MOTOR_NEURON_DIFFERENTIATION",
        "GOBP_SYMPATHETIC_NEURON_PROJECTION_EXTENSION",
        "GOBP_VENTRAL_SPINAL_CORD_INTERNEURON_DIFFERENTIATION",
        "GOBP_EXCITATORY_SYNAPSE_ASSEMBLY",
        "GOBP_IMMUNOLOGICAL_SYNAPSE_FORMATION",
        "GOBP_INHIBITORY_SYNAPSE_ASSEMBLY",
        "GOBP_MAINTENANCE_OF_SYNAPSE_STRUCTURE",
        "GOBP_NEGATIVE_REGULATION_OF_SYNAPSE_ASSEMBLY",
        "GOBP_NEGATIVE_REGULATION_OF_SYNAPSE_ORGANIZATION",
        "GOBP_POSITIVE_REGULATION_OF_PROTEIN_LOCALIZATION_TO_SYNAPSE",
        "GOBP_POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY",
        "GOBP_POSITIVE_REGULATION_OF_SYNAPSE_MATURATION",
        "GOBP_POSTSYNAPSE_ASSEMBLY",
        "GOBP_POSTSYNAPSE_ORGANIZATION",
        "GOBP_POSTSYNAPSE_TO_NUCLEUS_SIGNALING_PATHWAY",
        "GOBP_PRESYNAPSE_ASSEMBLY",
        "GOBP_PRESYNAPSE_ORGANIZATION",
        "GOBP_PROTEIN_LOCALIZATION_TO_POSTSYNAPSE",
        "GOBP_PROTEIN_LOCALIZATION_TO_PRESYNAPSE",
        "GOBP_PROTEIN_LOCALIZATION_TO_SYNAPSE",
        "GOBP_RECEPTOR_LOCALIZATION_TO_SYNAPSE",
        "GOBP_REGULATION_OF_EXCITATORY_SYNAPSE_ASSEMBLY",
        "GOBP_REGULATION_OF_INHIBITORY_SYNAPSE_ASSEMBLY",
        "GOBP_REGULATION_OF_POSTSYNAPSE_ASSEMBLY",
        "GOBP_REGULATION_OF_POSTSYNAPSE_ORGANIZATION",
        "GOBP_REGULATION_OF_PRESYNAPSE_ORGANIZATION",
        "GOBP_REGULATION_OF_PROTEIN_LOCALIZATION_TO_SYNAPSE",
        "GOBP_REGULATION_OF_RECEPTOR_LOCALIZATION_TO_SYNAPSE",
        "GOBP_REGULATION_OF_SYNAPSE_ASSEMBLY",
        "GOBP_REGULATION_OF_SYNAPSE_MATURATION",
        "GOBP_REGULATION_OF_SYNAPSE_PRUNING",
        "GOBP_REGULATION_OF_SYNAPSE_STRUCTURAL_PLASTICITY",
        "GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY",
        "GOBP_REGULATION_OF_TRANSLATION_AT_SYNAPSE",
        "GOBP_REGULATION_PROTEIN_CATABOLIC_PROCESS_AT_SYNAPSE",
        "GOBP_SYNAPSE_ASSEMBLY",
        "GOBP_SYNAPSE_MATURATION",
        "GOBP_SYNAPSE_ORGANIZATION",
        "GOBP_SYNAPSE_PRUNING",
        "GOBP_VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE",
    ]
    mainCategories = ["Neurotransmitter", "Myelin", "Neuron", "Synapse"]
    brainProcesses = {}
    with open("processes.gmt") as processFile:
        reader = csv.reader(processFile, delimiter="\t")
        for row in reader:
            processName = row[0]
            if processName in brainCategories:
                genes = row[2:] # Skip process name and URL
                brainProcesses[processName] = genes
    impactedBrainGenes = {}
    genes = set()
    for gene in brainGenes:
        for process in brainProcesses:
            if gene in brainProcesses[process]:
                mainCategory = None
                for category in mainCategories:
                    if process.find(category.upper()) != -1:
                        mainCategory = category
                        break
                if mainCategory is None:
                    continue
                currentGeneList = impactedBrainGenes.get(mainCategory, [])
                if gene not in currentGeneList:
                    genes.add(gene)
                    impactedBrainGenes[mainCategory] = impactedBrainGenes.get(
                        mainCategory, []
                    ) + [gene]
    print(f"Found {len(genes)} genes related to brain processes.")
    with open(f"output/brain-genes-{fileName}", "w") as outputFile:
        json.dump(impactedBrainGenes, outputFile)


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