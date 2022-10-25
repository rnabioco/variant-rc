import os
from types import NoneType

import typer
from rich import print

from collections import Counter
from email.mime import base

import pandas as pd
import pysam

app = typer.Typer()

if __name__ == "variantrc.variantrc" or __name__ == "__main__":
    from tqdm import tqdm
else:
    from tqdm.notebook import tqdm  # progress bar


import matplotlib
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure, savefig

# global matplotlib figure settings
matplotlib.rcParams["figure.figsize"] = [12, 6]
matplotlib.rcParams["figure.titlesize"] = "xx-large"
matplotlib.rcParams["axes.titlesize"] = "xx-large"
matplotlib.rcParams["font.weight"] = "bold"
matplotlib.rcParams["axes.titleweight"] = "bold"
matplotlib.rcParams["axes.labelweight"] = "bold"

# feature ideas
# choose the number of strains or which strains to use

### Read Classification Resolver Function
def read_resolver(
    list_read_results,
    read_result,
    list_variant_posqual,
    list_variant_confidence,
    mapq_threshold=20,
    base_quality_threshold=9,
):
    list_possible_results = list_read_results.copy()

    if (
        len(set(list_possible_results)) == 1
    ):  # if all reads are the same, return the first one
        read_result = list_possible_results[0]

    else:
        while (
            len(set(list_possible_results)) > 1
            and read_result != "mixed post filter"
            and read_result != "not determined"
            and read_result != "strain1"
            and read_result != "strain2"
            and read_result != "strain3"
        ):
            # print( "1", list_possible_results)

            # deal with low quality positions
            low_quality = "base at variant low quality " + str(base_quality_threshold)
            if low_quality in list_possible_results:
                if list_possible_results.count(low_quality) < len(
                    list_possible_results
                ):
                    list_possible_results = [
                        i for i in list_possible_results if i != low_quality
                    ]

            # deal with none as if it passed the initial scoring without getting a variant result
            if "incomplete indel match" in list_possible_results:
                if list_possible_results.count("incomplete indel match") < len(
                    list_possible_results
                ):
                    list_possible_results = [
                        i
                        for i in list_possible_results
                        if i != "incomplete indel match"
                    ]

            # deal with none as if it passed the initial scoring without getting a variant result
            if "none" in list_possible_results:
                if list_possible_results.count("none") < len(list_possible_results):
                    list_possible_results = [
                        i for i in list_possible_results if i != "none"
                    ]

            # deal with unknown
            if "unknown" in list_possible_results:
                if list_possible_results.count("unknown") < len(list_possible_results):
                    list_possible_results = [
                        i for i in list_possible_results if i != "unknown"
                    ]

            # deal with variant not in read
            if "variant not in read" in list_possible_results:
                if list_possible_results.count("variant not in read") < len(
                    list_possible_results
                ):
                    list_possible_results = [
                        i for i in list_possible_results if i != "variant not in read"
                    ]

            # deal with read base not in variant alleles
            if "base not variant option" in list_possible_results:
                if list_possible_results.count("base not variant option") < len(
                    list_possible_results
                ):
                    list_possible_results = [
                        i
                        for i in list_possible_results
                        if i != "base not variant option"
                    ]

            # deal with a variant that is not determinate of strain
            if "mixed snp locus" in list_possible_results:
                if list_possible_results.count("mixed locus") < len(
                    list_possible_results
                ):
                    list_possible_results = [
                        i for i in list_possible_results if i != "mixed locus"
                    ]

            # deal with a variant that is not determinate of strain
            if "mixed locus indel" in list_possible_results:
                if list_possible_results.count("mixed locus indel") < len(
                    list_possible_results
                ):
                    list_possible_results = [
                        i for i in list_possible_results if i != "mixed locus indel"
                    ]

            if "mixed locus indel 1+2" in list_possible_results:
                if list_possible_results.count("mixed locus indel 1+2") < len(
                    list_possible_results
                ):
                    list_possible_results = [
                        i for i in list_possible_results if i != "mixed locus indel 1+2"
                    ]

            if "mixed locus indel 1+3" in list_possible_results:
                if list_possible_results.count("mixed locus indel 1+3") < len(
                    list_possible_results
                ):
                    list_possible_results = [
                        i for i in list_possible_results if i != "mixed locus indel 1+3"
                    ]

            if "mixed locus indel 2+3" in list_possible_results:
                if list_possible_results.count("mixed locus indel 2+3") < len(
                    list_possible_results
                ):
                    list_possible_results = [
                        i for i in list_possible_results if i != "mixed locus indel 2+3"
                    ]

            # resolve multiple possible strain calls
            if len(set(list_possible_results)) >= 2:

                num_strain1 = list_possible_results.count("strain1")
                num_strain2 = list_possible_results.count("strain2")
                num_strain3 = list_possible_results.count("strain3")
                num_reference = list_possible_results.count("reference")
                per_strain1 = num_strain1 / (
                    num_strain1 + num_strain2 + num_strain3 + num_reference
                )
                per_strain2 = num_strain2 / (
                    num_strain1 + num_strain2 + num_strain3 + num_reference
                )
                per_strain3 = num_strain3 / (
                    num_strain1 + num_strain2 + num_strain3 + num_reference
                )
                per_reference = num_reference / (
                    num_strain1 + num_strain2 + num_strain3 + num_reference
                )

                # if any of the strains are at or over 50% of the variants
                if (
                    per_strain1 >= 0.5
                    or per_strain2 >= 0.5
                    or per_strain3 >= 0.5
                    or per_reference >= 0.5
                ):
                    if per_strain1 >= 0.5 and per_strain2 >= 0.5:
                        read_result = "mixed strain"
                        break
                    elif per_strain1 >= 0.5 and per_strain3 >= 0.5:
                        read_result = "mixed strain"
                        break
                    elif per_strain2 >= 0.5 and per_strain3 >= 0.5:
                        read_result = "mixed strain"
                        break
                    elif per_strain1 >= 0.5 and per_reference >= 0.5:
                        read_result = "mixed strain"
                        break
                    elif per_strain2 >= 0.5 and per_reference >= 0.5:
                        read_result = "mixed strain"
                        break
                    elif per_strain3 >= 0.5 and per_reference >= 0.5:
                        read_result = "mixed strain"
                        break
                    elif per_strain1 >= 0.5:
                        read_result = "strain1"
                        break
                    elif per_strain2 >= 0.5:
                        read_result = "strain2"
                        break
                    elif per_strain3 >= 0.5:
                        read_result = "strain3"
                        break
                    elif per_reference >= 0.5:
                        read_result = "reference"
                        break
                else:
                    # if set(list_possible_results):
                    read_result = "mixed strain"
                    break

            if len(set(list_possible_results)) == 0:
                read_result = "could not determine from variants"
                break

            if len(set(list_possible_results)) == 1:
                read_result = list_possible_results[0]
                break

            if len(set(list_possible_results)) > 1:
                read_result = "too many options"
                break

    return read_result


### Read Processing Function
@app.command()
def read_classification(
    bamfile_path: str,
    variantfile_path: str,
    # dict_gene_lookup: NoneType = None,
    mapq_threshold: int = 20,
    base_quality_threshold: int = 9,
    output_path: str = os.getcwd() + "/classification_results/",
    save_figures: bool = True,
    save_variants: bool = True,
):
    """Classify reads into strain specific bins based on the presence distinguishing variants."""

    ### Load Data
    # load reads
    print("\nLoading reads to classify..." + bamfile_path)
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    bamfile_name = bamfile_path.split("/")[-1].split(".")[0]

    # load variants
    print("Loading variants..." + variantfile_path)

    variantfile = pysam.VariantFile(variantfile_path)
    num_vcf_samples = len(str(variantfile.header).split("FORMAT\t")[1].split("\t"))
    vcf_sample_names = str(variantfile.header).split("FORMAT\t")[1].split("\t")

    print("\n")
    print("Number of samples in VCF: " + str(num_vcf_samples))
    for num, sample in enumerate(vcf_sample_names):
        print(f"Sample {num+1}: {sample}")

    ### Session level stats
    df_summary = pd.DataFrame()
    sorted_reads = []  # processed reads
    read_count = 0  # number of reads

    print("\nClassifying " + bamfile_name + "...")

    ### Main read parsing loop
    for read in tqdm(bamfile.fetch()):
        read_count += 1
        variant_num = 0  # reset variant count

        # reset read level stats
        read_result = ""  # reset read result
        result_confidence = -1  # confidence of read result
        read_sequence = ""  # reset read sequence
        genotype1 = ""  # reset genotype1
        genotype2 = ""  # reset genotype2
        genotype3 = ""  # reset genotype3

        ### Read level data
        list_variant_refseq = []
        list_read_results = []
        list_variant_alleles = []
        list_genotype1 = []
        list_genotype2 = []
        list_genotype3 = []
        list_isindel = []
        list_genes = []
        list_symbols = []
        list_variantpos = []
        list_variant_posqual = []
        list_variant_confidence = []

        if (
            int(read.mapping_quality) > mapq_threshold
        ):  # filter for overall read map quality # check presecense of read.qual or read.get_forward_qualities()?

            # generate the read level lookup dictionaries for base position and quality
            read_index = {}  # read base lookup dictionary
            mapping_quality = {}  # read base level mapping quality lookup dictionary

            for index, position in read.aligned_pairs:
                if index is not None and position is not None:
                    read_index[position] = read.seq[int(index)]
                    mapping_quality[position] = read.get_forward_qualities()[index]

            # fetch 1-base position variants found within the read sequence
            for variant in variantfile.fetch(
                contig=read.reference_name,
                start=read.reference_start,
                stop=read.reference_end,
            ):

                variant_position = variant.pos - 1  # fix the off by one error
                variant_result = "none"  # reset result
                variant_confidence = 0

                if "INDEL" in variant.info.keys():
                    is_indel = True
                else:
                    is_indel = False

                if "GENE" in variant.info.keys():
                    list_genes.append(variant.info["GENE"])

                    try:
                        list_symbols.append(
                            dict_gene_lookup[variant.info["GENE"]]
                        )  # translate gene name to symbol
                    except KeyError:
                        list_symbols.append(
                            variant.info["GENE"]
                        )  # if no symbol found, use gene name
                else:
                    list_genes.append("-")  # if no gene found in key(), use '-'
                    list_symbols.append("-")

                if variant_position in read_index.keys():
                    variant_num += 1
                    list_variant_posqual.append(
                        mapping_quality[variant_position]
                    )  # keep track of the read base quality

                    if (
                        mapping_quality[variant_position] > base_quality_threshold
                    ):  # check base quality at variant
                        # print(len(variant.samples.items()))
                        genotype1 = variant.samples.items()[0][1].values()[
                            0
                        ]  # get the genotype for the first sample
                        genotype2 = variant.samples.items()[1][1].values()[
                            0
                        ]  # get the genotype for the second sample

                        if num_vcf_samples == 3:
                            genotype3 = variant.samples.items()[2][1].values()[
                                0
                            ]  # get the genotype for the third sample
                        else:
                            genotype3 = (
                                9,
                                9,
                            )  # dummy int interator to parse through the genotype

                        # build allele/genotype lookup dictionaries
                        dict_alleles = {}  # [0:allele]
                        dict_genotypes = {}  # [allele:0]
                        allele_num = 0

                        for (
                            allele
                        ) in (
                            variant.alleles
                        ):  # iterate through alleles and build allele/genotype lookup dictionaries
                            dict_alleles[allele_num] = allele
                            dict_genotypes[allele] = allele_num
                            allele_num += 1

                        ## Variant to Read sequence matching and Variant Result Calling

                        ## handle INDELS first
                        list_indel_walk_positions = []
                        if is_indel == True:
                            list_index = list(
                                read_index.keys()
                            )  # get the read base position index
                            start_pos = list_index.index(variant_position)

                            list_possible_matches = []

                            for allele in dict_alleles.values():
                                try:
                                    list_indel_walk_positions = [
                                        list_index[start_pos + i]
                                        for i, x in enumerate(list(allele))
                                    ]  # construct new base positions for the indel based on allele length
                                except:
                                    pass

                                list_possible_matches.append(
                                    "".join(
                                        [
                                            read_index[i]
                                            for i in list_indel_walk_positions
                                        ]
                                    )
                                )  # add to list of possible indel matches

                            for indel_construct in sorted(
                                list_possible_matches, key=len, reverse=True
                            ):  # check in descending order of indel length
                                read_sequence = "-"  # reset read sequence

                                if indel_construct in dict_genotypes.keys():
                                    match_result = dict_genotypes[
                                        indel_construct
                                    ]  # returns the allele number for the match (0,1,2,3) etc
                                    variant_confidence = 0.25 * len(
                                        indel_construct
                                    )  # set variant confidence to 10x the length of the indel

                                    if (
                                        match_result in genotype1
                                        and match_result in genotype2
                                    ):
                                        variant_result = "mixed locus indel 1+2"
                                        read_sequence = indel_construct
                                        break

                                    if (
                                        match_result in genotype1
                                        and match_result in genotype3
                                    ):
                                        variant_result = "mixed locus indel 1+3"
                                        read_sequence = indel_construct
                                        break

                                    if (
                                        match_result in genotype2
                                        and match_result in genotype3
                                    ):
                                        variant_result = "mixed locus indel 2+3"
                                        read_sequence = indel_construct
                                        break

                                    elif match_result in genotype1:
                                        variant_result = "strain1"
                                        read_sequence = indel_construct
                                        break

                                    elif match_result in genotype2:
                                        variant_result = "strain2"
                                        read_sequence = indel_construct
                                        break

                                    elif match_result in genotype3:
                                        variant_result = "strain3"
                                        read_sequence = indel_construct
                                        break

                                    elif dict_alleles[match_result] in variant.ref:
                                        variant_result = "reference"
                                        read_sequence = indel_construct
                                        break

                                    else:
                                        variant_result = "unknown"
                                        read_sequence = indel_construct
                                        break
                                else:
                                    variant_result = "incomplete indel match"  # if no match found, set variant result to incomplete match

                        else:  # handle SNPs

                            read_sequence = read_index[variant_position]

                            if read_sequence in dict_genotypes.keys():
                                base_genotype_at_variant = dict_genotypes[
                                    read_sequence
                                ]  # get the genotype for the read base

                                if (
                                    base_genotype_at_variant in genotype1
                                    and base_genotype_at_variant in genotype2
                                ):
                                    variant_result = "mixed snp locus"
                                    variant_confidence += 0.25

                                elif base_genotype_at_variant in genotype1:
                                    variant_result = "strain1"
                                    variant_confidence += 0.25

                                elif base_genotype_at_variant in genotype2:
                                    variant_result = "strain2"
                                    variant_confidence += 0.25

                                elif base_genotype_at_variant in genotype3:
                                    variant_result = "strain3"
                                    variant_confidence += 0.25

                                elif (
                                    dict_alleles[base_genotype_at_variant]
                                    in variant.ref
                                ):
                                    variant_result = "reference"
                                    variant_confidence += 0.25

                                else:
                                    variant_result = "unknown"
                            else:
                                variant_result = "base not variant option"
                    else:
                        variant_result = "base at variant low quality " + str(
                            base_quality_threshold
                        )

                    # append variant level results
                    list_read_results.append(variant_result)
                    list_variant_confidence.append(variant_confidence)
                    list_variant_refseq.append(read_sequence)
                    list_variant_alleles.append(variant.alleles)
                    list_genotype1.append(genotype1)
                    list_genotype2.append(genotype2)
                    list_genotype3.append(genotype3)
                    list_isindel.append(is_indel)
                    list_variantpos.append(variant_position)

                    if variant_result == "none":
                        print(
                            variant_result,
                            read.reference_name,
                            variant_position,
                            read_sequence,
                            variant.alleles,
                            genotype1,
                            genotype2,
                            genotype3,
                            is_indel,
                            list_possible_matches,
                        )

                else:
                    variant_result = "variant not in read"

        else:
            read_result = "read failed MAPQ " + str(mapq_threshold)

        ### Resolve Read Result
        if variant_num > 0:
            read_result = read_resolver(
                list_read_results,
                read_result,
                list_variant_posqual,
                list_variant_confidence,
                mapq_threshold,
                base_quality_threshold,
            )

        elif variant_num == 0 and read_result == "":
            read_result = "no variants in read"

        elif variant_num > 0 and read_result == "":
            print(
                read_result,
                list_read_results,
                list_variant_refseq,
                list_variant_alleles,
            )

        # calculate overall read score
        result_confidence = sum(list_variant_confidence)

        ## Split the reads into multiple bam files based on the read result

        ### Append results to dataframe
        df2 = pd.DataFrame(
            {
                "read_count": read_count,
                "read_id": read.query_name,
                "chrom": read.reference_name,
                "read_pos_start": read.reference_start,
                "read_pos_end": read.reference_end,
                "read_length": read.infer_query_length(),
                "read_mapq": read.mapping_quality,
                "read_result": read_result,
                "result_confidence": result_confidence,
                "genes": str(set(list_genes)),
                "gene_symbols": str(set(list_symbols)),
                "variant_num": variant_num,
                "variant_results": str(list_read_results),
                "variant_positions": str(list_variantpos),
                "variant_posqual": str(list_variant_posqual),
                "read_seq": str(list_variant_refseq),
                "variant_alleles": str(list_variant_alleles),
                "variant_isindel": str(list_isindel),
                "geneotype1": str(list_genotype1),
                "genotype2": str(list_genotype2),
                "genotype3": str(list_genotype3),
            },
            index=[0],
        )

        df_summary = pd.concat((df_summary, df2), ignore_index=True)

    bamfile.close()
    variantfile.close()

    print("\nTotal reads classified:", sum(df_summary.read_result.value_counts()))
    print("\nBins:")
    print(df_summary.read_result.value_counts())

    ### Outputs
    try:
        os.mkdir(output_path)
    except:
        pass

    ### Save to CSV
    if save_variants:
        print("\nSaving..." + bamfile_name + ".csv")
        df_summary.to_csv(output_path + bamfile_name + ".csv")  # save to CSV

    if save_figures:
        print("\nSaving..." + bamfile_name + ".png")
        plot = df_summary.read_result.value_counts().plot(
            kind="bar", title=bamfile_name + " Read Results"
        )
        plot.figure.savefig(
            output_path + bamfile_name + ".png", bbox_inches="tight"
        )  # save to PNG


if __name__ == "__main__":

    variantfile_path = "vcfs/variant_calls_3strains_annotated_exons.vcf.gz"

    bamfile_path1 = "bams/cantons_sorted_subset.bam"
    bamfile_path2 = "bams/orer_sorted_subset.bam"
    bamfile_path3 = "bams/w1118_0-1replicate1_040121_subset.bam"

    bamfile = bamfile_path3
    bamfile_name = bamfile.split("/")[-1].split(".")[0]

    ## gene database lookup
    genedb = pd.read_csv("refs/pantherGeneList.txt", sep="\t", header=None)
    dict_gene_lookup = {}
    for gene in genedb.iterrows():
        dict_gene_lookup[gene[1][0].split("|")[1].split("=")[1]] = gene[1][1].strip(
            ";ortholog"
        )

    ### Read Sorting Function Testing
    read_sorting(bamfile, variantfile_path, dict_gene_lookup)

    df_readcalls = pd.read_csv(bamfile + ".csv")

    print("\nTotal reads sorted:", sum(df_readcalls.read_result.value_counts()))
    print("\nBins:")
    print(df_readcalls.read_result.value_counts())

    df_readcalls.read_result.value_counts().plot(
        kind="bar", title=bamfile_name + " Read Results"
    )
    plt.xticks(rotation=45)
    plt.savefig(bamfile_name + "_read_results.png")
    plt.show()
