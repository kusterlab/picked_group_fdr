import pandas as pd


"""
Create a peptides.txt from an evidence.txt file.
We sum the intensities of the PSMs per Modified sequence for all PSMs which were used in quantifications

USAGE: python create_peptides_txt  -e evidence_out_filtered.txt -p final_proteinGroups_filtered.txt -o peptides.txt
"""


def main(evidence_df_file: str, protein_group_file: str, peptide_df_file: str):
    """
    Executes the peptide-to-protein correlation analysis workflow.

    This function performs the following steps:
      1. Filters the evidence DataFrame to retain only entries referenced in the protein group file.
      2. Aggregates the filtered evidence data to the peptide level, calculating total intensities and PSM counts.
      3. Computes and prints the correlation between peptide-level summed intensities and protein-level iBAQ values.

    Args:
        evidence_df_file (str): Path to the evidence file (TSV format) containing proteomics evidence data.
        protein_group_file (str): Path to the protein group file (TSV format) containing protein-level data including iBAQ values.

    Returns:
        pandas.DataFrame: A peptide-level DataFrame with columns including:
            - 'Sequence'
            - 'Modified sequence'
            - 'Proteins'
            - 'Experiment'
            - 'Num_PSMs'
            - 'Sum_intensity'
    """

    filtered_evidence_df = make_filtered_evidence_df(
        evidence_df_file, protein_group_file
    )
    peptide_df = make_peptide_df(filtered_evidence_df)
    calculate_correlation_between_peptides_vs_protein_intensities(
        protein_group_file, peptide_df
    )
    peptide_df.to_csv(peptide_df_file)


def make_filtered_evidence_df(evidence_df_file: str, protein_group_file: str):
    """
    Filters an evidence DataFrame to include only rows that are referenced by a protein group file.

    This function reads an evidence DataFrame and a protein group DataFrame from tab-separated files.
    It extracts the unique evidence IDs from the 'Evidence IDs' column of the protein group file,
    flattens them into a single list, and filters the evidence DataFrame to include only rows
    whose 'id' values are present in that list.

    Args:
        evidence_df_file (str): Path to the evidence DataFrame file (TSV format) at 1% FDR level.
        protein_group_file (str): Path to the protein group DataFrame file (TSV format).

    Returns:
        pandas.DataFrame: A filtered version of the evidence DataFrame containing only rows
        that match the evidence IDs referenced in the protein group file.
    """

    evidence_df = read_tsv(evidence_df_file)  # evidence df at 1% FDR level
    protein_group_df = read_tsv(protein_group_file)
    list_evidence_ids = protein_group_df["Evidence IDs"].unique().tolist()
    nested_evidence_ids = [x.split(";") for x in list_evidence_ids]
    all_evidence_ids = list(set(sum(nested_evidence_ids, start=[])))
    evidence_df["ID"] = evidence_df["id"].astype(str)
    return evidence_df[evidence_df["ID"].isin(all_evidence_ids)]


def read_tsv(file: str):
    return pd.read_csv(file, sep="\t")


def make_peptide_df(filtered_evidence_df: pd.DataFrame):
    """
    Aggregates peptide-level data from a filtered evidence DataFrame.

    This function computes the number of peptide-spectrum matches (PSMs) and the total intensity
    for each modified peptide sequence. It returns a DataFrame with selected columns,
    removing duplicate entries based on 'Modified sequence', 'Proteins', and 'Experiment'.

    Args:
        filtered_evidence_df (pandas.DataFrame): A DataFrame containing evidence-level proteomics
            data, including 'Modified sequence', 'Intensity', 'Sequence', 'Proteins', and 'Experiment'.

    Returns:
        pandas.DataFrame: A peptide-level DataFrame with the following columns:
            - 'Sequence': Unmodified peptide sequence.
            - 'Modified sequence': Modified peptide sequence used for grouping.
            - 'Proteins': Associated protein(s) for the peptide.
            - 'Experiment': Experiment label or identifier.
            - 'Num_PSMs': Number of PSMs per modified sequence.
            - 'Sum_intensity': Total intensity summed over all PSMs per modified sequence.
    """

    filtered_evidence_df["Num_PSMs"] = filtered_evidence_df.groupby(
        "Modified sequence"
    ).transform("size")
    filtered_evidence_df["Sum_intensity"] = filtered_evidence_df.groupby(
        "Modified sequence"
    )["Intensity"].transform("sum")
    return filtered_evidence_df[
        [
            "Sequence",
            "Modified sequence",
            "Proteins",
            "Experiment",
            "Num_PSMs",
            "Sum_intensity",
        ]
    ].drop_duplicates(subset=["Modified sequence", "Proteins", "Experiment"])


def unnest_df(df: pd.DataFrame, delimiter=";") -> pd.DataFrame:
    """
    Unnests a DataFrame by splitting index values using a delimiter and expanding them into multiple rows.

    This function assumes that the DataFrame's index contains delimited strings (e.g., "A;B;C").
    It splits each index entry by the given delimiter, creates multiple rows (one per split value),
    and sets the resulting values as the new index.

    Args:
        df (pandas.DataFrame): The input DataFrame with a string-based index that includes delimited values.
        delimiter (str, optional): The delimiter used to split index values. Defaults to ';'.

    Returns:
        pandas.DataFrame: A new DataFrame where each split index value is unnested into its own row.

    Example:
        Input index: ['P1;P2', 'P3']
        Output index: ['P1', 'P2', 'P3']
    """

    temp_df = df
    temp_df["index"] = temp_df.index.str.split(delimiter)
    temp_df = temp_df.explode("index")
    temp_df = temp_df.set_index("index")
    return temp_df


def calculate_correlation_between_peptides_vs_protein_intensities(
    protein_group_file: str, peptide_df: pd.DataFrame
) -> None:
    """
    Calculates and prints the Pearson correlation coefficient between peptide-level summed intensities
    and protein-level iBAQ intensities.

    This function reads a protein group file, expands the protein group entries using `unnest_df`, and
    merges it with a peptide DataFrame on protein identifiers. It then computes the correlation between
    the `Sum_intensity` values from peptides and the corresponding protein `iBAQ` intensities.

    Args:
        protein_group_file (str): Path to the protein group TSV file containing protein-level iBAQ values.
        peptide_df (pandas.DataFrame): A DataFrame containing peptide-level data, including:
            - 'Proteins': protein IDs corresponding to each peptide,
            - 'Sum_intensity': summed intensity per modified sequence.

    Returns:
        None: The function prints the correlation coefficient to standard output.

    Prints:
        Correlation between peptide summed intensities and protein iBAQ intensities.
    """

    proteind_group_df = read_tsv(protein_group_file)
    merged_df = unnest_df(proteind_group_df.set_index("Protein IDs")).merge(
        peptide_df, right_on=["Proteins"], left_index=True
    )
    correlation_coef = merged_df["iBAQ"].corr(merged_df["Sum_intensity"])
    print(
        f"Correlation between peptides sum_intensities and protein iBAQ {correlation_coef}"
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate a peptides.txt file in MaxQuant format."
    )
    parser.add_argument(
        "-e", "--evidence", required=True, help="Path to the evidence TSV file"
    )
    parser.add_argument(
        "-p", "--protein", required=True, help="Path to the protein group TSV file"
    )

    parser.add_argument(
        "-o", "--output", required=True, help="Path to the peptides.txt TSV file"
    )
    args = parser.parse_args()

    # Run the main analysis
    main(args.evidence, args.protein, args.output)
