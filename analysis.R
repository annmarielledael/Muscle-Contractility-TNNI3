
# Molecular Mechanisms of Muscle Contractility
# Student Assignment 15: TNNI3

# Load required packages
library(rentrez)
library(ggplot2)

# Assignment details
gene_symbol <- "TNNI3"
transcript_accession <- "NM_000363.5"
clinvar_id <- "12424"
variant_protein <- "p.Arg192His"

# Create project folders if they do not exist
for (folder in c("data", "results", "figures")) {
  dir.create(folder, showWarnings = FALSE)
}

# Display assignment details
cat("Gene:", gene_symbol, "\n")
cat("Transcript:", transcript_accession, "\n")
cat("ClinVar Variation ID:", clinvar_id, "\n")
cat("Protein variant:", variant_protein, "\n")

# Step 6: Retrieve the human TNNI3 gene record from NCBI

hits <- entrez_search(
  db = "gene",
  term = paste0(gene_symbol, "[Gene Name] AND txid9606[Organism]"),
  retmax = 20
)

if (length(hits$ids) == 0) {
  stop("Gene not found in NCBI")
}

gene_id <- NA_character_

for (id in hits$ids) {
  rec <- entrez_summary(db = "gene", id = id)
  
  if (toupper(rec$name) == toupper(gene_symbol)) {
    gene_id <- id
    gene_info <- rec
    break
  }
}

if (is.na(gene_id)) {
  stop("Exact human TNNI3 gene not found")
}

# Display the retrieved record
print(gene_info)

# Save the gene summary
writeLines(
  c(
    paste("Gene symbol:", gene_info$name),
    paste("NCBI Gene ID:", gene_id),
    paste("Official name:", gene_info$description),
    paste("Chromosome:", gene_info$chromosome),
    "Normal role in muscle: TNNI3 encodes cardiac troponin I, an inhibitory subunit of the cardiac troponin complex that helps regulate heart muscle contraction and relaxation in response to calcium.",
    "Source: https://www.ncbi.nlm.nih.gov/gene/"
  ),
  "results/gene_summary.txt"
)

cat("NCBI Gene ID:", gene_id, "\n")
cat("Gene summary saved to results/gene_summary.txt\n")


# Step 7: Retrieve the matching reference protein sequence

gb <- entrez_fetch(
  db = "nuccore",
  id = transcript_accession,
  rettype = "gb",
  retmode = "text"
)

# Extract the protein accession from the transcript record
pat <- '/protein_id="NP_[0-9]+[.][0-9]+"'
fragment <- regmatches(gb, regexpr(pat, gb))

if (length(fragment) == 0 || !nzchar(fragment)) {
  stop("NP protein ID not found; inspect the transcript record")
}

protein_accession <- gsub(
  '^/protein_id="|"$',
  "",
  fragment
)

# Download the reference protein FASTA
protein_fasta <- entrez_fetch(
  db = "protein",
  id = protein_accession,
  rettype = "fasta",
  retmode = "text"
)

if (!startsWith(protein_fasta, ">")) {
  stop("Invalid FASTA response")
}

# Save the sequence
fasta_file <- file.path(
  "data",
  paste0(gene_symbol, "_protein.fasta")
)

writeLines(protein_fasta, fasta_file)

# Save the protein accession
writeLines(
  c(
    paste("Transcript accession:", transcript_accession),
    paste("NP accession:", protein_accession)
  ),
  "results/protein_accession.txt"
)

cat("Matched protein:", protein_accession, "\n")
cat("Protein FASTA saved to:", fasta_file, "\n")


# Step 3: Calculate amino acid composition

fasta_lines <- readLines("data/TNNI3_protein.fasta")

# Remove the FASTA header and combine the sequence
protein_sequence <- paste(
  fasta_lines[!startsWith(fasta_lines, ">")],
  collapse = ""
)

protein_sequence <- toupper(gsub("[[:space:]]", "", protein_sequence))

# Count the amino acids
amino_acids <- c(
  "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
  "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"
)

aa_counts <- table(
  factor(strsplit(protein_sequence, "")[[1]],
         levels = amino_acids)
)

aa_results <- data.frame(
  Amino_Acid = names(aa_counts),
  Count = as.integer(aa_counts),
  Frequency_Percent = round(
    as.integer(aa_counts) / nchar(protein_sequence) * 100,
    2
  )
)

protein_length <- nchar(protein_sequence)

write.csv(
  aa_results,
  "results/amino_acid_composition.csv",
  row.names = FALSE
)

cat("Protein length:", protein_length, "amino acids\n")
cat("Top 3 amino acids:\n")
print(head(aa_results[order(-aa_results$Count), ], 3))
cat("Amino acid composition saved to results/amino_acid_composition.csv\n")


# Step 4: Plot amino acid composition

aa_plot <- ggplot(
  aa_results,
  aes(x = reorder(Amino_Acid, -Count), y = Count)
) +
  geom_col(fill = "steelblue") +
  labs(
    title = "Amino Acid Composition of TNNI3",
    x = "Amino Acid (one-letter code)",
    y = "Count"
  ) +
  theme_minimal()

ggsave(
  "figures/01_amino_acid_composition.png",
  plot = aa_plot,
  width = 8,
  height = 5,
  dpi = 300
)

print(aa_plot)

cat("Graph saved to figures/01_amino_acid_composition.png\n")


# Step 5: Plot the variant position

variant_position <- 192

variant_data <- data.frame(
  Start = 1,
  End = protein_length,
  Variant = variant_position
)

variant_plot <- ggplot(variant_data) +
  geom_segment(
    aes(x = Start, xend = End, y = 1, yend = 1),
    linewidth = 6,
    color = "grey75"
  ) +
  geom_segment(
    aes(x = Variant, xend = Variant, y = 0.8, yend = 1.2),
    linewidth = 1.2,
    color = "red"
  ) +
  annotate(
    "text",
    x = variant_position,
    y = 1.35,
    label = "Arg192His",
    color = "red",
    hjust = 1
  ) +
  scale_x_continuous(
    limits = c(1, protein_length),
    breaks = seq(0, protein_length, by = 30)
  ) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(
    title = "TNNI3 Protein Variant Position",
    subtitle = paste("Protein length:", protein_length, "amino acids"),
    x = "Amino acid position"
  ) +
  theme_minimal()

ggsave(
  "figures/02_variant_position.png",
  plot = variant_plot,
  width = 8,
  height = 3.5,
  dpi = 300
)

print(variant_plot)

cat("Variant diagram saved to figures/02_variant_position.png\n")


# Step 6: Save the ClinVar variant summary

clinvar_notes <- c(
  "CLINVAR VARIANT SUMMARY",
  "=======================",
  "",
  "Gene: TNNI3",
  "Transcript: NM_000363.5",
  "DNA change: c.575G>A",
  "Protein change: p.Arg192His (p.R192H)",
  "Protein accession: NP_000354.4",
  "Variant type: Missense",
  "ClinVar Variation ID: 12424",
  "ClinVar accession: VCV000012424.31",
  "dbSNP: rs104894729",
  "Inheritance reported by expert panel: Autosomal dominant",
  "",
  "OVERALL AND EXPERT-PANEL INTERPRETATION",
  "Clinical significance: Pathogenic",
  "Condition: Hypertrophic cardiomyopathy",
  "Expert panel: ClinGen Cardiomyopathy Variant Curation Expert Panel",
  "Review status: 3 of 4 stars",
  "Last evaluated: November 13, 2025",
  "First submitted to ClinVar: April 4, 2013",
  "Last submission shown: September 22, 2026",
  "",
  "CONDITION-LEVEL SUBMISSIONS",
  "Familial restrictive cardiomyopathy 1: Pathogenic; 3 submissions;",
  "  2 stars; last evaluated February 14, 2023.",
  "Hypertrophic cardiomyopathy and restrictive cardiomyopathy:",
  "  Pathogenic; 1 submission; 1 star; last evaluated May 2, 2019.",
  "Condition not provided: Pathogenic; 4 submissions; 2 stars;",
  "  last evaluated June 12, 2025.",
  "Cardiovascular phenotype: Pathogenic; 1 submission; 1 star;",
  "  last evaluated December 18, 2017.",
  "Hypertrophic cardiomyopathy: Pathogenic; 2 submissions;",
  "  3 stars; last evaluated November 13, 2025.",
  "Restrictive cardiomyopathy: Pathogenic; 1 submission; 1 star;",
  "  last evaluated January 3, 2018.",
  "Sudden infant death syndrome: Pathogenic; 1 submission; 1 star;",
  "  last evaluated August 8, 2022.",
  "Dilated cardiomyopathy 2A: Pathogenic; 1 submission; 1 star;",
  "  evaluation date not provided in the supplied notes.",
  "TNNI3-related disorder: Pathogenic; 1 submission; 1 star;",
  "  last evaluated December 8, 2023.",
  "",
  "EXPERT-PANEL EVIDENCE SUMMARY",
  "The supplied summary reports clinical observations in affected",
  "individuals, de novo occurrence, rarity in population databases,",
  "functional studies, regional evidence, and computational predictions.",
  "Reported ACMG/AMP criteria: PS4_Supporting, PM2_Supporting,",
  "PS2, PM1, PS3_Moderate, and PP3.",
  "Functional research discussed changes in calcium sensitivity,",
  "thin-filament regulation, and cardiac muscle relaxation.",
  "",
  "SELECTED PUBLICATIONS REPORTED IN THE SUPPLIED SUMMARY",
  "Mogensen J et al. J Clin Invest. 2003. PMID: 12531876.",
  "Gomes AV et al. J Biol Chem. 2005. PMID: 15961398.",
  "Kobayashi T et al. J Biol Chem. 2006. PMID: 16531415.",
  "Du J et al. Arch Biochem Biophys. 2006. PMID: 17027633.",
  "Davis J et al. Circ Res. 2007. PMID: 17463320.",
  "Fujino M et al. J Cardiol Cases. 2018. PMID: 30279906.",
  "",
  "LIMITATIONS",
  "Condition-specific submissions differ in review status and evidence.",
  "Not every listed publication necessarily studies this exact variant;",
  "some discuss related variants or troponin biology more broadly.",
  "The supplied notes report little evidence for haploinsufficiency",
  "pathogenicity and no evidence available for triplosensitivity.",
  "",
  "Source: https://www.ncbi.nlm.nih.gov/clinvar/variation/12424/",
  "Access date: October 27, 2026",
  "Note: Summary transcribed from the ClinVar information supplied",
  "for this assignment; verify details against the live record."
)

writeLines(
  clinvar_notes,
  "results/clinvar_summary.txt"
)

cat("ClinVar summary saved to results/clinvar_summary.txt\n")
