# Ankitashree
# Central Dogma + GC Content + Pairwise Alignments (Global vs Local) + BLASTn (Top 3 alignments)
# Using only Nucleotide Accession ID

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align import PairwiseAligner
import urllib.error

# --- Step 1: Setup ---
Entrez.email = "ankitashree@example.com"   # your email
acc_id = "NM_001376310"                    # nucleotide accession number

# --- Step 2: Fetch Nucleotide Sequence ---
print("\n--- Fetching Nucleotide Sequence ---")
handle = Entrez.efetch(db="nucleotide", id=acc_id, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

seq = record.seq
print("Accession:", record.id)
print("Description:", record.description)
print("Sequence length:", len(seq))

# --- Step 3: GC Content ---
gc_content = round((seq.count("G") + seq.count("C")) / len(seq) * 100, 2)
print("\n--- GC Content ---")
print("GC Content:", gc_content, "%")

# --- Step 4: Central Dogma (CDS → Protein) ---
cds_seq, protein = None, None
for feature in record.features:
    if feature.type == "CDS":
        cds_seq = feature.extract(seq)
        protein = cds_seq.translate(to_stop=True)
        break

if cds_seq:
    print("\n--- Central Dogma ---")
    print("CDS length:", len(cds_seq))
    print("Protein length:", len(protein))
    print("Protein (first 30 aa):", protein[:30])
else:
    print("No CDS found in record.")

# --- Step 5: Global vs Local Alignment ---
def show_alignment(seq1, seq2):
    aligner = PairwiseAligner()

    # Global Alignment
    aligner.mode = "global"
    global_alignments = aligner.align(seq1, seq2)
    print("\n=== Global Alignment ===")
    print(global_alignments[0].format())
    print("Score:", global_alignments[0].score)

    # Local Alignment
    aligner.mode = "local"
    local_alignments = aligner.align(seq1, seq2)
    print("\n=== Local Alignment ===")
    print(local_alignments[0].format())
    print("Score:", local_alignments[0].score)

if cds_seq:
    # Use first 100 bases to keep output manageable
    show_alignment(cds_seq[:100], cds_seq.reverse_complement()[:100])
else:
    print("No CDS available for alignment.")

# --- Step 6: BLASTn (CDS vs nucleotide database, Top 3 Alignments) ---
try:
    if cds_seq:
        print("\n--- Running BLASTn (Top 3 Alignments) ---")
        blastn_handle = NCBIWWW.qblast("blastn", "nt", cds_seq[:200])  # shorter query for demo
        blastn_record = NCBIXML.read(blastn_handle)

        print("\nTop 3 BLASTn Hits:")
        for alignment in blastn_record.alignments[:3]:
            print(f"\nSequence: {alignment.title}")
            print(f"Length: {alignment.length}")
            for hsp in alignment.hsps:
                print(f"Score: {hsp.score}, E-value: {hsp.expect}")
                print(f"Query: {hsp.query}")
                print(f"Match: {hsp.match}")
                print(f"Subject: {hsp.sbjct}")
    else:
        print("No CDS available for BLASTn.")
except urllib.error.URLError as e:
    print(f"URL Error: {e.reason}")
except Exception as e:
    print(f"An error occurred: {e}")

