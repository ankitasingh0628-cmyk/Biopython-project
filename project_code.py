from Bio import Entrez, SeqIO  # type: ignore
from Bio.Blast import NCBIWWW, NCBIXML # type: ignore
from Bio.Align import PairwiseAligner # type: ignore



try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("⚠ matplotlib not installed. Skipping visualization step.")


Entrez.email = "uttamranjan123@gmail.com"   
nuc_accession = "NM_001371415"            


print("\n--- Fetching Nucleotide Sequence ---")
handle = Entrez.efetch(db="nucleotide", id=nuc_accession, rettype="gb", retmode="text")
record = SeqIO.read(handle, "genbank")
handle.close()

print("Accession:", record.id)
print("Description:", record.description)


seq = record.seq
gc_content = round((seq.count("G")+seq.count("C"))/len(seq)*100, 2)
print("Total length:", len(seq))
print("GC content:", gc_content, "%")

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
    print("Protein (first 60 aa):", protein[:60])
else:
    print("No CDS found.")
    exit()


aligner = PairwiseAligner()
aligner.mode = "global"
alns = aligner.align(protein, protein)
print("\n--- Alignment Demo ---")
print("Alignment score:", alns[0].score)
print(alns[0])


if MATPLOTLIB_AVAILABLE:
    bases = ["A","T","G","C"]
    counts = [seq.count(b) for b in bases]
    plt.bar(bases, counts, color="skyblue")
    plt.title("Nucleotide Composition")
    plt.xlabel("Base")
    plt.ylabel("Count")
    plt.savefig("nucleotide_composition.png")
    plt.close()
    print("Plot saved as nucleotide_composition.png")
else:
    print("Visualization skipped (matplotlib not available).")


print("\n--- Running BLASTn (CDS vs nt) ---")
blastn_handle = NCBIWWW.qblast("blastn", "nt", cds_seq)
blastn_record = NCBIXML.read(blastn_handle)
print("Top BLASTn hit:", blastn_record.alignments[0].hit_def)



print("\n--- Annotations (first 5 features) ---")
for i, feature in enumerate(record.features[:5]):
    print(i, feature.type, feature.location)
