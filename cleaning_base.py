from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
def correct_iupac_in_ab1(ab1_path, output_fasta):
    rec  = SeqIO.read(ab1_path, "abi")
    abif = rec.annotations["abif_raw"]

    # Instrument base calls & peak locations
    PBAS1 = abif["PBAS1"].decode("ascii")           # string of called bases
    PLOC2 = np.asarray(abif["PLOC2"], dtype=int)    # peak indices (one per base)

    # Channel order: maps DATA9..DATA12 to bases
    FWO_1 = abif["FWO_1"].decode("ascii")
    print("Channel order (FWO_1):", FWO_1)

    # Pull the four traces in the *correct* base order using FWO_1
    data_tags = ["DATA9", "DATA10", "DATA11", "DATA12"]
    traces_by_base = {b.upper(): np.asarray(abif[tag], dtype=int)
                    for b, tag in zip(FWO_1, data_tags)}


    PLOC2 = np.asarray(abif["PLOC2"], dtype=int)   # indices in the trace where peaks occur
    BASES = ("A","C","G","T")
    peak_heights = np.stack([traces_by_base[b][PLOC2] for b in BASES], axis=1)
    # shape: (n_bases, 4); columns correspond to A,C,G,T in that order
    peak_heights

    # Create base dict index
    base_to_index = {i: b for i, b in enumerate(BASES)}
    base_to_index

    # detected iupac 
    iupac_detected = {}
    for idx, base in  enumerate(PBAS1):
        if base != 'N' and base in iupac_list:
            if base not in iupac_detected:
                iupac_detected[base] = [idx]
            else:
                iupac_detected[base].append(idx)


    # replace iupac in PBAS1 with highest peak height base
    PBAS1_detected = PBAS1 # make a copy to modify
    for base, values in iupac_detected.items():
        for idx in values:
            peak = peak_heights[idx]
            #print(f"Index: {idx}, IUPAC: {base}, Peaks: {peak}")
            max_base = base_to_index[np.argmax(peak)]
            #print(f"Replacing {base} with {max_base}")
            # Replace in PBAS
            PBAS1_detected = PBAS1_detected[:idx] + max_base + PBAS1_detected[idx+1:]

    rec = SeqRecord(Seq(PBAS1_detected), id=f"{output_fasta}", description="")
    SeqIO.write(rec, output_fasta, "fasta")


# Get the path list of ab1 files to correct
from pathlib import Path
def define_file_lists(ab1_dir):
    root = Path(ab1_dir)          # use r"C:\path\to\dir" on Windows
    # Non-recursive: files directly inside root
    files = [p for p in root.iterdir() if p.is_file()]

    # Recursive: all files in root and subfolders
    files_recursive = [p for p in root.rglob("*") if p.is_file()]

    # Only certain extensions (e.g., .ab1 and .fasta), recursive
    wanted = {".ab1", ".fasta"}
    files_filtered = [p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in wanted]

    # Convert to absolute strings if you need plain paths
    paths = [str(p.resolve()) for p in files_filtered]

    # extract the name without suffix
    names = [p.stem for p in files_filtered]

    # Create output directory if it doesn't exist
    output_dir = Path(f"./{ab1_dir}_corrected")
    output_dir.mkdir(parents=True, exist_ok=True)
    return files_filtered, names, output_dir

if __name__ == "__main__":
    # execute correction for each file
    ab1_dir = "sequence_ab1"
    files_filtered, names, output_dir = define_file_lists(ab1_dir)
    for ab1_path, name in zip(files_filtered, names):
        output_fasta = output_dir / f"{name}.fasta"
        correct_iupac_in_ab1(ab1_path, output_fasta)