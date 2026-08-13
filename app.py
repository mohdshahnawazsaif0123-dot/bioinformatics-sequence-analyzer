from flask import Flask, render_template, request, jsonify
from collections import Counter
import re

app = Flask(__name__)


# =========================================================
# CODON TABLE
# =========================================================

CODONS = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGG": "W",

    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",

    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


# =========================================================
# AMINO ACID MOLECULAR WEIGHTS
# =========================================================

AA_MW = {
    "A": 89.09,
    "R": 174.20,
    "N": 132.12,
    "D": 133.10,
    "C": 121.15,
    "E": 147.13,
    "Q": 146.15,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "L": 131.17,
    "K": 146.19,
    "M": 149.21,
    "F": 165.19,
    "P": 115.13,
    "S": 105.09,
    "T": 119.12,
    "W": 204.23,
    "Y": 181.19,
    "V": 117.15
}


HYDROPHOBIC = set("AVILMFWY")
BASIC = set("KRH")
ACIDIC = set("DE")
AROMATIC = set("FWY")


# =========================================================
# CLEAN FASTA / DNA
# =========================================================

def clean(raw):

    sequence = "".join(
        line.strip()
        for line in raw.splitlines()
        if not line.strip().startswith(">")
    )

    sequence = re.sub(r"\s+", "", sequence).upper()

    if not sequence:
        raise ValueError("No DNA sequence found.")

    invalid = sorted(set(sequence) - set("ATGC"))

    if invalid:
        raise ValueError(
            "Invalid DNA character(s): "
            + ", ".join(invalid)
        )

    return sequence


# =========================================================
# DNA FUNCTIONS
# =========================================================

def complement(sequence):

    table = str.maketrans(
        "ATGC",
        "TACG"
    )

    return sequence.translate(table)


def reverse_complement(sequence):

    return complement(sequence)[::-1]


def transcribe(sequence):

    return sequence.replace("T", "U")


# =========================================================
# ORF
# =========================================================

def get_orf(rna):

    start = rna.find("AUG")

    if start == -1:
        return ""

    for i in range(
        start + 3,
        len(rna) - 2,
        3
    ):

        codon = rna[i:i + 3]

        if codon in ("UAA", "UAG", "UGA"):

            return rna[start:i + 3]

    return rna[start:]


# =========================================================
# TRANSLATION
# =========================================================

def translate(rna):

    start = rna.find("AUG")

    if start == -1:
        return ""

    protein = []

    for i in range(
        start,
        len(rna) - 2,
        3
    ):

        codon = rna[i:i + 3]

        amino_acid = CODONS.get(
            codon,
            ""
        )

        if not amino_acid:
            break

        if amino_acid == "*":
            break

        protein.append(amino_acid)

    return "".join(protein)


# =========================================================
# MOTIF POSITIONS
# =========================================================

def positions(sequence, motif):

    motif = motif.upper().strip()

    if not motif:
        return []

    found = []
    start = 0

    while True:

        index = sequence.find(
            motif,
            start
        )

        if index == -1:
            break

        found.append(index + 1)

        start = index + 1

    return found


# =========================================================
# GC WINDOW
# =========================================================

def gc_windows(sequence, window=10):

    if len(sequence) <= window:

        gc = (
            sequence.count("G")
            + sequence.count("C")
        )

        return [
            round(
                gc / len(sequence) * 100,
                2
            )
        ]

    values = []

    for i in range(
        0,
        len(sequence) - window + 1
    ):

        fragment = sequence[
            i:i + window
        ]

        gc = (
            fragment.count("G")
            + fragment.count("C")
        )

        values.append(
            round(
                gc / window * 100,
                2
            )
        )

    return values


# =========================================================
# PROTEIN PROPERTIES
# =========================================================

def protein_properties(protein):

    length = len(protein)

    if not protein:

        return {
            "molecular_weight": 0,
            "hydrophobic": 0,
            "hydrophilic": 0,
            "basic": 0,
            "acidic": 0,
            "charged": 0,
            "cys": 0,
            "aromatic": 0,
            "pI": None
        }

    # Molecular weight
    molecular_weight = sum(
        AA_MW.get(aa, 0)
        for aa in protein
    )

    # Peptide bond water loss
    molecular_weight -= (
        (length - 1) * 18.015
    )

    hydrophobic = sum(
        aa in HYDROPHOBIC
        for aa in protein
    )

    basic = sum(
        aa in BASIC
        for aa in protein
    )

    acidic = sum(
        aa in ACIDIC
        for aa in protein
    )

    charged = basic + acidic

    cysteine = protein.count("C")

    aromatic = sum(
        aa in AROMATIC
        for aa in protein
    )

    hydrophilic = length - hydrophobic

    # -----------------------------------------------------
    # Approximate theoretical pI
    # Educational estimate using Henderson-Hasselbalch
    # -----------------------------------------------------

    pka = {
        "Nterm": 8.0,
        "Cterm": 3.1,
        "K": 10.5,
        "R": 12.5,
        "H": 6.0,
        "D": 3.9,
        "E": 4.1,
        "C": 8.3,
        "Y": 10.1
    }

    def charge_at_pH(ph):

        charge = 0

        # N-terminus
        charge += (
            1 /
            (1 + 10 ** (ph - pka["Nterm"]))
        )

        # C-terminus
        charge -= (
            1 /
            (1 + 10 ** (pka["Cterm"] - ph))
        )

        # Basic residues
        for aa in "KRH":

            charge += (
                protein.count(aa)
                /
                (1 + 10 ** (ph - pka[aa]))
            )

        # Acidic residues
        for aa in "DECY":

            charge -= (
                protein.count(aa)
                /
                (1 + 10 ** (pka[aa] - ph))
            )

        return charge

    low = 0.0
    high = 14.0

    for _ in range(60):

        mid = (
            low + high
        ) / 2

        charge = charge_at_pH(mid)

        if charge > 0:
            low = mid
        else:
            high = mid

    theoretical_pI = round(
        (low + high) / 2,
        2
    )

    return {
        "molecular_weight": round(
            molecular_weight,
            2
        ),
        "hydrophobic": hydrophobic,
        "hydrophilic": hydrophilic,
        "basic": basic,
        "acidic": acidic,
        "charged": charged,
        "cys": cysteine,
        "aromatic": aromatic,
        "pI": theoretical_pI
    }


# =========================================================
# CUSTOM MOTIFS
# =========================================================

def parse_custom_motifs(
    raw,
    default
):

    if not raw:
        return default

    motifs = [
        item.strip().upper()
        for item in raw.split(",")
        if item.strip()
    ]

    return motifs or default


# =========================================================
# HOME
# =========================================================

@app.get("/")
def home():

    return render_template(
        "index.html"
    )


# =========================================================
# MAIN ANALYSIS API
# =========================================================

@app.post("/analyze")
def analyze():

    try:

        body = request.get_json() or {}

        sequence = clean(
            body.get(
                "sequence",
                ""
            )
        )

        counts = Counter(sequence)

        length = len(sequence)

        # -------------------------------------------------
        # GC / AT
        # -------------------------------------------------

        gc = round(
            (
                counts["G"]
                + counts["C"]
            )
            / length
            * 100,
            2
        )

        at = round(
            (
                counts["A"]
                + counts["T"]
            )
            / length
            * 100,
            2
        )

        # -------------------------------------------------
        # RNA / ORF / PROTEIN
        # -------------------------------------------------

        rna_sequence = transcribe(
            sequence
        )

        orf = get_orf(
            rna_sequence
        )

        protein = translate(
            rna_sequence
        )

        amino_acids = Counter(
            protein
        )

        # -------------------------------------------------
        # AMINO ACID %
        # -------------------------------------------------

        amino_acid_percent = {}

        if protein:

            for aa, count in sorted(
                amino_acids.items()
            ):

                amino_acid_percent[aa] = round(
                    count / len(protein) * 100,
                    2
                )

        # -------------------------------------------------
        # MOST COMMON AA
        # -------------------------------------------------

        if amino_acids:

            most_common = (
                amino_acids
                .most_common(1)[0]
            )

        else:

            most_common = (
                "—",
                0
            )

        # -------------------------------------------------
        # CUSTOM DNA MOTIFS
        # -------------------------------------------------

        dna_motifs = parse_custom_motifs(
            body.get("dna_motifs"),
            [
                "ATG",
                "TATA",
                "CAAT"
            ]
        )

        # -------------------------------------------------
        # CUSTOM PROTEIN MOTIFS
        # -------------------------------------------------

        protein_motifs = parse_custom_motifs(
            body.get("protein_motifs"),
            [
                "N-X-S/T",
                "S/T-P",
                "CXXC"
            ]
        )

        # -------------------------------------------------
        # DNA MOTIF RESULTS
        # -------------------------------------------------

        dna_motif_results = {}

        for motif in dna_motifs:

            dna_motif_results[
                motif
            ] = positions(
                sequence,
                motif
            )

        # -------------------------------------------------
        # PROTEIN MOTIF RESULTS
        # -------------------------------------------------

        protein_motif_results = {}

        for motif in protein_motifs:

            if motif == "N-X-S/T":

                matches = re.finditer(
                    r"N[^P][ST]",
                    protein
                )

                protein_motif_results[
                    motif
                ] = [
                    match.start() + 1
                    for match in matches
                ]

            elif motif == "S/T-P":

                matches = re.finditer(
                    r"[ST]P",
                    protein
                )

                protein_motif_results[
                    motif
                ] = [
                    match.start() + 1
                    for match in matches
                ]

            elif motif == "CXXC":

                matches = re.finditer(
                    r"C..C",
                    protein
                )

                protein_motif_results[
                    motif
                ] = [
                    match.start() + 1
                    for match in matches
                ]

            elif motif == "K/R-RICH":

                matches = re.finditer(
                    r"[KR]{2,}",
                    protein
                )

                protein_motif_results[
                    motif
                ] = [
                    match.start() + 1
                    for match in matches
                ]

            else:

                protein_motif_results[
                    motif
                ] = positions(
                    protein,
                    motif
                )

        # -------------------------------------------------
        # PROTEIN PROPERTIES
        # -------------------------------------------------

        properties = protein_properties(
            protein
        )

        # -------------------------------------------------
        # GC WINDOW
        # -------------------------------------------------

        gc_window_values = gc_windows(
            sequence,
            min(
                10,
                length
            )
        )

        # -------------------------------------------------
        # FINAL RESPONSE
        # -------------------------------------------------

        return jsonify({

            "success": True,

            "sequence": sequence,

            "length": length,

            "counts": {
                "A": counts["A"],
                "T": counts["T"],
                "G": counts["G"],
                "C": counts["C"]
            },

            "gc": gc,

            "at": at,

            "ratio": round(
                gc / at,
                2
            ) if at else None,

            "complement":
                complement(sequence),

            "reverse_complement":
                reverse_complement(sequence),

            "rna":
                rna_sequence,

            "orf":
                orf,

            "protein":
                protein,

            "protein_length":
                len(protein),

            "amino_acids":
                dict(
                    sorted(
                        amino_acids.items()
                    )
                ),

            "amino_acid_percent":
                amino_acid_percent,

            "most_aa":
                most_common[0],

            "most_freq":
                most_common[1],

            "dna_motifs":
                dna_motif_results,

            "protein_motifs":
                protein_motif_results,

            "protein_properties":
                properties,

            "gc_windows":
                gc_window_values
        })

    except Exception as error:

        return jsonify({

            "success": False,

            "error": str(error)

        }), 400


# =========================================================
# START SERVER
# =========================================================

if __name__ == "__main__":

    app.run(
        debug=True
    )