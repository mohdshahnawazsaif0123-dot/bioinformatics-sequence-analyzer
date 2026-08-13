let D = null;
let charts = [];

const $ = id => document.getElementById(id) || {
    set textContent(v) {},
    set innerHTML(v) {},
    set className(v) {},
    onclick: null
};

const demo = `>Sample_DNA_Sequence
ATGCGTACGCTAGCTAGCTAGCTA`;

$("seq").value = demo;


/* =========================
   FASTA FILE
========================= */

$("file").onchange = async e => {
    const f = e.target.files[0];
    if (!f) return;

    $("fname").textContent = f.name;
    $("seq").value = await f.text();
};


/* =========================
   ANALYZE
========================= */

async function run() {
    try {
        $("status").textContent = "● Analyzing...";
        $("status").className = "status busy";

        const response = await fetch("/analyze", {
            method: "POST",
            headers: {
                "Content-Type": "application/json"
            },
            body: JSON.stringify({
                sequence: $("seq").value,
                dna_motifs: $("dnaMotifs").value,
                protein_motifs: $("proteinMotifs").value
            })
        });

        D = await response.json();

        if (!response.ok) {
            throw new Error(D.error);
        }

        render(D);

        $("status").textContent = "● Analysis Complete";
        $("status").className = "status done";

    } catch (error) {

        $("status").textContent = "● Error";
        $("status").className = "status error";

        alert(error.message);
    }
}


/* Buttons */

$("analyze").onclick = run;

$("motifScan").onclick = run;

window.addEventListener("load", run);

$("clear").onclick = () => {

    $("seq").value = "";

    $("fname").textContent = "No FASTA file selected";

    D = null;
};


/* =========================
   MAIN RENDER
========================= */

function render(d) {

    /* DNA SUMMARY */

    $("len").textContent = d.length;
    $("gc").textContent = d.gc;
    $("at").textContent = d.at;
    $("pl").textContent = d.protein_length;


    /* DNA */

    $("dna").textContent = d.sequence;

    $("dlen").textContent = d.length + " bp";
    $("dgc").textContent = d.gc + "%";
    $("dat").textContent = d.at + "%";
    $("ratio").textContent = d.ratio;


    /* NUCLEOTIDES */

    ["A", "T", "G", "C"].forEach(base => {
        $(base).textContent = d.counts[base];
    });


    /* TRANSFORMATIONS */

    $("comp").textContent = d.complement;

    $("rcomp").textContent = d.reverse_complement;

    $("rna").textContent = d.rna;

    $("orf").textContent = d.orf || "Not found";


    /* =========================
       PROTEIN
    ========================= */

    $("protein").textContent =
        d.protein || "No protein translated";

    $("proteinLength").textContent =
        d.protein_length + " aa";

    $("maa").textContent =
        d.most_aa || "—";

    $("maf").textContent =
        d.most_freq || "0";


    /* =========================
       PROTEIN PROPERTIES
    ========================= */

    const p = d.protein_properties;

    $("mw").textContent =
        p.molecular_weight
            ? p.molecular_weight.toFixed(2)
            : "—";

    $("pi").textContent =
        p.pI ?? "—";

    $("hyd").textContent =
        p.hydrophobic ?? "0";

    $("hydrophilic").textContent =
        p.hydrophilic ?? "0";

    $("basic").textContent =
        p.basic ?? "0";

    $("acidic").textContent =
        p.acidic ?? "0";

    $("charged").textContent =
        p.charged ?? "0";

    $("cys").textContent =
        p.cys ?? "0";


    /* =========================
       AMINO ACID TABLE
    ========================= */

    const allAA =
        "ACDEFGHIKLMNPQRSTVWY".split("");

    $("aaGrid").innerHTML = allAA.map(aa => {

        const count =
            d.amino_acids[aa] || 0;

        const percent =
            d.amino_acid_percent[aa] || 0;

        return `
            <div class="aa-card">
                <b>${aa}</b>
                <span>
                    ${count}
                    <small>${percent}%</small>
                </span>
            </div>
        `;

    }).join("");


    /* =========================
       MOTIFS
    ========================= */

    renderMotifs(
        $("dnaMotifResults"),
        d.dna_motifs,
        "DNA"
    );

    renderMotifs(
        $("proteinMotifResults"),
        d.protein_motifs,
        "Protein"
    );


    /* =========================
       CHARTS
    ========================= */

    drawCharts(d);
}


/* =========================
   MOTIF DISPLAY
========================= */

function renderMotifs(element, data, type) {

    if (!data || Object.keys(data).length === 0) {

        element.innerHTML =
            `<div class="motif">
                No ${type} motifs analyzed.
            </div>`;

        return;
    }

    element.innerHTML =
        Object.entries(data).map(([motif, positions]) => {

            const found = positions.length > 0;

            return `
                <div class="motif">

                    <b>${motif}</b>

                    <span>
                        ${
                            found
                            ? `${positions.length}
                               occurrence(s)
                               • Positions:
                               ${positions.join(", ")}`
                            : "Not found"
                        }
                    </span>

                </div>
            `;

        }).join("");
}


/* =========================
   CHARTS
========================= */

function drawCharts(d) {

    charts.forEach(chart => chart.destroy());

    charts = [];


    const chartOptions = {

        responsive: true,

        maintainAspectRatio: false,

        plugins: {

            legend: {
                labels: {
                    color: "#cbd5e1"
                }
            }

        },

        scales: {

            x: {
                ticks: {
                    color: "#94a3b8"
                },

                grid: {
                    color: "#1e3a56"
                }
            },

            y: {

                beginAtZero: true,

                ticks: {
                    color: "#94a3b8"
                },

                grid: {
                    color: "#1e3a56"
                }
            }
        }
    };


    /* =========================
       NUCLEOTIDE CHART
    ========================= */

    charts.push(

        new Chart($("base"), {

            type: "bar",

            data: {

                labels: ["A", "T", "G", "C"],

                datasets: [{

                    label: "Count",

                    data: [
                        d.counts.A,
                        d.counts.T,
                        d.counts.G,
                        d.counts.C
                    ],

                    backgroundColor: [
                        "#38bdf8",
                        "#22c55e",
                        "#c084fc",
                        "#f59e0b"
                    ]

                }]
            },

            options: chartOptions
        })
    );


    /* =========================
       GC / AT
    ========================= */

    charts.push(

        new Chart($("gcchart"), {

            type: "doughnut",

            data: {

                labels: ["GC", "AT"],

                datasets: [{

                    data: [
                        d.gc,
                        d.at
                    ],

                    backgroundColor: [
                        "#22c55e",
                        "#38bdf8"
                    ]

                }]
            },

            options: {

                responsive: true,

                maintainAspectRatio: false,

                plugins: {

                    legend: {

                        labels: {
                            color: "#cbd5e1"
                        }
                    }
                }
            }
        })
    );


    /* =========================
       20 AMINO ACIDS
    ========================= */

    const aminoAcids =
        "ACDEFGHIKLMNPQRSTVWY".split("");

    charts.push(

        new Chart($("aa"), {

            type: "bar",

            data: {

                labels: aminoAcids,

                datasets: [{

                    label: "Amino Acid Count",

                    data: aminoAcids.map(
                        aa => d.amino_acids[aa] || 0
                    ),

                    backgroundColor: "#c084fc"

                }]
            },

            options: chartOptions
        })
    );


    /* =========================
       GC ALONG SEQUENCE
    ========================= */

    charts.push(

        new Chart($("gcwindow"), {

            type: "line",

            data: {

                labels:
                    d.gc_windows.map(
                        (_, i) => i + 1
                    ),

                datasets: [{

                    label: "GC %",

                    data: d.gc_windows,

                    borderColor: "#38bdf8",

                    backgroundColor: "#38bdf833",

                    fill: true,

                    tension: 0.25,

                    pointRadius: 2

                }]
            },

            options: chartOptions
        })
    );
}


/* =========================
   EXPORT REPORT
========================= */

$("export").onclick = () => {

    if (!D) {

        alert("Please analyze a sequence first.");

        return;
    }

    const p =
        D.protein_properties;

    const aa =
        Object.entries(D.amino_acids)
        .map(([a, c]) =>
            `${a}: ${c} (${D.amino_acid_percent[a]}%)`
        )
        .join("\n");


    const dnaMotifs =
        Object.entries(D.dna_motifs)
        .map(([m, pos]) =>
            `${m}: ${
                pos.length
                ? pos.join(", ")
                : "Not found"
            }`
        )
        .join("\n");


    const proteinMotifs =
        Object.entries(D.protein_motifs)
        .map(([m, pos]) =>
            `${m}: ${
                pos.length
                ? pos.join(", ")
                : "Not found"
            }`
        )
        .join("\n");


    const report = `

BIOINFORMATICS SEQUENCE ANALYZER
==============================================

Created by Mohd Shahnawaz
MSc Bioinformatics
Jamia Millia Islamia


DNA ANALYSIS
==============================================

Sequence:
${D.sequence}

DNA Length: ${D.length} bp

A: ${D.counts.A}
T: ${D.counts.T}
G: ${D.counts.G}
C: ${D.counts.C}

GC Content: ${D.gc}%
AT Content: ${D.at}%

GC/AT Ratio: ${D.ratio}


SEQUENCE TRANSFORMATIONS
==============================================

Complement:
${D.complement}

Reverse Complement:
${D.reverse_complement}

RNA:
${D.rna}

ORF:
${D.orf || "Not found"}


PROTEIN ANALYSIS
==============================================

Protein:
${D.protein || "Not found"}

Protein Length:
${D.protein_length} aa

Most Frequent Amino Acid:
${D.most_aa}

Frequency:
${D.most_freq}


PROTEIN PROPERTIES
==============================================

Molecular Weight:
${p.molecular_weight} Da

Theoretical pI:
${p.pI}

Hydrophobic:
${p.hydrophobic}

Hydrophilic:
${p.hydrophilic}

Basic:
${p.basic}

Acidic:
${p.acidic}

Charged:
${p.charged}

Cysteine:
${p.cys}


AMINO ACID COMPOSITION
==============================================

${aa}


DNA MOTIFS
==============================================

${dnaMotifs}


PROTEIN MOTIFS
==============================================

${proteinMotifs}


==============================================
END OF REPORT
==============================================
`;


    const blob =
        new Blob(
            [report],
            { type: "text/plain" }
        );
    const link =
        document.createElement("a");

    link.href =
        URL.createObjectURL(blob);

    link.download =
        "Bioinformatics_Sequence_Analysis_Report.txt";

    link.click();

    URL.revokeObjectURL(link.href);
};