# 🧬 DNA Hidden Repeats Detection Algorithm and Isochore Research

![Status](https://img.shields.io/badge/project-phase%20A-blue)
![Language](https://img.shields.io/badge/python-planned-yellow)
![License](https://img.shields.io/badge/license-TBD-lightgrey)
![Contributions](https://img.shields.io/badge/contributions-welcome-brightgreen)

---

## 🧠 Overview

This repository is part of a capstone research project focused on designing a novel algorithm for **DNA segmentation** and the **detection of hidden repeats**.  
The goal is to identify biologically relevant genomic regions—especially **isochores** and **hidden repetitive sequences**—by applying statistical models that separate **signal** from **noise**.

**📌 Project Title:**  
**_A new generation of DNA hidden repeats detection algorithm and its application for isochore research_**

---

## 🔍 Background & Motivation

**Isochores** are long, compositionally homogeneous DNA regions with relatively constant GC content. They play key roles in genome organization, gene regulation, and chromosomal structure.

However, a possible underlying cause of isochores—**hidden DNA repeats**—is not well understood or easily detected using traditional methods.

This project addresses the challenge by:

- Segmenting DNA sequences into fixed-length windows.
- Detecting **dominant k-mer patterns** using frequency analysis.
- Quantifying **statistical significance** using binomial models.
- Iteratively **merging segments** based on shared features and noise thresholds.

---

## ⚙️ Methodology Overview

### Step-by-Step Process:

1. **Initial Segmentation**  
   → Divide DNA into fixed-size segments (e.g., 12 bp).

2. **k-mer Frequency Analysis**  
   → Count all k-mers (default k=3) in each segment.

3. **Dominant k-mer Identification**  
   → Detect the most frequent k-mer per segment.

4. **Statistical Significance Testing**  
   → Compute:
   - **P-values** via binomial distribution.
   - **E-values** for multiple testing correction.

5. **Segment Merging**  
   → Merge adjacent segments with the same dominant k-mer.  
   → Account for noisy intermediate segments using heuristic rules.

6. **Classification**  
   → Label segments as **strong** (signal) or **weak** (noise).

---

## 📊 Visual Workflow

Here’s a simplified conceptual flowchart (can be replaced with an image later):

+---------------------+
| Initial Segmentation|
+---------------------+
↓
+---------------------+
| k-mer Frequency |
+---------------------+
↓
+---------------------+
| Dominant k-mer Found |
+---------------------+
↓
+-------------------------------+
| Statistical Significance Test|
+-------------------------------+
↓
+---------------------+
| Segment Merging |
+---------------------+
↓
+---------------------+
| Final Classification|
+---------------------+


## 📌 Current Status

- ✅ Theoretical model finalized.
- ✅ Segmentation and k-mer logic designed.
- 🧪 P-value & E-value computation tested.
- ⚙️ Algorithm implementation (in Python) is ongoing.
- 📖 Thesis writing in progress.
- 🔬 Future testing will involve real genomic datasets.

---

## 🚀 Future Plans

- [ ] Implement algorithm in Python.
- [ ] Add dataset loaders and parsers.
- [ ] Integrate data visualization tools (e.g., GC plots, repeat maps).
- [ ] Run experiments on human genome samples.
- [ ] Compare results with known annotations (GENCODE, UCSC).
- [ ] Publish a paper and release full code with documentation.

---

## 🤝 Contributors

- **Fatmeh Zoabi** — [fatmehzo3bi10@gmail.com](mailto:fatmehzo3bi10@gmail.com)  
- **Khalil Mansour** — [Khalel.Mnsor@e.braude.ac.il](mailto:Khalel.Mnsor@e.braude.ac.il)

GitHub: [github.com/f10zo](https://github.com/f10zo)

---

> ⚠️ *This repository is under active development and will be updated frequently as the research progresses.*

