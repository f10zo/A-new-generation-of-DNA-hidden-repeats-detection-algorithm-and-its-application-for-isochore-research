# 🧬 DNA Hidden Repeats Detection Algorithm and Isochore Research
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
- Quantifying **statistical significance** using P-values and a segment scoring strategy.
- Applying **Fisher's method** to combine multiple P-values for more robust signal detection.
- Iteratively **merging segments** based on shared features and noise thresholds.

---

## ⚙️ Methodology Overview
### Step-by-Step Process:

1. **Initial Segmentation**  
   → Divide DNA into fixed-size segments (e.g., 12 bp).

2. **k-mer Frequency Analysis**  
   → Count all k-mers (default k=3) in each segment.

3. **Dominant k-mer Identification**  
   → Identify the most frequent k-mer within each segment.

4. **Statistical Significance Testing**  
   → Evaluate how likely the dominant k-mer pattern is due to chance using:
   - **P-values**: Calculated from a binomial model based on expected vs. observed k-mer frequencies.
   - **Fisher’s Method**: Combine multiple P-values (from the occurrence matrix) to produce a single test statistic.
   - **Chi-squared Distribution**: Used to convert Fisher’s statistic into a combined P-value, enabling interpretation of overall segment significance.

5. **Segment Merging**  
   → Merge adjacent segments with the same dominant k-mer.  
   → Account for noisy or weak segments using conditional logic during merging.

6. **Final Classification**  
   → Classify segments as **strong** (signal) or **weak** (noise) based on:
   - The combined Fisher P-value,
   - Segment score,
   - And the overall structure and continuity of patterns.

---

## 📌 Current Status

- ✅ Theoretical model finalized.
- ✅ Segmentation and k-mer logic designed.
- 🧪 P-value, segment score, and Fisher combination logic under development.
- 📖 Thesis writing in progress.
- 🔬 Future testing will involve real genomic datasets.

---

## 🚀 Research Plan

### 🔬 Objectives

- Map **hidden repeats** within isochores using DNA sequences.
- Use reference data from genomic databases.
- Analyze genomes from multiple species (e.g., monkey, dog).
- Compare hidden repeat maps with known isochore structures to assess correlation.

### 🛠️ Planned Steps

1. **Construct Hidden Repeat Maps**
   - Generate repeat-based segmentations for each genome.
   - Represent each map as a data structure for comparison.

2. **Use Genomic Databases**
   - Source DNA sequences from public genome repositories.

3. **Analyze Multiple Species**
   - Perform segmentation on dozens of genomes across different species.
   - Focus on evolutionary comparison between organisms.

4. **Compare Maps Across Species**
   - Measure similarity between hidden repeat maps using metrics such as:
     - **Jaccard Index**
     - **Cosine Similarity**
   - Evaluate how well the structure of repeats aligns with isochore organization.

5. **Fragmentation Comparison**
   - Analyze and compare segmentations (fragment boundaries) between species.
   - Use existing segmentation comparison algorithms to assist with structural comparison.

---

*These steps will help test the hypothesis that hidden repeats influence or align with isochore structure across evolutionary scales.*

---

## 🤝 Contributors

- **Fatmeh Zoabi** — [fatmehzo3bi10@gmail.com](mailto:fatmehzo3bi10@gmail.com)  
- **Khalil Mansour** — [Khalel.Mnsor@e.braude.ac.il](mailto:Khalel.Mnsor@e.braude.ac.il)

GitHub: [github.com/f10zo](https://github.com/f10zo)

---

> ⚠️ *This repository is under active development and will be updated frequently as the research progresses.*
