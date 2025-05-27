# DNA Hidden Repeats Detection Algorithm and Isochore Research

## Overview

This repository contains ongoing research related to the development of a novel DNA segmentation and hidden repeats detection algorithm. 
The project aims to improve the identification of biologically relevant genomic regions—such as isochores and hidden repetitive sequences—by leveraging statistical models to distinguish signal from noise in DNA sequences.

The work is part of a broader capstone project titled:  
**"A new generation of DNA hidden repeats detection algorithm and its application for isochore research."**

---

## Project Motivation and Background

Isochores are large-scale, compositionally homogeneous regions of the genome that have important biological and evolutionary implications. 
Detecting hidden repeats and segmenting the genome into meaningful regions is challenging due to noise and random patterns in DNA sequences.

This project proposes a segmentation algorithm that:

- Divides the genome into fixed-length segments.
- Uses k-mer frequency analysis to identify dominant sequence patterns.
- Applies statistical significance testing (P-values and E-values from binomial models) to distinguish strong (biologically meaningful) segments from noise.
- Implements an iterative merging process resembling hierarchical clustering to refine segment boundaries.

---

## Current Status

- The algorithm and methodology are under active development.
- The theoretical framework, including statistical testing and merging logic, is being formalized.
- A written thesis/book chapter is in progress to detail the methods, results, and biological applications.
- Code implementations and example datasets will be added in the future as the project advances.

---

## Methodology Summary

1. **Initial Segmentation:** The genome is divided into fixed-size segments (e.g., 12 nucleotides).
2. **k-mer Frequency Analysis:** Each segment is analyzed for the frequency of k-mers (default k=3).
3. **Dominant k-mer Identification:** The most frequent k-mer in each segment is identified.
4. **Segment Merging:** Adjacent segments with the same dominant k-mer are merged; noisy intermediate segments are accounted for.
5. **Statistical Significance Testing:** Binomial distribution-based P-values and E-values quantify the likelihood that observed k-mer patterns are non-random.
6. **Ranking and Classification:** Segments are classified into strong (signal) or weak (noise) categories.

---

## Future Work

- Implementation of the algorithm in Python or other languages.
- Testing on real genomic datasets.
- Visualization of segmentation results.
- Integration with biological annotation databases for validation.
- Publication of research findings and source code.

---

## Contact

For questions or collaboration interests, please contact:

Your Name — your.email@example.com  
GitHub: https://github.com/f10zo

---

*This repository will be updated as the research progresses.*
