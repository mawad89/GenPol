# GENomePOLisher

**GENPOL** is a consensus module designed to polish contig assemblies, enhancing the accuracy of genome assemblies. It supports both short-read and long-read sequencing technologies.

### Key Features
- **Inputs**: Requires three files — contigs in FASTA format, reads in FASTA/FASTQ format, and mapping files in BED format with CIGAR strings.
- **Outputs**: Generates a set of polished contigs in FASTA format.

---

## Installation

To install GENPOL:

```bash
git clone https://github.com/mawad89/GenPol.git
cd GenPol
```
Add the `GenPol` directory to your `PATH` to enable command-line access.

## Usage Guide

1. **Read Mapping**  
   Choose any preferred mapping tool; however, the recommended tools are:
   - **Minimap2** for long reads.
   - **BWA** for short reads.

   > **Note**: Ensure the CIGAR string is in extended format (using `X=`).

2. **Convert Mapping Files**  
   Convert the mapped files into BED format, including CIGAR strings.

3. **Contig Separation**  
   To optimize processing time, GENPOL is designed to polish each contig individually. Separate the BED files by contig for parallel processing.

4. **Run Variant Calling**  
   Run the `vc.py` script to identify variants within each read, recording the variant positions in both the reads and contigs.

5. **Filter Variants**  
   Use the `filter.py` script to aggregate the variants identified by `vc.py`. This step filters out variants with less than 51% occurrence, retaining only those with sufficient support for accurate polishing.

6. **Generate Consensus**  
   Run the `polish.py` script to generate the consensus sequence for each contig.

7. **Final Assembly**  
   Concatenate all polished contigs into a single FASTA file.

## License

GENPOL is distributed under the MIT License. See the `LICENSE` file for details.



