#!/usr/bin/env python3
"""Convert PLINK .bed/.bim/.fam to VCF format.

Usage: python3 plink_to_vcf.py <plink_prefix> <output_vcf>

The output VCF uses the PLINK allele conventions:
  - BIM column 5 (allele 1, A1) = ALT (minor allele in PLINK 1.9)
  - BIM column 6 (allele 2, A2) = REF (major allele in PLINK 1.9)

Genotypes are encoded as 0/0, 0/1, 1/1, or ./. (missing).
"""

import sys
import struct
import gzip

def read_bim(bim_file):
    """Read .bim file and return list of (chr, rsid, cm, bp, a1, a2)."""
    markers = []
    with open(bim_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chr_val = parts[0]
            rsid = parts[1]
            cm = parts[2]
            bp = int(parts[3])
            a1 = parts[4].upper()  # allele 1 (ALT in PLINK convention)
            a2 = parts[5].upper()  # allele 2 (REF in PLINK convention)
            markers.append((chr_val, rsid, cm, bp, a1, a2))
    return markers

def read_fam(fam_file):
    """Read .fam file and return list of (fid, iid)."""
    samples = []
    with open(fam_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            fid = parts[0]
            iid = parts[1]
            samples.append((fid, iid))
    return samples

def plink_to_vcf(plink_prefix, output_vcf):
    bim_file = plink_prefix + '.bim'
    fam_file = plink_prefix + '.fam'
    bed_file = plink_prefix + '.bed'

    print(f"Reading {bim_file}...")
    markers = read_bim(bim_file)
    print(f"  {len(markers)} markers")

    print(f"Reading {fam_file}...")
    samples = read_fam(fam_file)
    n_samples = len(samples)
    print(f"  {n_samples} samples")

    n_markers = len(markers)
    bytes_per_marker = (n_samples + 3) // 4

    # PLINK genotype encoding (2 bits per sample):
    # 00 = homozygous A1 (ALT)  -> 1/1
    # 01 = missing               -> ./.
    # 10 = heterozygous          -> 0/1
    # 11 = homozygous A2 (REF)  -> 0/0
    geno_map = {0: '1/1', 1: './.', 2: '0/1', 3: '0/0'}

    print(f"Reading {bed_file} and writing {output_vcf}...")

    with open(bed_file, 'rb') as bed:
        # Read and verify magic bytes
        magic = bed.read(3)
        if magic[0] != 0x6c or magic[1] != 0x1b:
            raise ValueError("Not a valid PLINK .bed file (bad magic number)")
        if magic[2] != 0x01:
            raise ValueError("PLINK .bed file is not SNP-major mode")

        # Open output
        open_func = gzip.open if output_vcf.endswith('.gz') else open
        mode = 'wt' if output_vcf.endswith('.gz') else 'w'

        with open_func(output_vcf, mode) as out:
            # Write VCF header
            out.write("##fileformat=VCFv4.2\n")
            out.write("##source=plink_to_vcf.py\n")
            out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

            # Column header
            sample_ids = [s[1] for s in samples]  # Use IID
            out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
            for sid in sample_ids:
                out.write(f"\t{sid}")
            out.write("\n")

            # Process each marker
            for m_idx in range(n_markers):
                if (m_idx + 1) % 10000 == 0:
                    print(f"  Processed {m_idx + 1}/{n_markers} markers...")

                chr_val, rsid, cm, bp, a1, a2 = markers[m_idx]

                # Read genotype bytes for this marker
                raw = bed.read(bytes_per_marker)

                # Decode genotypes
                genotypes = []
                for i in range(n_samples):
                    byte_idx = i // 4
                    bit_idx = (i % 4) * 2
                    geno_code = (raw[byte_idx] >> bit_idx) & 0x03
                    genotypes.append(geno_map[geno_code])

                # In PLINK: A2 is REF, A1 is ALT
                ref = a2
                alt = a1

                # Write VCF line
                out.write(f"{chr_val}\t{bp}\t{rsid}\t{ref}\t{alt}\tPASS\tPASS\t.\tGT")
                for gt in genotypes:
                    out.write(f"\t{gt}")
                out.write("\n")

    print(f"Done. Wrote {n_markers} markers to {output_vcf}")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 plink_to_vcf.py <plink_prefix> <output_vcf>")
        sys.exit(1)

    plink_to_vcf(sys.argv[1], sys.argv[2])
