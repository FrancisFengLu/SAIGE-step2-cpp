#!/usr/bin/env python3
"""
Convert PLINK binary (.bed/.bim/.fam) to BGEN v1.2 format with zstd compression.
This is for testing the BGEN reader in the standalone C++ SAIGE Step 2.

BGEN v1.2 specification: https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html

Usage:
    python3 convert_plink_to_bgen.py <plink_prefix> <output_prefix>
"""

import sys
import struct
import numpy as np

# Try zstd first, fall back to zlib
try:
    import zstandard as zstd
    USE_ZSTD = True
except ImportError:
    import zlib
    USE_ZSTD = False

COMPRESSION_ZLIB = 1
COMPRESSION_ZSTD = 2


def read_fam(fam_path):
    """Read .fam file and return sample IDs (FID, IID)."""
    samples = []
    with open(fam_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                samples.append((parts[0], parts[1]))
    return samples


def read_bim(bim_path):
    """Read .bim file and return variant info."""
    variants = []
    with open(bim_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 6:
                chrom = parts[0]
                rsid = parts[1]
                gd = parts[2]
                pos = int(parts[3])
                allele1 = parts[4].upper()  # A1 (alt in PLINK)
                allele2 = parts[5].upper()  # A2 (ref in PLINK)
                variants.append((chrom, rsid, gd, pos, allele1, allele2))
    return variants


def read_bed_genotypes_vectorized(bed_path, n_samples, n_variants):
    """Read .bed file and return genotype matrix (n_variants x n_samples).

    Vectorized with numpy for speed.

    PLINK .bed encoding (SNP-major mode):
    - 3 magic bytes: 0x6c, 0x1b, 0x01
    - Then for each variant: ceil(n_samples/4) bytes
    - Each byte encodes 4 genotypes (2 bits each):
      00 = Homozygous A1/A1 -> 2 copies of A1
      01 = Missing
      10 = Heterozygous -> 1 copy of A1
      11 = Homozygous A2/A2 -> 0 copies of A1
    """
    bytes_per_variant = (n_samples + 3) // 4

    with open(bed_path, 'rb') as f:
        magic = f.read(3)
        if magic[0] != 0x6c or magic[1] != 0x1b:
            raise ValueError("Not a valid PLINK .bed file (bad magic bytes)")
        if magic[2] != 0x01:
            raise ValueError("Only SNP-major mode is supported")
        raw = np.frombuffer(f.read(), dtype=np.uint8)

    # Decode all 4 genotypes per byte using lookup table
    # PLINK encoding: 00=2, 01=-1(missing), 10=1, 11=0
    plink_to_geno = np.array([2, -1, 1, 0], dtype=np.int8)

    genotypes = np.zeros((n_variants, n_samples), dtype=np.int8)

    for v in range(n_variants):
        offset = v * bytes_per_variant
        variant_bytes = raw[offset : offset + bytes_per_variant]

        # Extract all 4 2-bit genotypes from each byte
        g0 = plink_to_geno[variant_bytes & 0x03]
        g1 = plink_to_geno[(variant_bytes >> 2) & 0x03]
        g2 = plink_to_geno[(variant_bytes >> 4) & 0x03]
        g3 = plink_to_geno[(variant_bytes >> 6) & 0x03]

        # Interleave: sample 0 is g0[0], sample 1 is g1[0], sample 2 is g2[0], ...
        expanded = np.empty(bytes_per_variant * 4, dtype=np.int8)
        expanded[0::4] = g0
        expanded[1::4] = g1
        expanded[2::4] = g2
        expanded[3::4] = g3

        genotypes[v, :] = expanded[:n_samples]

    return genotypes


def write_bgen(output_path, samples, variants, genotypes, use_zstd=True):
    """Write BGEN v1.2 file."""
    n_samples = len(samples)
    n_variants = len(variants)
    compression = COMPRESSION_ZSTD if use_zstd else COMPRESSION_ZLIB
    layout = 2

    # Precompute the probability lookup: geno -> (p11_byte, p10_byte)
    # geno: 0=hom_ref(BB), 1=het(AB), 2=hom_alt(AA), -1=missing
    # p11 = P(AA), p10 = P(AB)
    # For dosage = 2*p11 + p10 to equal the genotype:
    #   geno=2: p11=255, p10=0
    #   geno=1: p11=0, p10=255
    #   geno=0: p11=0, p10=0
    #   missing: p11=0, p10=0 (flagged via ploidy byte)
    prob_lookup_p11 = np.array([0, 0, 255, 0], dtype=np.uint8)   # index by geno (0,1,2); 3 unused
    prob_lookup_p10 = np.array([0, 255, 0, 0], dtype=np.uint8)

    with open(output_path, 'wb') as f:
        # Header
        L_H = 20
        flags = compression | (layout << 2)

        f.write(struct.pack('<I', L_H))    # offset
        f.write(struct.pack('<I', L_H))    # L_H
        f.write(struct.pack('<I', n_variants))
        f.write(struct.pack('<I', n_samples))
        f.write(b'bgen')
        f.write(struct.pack('<I', flags))

        if use_zstd:
            compressor = zstd.ZstdCompressor(level=3)

        # Pre-build the fixed part of probability data header
        prob_header = bytearray()
        prob_header += struct.pack('<I', n_samples)
        prob_header += struct.pack('<H', 2)   # K=2
        prob_header += struct.pack('<B', 2)   # Pmin=2
        prob_header += struct.pack('<B', 2)   # Pmax=2
        prob_header_bytes = bytes(prob_header)

        # Fixed trailer after ploidy bytes
        prob_trailer = struct.pack('<B', 0) + struct.pack('<B', 8)  # Phased=0, B=8

        for v_idx in range(n_variants):
            if v_idx % 10000 == 0:
                print(f"  Writing variant {v_idx}/{n_variants}...")

            chrom, rsid, gd, pos, allele1, allele2 = variants[v_idx]
            genos = genotypes[v_idx]

            # Build variant identifying data
            snpid_b = rsid.encode('ascii')
            rsid_b = rsid.encode('ascii')
            chrom_b = chrom.encode('ascii')
            a1_b = allele1.encode('ascii')
            a2_b = allele2.encode('ascii')

            var_id_data = b''
            var_id_data += struct.pack('<H', len(snpid_b)) + snpid_b
            var_id_data += struct.pack('<H', len(rsid_b)) + rsid_b
            var_id_data += struct.pack('<H', len(chrom_b)) + chrom_b
            var_id_data += struct.pack('<I', pos)
            var_id_data += struct.pack('<H', 2)  # K=2
            var_id_data += struct.pack('<I', len(a1_b)) + a1_b
            var_id_data += struct.pack('<I', len(a2_b)) + a2_b

            # Build probability data (vectorized)
            # Ploidy/missingness bytes: 2 for normal, 130 for missing
            ploidy_bytes = np.where(genos == -1,
                                    np.uint8(130),
                                    np.uint8(2))

            # Genotype probabilities (2 bytes per sample)
            # Map genos: -1->index 3, 0->0, 1->1, 2->2
            geno_idx = np.where(genos == -1, 3, genos).astype(np.intp)
            p11_bytes = prob_lookup_p11[geno_idx]
            p10_bytes = prob_lookup_p10[geno_idx]

            # Interleave p11 and p10
            prob_pairs = np.empty(n_samples * 2, dtype=np.uint8)
            prob_pairs[0::2] = p11_bytes
            prob_pairs[1::2] = p10_bytes

            # Assemble full probability data
            prob_data = (prob_header_bytes
                        + ploidy_bytes.tobytes()
                        + prob_trailer
                        + prob_pairs.tobytes())

            # Compress
            if use_zstd:
                compressed = compressor.compress(prob_data)
            else:
                compressed = zlib.compress(prob_data)

            C = len(compressed) + 4
            D = len(prob_data)

            # Write variant block
            f.write(var_id_data)
            f.write(struct.pack('<I', C))
            f.write(struct.pack('<I', D))
            f.write(compressed)

    print(f"Wrote {n_variants} variants for {n_samples} samples to {output_path}")


def write_sample_file(output_path, samples):
    """Write BGEN .sample file."""
    with open(output_path, 'w') as f:
        f.write("ID_1 ID_2 missing\n")
        f.write("0 0 0\n")
        for fid, iid in samples:
            f.write(f"{fid} {iid} 0\n")
    print(f"Wrote {len(samples)} samples to {output_path}")


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 convert_plink_to_bgen.py <plink_prefix> <output_prefix>")
        sys.exit(1)

    plink_prefix = sys.argv[1]
    output_prefix = sys.argv[2]

    fam_path = plink_prefix + ".fam"
    bim_path = plink_prefix + ".bim"
    bed_path = plink_prefix + ".bed"

    print(f"Reading PLINK files with prefix: {plink_prefix}")

    samples = read_fam(fam_path)
    variants = read_bim(bim_path)
    n_samples = len(samples)
    n_variants = len(variants)
    print(f"  {n_samples} samples, {n_variants} variants")

    print("Reading genotypes from .bed file...")
    genotypes = read_bed_genotypes_vectorized(bed_path, n_samples, n_variants)
    print(f"  Genotype matrix shape: {genotypes.shape}")

    n_missing = np.sum(genotypes == -1)
    print(f"  Missing values: {n_missing}")

    use_zstd = USE_ZSTD
    print(f"Writing BGEN file (compression: {'zstd' if use_zstd else 'zlib'})...")
    write_bgen(output_prefix + ".bgen", samples, variants, genotypes, use_zstd=use_zstd)

    print("Writing sample file...")
    write_sample_file(output_prefix + ".sample", samples)

    print("Done!")


if __name__ == "__main__":
    main()
