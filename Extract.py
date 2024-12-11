import argparse
import os
import pysam
#import pyBigWig

def bam_to_bed(bam_path, bed_path):
    """
    Converts a BAM file to a BED file.

    Parameters:
        bam_path (str): Path to the input BAM file.
        bed_path (str): Path to the output BED file.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam_file, open(bed_path, "w") as bed_file:
        for read in bam_file.fetch():
            if not read.is_unmapped: 
                chrom = bam_file.get_reference_name(read.reference_id)
                start = read.reference_start
                end = read.reference_end
                bed_file.write(f"{chrom}\t{start}\t{end}\n")

def bed_format(bed_path,chr_label):
    """
    custom bed files

    Parameters:
        bed_path (str): Path to the input BED file.
    """
    prefix=os.path.split(args.bed_file)[-1].split('.')[0]
    with open('tmp_'+prefix+'_bed.txt','w') as am:
        with open(bed_path, "r") as bed_file:
            for line in bed_file:
                aa=line.rstrip().split()
                if len(aa)>=5:
                    chrom = aa[0]
                    start=aa[3]
                    end=aa[4]
                    if chrom!='chr':
                        if int(start)<=int(end):
                            aa[0]=chr_label+aa[0]
                            am.write('{0}\n'.format('\t'.join(aa[:])))
                        else:
                            aa[0]=chr_label+aa[0]
                            aa[3]=end
                            aa[4]=start
                            am.write('{0}\n'.format('\t'.join(aa[:])))
    return 'tmp_'+prefix+'_bed.txt'
def extract_from_bam(bam_path, bed_path, out_path):
    """
    Count reads from a BAM file for specified regions in a BED file.

    Parameters:
        bam_path (str): Path to the input BAM file.
        bed_path (str): Path to the input BED file.
        out_path (str): Path to the output file.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam_file, open(bed_path, "r") as bed_file, open(out_path, 'w') as out_file:
        for line in bed_file:
            chrom, lstart,lend,start, end = line.strip().split()[:5]
            if int(start) <= int(end):
                count = bam_file.count(chrom, int(start), int(end))
            elif int(end) < int(start):
                count = bam_file.count(chrom, int(end), int(start))
            out_file.write(f"{line.strip()}\t{count}\n")

def extract_from_bigwig(bigwig_path, bed_path, out_path):
    """
    Extract and count reads from a bigWig file for specified regions in a BED file.

    Parameters:
        bigwig_path (str): Path to the input bigWig file.
        bed_path (str): Path to the input BED file.
        out_path (str): Path to the output file.
    """
    with pyBigWig.open(bigwig_path) as bw_file, open(bed_path, "r") as bed_file, open(out_path, 'w') as out_file:
        for line in bed_file:
            chrom, lstart,lend, start, end = line.strip().split()[:5]
            if int(start) <= int(end):
                count = bw_file.stats(chrom, int(start), int(end), type='mean')[0] if bw_file.stats(chrom, int(start), int(end), type='mean') else 0
            elif int(end) < int(start):
                count = bw_file.stats(chrom, int(end), int(start), type='mean')[0] if bw_file.stats(chrom, int(start), int(end), type='mean') else 0
            out_file.write(f"{line.strip()}\t{count}\n")

def main(args):
    tmp_file=bed_format(args.bed_file,args.chr_label)
    if args.bam_file:
        extract_from_bam(args.bam_file, tmp_file, args.out_file)
    elif args.bigwig_file:
        extract_from_bigwig(args.bigwig_file, tmp_file, args.out_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process genomic data files.")
    parser.add_argument("--bam_file", type=str, help="Path to the input BAM file", required=False)
    parser.add_argument("--bigwig_file", type=str, help="Path to the input bigWig file", required=False)
    parser.add_argument("--bed_file", type=str, help="Path to the input BED file", required=True)
    parser.add_argument("--chr_label", type=str, help="add CC001 or CC074 in the chr", required=True)
    parser.add_argument("--out_file", type=str, help="Path to the output file", required=True)
    
    args = parser.parse_args()
    main(args)
