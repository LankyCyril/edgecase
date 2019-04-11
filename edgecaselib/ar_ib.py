from sys import stdout
from pysam import AlignmentFile
from edgecaselib.util import ReadFileChain


def get_readname_set(bams):
    """Get set of read names in BAM files"""
    with ReadFileChain(bams, AlignmentFile) as bam_data:
        return {entry.qname for entry in bam_data}


def main(ar, acs, file=stdout, **kwargs):
    # blacklist names of reads in AC files:
    blacklist = get_readname_set(acs)
    # print out mapped non-blacklisted reads:
    with AlignmentFile(ar) as bam_data:
        for entry in bam_data:
            if not (entry.is_unmapped or entry.is_secondary):
                if not (entry.is_supplementary or (entry.qname in blacklist)):
                    print(">{}\n{}".format(entry.qname, entry.seq), file=file)
