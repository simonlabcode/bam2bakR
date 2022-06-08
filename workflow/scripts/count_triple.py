import sys
import argparse
import operator
import itertools
import warnings
import traceback
import os.path
import multiprocessing
import pysam
import random

import HTSeq

'''
    Note: 
        This is a modified version of htseq-count script from HTSeq/0.13.5 
        that allows running multiple different counting jobs on a single pass
        of .sam and GFF/GTF files.

        e.g. three counting jobs: #1 gene, gene_id, union;  #2 exon, gene_id, union;  #3 exon, gene_id, intersection-strict 
        python count_triple.py \
            -f sam \
            --samout triple_htseq.sam \
            -t gene,exon,exon \
            -i gene_id,gene_id,gene_id \
            -m union,union,intersection-strict \
            -c EF_htseq.txt,GF_htseq.txt,XF_htseq.txt \
            input.sam  \
            annotation.gtf

        Counts are returned in separate files.
        Output is a single file with each read containing new flags: GF:Z, EF:Z, XF:Z, [XG:Z, XH:Z, ...]
    Modified_by: 
        Martin Machyna, 8/27/21
'''

class UnknownChrom(Exception):
    pass


def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2


def make_feature_genomicarrayofsets_triple(
        feature_sequence,
        id_attribute,
        feature_type=None,
        feature_query=None,
        additional_attributes=None,
        stranded=False,
        verbose=False,
        ):
    """Organize a sequence of Feature objects into a GenomicArrayOfSets.
    Note: 
        This is a modified function from HTSeq/0.13.5 to accomocate creating
            multiple genomicarrayofsets on a single pass of GFF/GTF file.
        id_attribute: is a list of strings
        feature_type: is a list of strings
        return: is a list of dictionaries with features, attriburtes keys
    Modified_by: 
        Martin Machyna, 8/27/21

    Args:
        feature_sequence (iterable of Feature): A sequence of features, e.g. as
            obtained from GFF_reader('myfile.gtf')
        id_attribute (string): An attribute to use to identify the feature in
            the output data structures (e.g. 'gene_id')
        feature_type (string or None): If None, collect all features. If a
            string, restrict to only one type of features, e.g. 'exon'.
        feature_query (string or None): If None, all features of the selected
            types will be collected. If a string, it has to be in the format:

        <feature_attribute> == <attr_value>

        e.g.

        'gene_id == "Fn1"'

        (note the double quotes inside).

        Then only that feature will be collected. Using this argument is more
        efficient than collecting all features and then pruning it down to a
        single one.

        additional_attributes (list or None): A list of additional attributes
            to be collected into a separate dict for the same features, for
            instance ['gene_name']
        stranded (bool): Whether to keep strandedness information
        verbose (bool): Whether to output progress and error messages

    Returns:
        dict with two keys, 'features' with the GenomicArrayOfSets populated
        with the features, and 'attributes' which is itself a dict with
        the id_attribute as keys and the additional attributes as values.   

    Example: Let's say you load the C. elegans GTF file from Ensembl and make a
    feature dict:

    >>> gff = HTSeq.GFF_Reader("Caenorhabditis_elegans.WS200.55.gtf.gz")
    >>> worm_features = HTSeq.make_feature_genomicarrayofsets(gff)

    (This command may take a few minutes to deal with the 430,000 features
    in the GTF file. Note that you may need a lot of RAM if you have millions
    of features.)

    This function is related but distinct from HTSeq.make_feature_dict. This
    function is used in htseq-count and its barcoded twin to count gene
    expression because the output GenomicArrayofSets is very efficient. You
    can use it in performance-critical scans of GFF files.
    """

    if additional_attributes is None:
        additional_attributes = []

    if feature_query is not None:
        feature_qdic = _parse_feature_query(feature_query)

    features = [HTSeq.GenomicArrayOfSets("auto", stranded) for feature in feature_type]
    attributes = [{} for feature in feature_type]
    i = 0
    try:
        for f in feature_sequence:
            for p in range(0, len(feature_type)):
                if feature_type[p] in (None, f.type):
                    try:
                        feature_id = f.attr[id_attribute[p]]
                    except KeyError:
                        raise ValueError(
                                "Feature %s does not contain a '%s' attribute" %
                                (f.name, id_attribute[p]))
                    if stranded and f.iv.strand == ".":
                        raise ValueError(
                                "Feature %s at %s does not have strand information but you are "
                                "using stranded mode. Try with unstrnded mode." %
                                (f.name, f.iv))

                    if feature_query is not None:
                        # Skip the features that don't even have the right attr
                        if feature_qdic['attr_cat'] not in f.attr:
                            continue
                        # Skip the ones with an attribute with a different name
                        # from the query (e.g. other genes)
                        if f.attr[feature_qdic['attr_cat']] != feature_qdic['attr_name']:
                            continue

                    features[p][f.iv] += feature_id
                    attributes[p][feature_id] = [
                            f.attr[attr] if attr in f.attr else ''
                            for attr in additional_attributes]
            i += 1
            if i % 100000 == 0 and verbose:
                if hasattr(feature_sequence, 'get_line_number_string'):
                    msg = "{:d} GFF lines processed.".format(i)
                else:
                    msg = "{:d} features processed.".format(i)
                sys.stderr.write(msg+'\n')
                sys.stderr.flush()
    except(KeyError, ValueError):
        if verbose:
            if hasattr(feature_sequence, 'get_line_number_string'):
                msg = "Error processing GFF file ({:}):".format(
                    feature_sequence.get_line_number_string())
            else:
                msg = "Error processing feature sequence ({:}):".format(
                    str(i+1))
            sys.stderr.write(msg+'\n')
        raise

    if verbose:
        if hasattr(feature_sequence, 'get_line_number_string'):
            msg = "{:d} GFF lines processed.".format(i)
        else:
            msg = "{:d} features processed.".format(i)
        sys.stderr.write(msg+"\n")
        sys.stderr.flush()

    return [{
        'features': features[p],
        'attributes': attributes[p],
        } for p in range(0, len(feature_type))]



def count_reads_single_file(
        isam,
        sam_filename,
        features,
        feature_attr,
        order,
        max_buffer_size,
        stranded,
        overlap_mode,
        multimapped_mode,
        secondary_alignment_mode,
        supplementary_alignment_mode,
        feature_type,
        id_attribute,
        additional_attributes,
        quiet,
        minaqual,
        samout_format,
        samout_filename,
        ):
    '''
    Note: 
        This is a modified function from HTSeq/0.13.5 to accomocate creating
            multiple geno`micarrayofsets on a single pass of GFF/GTF file.
        features:     is a list of strings
        feature_attr: is a list of strings
        overlap_mode: is a list of strings
        feature_type: is a list of strings
        id_attribute: is a list of strings
        return: list of dictionatires
    Modified_by: 
        Martin Machyna, 8/27/21

    '''

    def write_to_samout(r, assignment, samoutfile, template=None):
        '''
        assignment: is a list of strings
        '''
        if samoutfile is None:
            return
        if not pe_mode:
            r = (r,)
        tags = ['GF','EF','XF','XG','XH','XI','XJ','XK','XL','XN','XO','XP','XQ','XR','XS','XT','XU']
        for read in r:
            if read is not None:
                for p in range(0, len(assignment)):
                    read.optional_fields.append((tags[p], assignment[p]))
                if samout_format in ('SAM', 'sam'):
                    samoutfile.write(read.get_sam_line() + "\n")
                else:
                    samoutfile.write(read.to_pysam_AlignedSegment(template))

    try:
        if sam_filename == "-":
            read_seq_file = HTSeq.BAM_Reader(sys.stdin)
        else:
            read_seq_file = HTSeq.BAM_Reader(sam_filename)

        # Get template for output BAM
        if samout_filename is None:
            template = None
            samoutfile = None
        elif samout_format in ('bam', 'BAM'):
            template = read_seq_file.get_template()
            samoutfile = pysam.AlignmentFile(
                    samout_filename, 'wb',
                    template=template,
                    )
        else:
            template = None
            samoutfile = open(samout_filename, 'w')

        read_seq_iter = iter(read_seq_file)
        # Catch empty BAM files
        try:
            first_read = next(read_seq_iter)
            pe_mode = first_read.paired_end
        # FIXME: catchall can hide subtle bugs
        except:
            first_read = None
            pe_mode = False
        if first_read is not None:
            read_seq = itertools.chain([first_read], read_seq_iter)
        else:
            read_seq = []
    except:
        sys.stderr.write(
            "Error occured when reading beginning of SAM/BAM file.\n")
        raise

    # CIGAR match characters (including alignment match, sequence match, and
    # sequence mismatch
    com = ('M', '=', 'X')
    counts = [{key: 0 for key in feature_attr[p]} for p in range(0, len(feature_type))]

    try:
        if pe_mode:
            if ((supplementary_alignment_mode == 'ignore') and
               (secondary_alignment_mode == 'ignore')):
                primary_only = True
            else:
                primary_only = False
            if order == "name":
                read_seq = HTSeq.pair_SAM_alignments(
                        read_seq,
                        primary_only=primary_only)
            elif order == "pos":
                read_seq = HTSeq.pair_SAM_alignments_with_buffer(
                        read_seq,
                        max_buffer_size=max_buffer_size,
                        primary_only=primary_only)
            else:
                raise ValueError("Illegal order specified.")
        empty = [0] * len(feature_type)
        ambiguous = [0] * len(feature_type)
        notaligned = 0
        lowqual = 0
        nonunique = 0
        i = 0
        for r in read_seq:
            if i > 0 and i % 100000 == 0 and not quiet:
                sys.stderr.write(
                    "%d alignment record%s processed.\n" %
                    (i, "s" if not pe_mode else " pairs"))
                sys.stderr.flush()

            i += 1
            if not pe_mode:
                if not r.aligned:
                    notaligned += 1
                    write_to_samout(
                            r, ["__not_aligned"] * len(feature_type), samoutfile,
                            template)
                    continue
                if ((secondary_alignment_mode == 'ignore') and
                   r.not_primary_alignment):
                    continue
                if ((supplementary_alignment_mode == 'ignore') and
                   r.supplementary):
                    continue
                try:
                    if r.optional_field("NH") > 1:
                        nonunique += 1
                        write_to_samout(
                                r, ["__alignment_not_unique"] * len(feature_type), samoutfile,
                                template)
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if r.aQual < minaqual:
                    lowqual += 1
                    write_to_samout(
                            r, ["__too_low_aQual"] * len(feature_type), samoutfile,
                            template)
                    continue
                if stranded != "reverse":
                    iv_seq = (co.ref_iv for co in r.cigar if co.type in com
                              and co.size > 0)
                else:
                    iv_seq = (invert_strand(co.ref_iv)
                              for co in r.cigar if (co.type in com and
                                                    co.size > 0))
            else:
                if r[0] is not None and r[0].aligned:
                    if stranded != "reverse":
                        iv_seq = (co.ref_iv for co in r[0].cigar
                                  if co.type in com and co.size > 0)
                    else:
                        iv_seq = (invert_strand(co.ref_iv) for co in r[0].cigar
                                  if co.type in com and co.size > 0)
                else:
                    iv_seq = tuple()
                if r[1] is not None and r[1].aligned:
                    if stranded != "reverse":
                        iv_seq = itertools.chain(
                                iv_seq,
                                (invert_strand(co.ref_iv) for co in r[1].cigar
                                if co.type in com and co.size > 0))
                    else:
                        iv_seq = itertools.chain(
                                iv_seq,
                                (co.ref_iv for co in r[1].cigar
                                 if co.type in com and co.size > 0))
                else:
                    if (r[0] is None) or not (r[0].aligned):
                        write_to_samout(
                                r, ["__not_aligned"] * len(feature_type), samoutfile,
                                template)
                        notaligned += 1
                        continue
                if secondary_alignment_mode == 'ignore':
                    if (r[0] is not None) and r[0].not_primary_alignment:
                        continue
                    elif (r[1] is not None) and r[1].not_primary_alignment:
                        continue
                if supplementary_alignment_mode == 'ignore':
                    if (r[0] is not None) and r[0].supplementary:
                        continue
                    elif (r[1] is not None) and r[1].supplementary:
                        continue
                try:
                    if ((r[0] is not None and r[0].optional_field("NH") > 1) or
                       (r[1] is not None and r[1].optional_field("NH") > 1)):
                        nonunique += 1
                        write_to_samout(
                                r, ["__alignment_not_unique"] * len(feature_type), samoutfile,
                                template)
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if ((r[0] and r[0].aQual < minaqual) or
                   (r[1] and r[1].aQual < minaqual)):
                    lowqual += 1
                    write_to_samout(
                            r, ["__too_low_aQual"] * len(feature_type), samoutfile,
                            template)
                    continue

            iv_seq = list(iv_seq)
            fs = [''] * len(feature_type)
            samtagvalue = [''] * len(feature_type)
            for p in range(0, len(feature_type)):
                try:
                    if overlap_mode[p] == "union":
                        fs[p] = set()
                        for iv in iv_seq:
                            if iv.chrom not in features[p].chrom_vectors:
                                raise UnknownChrom
                            for iv2, fs2 in features[p][iv].steps():
                                fs[p] = fs[p].union(fs2)
                    elif overlap_mode[p] in ("intersection-strict",
                                          "intersection-nonempty"):
                        fs[p] = None
                        for iv in iv_seq:
                            if iv.chrom not in features[p].chrom_vectors:
                                raise UnknownChrom
                            for iv2, fs2 in features[p][iv].steps():
                                if ((len(fs2) > 0) or
                                   (overlap_mode[p] == "intersection-strict")):
                                    if fs[p] is None:
                                        fs[p] = fs2.copy()
                                    else:
                                        fs[p] = fs[p].intersection(fs2)
                    else:
                        sys.exit("Illegal overlap mode.")

                    if fs[p] is None or len(fs[p]) == 0:
                        samtagvalue[p] = "__no_feature"
                        empty[p] += 1
                    elif len(fs[p]) > 1:
                        samtagvalue[p] = "__ambiguous[" + '+'.join(sorted(fs[p])) + "]"
                        ambiguous[p] += 1
                    else:
                        samtagvalue[p] = list(fs[p])[0]

                    if fs[p] is not None and len(fs[p]) > 0:
                        if multimapped_mode == 'none':
                            if len(fs[p]) == 1:
                                counts[p][list(fs[p])[0]] += 1
                        elif multimapped_mode == 'all':
                            for fsi in list(fs[p]):
                                counts[p][fsi] += 1
                        elif multimapped_mode == 'fraction':
                            for fsi in list(fs[p]):
                                counts[p][fsi] += 1.0 / len(fs[p])
                        elif multimapped_mode == 'random':
                            fsi = random.choice(list(fs[p]))
                            counts[p][fsi] += 1
                        else:
                            sys.exit("Illegal multimap mode.")

                except UnknownChrom:
                    samtagvalue[p] = "__no_feature"
                    empty[p] += 1
            write_to_samout(
                        r, samtagvalue, samoutfile,
                        template)
    except:
        sys.stderr.write(
            "Error occured when processing input (%s):\n" %
            (read_seq_file.get_line_number_string()))
        raise

    if not quiet:
        sys.stderr.write(
            "%d %s processed.\n" %
            (i, "alignments " if not pe_mode else "alignment pairs"))
        sys.stderr.flush()

    if samoutfile is not None:
        samoutfile.close()

    return [{
        'isam': isam,
        'counts': counts[p],
        'empty': empty[p],
        'ambiguous': ambiguous[p],
        'lowqual': lowqual,
        'notaligned': notaligned,
        'nonunique': nonunique,
    } for p in range(0, len(feature_type))]


def count_reads_in_features(
        sam_filenames,
        gff_filename,
        order,
        max_buffer_size,
        stranded,
        overlap_mode,
        multimapped_mode,
        secondary_alignment_mode,
        supplementary_alignment_mode,
        feature_type,
        id_attribute,
        additional_attributes,
        quiet,
        minaqual,
        samouts,
        samout_format,
        output_delimiter,
        output_filename,
        output_append,
        nprocesses,
        feature_query,
        ):
    '''Count reads in features, parallelizing by file
    Note: 
        This is a modified function from HTSeq/0.13.5 to accomocate creating
            multiple geno`micarrayofsets on a single pass of GFF/GTF file.
        overlap_mode: is a list of strings
        feature_type: is a list of strings
        id_attribute: is a list of strings
    Modified_by: 
        Martin Machyna, 8/27/21
    '''
    # Never use more CPUs than files
    nprocesses = min(nprocesses, len(sam_filenames))

    if samouts != []:
        if len(samouts) != len(sam_filenames):
            raise ValueError(
                    'Select the same number of input and output files')
        # Try to open samout files early in case any of them has issues
        if samout_format in ('SAM', 'sam'):
            for samout in samouts:
                with open(samout, 'w'):
                    pass
        else:
            # We don't have a template if the input is stdin
            if (len(sam_filenames) != 1) or (sam_filenames[0] != '-'):
                for sam_filename, samout in zip(sam_filenames, samouts):
                    with pysam.AlignmentFile(sam_filename, 'r') as sf:
                        with pysam.AlignmentFile(samout, 'w', template=sf):
                            pass
    else:
        samouts = [None for x in sam_filenames]

    # Try to open samfiles to fail early in case any of them is not there
    if (len(sam_filenames) != 1) or (sam_filenames[0] != '-'):
        for sam_filename in sam_filenames:
            with pysam.AlignmentFile(sam_filename, 'r') as sf:
                pass

    # Prepare features
    gff = HTSeq.GFF_Reader(gff_filename)
    feature_scan = make_feature_genomicarrayofsets_triple(
        gff,
        id_attribute,
        feature_type=feature_type,
        feature_query=feature_query,
        additional_attributes=additional_attributes,
        stranded=stranded != 'no',
        verbose=not quiet,
        )
    features = [feature_scan[p]['features'] for p in range(0, len(feature_type))]
    attributes = [feature_scan[p]['attributes'] for p in range(0, len(feature_type))]
    feature_attr = [sorted(attributes[p].keys()) for p in range(0, len(feature_type))]

    for p in range(0, len(feature_type)):
        if len(feature_attr[p]) == 0:
            sys.stderr.write(
                "Warning: No features of type '%s' found.\n" % feature_type[p])

    # Prepare arguments for counting function
    args = []
    for isam, (sam_filename, samout_filename) in enumerate(zip(sam_filenames, samouts)):
        args.append((
            isam,
            sam_filename,
            features,
            feature_attr,
            order,
            max_buffer_size,
            stranded,
            overlap_mode,
            multimapped_mode,
            secondary_alignment_mode,
            supplementary_alignment_mode,
            feature_type,
            id_attribute,
            additional_attributes,
            quiet,
            minaqual,
            samout_format,
            samout_filename,
            ))

    # Count reads
    if nprocesses > 1:
        with multiprocessing.Pool(nprocesses) as pool:
            results = pool.starmap(count_reads_single_file, args)
        results.sort(key=operator.itemgetter('isam'))
    else:
        results = list(itertools.starmap(count_reads_single_file, args))

    # Write output
    other_features = [
        ('__no_feature', 'empty'),
        ('__ambiguous', 'ambiguous'),
        ('__too_low_aQual', 'lowqual'),
        ('__not_aligned', 'notaligned'),
        ('__alignment_not_unique', 'nonunique'),
        ]
    pad = ['' for attr in additional_attributes]
    for p in range(0, len(feature_type)):
        for ifn, fn in enumerate(feature_attr[p]):
            fields = [fn] + attributes[p][fn] + [str(r[p]['counts'][fn]) for r in results]
            line = output_delimiter.join(fields)
            if output_filename == '':
                print(line)
            else:
                omode = 'a' if output_append or (ifn > 0) else 'w'
                with open(output_filename[p], omode) as f:
                    f.write(line)
                    f.write('\n')

        for title, fn in other_features:
            fields = [title] + pad + [str(r[p][fn]) for r in results]
            line = output_delimiter.join(fields)
            if output_filename == '':
                print(line)
            else:
                with open(output_filename[p], 'a') as f:
                    f.write(line)
                    f.write('\n')


def my_showwarning(message, category, filename, lineno=None, file=None,
                   line=None):
    sys.stderr.write("Warning: %s\n" % message)


def main():

    pa = argparse.ArgumentParser(
        usage="%(prog)s [options] alignment_file gff_file",
        description="This script takes one or more alignment files in SAM/BAM " +
        "format and a feature file in GFF format and calculates for each feature " +
        "the number of reads mapping to it. See " +
        "http://htseq.readthedocs.io/en/master/count.html for details.",
        epilog="Written by Simon Anders (sanders@fs.tum.de), " +
        "European Molecular Biology Laboratory (EMBL) and Fabio Zanini " +
        "(fabio.zanini@unsw.edu.au), UNSW Sydney. (c) 2010-2020. " +
        "Released under the terms of the GNU General Public License v3. " +
        "Part of the 'HTSeq' framework, version %s." % HTSeq.__version__)

    pa.add_argument(
            "--version", action="store_true",
            help='Show software version and exit')
    args, argv = pa.parse_known_args()
    # Version is the only case where the BAM and GTF files are optional
    if args.version:
        print(HTSeq.__version__)
        sys.exit()

    pa.add_argument(
            "samfilenames", nargs='+', type=str,
            help="Path to the SAM/BAM files containing the mapped reads. " +
            "If '-' is selected, read from standard input")

    pa.add_argument(
            "featuresfilename", type=str,
            help="Path to the GTF file containing the features")

    pa.add_argument(
            "-f", "--format", dest="samtype",
            choices=("sam", "bam", "auto"), default="auto",
            help="Type of <alignment_file> data. DEPRECATED: " +
            "file format is detected automatically. This option is ignored.")

    pa.add_argument(
            "-r", "--order", dest="order",
            choices=("pos", "name"), default="name",
            help="'pos' or 'name'. Sorting order of <alignment_file> (default: name). Paired-end sequencing " +
            "data must be sorted either by position or by read name, and the sorting order " +
            "must be specified. Ignored for single-end data.")

    pa.add_argument(
            "--max-reads-in-buffer", dest="max_buffer_size", type=int,
            default=30000000,
            help="When <alignment_file> is paired end sorted by position, " +
            "allow only so many reads to stay in memory until the mates are " +
            "found (raising this number will use more memory). Has no effect " +
            "for single end or paired end sorted by name")

    pa.add_argument(
            "-s", "--stranded", dest="stranded",
            choices=("yes", "no", "reverse"), default="yes",
            help="Whether the data is from a strand-specific assay. Specify 'yes', " +
            "'no', or 'reverse' (default: yes). " +
            "'reverse' means 'yes' with reversed strand interpretation")

    pa.add_argument(
            "-a", "--minaqual", type=int, dest="minaqual",
            default=10,
            help="Skip all reads with MAPQ alignment quality lower than the given " +
            "minimum value (default: 10). MAPQ is the 5th column of a SAM/BAM " +
            "file and its usage depends on the software used to map the reads.")

    pa.add_argument(
            "-t", "--type", type=str, dest="featuretype",
            default="exon",
            help="Feature type (3rd column in GTF file) to be used, " +
            "all features of other type are ignored (default, suitable for Ensembl " +
            "GTF files: exon)")

    pa.add_argument(
            "-i", "--idattr", type=str, dest="idattr",
            default="gene_id",
            help="GTF attribute to be used as feature ID (default, " +
            "suitable for Ensembl GTF files: gene_id). All feature of the " +
            "right type (see -t option) within the same GTF attribute will " +
            "be added together. The typical way of using this option is to " +
            "count all exonic reads from each gene and add the exons " +
            "but other uses are possible as well.")

    pa.add_argument(
            "--additional-attr", type=str,
            action='append',
            default=[],
            help="Additional feature attributes (default: none, " +
            "suitable for Ensembl GTF files: gene_name). Use multiple times " +
            "for more than one additional attribute. These attributes are " +
            "only used as annotations in the output, while the determination " +
            "of how the counts are added together is done based on option -i.")

    pa.add_argument(
            "-m", "--mode", dest="mode",
            # choices=("union", "intersection-strict", "intersection-nonempty"),
            default="union",
            help="Mode to handle reads overlapping more than one feature " +
            "(choices: union, intersection-strict, intersection-nonempty; default: union)")

    pa.add_argument(
            "--nonunique", dest="nonunique", type=str,
            choices=("none", "all", "fraction", "random"), default="none",
            help="Whether and how to score reads that are not uniquely aligned " +
            "or ambiguously assigned to features " +
            "(choices: none, all, fraction, random; default: none)")

    pa.add_argument(
            "--secondary-alignments", dest="secondary_alignments", type=str,
            choices=("score", "ignore"), default="ignore",
            help="Whether to score secondary alignments (0x100 flag)")

    pa.add_argument(
            "--supplementary-alignments", dest="supplementary_alignments", type=str,
            choices=("score", "ignore"), default="ignore",
            help="Whether to score supplementary alignments (0x800 flag)")

    pa.add_argument(
            "-o", "--samout", type=str, dest="samouts",
            action='append',
            default=[],
            help="Write out all SAM alignment records into " +
            "SAM/BAM files (one per input file needed), annotating each line " +
            "with its feature assignment (as an optional field with tag 'XF')" +
            ". See the -p option to use BAM instead of SAM.")

    pa.add_argument(
            "-p", '--samout-format', type=str, dest='samout_format',
            choices=('SAM', 'BAM', 'sam', 'bam'), default='SAM',
            help="Format to use with the --samout option."
            )

    pa.add_argument(
            "-d", '--delimiter', type=str, dest='output_delimiter',
            default='\t',
            help="Column delimiter in output (default: TAB)."
            )
    pa.add_argument(
            "-c", '--counts_output', type=str, dest='output_filename',
            default='',
            help="Filename to output the counts to instead of stdout."
            )

    pa.add_argument(
            '--append-output', action='store_true', dest='output_append',
            help='Append counts output to an existing file instead of ' +
            'creating a new one. This option is useful if you have ' +
            'already creates a TSV/CSV/similar file with a header for your ' +
            'samples (with additional columns for the feature name and any ' +
            'additionl attributes) and want to fill in the rest of the file.'
            )

    pa.add_argument(
            "-n", '--nprocesses', type=int, dest='nprocesses',
            default=1,
            help="Number of parallel CPU processes to use (default: 1). " +
            "This option is useful to process several input files at once. " +
            "Each file will use only 1 CPU. It is possible, of course, to " +
            "split a very large input SAM/BAM files into smaller chunks " +
            "upstream to make use of this option."
            )

    pa.add_argument(
            '--feature-query', type=str, dest='feature_query',
            default=None,
            help='Restrict to features descibed in this expression. Currently ' +
            'supports a single kind of expression: attribute == "one attr" to ' +
            'restrict the GFF to a single gene or transcript, e.g. ' +
            '--feature-query \'gene_name == "ACTB"\' - notice the single ' +
            'quotes around the argument of this option and the double ' +
            'quotes around the gene name. Broader queries might become ' +
            'available in the future.',
            )

    pa.add_argument(
            "-q", "--quiet", action="store_true", dest="quiet",
            help="Suppress progress report")  # and warnings" )

    args = pa.parse_args()

    # Parse multiple conditions into list
    args.featuretype = args.featuretype.split(",")
    args.idattr = args.idattr.split(",")
    args.mode = args.mode.split(",")
    args.output_filename = args.output_filename.split(",")

    warnings.showwarning = my_showwarning
    try:
        count_reads_in_features(
                args.samfilenames,
                args.featuresfilename,
                args.order,
                args.max_buffer_size,
                args.stranded,
                args.mode,
                args.nonunique,
                args.secondary_alignments,
                args.supplementary_alignments,
                args.featuretype,
                args.idattr,
                args.additional_attr,
                args.quiet,
                args.minaqual,
                args.samouts,
                args.samout_format,
                args.output_delimiter,
                args.output_filename,
                args.output_append,
                args.nprocesses,
                args.feature_query,
                )
    except:
        sys.stderr.write("  %s\n" % str(sys.exc_info()[1]))
        sys.stderr.write("  [Exception type: %s, raised in %s:%d]\n" %
                         (sys.exc_info()[1].__class__.__name__,
                          os.path.basename(traceback.extract_tb(
                              sys.exc_info()[2])[-1][0]),
                          traceback.extract_tb(sys.exc_info()[2])[-1][1]))
        sys.exit(1)


if __name__ == "__main__":
    main()