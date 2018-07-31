/**
Simply remove Casava 1.8 vendor quality failed reads from
fastq files. Both raw and gzip files are ok.

Kyle Hernandez
**/
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>

#include "seqtk-1.3/kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef int bool;
#define true 1
#define false 0

bool is_failed(kseq_t *seq)
{
    bool ret = false;

    if (seq->comment.l) {
        char *s = index(seq->comment.s, ':');
        if (s && *(++s) == 'Y') {
            ret = true;
        }
    }
    return ret;
}

void write_record(kseq_t *seq, gzFile op )
{
    // seqid
    gzputc(op, '@');
    gzputs(op, seq->name.s);
    if (seq->comment.l) {
        gzputc(op, ' ');
        gzputs(op, seq->comment.s);
    }
    gzputc(op, '\n');
    // sequence
    gzputs(op, seq->seq.s);
    gzputc(op, '\n');
    // comment
    gzputs(op, "+\n");
    // qual
    gzputs(op, seq->qual.s);
    gzputc(op, '\n');
}

int next_pe(kseq_t *seqf, kseq_t *seqr)
{
    int rf = kseq_read(seqf);
    int rr = kseq_read(seqr);
    if ( rf >= 0 && rr >= 0 ) {
        return 1;
    }
    else if ( rf >= 0 && rr < 0 ) {
        return -2;
    }
    else if ( rf < 0 && rr >= 0 ) {
        return -3;
    }
    else {
        return rf;
    }
}

int process_se(const char *f, const char *opfx)
{
    gzFile fp;
    gzFile op;
    kseq_t *seq;
    int total = 0;
    int removed = 0;

    fp = gzopen(f, "r");

    if (fp == 0) {
        fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
        return 1;
    }

    const char * suffix = "_R1.fq.gz";
    char * outfqname = malloc(strlen(opfx) + strlen(suffix) + 1);
    strcpy(outfqname, opfx); 
    strcat(outfqname, suffix); 

    fprintf(stderr, "Writing SE Fastq to: %s\n", outfqname);
    op = gzopen(outfqname, "wb");
    free(outfqname);

    if (op == 0) {
        fprintf(stderr, "[E::%s] failed to open the output file/stream.\n", __func__);
        return 1;
    }
    seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        total++;

        if ( total % 1000000 == 0 ) {
            fprintf(stderr, "Processed %d records...\n", total);
        }

        if (is_failed(seq)) {
            removed++;
            continue;
        } 

        // Write
        write_record( seq, op );
    }

    // Cleanup
    kseq_destroy(seq);
    gzclose(fp);
    gzclose(op);

    fprintf(stderr, "\n");
    fprintf(stderr, "Processed %d records.\n", total);
    fprintf(stderr, "Removed %d records.\n", removed);
    return 0;
}

int process_pe(const char *ff, const char *fr, const char *opfx)
{
    gzFile ffp;
    gzFile frp;
    gzFile ofp;
    gzFile orp;

    kseq_t *seqf;
    kseq_t *seqr;
    int total = 0;
    int removed = 0;


    ffp = gzopen(ff, "r");
    if (ffp == 0) {
        fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
        return 1;
    }

    frp = gzopen(fr, "r");
    if (frp == 0) {
        fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
        return 1;
    }

    const char * suffixf = "_R1.fq.gz";
    char * outfq1name = malloc(strlen(opfx) + strlen(suffixf) + 1);
    strcpy(outfq1name, opfx); 
    strcat(outfq1name, suffixf); 
    fprintf(stderr, "Writing PE Fastq R1 to: %s\n", outfq1name);

    ofp = gzopen(outfq1name, "wb");
    if (ofp == 0) {
        fprintf(stderr, "[E::%s] failed to open the output file/stream.\n", __func__);
        return 1;
    }
    free(outfq1name);

    const char * suffixr = "_R2.fq.gz";
    char * outfq2name = malloc(strlen(opfx) + strlen(suffixr) + 1);
    strcpy(outfq2name, opfx); 
    strcat(outfq2name, suffixr); 
    fprintf(stderr, "Writing PE Fastq R2 to: %s\n", outfq2name);

    orp = gzopen(outfq2name, "wb");
    if (orp == 0) {
        fprintf(stderr, "[E::%s] failed to open the output file/stream.\n", __func__);
        return 1;
    }
    free(outfq2name);

    seqf = kseq_init(ffp);
    seqr = kseq_init(frp);

    int state = 0;
    while ((state = next_pe(seqf, seqr)) >= 0)
    {
        total++;

        if (is_failed(seqf) || is_failed(seqr)) {
            removed++;
            continue;
        }

        else {
            write_record(seqf, ofp);
            write_record(seqr, orp);
        }
    }

    // Cleanup
    kseq_destroy(seqf);
    kseq_destroy(seqr);
    gzclose(ffp);
    gzclose(frp);
    gzclose(ofp);
    gzclose(orp);

    if( state != -1 ) {
        fprintf(stderr, "[E::] Uneven R1 and R2 files!\n");
        return 1;
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "Processed %d records.\n", total);
    fprintf(stderr, "Removed %d records.\n", removed);
    return 0;
}

int main(int argc, char *argv[])
{

    const char *outpfx;
    int c;
    while ((c = getopt(argc, argv, "o:")) >= 0) {
        switch (c) {
            case 'o': outpfx = optarg; break; 
        }
    }
    if (optind == argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage: fqvendorfail -o <output prefix> <R1.fq> [<R2.fq>]\n");
        fprintf(stderr, "    -o STRING output prefix for filtered fastq files.\n");
        fprintf(stderr, "\n");
        return 1;
    }

    if(argc-optind == 1) {
        return process_se(argv[optind], outpfx);
    }
    else if(argc-optind == 2) {
        return process_pe(argv[optind], argv[optind+1], outpfx);
    }

    return 0;
}
