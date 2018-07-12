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

int process_se(const char *f, const char *opfx)
{
    gzFile fp;
    gzFile op;
    kseq_t *seq;
    int total = 0;
    int removed = 0;
    const char * suffix = "_R1.fq.gz";
    char * outfqname = malloc(strlen(opfx) + strlen(suffix) + 1);
    strcpy(outfqname, opfx); 
    strcat(outfqname, suffix); 
    fprintf(stderr, "%s\n", outfqname);
    fp = strcmp(f, "-")? gzopen(f, "r") : gzdopen(fileno(stdin), "r");

    if (fp == 0) {
        fprintf(stderr, "[E::%s] failed to open the input file/stream.\n", __func__);
        return 1;
    }

    op = gzopen(outfqname, "wb");
    if (op == 0) {
        fprintf(stderr, "[E::%s] failed to open the output file/stream.\n", __func__);
        return 1;
    }
    seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        total++;
         
        if (seq->comment.l) {
            char *s = index(seq->comment.s, ':'); 
            if (s && *(++s) == 'Y') {
                removed++;
                continue;
            }
        }

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

    free(outfqname);
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
    const char * suffixf = "_R1.fq.gz";
    const char * suffixr = "_R2.fq.gz";
    char * outfq1name = malloc(strlen(opfx) + strlen(suffixf) + 1);
    strcpy(outfq1name, opfx); 
    strcat(outfq1name, suffixf); 
    char * outfq2name = malloc(strlen(opfx) + strlen(suffixr) + 1);
    strcpy(outfq2name, opfx); 
    strcat(outfq2name, suffixr); 

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

    ofp = gzopen(outfq1name, "wb");
    if (ofp == 0) {
        fprintf(stderr, "[E::%s] failed to open the output file/stream.\n", __func__);
        return 1;
    }
    free(outfq1name);

    orp = gzopen(outfq2name, "wb");
    if (orp == 0) {
        fprintf(stderr, "[E::%s] failed to open the output file/stream.\n", __func__);
        return 1;
    }
    free(outfq2name);

    seqf = kseq_init(ffp);
    seqr = kseq_init(frp);

    while (kseq_read(seqf) >= 0) {
        if (kseq_read(seqr) < 0) {
            fprintf(stderr, "[E::] different file lengths.\n");
            gzclose(ffp);
            gzclose(frp);
            gzclose(ofp);
            gzclose(orp);
            return 1;
        }

        total++;

        if (strcmp(seqf->name.s, seqr->name.s) != 0) {
            char *sf = index(seqf->name.s, '/'); 
            char *sr = index(seqr->name.s, '/'); 
            if (sf && sr) {
              char * fname = malloc(sf - seqf->name.s);
              for (int i = 0; i<sf-seqf->name.s; i++) {
                  fname[i] = seqf->name.s[i];
              }

              char * rname = malloc(sr - seqr->name.s);
              for (int i = 0; i<sr-seqr->name.s; i++) {
                  rname[i] = seqr->name.s[i];
              }
              fprintf(stderr, "%s %s\n", fname, rname);
              if(strcmp(fname, rname) != 0) {
                  free(fname);
                  free(rname);
                  fprintf(stderr, "[E::] Seqname mismatch! %s %s.\n", seqf->name.s, seqr->name.s);
                  gzclose(ffp);
                  gzclose(frp);
                  gzclose(ofp);
                  gzclose(orp);
                  return 1;
              } 
              free(fname);
              free(rname);
            }
            else { 
                fprintf(stderr, "[E::] Seqname mismatch! %s %s.\n", seqf->name.s, seqr->name.s);
                gzclose(ffp);
                gzclose(frp);
                gzclose(ofp);
                gzclose(orp);
                return 1;
            }
        }

        if (seqf->comment.l) {
            char *s = index(seqf->comment.s, ':'); 
            if (s && *(++s) == 'Y') {
                removed++;
                continue;
            }
        }

        else if (seqr->comment.l) {
            char *s = index(seqr->comment.s, ':'); 
            if (s && *(++s) == 'Y') {
                removed++;
                continue;
            }
        }

        // Read 1
        // seqid
        gzputc(ofp, '@');
        gzputs(ofp, seqf->name.s);
        if (seqf->comment.l) {
            gzputc(ofp, ' ');
            gzputs(ofp, seqf->comment.s);
        }
        gzputc(ofp, '\n');
        // sequence
        gzputs(ofp, seqf->seq.s);
        gzputc(ofp, '\n');
        // comment
        gzputs(ofp, "+\n");
        // qual
        gzputs(ofp, seqf->qual.s);
        gzputc(ofp, '\n');

        // Read 2
        // seqid
        gzputc(orp, '@');
        gzputs(orp, seqr->name.s);
        if (seqr->comment.l) {
            gzputc(orp, ' ');
            gzputs(orp, seqr->comment.s);
        }
        gzputc(orp, '\n');
        // sequence
        gzputs(orp, seqr->seq.s);
        gzputc(orp, '\n');
        // comment
        gzputs(orp, "+\n");
        // qual
        gzputs(orp, seqr->qual.s);
        gzputc(orp, '\n');
    }

    kseq_destroy(seqf);
    kseq_destroy(seqr);
    gzclose(ffp);
    gzclose(frp);
    gzclose(ofp);
    gzclose(orp);

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
