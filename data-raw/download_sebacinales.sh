#!/usr/bin/env bash
# download reannotated RDP training set sequences from Sebacinales
# the ones that are annotated to genus get kept as references.
# the rest have their names removed and are used for query sequences.

[ -d ../inst/extdata ] || mkdir -p ../inst/extdata
wget -O - "https://github.com/brendanf/reannotate/releases/download/v0.1/rdp_train.LSU.sintax.fasta.gz" |
zcat |
sed '/^>/!{{H; $!d}}; /^>/ x; $x; /^>.*o:Sebacinales/ !d; /unidentified/ d; 1 d' |
gzip -c > ../inst/extdata/sebacinales.sintax.fasta.gz

wget -O - "https://github.com/brendanf/reannotate/releases/download/v0.1/rdp_train.LSU.sintax.fasta.gz" |
zcat |
sed '/^>/!{{H; $!d}}; /^>/ x; $x; /^>.*o:Sebacinales/ !d; /unidentified/ !d; 1 d' |
awk '/^>/{gsub(/^>.*/,">Seq"i++);}1' i=1 |
gzip -c > ../inst/extdata/unknowns.fasta.gz

wget -O - "https://github.com/brendanf/reannotate/releases/download/v0.1/rdp_train.LSU.dada2.fasta.gz" |
zcat |
sed '/^>/!{{H; $!d}}; /^>/ x; $x; /^>.*Sebacinales/ !d; /unidentified/ d; 1 d' |
gzip -c > ../inst/extdata/sebacinales.dada2.fasta.gz