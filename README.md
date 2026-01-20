**eRNAFinder** is a simple R package to identify putative eRNAs from StringTie output (work in progress). 

### Transcript Filtering Strategy

Transcripts overlapping with GENCODE v49–annotated transcripts on the same strand are removed, except for lncRNAs. Remaining transcripts overlapping with upstream antisense RNAs (uaRNAs)—defined as RNA transcribed from the −300 bp to +1000 bp region of protein-coding gene TSSs on the opposite strand—are further excluded.

This strategy was adapted from Lee, Joo-Hyung et al., *Molecular Cell* (2021).
