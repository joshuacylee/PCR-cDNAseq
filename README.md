# PCR-cDNAseq
Code used to map raw reads generated using ONT long-read PCR-cDNAseq

Raw reads (Q > 7, filtered by MinKNOW/Guppy basecalling) were first trimmed (cutadapt v3.4) and orientated (pychopper v2.7) before mapping to the reference genome (GRCh38) using minimap2 (v 2.24). Samtools (v 1.16) was used to filter out alignments with low mapping quality (MAPQ < 10). Reads per gene were subsequently counted using featureCounts (v 2.0.0) using long-read counting mode (-L). 
