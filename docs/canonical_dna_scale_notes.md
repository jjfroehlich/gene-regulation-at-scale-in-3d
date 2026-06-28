# Canonical DNA Scale Notes

The current V5 scene uses the practical full-gene ACTB default: the canonical ACTB gene span plus an upstream promoter segment. It is not a full broad genomic ACTB locus.

Current V5 values:

- DNA: `3,954 bp` total, made from `3,454 bp` ACTB canonical gene span plus `500 bp` upstream promoter.
- DNA contour length: `537.744 mm` at `0.136 mm/bp`.
- mRNA: `1,852 nt`, `222.24 mm` contour length at `0.12 mm/nt`.
- Shared scene scale: `1 nm = 0.4 mm`.
- Working interpretation: practical promoter-plus-gene ACTB scene with exon/intron annotation and a full-length actin mRNA product.

Human ACTB comparison:

- Practical full-gene default: NCBI Gene ACTB / Gene ID 60 and the Ensembl canonical transcript span are about `3,454 bp`.
- Including a `500 bp` upstream promoter region gives `3,954 bp`.
- Ensembl canonical ACTB has 6 exons totaling about `1,812 bp`.
- The broad Ensembl ACTB gene object spans about `37,494 bp`, which would be about `5.10 m` at the current scene scale or `510 mm` at `1:10`; this is not the recommended default for the sculpture.

Legacy note:

- Older notes and preserved V3 output used a shorter `1,900 bp` DNA interpretation, close to an exon-scale actin-region representation. That is no longer the current scene interpretation.

Sources:

- https://www.ncbi.nlm.nih.gov/gene/60
- https://rest.ensembl.org/lookup/symbol/homo_sapiens/ACTB?expand=1;content-type=application/json
