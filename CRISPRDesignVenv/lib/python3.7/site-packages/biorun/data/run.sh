# Default alignment.
bio align GATTACA GATCA > gattaca1.txt
# Default alignment.
bio align GATTACA GATCA --global > gattaca2.txt
# Default alignment.
bio align GATTACA GATCA --local > gattaca3.txt
bio align GATTACA GATCA --vcf > gattaca.vcf
bio align GATTACA GATCA --diff  > gattaca.diff
bio align align_input.fa --vcf > align_input.vcf
bio align GATTACA GTTAACA GTTTATA GTTT --pile > pile_1.txt
bio align s.fa > align-s-pairwise.txt
bio align s.fa --table > align-s-table.txt
bio align s.fa --vcf > align-s.vcf
