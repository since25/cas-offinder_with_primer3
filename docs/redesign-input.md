# Existing Off-Target Table Input

`otp.redesign --input` accepts Excel and CSV off-target result tables. Every non-blank row must contain `crRNA`, `DNA`, `Chromosome`, `Position`, `Direction`, `Mismatches`, and `Bulge Size`.

The first completely blank row marks the end of the off-target table; content after that row is ignored. A partially filled row before the blank separator is rejected with its data-row number, Excel row number, and the missing required columns so the source table can be corrected.
