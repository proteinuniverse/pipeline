# `stitch`

Input is sorted `m8` output of legacy `blast` or when using `blast+` or
`diamond`with:

```
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
```

Followed by:

```
sort -k1,1 -k2,2 -k7,7n blast_output.tsv > blast_output.sorted.tsv
```

## Non-overlapping Hits
Multiple, non-overlapping hits of the same query-subject pair are stitched
together.


```
node3	target3	90	100	10	1	11	110	901	1000	5	10
node3	target3	90	100	10	1	151	250	1101	1200	5	10
```

Becomes:

```
node3	target3	90.00	200	20	0	11	250	901	1200	5	20.0
```

## Overlapping Hits
Overlapping hits are merged with respect the query's overlap.

```
node1	target1	90	10	1	1	10	21	20	31	5	10
node1	target1	90	10	1	1	15	26	25	36	5	10
node1	target1	90	10	1	1	30	41	40	51	5	10
```

Becomes:

```
node1	target1	84.62	26	2	0	10	41	20	51	5	25.0
```

For this work as intended, input must be sorted by query sequence, subject
sequence, then by query start. To do so, you would:

```
sort -k1,1 -k2,2 -k7,7n blast_output.tsv > blast_output.sorted.tsv
./stitch blast_output.sorted.tsv > blast_output.sorted.stitched.tsv
```

## Columns
+ For overlapping hits: pident, mismatch, gapopen, and bitscore are scaled
with respect to the amount of overlap occurring on the query.
+ The stitched alignment length is the sum of the individual alignment lengths
+ The stitched pident is calculated by determining the number of
identical positions in each separate hit `(%id * length)`, then summing the
total number of identical positions, and dividing by the stitched alignment
length to determine the pident for the stitched length
+ bitscore is the sum of the individual bit scores
+ evalue is not recalculated; the evalue reported is the first evalue found
for a given hsp
