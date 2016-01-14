# `stitch`

Input is `m8` output of legacy `blast` or when using `blast+` or `diamond`
with:

```
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
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
node3	target3	90.00	200	10	0	11	250	901	1200	5	20.0
```

## Overlapping Hits
Overlapping hits after the initial observation are removed and can
have weird effects in the output.

```
node1	target1	90	100	10	1	100	1	901	1000	5	10
node1	target1	90	100	10	1	51	150	1050	951	5	10
```

Becomes:

```
node1	target1	0.00	18446744073709551518	10	0	100	1	901	1000	5	10.0
```

## Columns
+ The stitched alignment length is the sum of the individual alignment lengths
+ The stitched percent-identity is calculated by determining the number of
identical positions in each separate hit `(%id * length)`, then summing the
total number of identical positions, and dividing by the stitched alignment
length to determine the %id for the stitched length
+ The stitched bit score is the sum of the individual bit scores
+ The e-value is not recalculated; the e-value reported in the stitched output
is the first e-value found for a given query-subject pair
