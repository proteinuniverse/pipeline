# pipeline

# stitch
the stitch code will take blast or diamond m8 output either as a file or from STDIN, and stitch together multiple lines for the same query-hit pair, provided the hits are non-overlapping

the stitched alignment length is the sum of the individual alignment lengths

the stitched %identity is calculated by determining the number of identical positions in each separate hit (%id*length), then summing the total number of identical positions, and dividing by the stitched alignment length to determine the %id for the stitched length

the stitched bit score is the sum of the individual bit scores

the evalue is not recalculated; the evalue reported in the stitched m8 output is the first evalue found for a given query-hit pair
