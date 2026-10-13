# User-reported scored-run implementation repair

The user reports all16 v2 GPU shards completed after correcting chromosome-boundary
padding. The input code identities are retained in the new audit provenance.
This CPU package does not apply the patch, edit old hashes, or rescore variants.

The reported correction was:

```python
# With a=max(0,start), b=min(chromosome_size,start+width):
mid = seq if b > a else ''
left = min(width, max(0, -start))
right = width - left - len(mid)
result = 'N' * left + mid + 'N' * right
assert len(result) == width
```

Keep the actual patched scoring source with the final reproducibility archive.
The snippet above records the user's explanation, not a separate verification
of source files that remain on Wahab.
