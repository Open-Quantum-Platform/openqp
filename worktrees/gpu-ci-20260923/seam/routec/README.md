# OpenQP Route-C host integration

`int2_routec.F90` is the committed bridge from OpenQP tip
`d9a752d4002a` (full SHA and source hash in `docs/integration/source-files.json`).
The `patches/` directory preserves the six source-only integration commits in
this order: f3c10399b6, 4ca6cc97e8, 8b3b5ab481, 438a4febed, a9f66e9eaa, 4d8afd6e12.
They target historical OpenQP, not this standalone repository. Check the current
internal OpenQP host before applying any patch. The relevant exported GPU
implementations already live in this repository's `src/` and public header.
