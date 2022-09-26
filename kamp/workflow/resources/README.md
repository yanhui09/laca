This folder contains necessary resource which can't be easily integrated in a `conda` environment or a `singularity` container.

Due to our barcode design, we stick to `Guppy` `Version 3.3.3+fa743a6`, built in a `singularity` container.
If using custom barcode sets, consider preparing the barcode sequence file, e.g., [`barcodes_custome2.fasta`](guppy_barcoding/barcodes_custom2.fasta)
and corresponding config file, e.g., [`barcode_arrs_16S-GXO.cfg`](guppy_barcoding/barcode_arrs_16S-GXO.cfg), accordingly.
Remeber to choose the correct `barcode_kits` in the config file before you start analysis.
