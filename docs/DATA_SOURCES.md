# Dataset provenance

This repository includes the original meeting-2 snapshot so team members share
identical input data and split assignments. Metadata sidecars beside each raw
file record download time, source URL/query, available release information, file
size, and SHA256. `source_manifest.json` is a second copy of the source metadata.

- GO GAF files: https://geneontology.org/docs/download-go-annotations/downloads/
- GAF 2.2 format: https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
- Evidence codes: https://geneontology.org/docs/guide-go-evidence-codes/
- UniProt API: https://www.uniprot.org/api-documentation/uniprotkb
- UniProt terms: https://www.uniprot.org/help/license
- GO citation policy: https://geneontology.org/docs/go-citation-policy/
- MMseqs2 source, documentation, and license: https://github.com/soedinglab/MMseqs2

The recorded UniProt release is 2026_03; GAF headers reference GO 2026-07-26.
GO identifiers in the masters are direct annotations, without ancestor propagation.
This dataset was prepared by the team from these public sources, not released
as a new experimentally validated annotation database.

`reports/clusters_run.json` and `homology_audit_run.json` preserve tool versions,
commands and input/output hashes. Local paths in original records are historical;
they do not need to exist on a teammate's computer. New runs update the records.

Use the frozen source files instead of re-downloading current endpoints when
comparing model results. Keep source attribution when sharing derived data.
