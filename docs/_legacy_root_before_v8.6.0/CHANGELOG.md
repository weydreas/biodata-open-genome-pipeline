# Dante WGS pipeline — CHANGELOG

Единая история релизов, технических изменений и зафиксированных performance-настроек проекта.

## Текущий релиз: v8.5.2
- Аналитические контракты Stage 1–5, схемы Stage 15/16, пороги calling и форматы результатов **не изменены**.
- Удалены разрозненные `PATCH_NOTES_*.txt`, `RELEASE_MANIFEST_*.txt` и `HOME_COMPUTER_WSL_PERFORMANCE.txt`.
- Их содержимое объединено в этот `CHANGELOG.md`.
- Текущая N-safe логика Stage 16 из v8.2.1 сохранена без изменений.


## v8.2.3

- В корне оставлен один пользовательский launcher: `start.cmd`. Удалены `RUN_DANTE_MINIMAL.cmd`, `RUN_DANTE_WGS.cmd` и `RUN_TESTS.cmd`; тесты разработчика по-прежнему доступны через `Run-Tests.ps1`. При первом запуске `start.cmd` также удаляет эти legacy `.cmd` и старые раздельные `PATCH_NOTES_*.txt` / `RELEASE_MANIFEST_*.txt`, если они остались в рабочей папке после overlay-обновлений.
- Исправлена упаковка релиза: ZIP теперь предназначен для распаковки **непосредственно в корень проекта**, без дополнительной вложенной папки версии. Это исключает ситуацию, когда новый ZIP лежит рядом, а запускаются старые `.py/.ps1`.
- Stage 16 получил fail-fast software-contract probe: PowerShell v8.2.3 не запустит старый `reference_calling.py` и явно сообщит о смешанной установке.
- N-safe GRCh37 логика v8.2.1/v8.2.2 подтверждена и сохранена: loci с `REF=N` исключаются из constrained calling, записываются в `Dante_reference_calls_uncallable.tsv`, а pipeline продолжает работу.
- Аналитические схемы Stage 16 остаются `reference-calls-v3-known-alleles-ploidy-nsafe` / `reference-call-checkpoints-v3-known-alleles-ploidy-nsafe`; готовый Stage 15 BAM не инвалидируется.

## Текущий измеренный performance baseline — «Домашний компьютер»

Домашний компьютер — рекомендуемый WSL/Docker memory ceiling

Текущая конфигурация: Ryzen 5 7500F, 12 logical CPUs, 32 GB RAM.
Docker уже видит все 12 CPU, поэтому CPU-limit менять не нужно.
WSL по умолчанию ограничен примерно половиной RAM; для genomics workload рекомендуется 24 GB.

Файл: %UserProfile%\.wslconfig

[wsl2]
memory=24GB
processors=12

После изменения выполнить в PowerShell:
  wsl --shutdown
и заново запустить Docker Desktop.

При переходе на 64 GB RAM рекомендуемый memory ceiling: 48GB.
Измеренный baseline Stage 15 для текущей станции: BWA 12 threads, samtools sort 4 threads, 1536M/thread, 10M-read checkpoints, double buffering.
Повышение samtools sort до 6 threads ухудшило wall-clock и не рекомендуется.
Pipeline v8.2.1 читает фактически доступные Docker CPU/RAM и использует этот MAX-профиль на 12-thread Home Computer.

## История релизов


## v8.5.2

- Fixed PanelApp HTTP 429 failure in research bundle installation.
- Added proactive ~1 request/sec PanelApp pacing.
- Added `Retry-After` handling and conservative 60–300 s cooldown.
- Increased PanelApp retry budget to 12 attempts.
- PanelApp SQLite snapshot commits every successful detail request.
- Partial snapshots resume by panel ID/version without restarting from zero.
- Patient layout and genomic analytical contracts are unchanged.



## v8.5.1

- Fixed Windows `WinError 123` before patient migration/selection.
- Root cause: `%~dp0` ends in `\`; passing that directly as a quoted native-process
  argument could escape the closing quote and produce a root ending in a literal `"`.
- `start.cmd` now normalizes `DANTE_ROOT` by removing its trailing backslash before
  invoking `patient_manager.py`.
- No patient data migration had occurred before this failure; no rollback is needed.
- Genomic analytical contracts are unchanged.



## v8.5.0

- Added multi-patient workspace isolation under `patients/<Patient Name>/`.
- `start.cmd` now runs the patient selector before any stage status is calculated.
- One-time legacy migration moves current source VCF/FASTQ, patient reference-panel inputs,
  `Dante_work`, `Dante_qc`, and `Dante_results` into `patients/Demo 1`.
- Shared `annotation_data` and `reference_data` remain at root and are reused by all patients.
- Patient-specific sample metadata (`DANTE_SAMPLE_DATE/YEAR/MATERIAL/SEX`) migrate from
  root `.env` to the Demo 1 patient `.env`; root `.env` is backed up first.
- FASTQ validation identity upgraded from path-dependent v2 to relocation-stable v3.
  A migrated v2 cache is adopted without rereading/hashing the full FASTQ stream.
- Input discovery, array-panel discovery, telomere metadata scan, results, QC, work state,
  and local genome browser are patient-scoped.
- New patients start with clean `input`, `Dante_work`, `Dante_qc`, and `Dante_results` spaces.
- Stage analytical schemas and existing completed genomic outputs are not recalculated by migration.



## v8.4.2

- Fixed research-bundle HTTPS robustness without disabling certificate verification.
- On Python `CERTIFICATE_VERIFY_FAILED`, downloader/API may fall back to system curl using the OS trust store; no insecure TLS mode is used.
- Removed UCSC Chain Files from automatic research installation after re-checking current UCSC terms: free chain-file use is limited to Independent Researchers and Nonprofit Organizations; for-profit-directed/benefiting use requires a commercial license.
- Research bundle READY contract now requires PanelApp + CADD only.
- UCSC and OMIM are represented as gated/not-installed sources in the manifest.
- Planned target GRCh37/GRCh38 reconciliation remains covered by the open dbSNP dual-build layer.
- Stage 1–19 analytical contracts and genome-browser behavior are unchanged.



## v8.4.1

- Fixed Stage 19 immediate argparse failure when sample date/year/material metadata are empty.
- PowerShell now appends `--sample-date`, `--sample-year`, and `--sample-material` only when non-empty.
- Python adds defensive `nargs='?'` handling for the same optional metadata flags.
- Audited the same optional-argument class in the current PowerShell launcher; other dynamic optional CLI arguments are already guarded before append.
- Telomere algorithm, FASTQ validation, Stage 15–17 analytical contracts, and browser behavior are unchanged.



## v8.4.0

- Added menu item 24: local IGV.js genome browser on `127.0.0.1`.
- Uses the exact pipeline GRCh37 FASTA, Stage 15 target BAM/BAI, Stage 16 calls VCF/TBI when available, and ClinVar GRCh37 when available.
- Adds local rsID lookup via a background SQLite index built from `Dante_reference_targets.tsv`.
- Genomic file HTTP routes are allowlisted; arbitrary local filesystem paths are never accepted.
- Flask 3.1.3 is installed on first browser launch if absent.
- IGV.js 3.8.5 is fetched once from the pinned official jsDelivr npm package and cached locally; genomic data are served only from localhost.
- Current browser BAM is target-only; empty inter-target regions must not be interpreted as absent original WGS coverage.
- Stage 1–17 analytical contracts are unchanged.



## v8.3.3

- Fixed Stage 17 WinError 32 at 97% while publishing `Dante_reference_FASTQ_evidence.tsv.gz.part`.
- Root cause: `gzip.GzipFile(fileobj=...)` does not close the externally supplied raw file object; the previous `elif` cleanup left that raw handle open in the child process.
- Stage 17 now explicitly closes all three I/O layers (TextIOWrapper, GzipFile, raw file object) before publication.
- All atomic publications now use a short Windows retry for transient WinError 5/32 locks from Defender/indexers.
- No Stage 15/16 recalculation and no A/B/C analytical-contract change.



## v8.3.2

- Fixed the v8.3.1 Stage 16 false "mixed/incomplete installation" error.
- Root cause: the software-release string, expected analytical contract, and regex used by the probe had drifted to three different versions.
- Stage 16 compatibility is now based only on the actual analytical schemas (`reference-calls-v5...` / `reference-call-checkpoints-v4...`), not on the package release number.
- Fixed a second false-stale condition: overlay extraction rewrote the 149-byte GRCh37 ploidy file and changed its mtime, which could mark a completed Stage 16 INCOMPLETE despite identical content. This file is now validated by name + length + SHA256.
- Stage 5 menu is now explicitly `REFERENCE KNOWLEDGE INSTALLATION (DOWNLOAD / SNAPSHOT ONLY)`.
- The menu warns that installed bundles do not create findings by themselves. A downstream enrichment stage must scan the full genomic outputs before AI/report publication.
- Stage 16/17 analytical thresholds and schemas are unchanged.



## v8.3.1

- Stage 6 split into two legally isolated bundles: item 22 commercial/open-safe and item 23 research-only.
- New research bundle is human-only (`Homo sapiens`), `RESEARCH_ONLY`, `commercial_allowed=false`.
- Includes Genomics England PanelApp current-panel SQLite snapshot, CADD v1.7 GRCh37 tabix index + remote-score metadata, and UCSC hg19↔hg38 liftOver chain files.
- PanelApp restrictions are encoded separately: absent a separate GEL agreement, no commercial use, commercial research, diagnostic use, medical decision-making, or healthcare-service use.
- CADD full 79 GiB score file is not downloaded; future enrichment can query only requested loci.
- OMIM remains optional/gated and is not redistributed or auto-downloaded without user credentials/acceptance of its separate terms.
- Added `report_publication_guard.py` contract for future PDF publication: COMMERCIAL blocks actually-used research-only sources; RESEARCH permits them with mandatory restrictions/attributions.
- Stage 1–17 analytical contracts are unchanged.


## v8.3.0

- Added Stage 6 / item 22: resumable installer for open reference knowledge.
- Sources: HGNC, GenCC, ClinGen, selected GIAB v3.6 GRCh37 stratifications, AlphaMissense hg19, Mondo, HPO ontology-only, and dbSNP GRCh37/GRCh38 remote-index metadata.
- No new Python runtime dependencies: installer uses the standard library plus the already-approved bcftools image for AlphaMissense tabix indexing.
- Adds a machine-readable license/attribution manifest and `reference_data/open_reference/ATTRIBUTIONS.md` for later PDF legal attribution.
- Full dbSNP VCFs (26 GiB GRCh37 + 28 GiB GRCh38) are intentionally not downloaded; only indexes, assembly reports, and remote VCF metadata are retained for later target-sliced normalization.
- HPO is ontology-only; third-party-derived disease/gene annotation files are excluded.
- Stage 16/17 analytical contracts are unchanged from v8.2.6.

## v8.2.6

- Stage 17 switched from one strict pass/no-call threshold to explicit A/B/C FASTQ confidence tiers.
- Tier A HIGH: DP>=10, GQ/REFQ>=30; homozygous AF>=0.90; heterozygous >=3 reads per called allele and AF 0.20-0.80.
- Tier B SUPPORTIVE: DP>=5, GQ/REFQ>=10; homozygous >=4 supporting reads and AF>=0.80; heterozygous >=2 reads per called allele and AF 0.20-0.80.
- Tier C HYPOTHESIS: DP>=5, GQ/REFQ>=5 and >=1 observed read for every called allele; retained for research review only.
- Default reference-complete genotype export publishes Tier A+B. `DANTE_REFERENCE_EXPORT_MIN_TIER=A` can restrict the main export to Tier A without discarding B/C evidence.
- New user-facing research artifact: `Dante_results/Dante_reference_FASTQ_evidence.tsv.gz`, containing genotype, tier, DP, GQ/REFQ metric/value, AD and allele fractions.
- Provenance now records the same tier/quality metrics.
- Stage 17 state schema becomes `reference-export-v4-tiered`; Stage 16 v5 calls and v4 checkpoints remain fully compatible and are not recalculated.
- Fixed the Stage 16 completion summary to report usable GQ/REFQ coverage rather than the removed `hom_ref_gq_fraction` field.

## v8.2.5

- Fixed the Stage 16 semantic-QC failure at the end of an otherwise successful run.
- Root cause: bcftools multiallelic caller does not calculate FORMAT/GQ for non-variant homozygous-reference records; with `-A` it retains constrained PL likelihoods instead.
- Homozygous-reference confidence is now evaluated from native GQ when present, otherwise from a PL-derived reference-quality: likelihood gap between the reference genotype and the best non-reference genotype.
- ALT-carried genotypes continue to require native bcftools GQ.
- Stage 17 uses the same confidence rule, so strict `min_gq=30` remains unchanged in meaning.
- The 13 completed v8.2.4 Stage 16 checkpoints remain valid and are intentionally reused; only final concat/QC/publication is repeated.
- Stage 16 checkpoint schema remains v4. Final call/state schema becomes `reference-calls-v5-known-alleles-ploidy-refq`.

## v8.2.4

- Stage 16 reference audit is now coordinate-safe, not only N-safe.
- `NON_ACGT_REFERENCE`, `OUTSIDE_FASTA_CONTIG` and `MISSING_FASTA_CONTIG` loci are reported and skipped instead of aborting a multi-million-locus run.
- No REF allele is invented for an invalid GRCh37 coordinate.
- The Stage 15 BAM remains valid; only the Stage 16 call/checkpoint contract is bumped.

## v8.2.0

- Stage 16 switched to explicit known REF/ALT hypotheses using constrained bcftools calling.
- Added homozygous-reference GQ validation and sex-aware GRCh37 ploidy handling.
- Kept resumable Stage 16 checkpoint architecture introduced in v8.1.0.
- Superseded immediately by v8.2.1 after production revealed legitimate GRCh37 `N` positions that required N-safe handling.

## v8.2.1

Stage 16 GRCh37 ambiguous-reference hotfix.

Why:
- GRCh37 legitimately contains N bases in assembly gaps.
- v8.2.0 treated any non-ACGT FASTA symbol at a target coordinate as a fatal error.
- Example observed in production: chr1:3847570 -> N.

New contract:
- every unique target locus is audited against the exact GRCh37 FASTA before calling;
- A/C/G/T reference loci continue through known-alleles + GQ + sex-aware ploidy calling;
- non-ACGT reference loci are NOT assigned an invented REF allele and are NOT sent to bcftools;
- they are written to Dante_qc/Dante_reference_calls_uncallable.tsv;
- Stage 17 still sees the original master target panel: an agreed external/reference-source genotype may be used, otherwise the locus remains --;
- Stage 16 progress/ETA denominator is the actually callable A/C/G/T target set;
- sex/ploidy inference remains normalized to the same full target BED used to build the target-filtered BAM.

Schema bump:
- checkpoints: reference-call-checkpoints-v3-known-alleles-ploidy-nsafe
- final calls: reference-calls-v3-known-alleles-ploidy-nsafe

v8.2.0 Stage 16 checkpoints are intentionally stale. Stage 15 BAM remains valid and is not recalculated.

### Historical release manifest

Base: v8.2.0
Stage 15 analytical contract: unchanged.
Stage 16:
- explicit known REF/ALT hypotheses;
- homozygous-reference GQ contract;
- GRCh37 sex-aware ploidy;
- resumable 250k-locus checkpoints;
- Docker transport recovery;
- non-ACGT GRCh37 reference loci skipped safely and reported, never invented.

Stage 16 checkpoints: reference-call-checkpoints-v3-known-alleles-ploidy-nsafe
Stage 16 final state: reference-calls-v3-known-alleles-ploidy-nsafe

### Qc:
- Dante_reference_calls.metadata.json
- Dante_reference_calls_uncallable.tsv

## v8.1.0

Stage 16 resumable reference-aware calling
- Target BED is processed in contiguous checkpoints, 250,000 target loci by default.
- Each completed checkpoint is a validated bgzip VCF + tabix index and is recorded atomically in checkpoint_manifest.json.
- Crashes/restarts reuse all completed checkpoints; only the current unfinished block is recomputed.
- Progress and ETA are based on completed target loci, not output-file size.
- Transient Docker/WSL exit 125 is retried up to 3 times for the current chunk.
- Disk guard estimates remaining checkpoint storage and safe-stops before exhausting the volume.
- Final VCF is published only after all checkpoint VCFs are concatenated, indexed and fully readable.
- Stage 16 checkpoints are deleted only after outer PowerShell provenance/state validation passes.
- Existing complete reference-calls-v1 outputs remain compatible; an unfinished legacy monolithic .part is intentionally not reused because it has no exact completed-target boundary.

### Stage 15 baseline consolidated
- BWA 12 threads.
- samtools sort 4 threads, 1536M per thread.
- 10M-read resumable checkpoints.
- double-buffered single-end FASTQ staging.
- This is the measured wall-clock winner on the Home Computer; sort=6 experiments were slower.

### Historical release manifest

Base: v8.0.1 + measured Stage 15 performance hotfixes.
Primary change: resumable Stage 16 target-genotype calling with validated VCF/TBI checkpoints.
Stage 15 baseline: BWA 12 / samtools sort 4 / 1536M per sort thread / 10M checkpoints / double buffer.
All v8.0 genome-wide uncertain-coding/CHESS and Stage 4 telomere features retained.

## v8.0.1

Stage 15 Docker resilience hotfix
- Detect Docker exit code 125 / BrokenPipe as infrastructure transport failures.
- Preserve all completed BAM checkpoints.
- Wait up to 120 seconds for Docker Engine recovery.
- Retry the current single-end chunk automatically, up to 3 attempts total.
- Rewind the FASTQ stream only to the current chunk boundary; previously completed alignment checkpoints are not recomputed.
- Apply the same transient-error classification to paired-end chunk containers.
- Existing 5M/20M v8 checkpoint manifests remain compatible.

Analytical alignment/calling thresholds and reference contracts are unchanged.

### Historical release manifest

Base: v8.0
Primary change: Stage 15 transient Docker/WSL engine recovery and automatic current-chunk retry.
All v8.0 genome-wide uncertain-coding/CHESS and telomere features retained.

## v8.0

### Full Release / Consolidation
- Stage 2 item 11: genome-wide VUS/unclassified coding+splice sweep; no priority gene list.
- Embedded CHESS 2.2 full reference: 42,611 rows, 38,173 unique named symbols.
- Dedicated result: Dante_results\Dante_genomewide_uncertain_coding.tsv.
- CHESS metadata appended to shortlist candidates; CHESS coordinates are not used for GRCh37 mapping.
- Explicit ClinVar VUS functional variants are included without SIFT/PolyPhen gating.
- Unclassified variants use the same functional classes with max AF <= 1% or unknown frequency.

### Stage 15
- Resumable coordinate-BAM checkpoints retained.
- Single-end checkpoint cadence: 20,000,000 records. Existing v7.5.11 5M checkpoints remain resumable.
- Single-end input streams from FASTQ(.gz) directly to BWA stdin; no uncompressed input_chunk FASTQ.
- Home Computer MAX profile: BWA 12 threads, samtools sort 4, index/VEP 12, bcftools 10.
- Disk guard and checkpoint growth projection retained.
- RESUME SCAN reports progress while skipping already-completed FASTQ records.
- Long-running Stage 15 Docker containers carry private labels/names; startup/finally cleanup removes only v8 Stage 15 orphans.
- Generic watched-output heartbeat prints once and only repeats when output state changes.

### Analytical Contract
- GRCh37 reference, BWA/samtools/bcftools versions and strict Stage 16 genotype thresholds are unchanged.
- Chunk size/thread count/streaming change execution and recovery only, not which FASTQ reads are aligned.
- Whole FASTQ is still aligned; target BED filtering remains after global merge/markdup.

### Historical release manifest

```text
ee0e5571e89fc041689bb813c586012aa0b9086ddde9f8b2204003df7a6dcb2c  dante_wgs_pipeline.ps1
1fab2a25bf9bc00d528bdc4607b50496570ee400af5651b188972d9e17ea5b23  reference_alignment.py
c8360e54f2c4576596d48efcfc6c733133886686169675e33a2bde5924da139e  dante_activity.py
f01884618ccc5c12e463a505fee277e689ac84219b9ae43c15d2c67ecb9d7519  build_shortlist.py
baf8bafce2a5bf2039e7262667edfdfabae5f6131f487be770a8304bdd8b973a  reference_data/CHESS_2.2_genes.tsv
7e40ee0663e6c6a5339d023efe30a71602769748d8d302867fed32714a7896f1  reference_data/CHESS_2.2_genes.xlsx

### Chess 2.2:
  rows excluding header: 42,611
  unique named symbols: 38,173

Validation:
  Python unittest suite: 48/48 passed
  Python syntax/static checks: passed
  Windows PowerShell parser: run automatically by RUN_DANTE_WGS.cmd on the target Windows host
```

## v7.5.11

### Stage 15 reliability/performance redesign
- Removes the whole-genome query-name `samtools sort -n` that exhausted C: in v7.5.10.
- Single-end R1: BWA -> coordinate sort per resumable chunk; no name collation/fixmate is needed for unpaired duplicate criteria.
- Paired-end R1+R2: BWA -> samtools collate -> fixmate -m -> coordinate sort per resumable chunk.
- Default checkpoint size: 5,000,000 FASTQ records per mate.
- Every completed chunk is atomically published as a full-genome coordinate-sorted BAM.
- Restart reuses all contiguous completed chunks; at most the currently incomplete chunk is redone.
- Resume may rescan/decompress the source FASTQ to reach the checkpoint boundary, but completed chunks are not realigned.
- Finalization streams `samtools merge` of all full-genome checkpoints through one global `samtools markdup`, then applies the target BED. Thus duplicate marking still sees reads across chunk boundaries.
- Heavy checkpoint BAMs are deleted only after the final target BAM+BAI are atomically published.
- Adds disk-space guard with projected checkpoint growth and a 15 GiB safety reserve. Low disk causes a SAFE STOP between chunks while preserving progress instead of allowing filesystem exhaustion.
- Optional `.env`: DANTE_REFERENCE_CHECKPOINT_DIR can place checkpoint storage on another SSD.

### Performance
- MAX profile on a 12-logical-CPU machine now selects BWA=10 threads and samtools=2 threads.
- 24 GiB WSL/Docker memory continues to select 768M per samtools sort thread.

### Analytical contract
- The entire FASTQ is still read and aligned to GRCh37.
- Chunking changes execution/restart behavior, not which reads are aligned.
- Full-genome checkpoint BAMs are merged before one global duplicate-marking pass.
- Only after global markdup is the final BAM restricted to the reference target BED, as before.
- Stage 16 thresholds and genotype/export logic are unchanged.

## v7.5.10

Stage 15 progress correction
- BWA MEM emits two messages for each batch: `[M::process] read N sequences` when a batch is read and `[M::mem_process_seqs] Processed N reads` when that same batch completes.
- v7.5.7-v7.5.9 counted both messages, approximately doubling the displayed alignment percentage and reads/s and shortening the displayed ETA.
- v7.5.10 advances progress only on the completed `[M::mem_process_seqs]` message.
- Analytical alignment/calling behavior is unchanged; this is a UI/progress-only correction.
- Existing Stage 15 execution is not resumable. FASTQ validation, Docker images and the persistent BWA index are reusable, but interrupted BWA/samtools stream work restarts from read 1.

## v7.5.9

### New Stage 4 - Telomere Calculation
- Menu item 19 estimates mean telomere length directly from original WGS FASTQ.
- R1-only is supported now; when R2 appears the same stage automatically switches to R1+R2 and marks the previous result stale.
- Primary method is a transparent TelSeq-compatible canonical-repeat estimator:
  TTAGGG / CCCTAA, k >= 7, 48-52% GC normalizer, human TelSeq constants 332,720,800 and 46 telomere ends.
- Primary published artifact: Dante_results\Dante_telomere_report.txt.
- Reusable raw metrics cache: Dante_qc\Dante_telomere_metrics.json.
- The report contains input fingerprints, read counts, read lengths, telomeric-read counts, GC-normalizer counts, k=5..10 sensitivity QC, method/limitations and metadata clues.
- Raw FASTQ lacks duplicate flags used by BAM-based TelSeq. The report explicitly labels the output as a TelSeq-compatible historical WGS estimate, not a laboratory assay.

### Date / Sample Provenance
- Optional .env fields: DANTE_SAMPLE_DATE, DANTE_SAMPLE_YEAR, DANTE_SAMPLE_MATERIAL.
- No collection date is guessed from a standard Illumina FASTQ header.
- VCF date-like header fields and filesystem mtimes are printed only as provenance clues.
- Changing date/material later regenerates the TXT report from cached telomere metrics without rereading the 100+ GiB FASTQ.

### Menu
- Stage 4 is now telomere calculation (item 19).
- AI evidence research moved to Stage 5 (item 20).
- Local stages 1+2 automation moved to item 21.

### Prospective / Regression Checks
- Single-end and paired-end telomere modes have separate stale/ready identity.
- R2 arrival invalidates only the telomere report/metrics identity that depends on FASTQ mode; other completed deterministic outputs retain their existing rules.
- Telomere scan has record-based progress and ETA using the exact FASTQ record counts cached by Stage 15 validation.
- Existing Stage 1-3 analytical contracts and performance profile are unchanged.

## v7.5.8

### Performance
- RUN_DANTE_WGS.cmd defaults to -PerformanceProfile max.
- Docker/WSL CPU and RAM are detected with docker info.
- On 12 CPUs the max profile allocates approximately 8 BWA threads + 3 samtools sort threads (one logical CPU reserved for OS/WSL scheduling), rather than the old hard cap of 4.
- samtools sort receives an explicit per-thread memory budget: 512M below 18 GiB Docker RAM, 768M at 18-32 GiB, 1024M at >=32 GiB.
- bcftools gets a separate CPU budget to avoid blindly giving both piped processes the full BWA count.
- Resource profile is printed once at startup.

### Analytical Contract
- No genotype thresholds, GRCh37 reference contract, tool versions, target panel, calling filters, or output schemas changed.
- Resource settings are performance-only and do not invalidate completed analytical outputs.

### Home Computer
- With current 32 GB RAM, recommended WSL memory ceiling is 24 GB.
- Docker already sees 12 CPUs; no additional Docker CPU setting is required.

## v7.5.7

### ### Root cause fixed
v7.5.6 removed implicit Docker ENTRYPOINT assumptions but used complex `sh -c` command-substitution strings in PowerShell preflight. On Windows this could be re-quoted by the native-command boundary and reach the container as a broken shell program (`unexpected EOF while looking for matching )`).

### ### Changes
1. BWA preflight now invokes `docker run --entrypoint bwa IMAGE` directly and parses the combined output in PowerShell; no shell and no command substitution.
2. samtools preflight uses direct `samtools --version`; fragile help-text parsing was removed because the binary is pinned to exact 1.22.1.
3. bcftools Stage 2 and Stage 16 preflights use direct `bcftools --version`; fragile option/help parsing was removed because the binary is pinned to exact 1.24.
4. bash availability is probed directly with `--entrypoint bash --version`.
5. Diagnostic errors report the specific tool, exit code and a bounded tail of actual probe output.
6. Prospective UI audit fixed BWA-MEM alignment progress parsing: real BWA messages such as `[M::mem_process_seqs] Processed 196608 reads ...` now advance ALIGNMENT %, speed and ETA.
7. Runtime analytical contracts, thresholds, target panel, BAM/VCF formats and state fingerprints are unchanged.

### ### Prospective audit scope
- Stage 15 tool resolution/version/capability probes.
- Stage 15 BWA index/MEM and samtools streaming/QC/index invocation shapes.
- Stage 16 bcftools version/mpileup/call checks and runtime calling/index/readability.
- Stage 17 export path (no Docker shell preflight; streaming merge retained).
- Stage 18 orchestration inherits 15-17 fixes.
- Stage 19 is API-only and unrelated to Docker shell quoting.

## v7.5.6

### ### Root cause fixed
v7.5.5 could successfully pull the approved BioContainer images and still fail its own preflight because container commands were invoked as `IMAGE sh -c ...` / `IMAGE bwa ...` / `IMAGE samtools ...`. That assumes a neutral image ENTRYPOINT. It also chained samtools `--help` calls with `&&`, so a valid command whose help path returned non-zero could be misclassified as missing.

### ### Changes
1. Explicit Docker --entrypoint for BWA, samtools and bcftools in preflight and runtime.
2. Robust version checks: BWA 0.7.19, samtools 1.22.1, bcftools 1.24.
3. Capability help output is captured with `|| true`; its text, not its exit code, determines option availability.
4. Samtools required subcommands are checked from the command listing; bash is checked because the streaming pipelines require pipefail.
5. Prospective audit fixed the same ENTRYPOINT assumption in:
   - Stage 2 bcftools normalization/indexing and reference-contig extraction;
   - Stage 15 BWA indexing/MEM, samtools streaming pipeline/QC/index/readability;
   - Stage 16 bcftools mpileup/call/index/readability.
6. Analytical tool contracts and output fingerprints are unchanged: this patch changes invocation robustness, not the analysis algorithm.

## v7.5.5

### ### Observed failure
Stage 15 completed the full R1 validation (326,330,259 records; 111.70 GiB decompressed) and then failed because Docker could not reach quay.io through the current CONNECT proxy/VPN path. This is a container-registry transport failure, not FASTQ corruption.

### ### Direct fix
1. BWA 0.7.19 now resolves from approved sources in order: existing local image -> Quay BioContainers -> Docker Hub StaPH-B (`staphb/bwa:0.7.19`).
2. samtools 1.22.1 uses the same model with `staphb/samtools:1.22.1`.
3. bcftools 1.24 uses the same model with `staphb/bcftools:1.24`.
4. Each selected image must pass an exact-version and required-capability check before analytical processing starts.
5. Registry source is separated from the analytical tool contract. Switching from Quay to the approved same-version Docker Hub image does not itself invalidate compatible state.
6. Stage 15 tool preflight now occurs before a new full FASTQ validation, so a missing/unreachable tool image fails fast instead of wasting an 8-minute 100+ GiB pass. Existing completed validation remains reusable.

### ### Prospective audit after the failing point
7. Stage 16 had the same single-Quay dependency for bcftools. Fixed prospectively with the exact 1.24 fallback and exact-version capability check.
8. Stage 8 also had direct unconditional bcftools pulls and would contact Quay even with completed local work. Replaced with the same local-first/multi-registry resolver and a stable `bcftools-1.24` tool contract fingerprint.
9. VEP and AnnotSV installer paths contained unconditional `docker pull` calls. They are now local-first; they do not contact a registry when the already-approved image exists locally.
10. Stage 17 loaded the complete FASTQ-call VCF into a Python dictionary. With a multi-million target panel this could consume several GiB RAM and fail after expensive alignment/calling. It now performs a sorted streaming merge with memory bounded to the current coordinate group.
11. Stages 14/17 now report merge progress from the exact number of target-panel rows, including rows/s and ETA.
12. Stage 16 now reports phase progress around mpileup/call, tabix indexing, full readability check and atomic publication. During the long bcftools call the existing heartbeat/file-growth monitor remains the truthful activity signal; no fabricated target percentage is shown.
13. Stage 18 inherits all fixes because it invokes items 13-17.
14. Stage 19 was reviewed: the OpenAI SDK is already local-first and installation uses retry only if missing; failure affects only optional AI research and does not invalidate local genome outputs.
15. Stage 20 was reviewed: its VEP/AnnotSV/bcftools registry paths now use the local-first resolver changes above.

### ### Restart behavior
The completed FASTQ validation state from v7.5.4 contains exact record counts and hashes and is intentionally compatible with v7.5.5. After replacing program files, item 15 should reuse validation rather than reread R1, then resolve BWA/samtools images and continue.

## v7.5.4

Stage 15 progress and Docker-registry resilience:
1. Full FASTQ validation now persists exact FASTQ record counts and decompressed byte counts together with content SHA256.
2. Stage 15 reports an overall progress scale: validation 0-5%, tool/reference preparation 5-8%, BWA read alignment 8-90%, samtools drain/target extraction 90-95%, BAM QC/index/finalization 95-100%.
3. During BWA MEM, progress is driven by BWA's processed-sequence batches against the exact validated FASTQ record count; processed/total reads, speed and ETA are printed.
4. Existing v7.5.3 FASTQ validation metadata without record counts is upgraded once by re-reading the FASTQ. Thereafter the exact counts are cached.
5. BWA and samtools pinned images are now local-first. If the exact pinned image exists in Docker, Stage 15 does not contact quay.io at all.
6. If a required pinned image is absent locally, Docker pull is retried with backoff. Failure produces a registry/VPN/proxy-specific error while preserving the completed FASTQ validation for restart.
7. The same local-first rule is applied to the pinned bcftools image used by Stage 16, preventing the analogous failure there.
8. The reference-alignment state schema remains v2 because the biological/alignment contract did not change; a successful existing BAM is not invalidated merely by this progress/UI patch.

Observed failure addressed:
  failed to resolve reference quay.io/...; reading HTTP CONNECT ... wsarecv ... connection aborted
This occurs after FASTQ validation and is a Docker registry/proxy connectivity failure, not FASTQ corruption.
