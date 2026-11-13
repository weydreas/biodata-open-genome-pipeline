# Dante WGS v9.0.0 — технический контракт Stage 24 Knowledge Enrichment

Статус: реализовано в `src/knowledge_enrichment.py`.

## 1. Назначение

Stage 24 — первый этап, который применяет установленные в пунктах 22 и 23 knowledge-bundles к каноническим геномным результатам конкретного пациента. Этап детерминированный: OpenAI не используется.

Главная цепочка v9.0.0:

```text
Stage 1 CNV/SV + Stage 2 full small-variant catalog + optional Stage 17 FASTQ evidence
        ↓
22 open/commercial-safe bundle + 23 research-only bundle
        ↓
24 FULL-GENOME KNOWLEDGE ENRICHMENT
        ↓
open findings + research-only overlays + provenance/QC
        ↓
25 AI evidence research
```

Stage 24 обязан просканировать весь `Dante_small_variants_annotated.tsv.gz`. Однако полный каталог не копируется второй раз в SQLite: в patient knowledge DB материализуются только события, которые пересекли хотя бы одно knowledge-rule. Число реально просмотренных строк сохраняется в `scan_stats` и текстовом отчёте.

## 2. Preconditions / READY dependencies

Пункт 24 доступен только если:

1. Stage 1 AnnotSV CNV/SV outputs имеют `READY`;
2. Stage 2 unified catalog имеет `READY`;
3. пункт 22 open reference bundle имеет `READY`;
4. пункт 23 research-only bundle имеет `READY`.

Stage 17 не блокирует Stage 24. Если Stage 17 уже существует, его Tier A/B ALT evidence включается. Если Stage 17 создаётся или обновляется позже, Stage 24 становится `STALE` и должен быть пересчитан.

## 3. Канонические входы

Обязательные patient inputs:

- `Dante_results/Dante_small_variants_annotated.tsv.gz`;
- `Dante_results/Dante_CNV_AnnotSV.tsv`;
- `Dante_results/Dante_SV_AnnotSV.tsv`.

Условные inputs:

- `Dante_qc/Dante_reference_targets.tsv`;
- `Dante_results/Dante_reference_FASTQ_evidence.tsv.gz`;
- актуальный `annotation_data/clinvar_grch37/clinvar.vcf.gz` для exact lookup Stage-17-discovered alleles.

Reference inputs:

- `reference_data/open_reference/open_reference_manifest.json`;
- `reference_data/research_reference/research_reference_manifest.json`.

Assembly contract: GRCh37.

## 4. Shared knowledge index

При первом запуске Stage 24 создаётся общий для всех пациентов индекс:

`reference_data/knowledge_index/knowledge_index.sqlite3`

Он строится из точных snapshots пунктов 22/23 и пересобирается, если меняется hash любого bundle manifest или `RULES_SCHEMA`.

Таблицы:

- `genes` — канонические HGNC symbols;
- `gene_alias` — текущие/alias/previous symbols;
- `ontology_terms` — Mondo и HPO vocabulary;
- `gene_disease` — GenCC + ClinGen gene-disease validity;
- `actionability` — ClinGen Actionability;
- `dosage` — ClinGen dosage sensitivity;
- `panelapp_gene` — research-only PanelApp snapshot;
- `source_snapshot` / `meta` — provenance и schema state.

HPO в текущем bundle является только ontology vocabulary. Stage 24 не строит phenotype→gene inference и не изображает ontology как disease-gene evidence.

## 5. Источники и их роль

### OPEN / commercial-safe layer

- HGNC — canonical gene symbol / aliases;
- GenCC — gene-disease validity;
- ClinGen gene validity — gene-disease validity;
- ClinGen dosage — CNV/SV dosage evidence;
- ClinGen actionability — actionability context;
- ClinVar — exact allele clinical assertions (частично уже в Stage 2, дополнительно target lookup для новых Stage 17 alleles);
- AlphaMissense hg19 — research-only supporting missense predictor (CC BY-NC-SA 4.0);
- GIAB GRCh37 stratifications — technical-confidence flags;
- dbSNP GRCh37 target slice — актуализация rsID для knowledge-relevant alleles;
- Mondo/HPO — ontology vocabulary.

### RESEARCH ONLY layer

- PanelApp — research gene/panel overlay;
- CADD v1.7 GRCh37 — target-only research predictor overlay.

Research sources физически не записываются в open findings TSV и не должны загрязнять commercial publication manifest.

## 6. Полный scan малых вариантов

Каждая строка unified catalog учитывается в `catalog_rows`.

Базовый gate:

- `alt_carried=yes`;
- `vcf_filter=PASS`;
- canonical human chromosome 1–22/X/Y/MT.

По умолчанию rare threshold: `max_af <= 0.01`; неизвестная AF считается допустимой для rare-candidate screening.

Варианты, не пересёкшие knowledge-rules, не дублируются в `Dante_knowledge.sqlite3`, но факт полного прохода фиксируется счётчиками.

## 7. Canonical event identity

Small variants используют Stage-2 `variant_key` либо fallback:

`GRCh37:chrom:position:REF:ALT`

Structural events:

`GRCh37:CNV|SV:chrom:start:end:sv_type:AnnotSV_ID`

Sources объединяются в `source_layers`. Stage 17 может добавить новый allele-level event, отсутствовавший в исходном variant-only Dante VCF.

## 8. Stage 17 integration

При наличии Stage 17 используются только ALT genotypes Tier A/B. Tier C не создаёт open finding.

Приоритет данных остаётся совместимым с существующим pipeline:

`Dante VCF > FASTQ Tier A/B > согласованный array source > no-call`.

Новый Stage-17 allele проходит targeted current ClinVar lookup перед построением findings.

## 9. Finding reason codes

Stage 24 не создаёт искусственный «pathogenicity score». Используются дискретные reason codes.

### OPEN

- `K1_EXACT_CLINVAR_PLP` — exact allele с чистой ClinVar Pathogenic/Likely pathogenic интерпретацией;
- `K1_EXACT_CLINVAR_PLP_COMMON` — exact P/LP ClinVar assertion для частого allele при недостаточной review strength; сохраняется для проверки, но не получает высокий приоритет;
- `K1_EXACT_CLINVAR_MIXED` — aggregate ClinVar assertion одновременно содержит P/LP и association/risk/protective/drug-response смысл; не считается чистым P/LP;
- `K1_EXACT_CLINVAR_CONFLICT` — exact allele conflicting ClinVar assertion;
- `K1_EXACT_CLINVAR_ASSOCIATION` — exact allele association/risk/drug-response class, если он поддержан текущим ClinVar mapping;
- `K2_RARE_LOF_VALIDATED_GENE` — LoF с наблюдаемой population AF <= threshold в gene с Strong/Definitive gene-disease validity;
- `K2_LOF_VALIDATED_GENE_AF_UNKNOWN` — LoF в таком gene, но population AF отсутствует; это candidate с явно зафиксированной неопределённостью, а не «rare» claim;
- `K3_RARE_MISSENSE_VALIDATED_GENE` — missense с наблюдаемой population AF <= threshold в gene с Strong/Definitive validity;
- `K3_MISSENSE_VALIDATED_GENE_AF_UNKNOWN` — missense в таком gene без population AF; низкоприоритетный candidate;
- `K4_CLINGEN_ACTIONABLE` — rare plausible functional allele в ClinGen-actionable gene; наличие gene в списке само по себе не является pathogenicity claim;
- `K4_CLINGEN_ACTIONABLE_AF_UNKNOWN` — plausible functional allele в actionable gene без доступной population AF; не трактуется как rare;
- `K5_DOSAGE_SENSITIVE_CNV` — DEL/DUP с AnnotSV class >=3 и direction-specific ClinGen dosage evidence; gene-level evidence применяется только к соответствующему gene, а recurrent region — только при >=50% покрытии curated GRCh37 region;
- `K6_HIGH_EVIDENCE_SV_GENE` — AnnotSV class 4/5 structural event затрагивает gene с Strong/Definitive gene-disease validity;
- `K6_STRUCTURAL_GENE_CANDIDATE` — AnnotSV class 3 или unclassified structural event в таком gene; это candidate, а не high-evidence claim.

AlphaMissense не может попадать в open/commercial-safe finding. После exact target match создаётся отдельный research-only `K8_ALPHAMISSENSE_SUPPORT` overlay.

### RESEARCH ONLY

- `K7_PANELAPP_RESEARCH` — candidate allele/gene пересекается с PanelApp snapshot;
- `K8_CADD_SUPPORT_TARGET` — candidate small variant получил CADD target score.

PanelApp и CADD не повышают open clinical class и не превращают research overlay в commercial finding.

## 10. Причинность и inheritance boundaries

Gene-disease relations имеют наследование из исходных knowledge sources.

Для AR:

- одна HET ALT → `AR_SINGLE_HET_CARRIER_OR_CANDIDATE`;
- HOM ALT → `AR_HOMOZYGOUS_ALT`;
- >=2 HET alleles с численно подтверждённой редкостью и non-benign ClinVar status в одном AR gene → `possible_biallelic_unphased=yes`.

Последний статус означает только «возможны два allele candidates». Без phasing нельзя писать `compound heterozygous` и нельзя утверждать cis/trans.

## 11. GIAB technical confidence

GIAB difficult-region BEDs используются как QC, а не как биологический фильтр.

- small variant: проверяется locus;
- CNV/SV: проверяются оба breakpoint, а не вся протяжённость события.

При пересечении:

`technical_confidence=CAUTION`

и `giab_flags` содержит соответствующие strata. Finding не удаляется.

## 12. Remote indexed target queries

Stage 24 не скачивает гигантские базы целиком только ради небольшого candidate set.

После полного локального scan формируются target slices:

- AlphaMissense — только research-only relevant rare missense targets;
- dbSNP — только candidate variants без resolved rsID;
- CADD — только research candidate targets.

Запросы выполняются batched через pinned bcftools/htslib Docker runtime.

По умолчанию сетевой сбой этих lookup не уничтожает локальный full-genome result: Stage 24 публикуется как `READY / REMOTE WARNINGS`, а deferred sources записываются в state/report. Для fail-closed режима:

`DANTE_KNOWLEDGE_STRICT_REMOTE=1`

## 13. Patient knowledge database

Canonical Stage-24 DB:

`Dante_results/knowledge_enrichment/Dante_knowledge.sqlite3`

Таблицы:

- `event_candidates`;
- `open_findings`;
- `research_overlays`;
- `scan_stats`;
- `meta`.

DB содержит только knowledge-relevant events; lossless full genome остаётся в upstream catalog/VCF/BAM.

## 14. Published outputs

Каталог:

`Dante_results/knowledge_enrichment/`

Файлы:

1. `Dante_knowledge_findings_open.tsv`
   - commercial-safe/open findings;
   - без PanelApp/CADD fields;
   - содержит provenance, inheritance, ClinVar, actionability, dosage и GIAB; AlphaMissense исключён из open output.

2. `Dante_knowledge_findings_research.tsv`
   - research-only overlays;
   - явный banner `RESEARCH ONLY — NOT FOR COMMERCIAL USE`;
   - PanelApp и CADD сохраняются reason-specific, без взаимного «протекания» полей.

3. `Dante_knowledge_qc.tsv`
   - GIAB/difficult-region candidate events.

4. `Dante_knowledge_enrichment_report.txt`
   - full scan counts;
   - finding counts by reason;
   - remote lookup status;
   - interpretation boundaries.

5. `Dante_knowledge.sqlite3`
   - patient knowledge graph/index.

6. `Dante_knowledge_usage_manifest.json`
   - только реально использованные OPEN sources;
   - предназначен для commercial publication guard.

7. `Dante_knowledge_usage_research_manifest.json`
   - OPEN + research sources;
   - предназначен только для research publication mode.

8. `Dante_knowledge_enrichment.state.json`
   - READY/reproducibility contract.

## 15. READY/state contract

Schemas:

- enrichment: `dante-knowledge-enrichment-v1`;
- shared index: `dante-knowledge-reference-index-v1`;
- rules: `dante-knowledge-rules-v7`;
- usage manifest: `knowledge-usage-v1`.

State сохраняет:

- SHA256 genomic inputs;
- SHA256 open/research manifests;
- SHA256 shared index state;
- `max_af`;
- output name/size/SHA256;
- scan stats;
- remote deferred sources;
- runtime duration.

Повторный запуск переиспользует result только при совпадении signatures и при проверке size + SHA256 каждого output.

Меняется ClinGen/PanelApp/другая bundle snapshot → пересчитывается knowledge index и Stage 24, но Stage 1–17 не пересчитываются.

Меняется Stage 17/AnnotSV/ClinVar output → Stage 24 становится `STALE`.

## 16. Stage 25 contract

AI research перемещён в пункт 25 и блокируется до READY Stage 24.

Input:

`Dante_knowledge_findings_open.tsv`

Stage 25 сортирует findings по `priority_rank` и дедуплицирует строки по `event_key`, чтобы один allele с несколькими disease-relations не съедал несколько API slots.

Research-only overlays не отправляются автоматически в commercial-oriented AI pass. Их можно анализировать отдельным research workflow с соответствующим usage manifest.

## 17. Нумерация меню v9.0.0

```text
STAGE 5 — REFERENCE KNOWLEDGE + FULL-GENOME ENRICHMENT
22  Install/update open reference knowledge bundle
23  Install research-only knowledge bundle
24  Enrich full genomic outputs with reference knowledge

STAGE 6 — AI EVIDENCE RESEARCH
25  Research prioritized Stage 24 findings

STAGE 7 — LOCAL GENOME BROWSER
26  Launch local IGV.js browser

AUTOMATION
21  Install and run stages 1 + 2
```

## 18. Осознанно вне Stage 24

Stage 24 не притворяется универсальным caller для специальных генетических задач. Отдельные будущие модули нужны для:

- PRS/GWAS polygenic scores;
- HLA typing;
- CYP/PGx star alleles/haplotypes;
- repeat expansions / STR;
- phasing;
- mitochondrial heteroplasmy.

Отсутствие этих анализов не является ошибкой текущего knowledge-enrichment contract.

## 19. Performance contract

Целевая метрика — wall-clock до валидного `Dante_knowledge_enrichment.state.json`, а не CPU/RAM utilization.

Shared knowledge bundle компилируется один раз для всех пациентов. Patient scan потоковый; полный каталог не загружается в RAM целиком. Remote resources применяются только к заранее выбранным target loci.

Любая будущая performance-оптимизация Stage 24 должна сравниваться A/B на одинаковом полном patient workload до READY state.


## Structural call confidence (rules-v7)

AnnotSV `ACMG_class` оценивает потенциальную клиническую значимость варианта при условии,
что structural call реален. Он не является мерой технической достоверности самого call.

Поэтому Stage 24 сначала нормализует representation:
- использует только `Annotation_mode=full`;
- `split` gene rows не становятся отдельными событиями;
- reciprocal Manta BND mates объединяются в один source cluster;
- DEL/DUP Manta и DRAGEN CNV при reciprocal overlap >= 0.70 объединяются в одно событие
  с `CROSS_CALLER_CONFIRMED` provenance.

Technical call confidence:
- `CROSS_CALLER_CONFIRMED` — независимые Manta + DRAGEN DEL/DUP;
- `SINGLE_CALLER_STRONG` — precise PASS Manta, GQ >= 20 и достаточная PR/SR support;
- `CALLER_SUPPORTED_DRAGEN_CNV` — PASS DRAGEN CNV;
- `CAUTION_IMPRECISE`, `CAUTION_LOW_SUPPORT`, `CAUTION_LARGE_UNCONFIRMED` — технические ограничения;
- GIAB difficult-region overlap добавляет `GIAB_CAUTION`.

Multi-megabase DEL/DUP от одного Manta caller без независимого CNV подтверждения не может
оставаться `K6_HIGH_EVIDENCE_SV_GENE`: biological evidence сохраняется, но finding
понижается до structural candidate и patient priority 5.
