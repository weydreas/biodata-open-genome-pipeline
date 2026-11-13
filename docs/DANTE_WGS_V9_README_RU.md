# Dante WGS local pipeline v9.0.0

Текущий major release добавляет Stage 24 — детерминированное применение установленных open + research knowledge bundles ко всему каноническому геномному результату пациента.

История предыдущих релизов сохранена в `CHANGELOG.md` и `DANTE_WGS_V8_README_RU.md`. Полный технический контракт Stage 24: `STAGE24_KNOWLEDGE_ENRICHMENT_CONTRACT_RU.md`. Контракт финальных PDF: `FINAL_REPORTS_CONTRACT_RU.md`.

## Главная архитектура v9.0.0-hotfix8

```text
Stage 1  CNV/SV AnnotSV
Stage 2  full deterministic small-variant catalog
Stage 3  optional array + FASTQ reference-genome completion
Stage 4  telomere estimate

Stage 5
  22 open/commercial-safe reference knowledge
  23 research-only reference knowledge
  24 full-genome knowledge enrichment (NO OPENAI)
  25 CLIENT PDF: complete open/commercial-safe Stage-24 register (NO AI)
  26 RESEARCH ONLY PDF: complete research overlay register (NO AI)

Stage 6
  27 local IGV.js genome browser

Automation
  21 stages 1+2
```

## Что изменилось по сравнению с v8.6.0

### Новый пункт 24

Пункт 24 сканирует:

- `Dante_small_variants_annotated.tsv.gz` — весь unified catalog;
- `Dante_CNV_AnnotSV.tsv`;
- `Dante_SV_AnnotSV.tsv`;
- при наличии — Stage 17 Tier A/B FASTQ ALT evidence.

Stage 24 не использует OpenAI. Open/commercial-safe слой объединяет genomic events с HGNC, GenCC, ClinGen, ClinVar, GIAB, Mondo/HPO ontology vocabulary и target-only dbSNP. AlphaMissense, PanelApp и CADD физически публикуются только в research-only overlay.

Open и research-only knowledge физически разделены. Research-only predictors/panel membership никогда не записываются в commercial-safe findings TSV.

### Пункт 25 — коммерческий PDF

Пункт 25 требует READY Stage 24 и читает только:

- `Dante_knowledge_enrichment.state.json`;
- `Dante_knowledge_findings_open.tsv`;
- `Dante_knowledge_usage_manifest.json`.

Он не открывает research-only TSV/manifest. В PDF включаются **все** open finding rows, сгруппированные по `reason_code` и `event_key`. Один variant event показывается один раз внутри категории, а его несколько gene/disease relations перечисляются внутри карточки.

Никаких OpenAI/API вызовов, top-N shortlist или `OPENAI_MAX_VARIANTS` в пункте 25 нет.

### Пункт 26 — research-only PDF

Пункт 26 является физически отдельной кнопкой и требует явного подтверждения `RESEARCH` перед открытием research-only layer.

В PDF включаются **все** research overlays K7/K8, сгруппированные по `reason_code` и `event_key`. Для каждого research event добавляется только локальный open context того же `event_key`; полный open register остается в пункте 25.

Никаких OpenAI/API вызовов в пункте 26 нет.

### Пункт 27

Local IGV.js browser перенумерован в пункт 27. Его геномный контракт не изменен.

## Knowledge bundles

### 22 — open/commercial-safe

Open bundle содержит HGNC, GenCC, ClinGen, GIAB GRCh37 stratifications, Mondo, HPO ontology-only и dbSNP dual-build remote indexes/assembly reports. AlphaMissense из-за CC BY-NC-SA 4.0 находится в research-only bundle вместе с PanelApp/CADD. Полные dbSNP VCF намеренно не скачиваются локально.

HPO в текущем bundle — ontology vocabulary. Disease/gene annotation products, содержащие внешние лицензируемые mappings, автоматически не добавляются.

### 23 — research only

Автоматически устанавливаются:

- PanelApp snapshot;
- AlphaMissense;
- CADD v1.7 GRCh37 index/remote-score metadata.

UCSC chain files намеренно **не устанавливаются** автоматически; target build reconciliation использует open dbSNP dual-build layer. OMIM остается gated/not-installed.

## Требования к пункту 24

Обязательны READY:

- Stage 1 AnnotSV outputs;
- Stage 2 unified catalog;
- item 22 open bundle;
- item 23 research bundle.

Stage 17 опционален. Если он рассчитан после Stage 24, knowledge enrichment становится stale и должен быть повторен.

## Результаты Stage 24

Каталог:

`patients/<Patient>/Dante_results/knowledge_enrichment/`

Основные outputs:

| Файл | Назначение |
| --- | --- |
| `Dante_knowledge_findings_open.tsv` | open/commercial-safe deterministic findings |
| `Dante_knowledge_findings_research.tsv` | research-only PanelApp/AlphaMissense/CADD overlays |
| `Dante_knowledge_qc.tsv` | GIAB difficult-region/breakpoint flags |
| `Dante_knowledge.sqlite3` | patient knowledge DB только для релевантных events |
| `Dante_knowledge_enrichment_report.txt` | scan counts, reason counts, remote lookup status, boundaries |
| `Dante_knowledge_usage_manifest.json` | реально использованные open sources для commercial publication guard |
| `Dante_knowledge_usage_research_manifest.json` | open + research source usage для research mode |
| `Dante_knowledge_enrichment.state.json` | SHA256/signature/READY contract |

Shared index:

`reference_data/knowledge_index/knowledge_index.sqlite3`

Он строится один раз на bundle snapshots и переиспользуется между пациентами.

## Finding rules

Open reason codes:

- `K1_EXACT_CLINVAR_PLP`;
- `K1_EXACT_CLINVAR_PLP_COMMON`;
- `K1_EXACT_CLINVAR_MIXED`;
- `K1_EXACT_CLINVAR_CONFLICT`;
- `K1_EXACT_CLINVAR_ASSOCIATION`;
- `K2_RARE_LOF_VALIDATED_GENE`;
- `K2_LOF_VALIDATED_GENE_AF_UNKNOWN`;
- `K3_RARE_MISSENSE_VALIDATED_GENE`;
- `K3_MISSENSE_VALIDATED_GENE_AF_UNKNOWN`;
- `K4_CLINGEN_ACTIONABLE`;
- `K4_CLINGEN_ACTIONABLE_AF_UNKNOWN`;
- `K5_DOSAGE_SENSITIVE_CNV`;
- `K6_HIGH_EVIDENCE_SV_GENE`;
- `K6_STRUCTURAL_GENE_CANDIDATE`.

Research-only:

- `K7_PANELAPP_RESEARCH`;
- `K8_ALPHAMISSENSE_SUPPORT`;
- `K8_CADD_SUPPORT_TARGET`.

AlphaMissense используется только как research-only supporting predictor и никогда не попадает в commercial-safe/open TSV.

Для AR двух heterozygous candidate alleles Stage 24 может поставить `possible_biallelic_unphased=yes`, но не называет их compound heterozygous без phasing.

GIAB overlap меняет `technical_confidence`, но не удаляет finding. Для structural events применяются breakpoint-safe проверки.

## Remote target lookups Stage 24

После полного локального scan remote/indexed базы вызываются только для candidate targets, а не для миллионов строк:

- AlphaMissense — research-only relevant missense targets;
- dbSNP — candidate alleles без resolved rsID;
- CADD — research targets.

По умолчанию сбой optional remote lookup приводит к `READY / REMOTE WARNINGS`, сохраняя результат полного локального scan. Для fail-closed режима в `.env`:

```dotenv
DANTE_KNOWLEDGE_STRICT_REMOTE=1
```

## Финальные PDF

Каталог:

`patients/<Patient>/Dante_results/reports/`

Выходы:

- `Dante_genome_report_COMMERCIAL.pdf`;
- `Dante_genome_report_COMMERCIAL.manifest.json`;
- `Dante_genome_report_RESEARCH_ONLY.pdf`;
- `Dante_genome_report_RESEARCH_ONLY.manifest.json`.

Manifest schema hotfix8: `dante-genome-report-v2`.

Оба manifest фиксируют:

- `generation_mode=DETERMINISTIC_STAGE24_ONLY`;
- `ai_used=false`;
- `network_used_by_report_generator=false`;
- `all_stage24_findings_included=true`;
- category counts;
- source SHA256;
- output SHA256.

Каждый PDF содержит:

1. юридический статус;
2. Stage-24 summary;
3. counts по категориям;
4. фиксированное описание смысла каждой категории;
5. полный реестр соответствующих событий и evidence rows;
6. database-derived поля и кликабельные ссылки на базы;
7. knowledge usage/licensing;
8. карту pipeline items 1–24;
9. ограничения текущего analytical coverage.

## `.env`

Для пунктов 25/26 OpenAI API key **не нужен**.

```dotenv
# Stage 24: optional fail-closed indexed target lookups
DANTE_KNOWLEDGE_STRICT_REMOTE=0
```

Legacy/explicit AI-инструменты, если они запускаются отдельно от финальных PDF, могут иметь собственные требования, но клиентские и research PDF от них не зависят.

## Запуск

Использовать один launcher:

```text
start.cmd
```

Если пункты 22 и 23 уже имеют `[READY]`, следующий аналитический шаг — пункт `24`.

После успешного 24:

- `25` — полный коммерческий/open PDF;
- `26` — отдельный research-only PDF;
- `27` — локальный browser.

## GitHub-first обновления кода

Начиная с hotfix 5 изменяемая часть pipeline хранится в приватном репозитории
`weydreas/biodata-open-genome-pipeline`. `start.cmd` запускает `src/bootstrap_update.ps1`,
который использует скрытый локальный cache `.pipeline_repo`, делает `git fetch/reset`
до `origin/main`, копирует только пути из `release/pipeline_manifest.json` и не касается
`patients/`, больших reference datasets или результатов.

После синхронизации кода bootstrap вызывает `status` для open/research reference managers.
Если схема bundle устарела или отсутствует новый обязательный artifact, manager запускается
в resumable install mode. Поэтому изменение состава научных баз в Git автоматически
превращается в локальную догрузку на следующем старте.

Отключение при необходимости:
- `DANTE_SKIP_GIT_UPDATE=1` — не проверять Git;
- `DANTE_SKIP_REFERENCE_CHECK=1` — не проверять reference bundles;
- `DANTE_SKIP_REFERENCE_UPDATE=1` — проверить, но не устанавливать недостающее.
