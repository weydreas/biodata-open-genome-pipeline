# Dante WGS v9.0.0-hotfix9 — контракт финальных отчетов

## Цель

Пункты 25 и 26 являются двумя физически раздельными публикационными путями поверх READY Stage 24. Они не являются двумя режимами одной кнопки.

Оба отчета формируются **детерминированно и без генеративного ИИ/API**. Пункты 25/26 не отправляют варианты, геном, VCF/BAM/FASTQ или другие patient artifacts во внешние AI-сервисы.

Stage 24 отвечает за выявление и структурированную аннотацию находок. Пункты 25/26 превращают готовый Stage-24 слой в читаемые PDF с фиксированными пояснениями категорий, фактическими полями из баз знаний, ссылками на базы, provenance IDs и картой pipeline items 1–24.

Hotfix9 добавляет **publication filter**, не изменяя Stage 24: события с доступным `Max AF > 0.70` не подсвечиваются в основном отчете и автоматически переносятся в отдельное справочное PDF-приложение. Исходные Stage-24 TSV/SQLite при этом не изменяются.

Термин `common/major allele` в этом контракте означает только высокую наблюдаемую популяционную частоту. Это не отдельная научная классификация варианта и не абсолютное доказательство функциональной нейтральности.

## Пункт 25 — CLIENT PDF / OPEN / COMMERCIAL-SAFE / NO AI

Пункт 25 читает только:

- `Dante_knowledge_enrichment.state.json`;
- `Dante_knowledge_findings_open.tsv`;
- `Dante_knowledge_usage_manifest.json`.

Он не требует наличия и не открывает:

- `Dante_knowledge_findings_research.tsv`;
- `Dante_knowledge_usage_research_manifest.json`.

Перед публикацией generator требует `scope=OPEN_FINDINGS_ONLY` и блокирует отчет, если в usage manifest присутствует источник без `publication_scope=COMMERCIAL_ALLOWED`.

Все Stage-24 open finding rows сначала группируются по `reason_code + event_key`, после чего publication filter разделяет события:

- `Max AF <= 0.70` или AF неизвестна → основной коммерческий PDF;
- `Max AF > 0.70` → справочное common-major приложение.

Один геномный вариант показывается как одно событие внутри категории, а его несколько gene/disease relations сохраняются внутри карточки события.

Основные разделы:

- K1 exact ClinVar: P/LP, P/LP-common, conflicts, mixed, association/risk/protective/drug-response context;
- K2 predicted loss-of-function в validated disease genes;
- K3 missense в validated disease genes;
- K4 ClinGen Actionability context;
- K5 dosage-sensitive CNV;
- K6 structural gene findings/candidates.

Каждый раздел начинается с фиксированного описания того, что означает категория и чего она **не** доказывает. Карточка конкретного события заполняется значениями Stage 24: генотип/зиготность, VEP consequence/impact, AF, ClinVar significance/conditions/review stars, gene-disease classification, inheritance, actionability/dosage, technical confidence, GIAB flags, evidence sources, Stage-24 IDs и pipeline items. Для доступных источников добавляются кликабельные ссылки на базы.

### Заголовок события

Hotfix9 делает заголовок человекочитаемым. После rsID/координаты выводятся:

1. ген(ы) из Stage 24 в скобках;
2. фактический аллель пациента, восстановленный из `REF`, `ALT` и `genotype`;
3. reason code.

Пример:

```text
38. rs4148632 (ABCC8, LOC110121471) | аллель пациента G/G | K1_EXACT_CLINVAR_CONFLICT
```

Для `0/1` при `REF=A`, `ALT=G` будет показано `A/G`; для `1/1` — `G/G`.

### Выходы пункта 25

- `Dante_results/reports/Dante_genome_report_COMMERCIAL.pdf` — основной отчет;
- `Dante_results/reports/Dante_genome_report_COMMERCIAL_COMMON_MAJOR_REFERENCE.pdf` — справочное приложение для событий с `Max AF > 0.70`;
- `Dante_results/reports/Dante_genome_report_COMMERCIAL.manifest.json`.

Manifest schema остается `dante-genome-report-v2` для backward-compatible READY validator и дополнительно содержит:

- `generation_mode=DETERMINISTIC_STAGE24_ONLY`;
- `ai_used=false`;
- `network_used_by_report_generator=false`;
- `all_stage24_findings_included=true` в package-level смысле: все строки учтены между primary PDF и reference appendix;
- `all_stage24_findings_accounted_for=true`;
- `primary_pdf_filtered=true`;
- `research_layer_opened=false`;
- `common_major_filter.rule=max_af > 0.70`;
- отдельные counts primary/reference events и rows;
- SHA256 Stage-24 inputs, основного PDF и reference appendix.

## Пункт 26 — RESEARCH ONLY PDF / NON-COMMERCIAL / NO AI

Пункт 26 является отдельным действием меню. До открытия research-only findings требуется явное подтверждение `RESEARCH` в интерактивном режиме либо отдельный CLI-флаг `--accept-research-restrictions`.

Он читает:

- Stage-24 open findings — только для локального context того же `event_key`;
- `Dante_knowledge_findings_research.tsv`;
- `Dante_knowledge_usage_research_manifest.json`.

Research overlays группируются по `reason_code + event_key`:

- K7 PanelApp;
- K8 AlphaMissense;
- K8 CADD;
- любые будущие research-only reason codes показываются отдельной категорией с явной маркировкой.

Hotfix9 применяет тот же common-major publication filter к research overlay. Если AF отсутствует в overlay, используется `Max AF` соответствующего open event того же `event_key`.

- `Max AF <= 0.70` или AF неизвестна → основной research PDF;
- `Max AF > 0.70` → отдельное RESEARCH ONLY common-major приложение.

Research PDF не дублирует весь коммерческий реестр. Для каждого research event он показывает необходимый open context того же `event_key`: genes, open reason codes, genotype/zygosity, VEP consequence, AF, ClinVar и disease relations, если они присутствуют в open layer.

Research layer может содержать PanelApp, AlphaMissense и CADD. Их ограничения берутся из фактического Stage-24 usage manifest и выводятся в PDF.

### Выходы пункта 26

- `Dante_results/reports/Dante_genome_report_RESEARCH_ONLY.pdf`;
- `Dante_results/reports/Dante_genome_report_RESEARCH_ONLY_COMMON_MAJOR_REFERENCE.pdf`;
- `Dante_results/reports/Dante_genome_report_RESEARCH_ONLY.manifest.json`.

## Полнота без top-N

Финальные PDF не выбирают первые 20 или иной top-N shortlist.

Hotfix9 меняет только место отображения по популяционной частоте:

```text
Stage 24 full layer
        ↓
reason_code + event_key grouping
        ↓
Max AF <= 0.70 / unknown     Max AF > 0.70
        ↓                         ↓
PRIMARY PDF                COMMON-MAJOR REFERENCE PDF
```

Ни одна Stage-24 строка не удаляется. Manifest сохраняет counts обеих частей и SHA256 обоих PDF.

## Состав PDF

1. юридический статус и границы применения;
2. явная маркировка `DETERMINISTIC / NO AI`;
3. assembly, pipeline/rules/report schema и Stage-24 summary;
4. publication filter и counts primary/reference;
5. counts по reason codes: unique genomic events и evidence rows;
6. фиксированное описание каждой категории находок;
7. детальный реестр соответствующих Stage-24 событий;
8. gene + patient allele в заголовке каждого события;
9. database-derived annotations и ссылки на базы;
10. knowledge-source / license section из usage manifest;
11. таблица воспроизводимости pipeline items 1–24;
12. ограничения текущего аналитического покрытия.

## Что hotfix9 не меняет

- Stage 24 не пересчитывается;
- `dante-knowledge-rules-v8` не меняется;
- исходные open/research TSV и SQLite не фильтруются и не перезаписываются;
- K1/K2/K3/K4/K5/K6/K7/K8 scientific semantics не меняются;
- ClinVar conflicts пока остаются агрегированным `Conflicting_classifications_of_pathogenicity`: для разложения на `P/LP ↔ B/LB`, `P/LP ↔ VUS`, `VUS ↔ B/LB` и другие подтипы нужен отдельный submission-level ClinVar слой.

## Ограничения текущего анализа

В текущий детерминированный contract не входят как завершенные валидированные модули:

- PRS/GWAS;
- полные CYP/PGx haplotypes;
- HLA typing;
- repeat expansions;
- global phasing;
- mtDNA heteroplasmy.

Это явно отражается в финальном отчете, чтобы отсутствие находки не трактовалось как доказанное отсутствие соответствующего риска.

## Пункт 27

Local IGV.js browser остается пунктом 27. Это сервис ручной технической проверки локусов, а не обязательная стадия публикации PDF.

## Юридический статус текста

Встроенные disclaimer/usage тексты являются инженерным fail-closed контрактом pipeline и отражают ограничения используемых data sources. Они не заменяют отдельную юридическую экспертизу договоров, оферты и формулировок, используемых организацией в конкретной юрисдикции.