# Clinical Risk Evidence — контракт публикационного слоя

## Назначение

`Clinical Risk Evidence` — отдельный детерминированный слой между Stage 24 и финальными PDF 25/26.

Он **не меняет** правила Stage 24, reason codes, pathogenicity, disease-chapter assignment или набор найденных вариантов. Его задача — дополнить уже найденные Stage-24 события полями, необходимыми для кратких клинических таблиц:

- аллель и читаемый генотип пациента;
- нормализованная ClinVar-классификация;
- количественный риск / пенетрантность, только если локальный источник содержит явную численную оценку;
- provenance и явное указание уровня оценки (`variant-specific` или `gene-condition`).

## Источники v1

Первая версия использует только локальные open/commercial-safe данные:

1. `Dante_results/knowledge_enrichment/Dante_knowledge_findings_open.tsv`;
2. `Dante_knowledge_enrichment.state.json`;
3. локальный `reference_data/open_reference/open_reference_manifest.json`;
4. уже установленные пунктом 22 ClinGen Actionability Adult/Pediatric JSON snapshots.

При построении слоя сетевые запросы не выполняются. rsID, гены, генотипы и иные пациентские геномные идентификаторы не передаются во внешние API.

## Допустимые количественные метрики

Из источника принимается только явная количественная формулировка:

- penetrance;
- lifetime / cumulative / absolute risk;
- probability;
- incidence / prevalence / frequency;
- RR;
- HR;
- OR.

Допускаются явные проценты, диапазоны процентов, `1 in N` и явно подписанные RR/HR/OR.

Числа без семантики риска не преобразуются в риск. В частности, semi-quantitative ClinGen likelihood/actionability scores не являются процентом пенетрантности и не пересчитываются.

## Сопоставление disease context

Совпадение только по гену недостаточно для переноса количественной оценки. Gene-condition estimate используется только при согласованном disease/condition context между Stage 24 и ClinGen snapshot.

Если конкретная запись источника содержит rsID пациента, оценка может быть помечена `variant-specific`. В остальных случаях она должна оставаться `gene-condition` и не представляться как пенетрантность конкретного варианта.

Если подходящей численной оценки нет, результатом является явное:

`Количественная оценка риска не установлена по локальному ClinGen Actionability snapshot`.

Отсутствие числа не является ошибкой и не должно заменяться генеративной оценкой.

## Артефакты

Слой записывается в:

`Dante_results/clinical_risk_evidence/`

- `Dante_clinical_risk_evidence.json` — event-oriented представление для PDF;
- `Dante_clinical_risk_evidence.tsv` — плоский аудит;
- `Dante_clinical_risk_evidence.state.json` — identities, hashes, source policy и counts.

READY зависит от хэшей Stage 24, open-reference manifest и обоих ClinGen Actionability snapshots. Изменение любого входа делает слой STALE, но не изменяет и не инвалидирует сам Stage 24.

## Интеграция с отчетами

Пункты 25/26 автоматически:

1. проверяют READY Stage 24;
2. проверяют Clinical Risk Evidence;
3. если слой отсутствует/STALE — пересобирают только его;
4. создают новую versioned PDF package без перезаписи предыдущих отчетов.

PDF renderer сам не выполняет сетевые запросы и не использует AI/API.

Краткая таблица раздела содержит:

1. `rsID / вариант`;
2. `Аллель пациента`;
3. `Ген(ы)`;
4. `Классификация ClinVar`;
5. `Риск / пенетрантность`;
6. `Связанные состояния`.

Та же информация добавляется в детальную карточку события.

## Интерпретационная граница

Количественный показатель — это контекст опубликованного источника, а не индивидуальный прогноз, диагноз или назначение. Gene-condition penetrance не должна автоматически интерпретироваться как variant-specific penetrance. Клинический смысл зависит от возраста, пола, наследования, семейного анамнеза, фенотипа, качества исходного генотипа и других факторов.
