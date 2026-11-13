# VCF -> 23andMe v5-compatible TXT — контракт v9.1.0

## Назначение

Начиная с Biodata Genome Pipeline `9.1.0` одиночный single-sample VCF можно преобразовать в четырехколоночный TXT:

```text
# rsid    chromosome    position    genotype
```

Формат совместим по структуре с raw-data TXT 23andMe v5. Это **не означает**, что результат искусственно урезается до набора SNP, реально находившихся на чипе 23andMe v5.

Цель обратная: сохранить максимально полный набор значимых SNP, присутствующих в исходном VCF, не теряя строки только потому, что в поле VCF `ID` отсутствует rsID.

## Главный invariant

Для каждого валидного CALLED SNP на канонических хромосомах `1..22/X/Y/MT`, который можно представить одной или двумя буквами `A/C/G/T`:

```text
called canonical SNV record -> exported record
```

В отчете обязательно выполняется:

```text
called_canonical_snv_records == exported_called_snv_records
meaningful_called_snv_loss_records == 0
loss_accounting_invariant == PASS
```

Если rsID отсутствует, SNP не удаляется. Он получает детерминированный marker ID вида:

```text
DANTE_GRCh37_1_123456_A_G
```

и одновременно фиксируется в audit TSV.

## Выходные файлы

Для generic single-sample VCF создаются:

```text
patients/<Patient>/Dante_results/Dante_23andMe_v5.txt
patients/<Patient>/Dante_results/Dante_23andMe_v5.state.json
patients/<Patient>/Dante_qc/Dante_23andMe_v5_unresolved.tsv
patients/<Patient>/Dante_qc/Dante_23andMe_v5_report.txt
```

`Dante_23andMe_v5_unresolved.tsv` — не только список ошибок. Он является audit ledger и включает также успешно экспортированные SNP без rsID, для которых использован `DANTE_*` marker.

## Поддерживаемые генотипы

Поддерживаются:

- `0/0` — если он реально присутствует в исходном VCF;
- `0/1`, `1/0`;
- `1/1`;
- haploid `0`/`1` для X/Y/MT и иных корректных гаплоидных calls;
- multiallelic SNP, например `REF=A ALT=C,G GT=1/2` -> `CG`.

Фазирование в четырехколоночном формате не представимо. Для `0|1` сохраняются сами аллели, но разделитель фазы исчезает; такие записи отдельно считаются в отчете.

## Что не попадает в SNP TXT

Не смешиваются с SNP-экспортом:

- INDEL;
- symbolic ALT;
- CNV/SV;
- некорректные/неизвестные хромосомы;
- missing genotype;
- человечески нетипичные ploidy > 2, которые нельзя корректно представить колонкой genotype 23andMe.

Каждый такой случай фиксируется в audit TSV с `reason` и не исчезает молча.

## rsID resolution

При обычном запуске конвертер использует публичный NCBI Variation Services endpoint:

```text
/variation/v0/vcf/file/set_rsids
```

Запросы отправляются batch-ами. Это соответствует публичному dbSNP/NCBI workflow для VCF -> RefSNP annotation.

Сбой сети не делает conversion fail-open по потере данных: для записи используется `DANTE_*` marker, поэтому genotype остается в выходном TXT.

Для полностью offline запуска:

```powershell
$env:DANTE_VCF_23ME_NO_RSID_RESOLUTION="1"
```

После этого внешние запросы не выполняются, а все отсутствующие rsID заменяются детерминированными локальными marker IDs.

## Assembly policy

### GRCh37

`GRCh37/hg19` является канонически совместимым по координатам с raw-data convention 23andMe v5.

### GRCh38

`GRCh38/hg38` также может быть экспортирован в ту же четырехколоночную структуру, но координаты остаются GRCh38.

Pipeline **не выполняет фиктивный или эвристический liftover** и не меняет strand по догадке.

Поэтому GRCh38-файл маркируется в header как:

```text
reference_build=GRCh38
coordinates=GRCh38; column-compatible only
```

Для downstream-инструмента, требующего именно GRCh37 coordinates, такой файл нельзя выдавать за оригинальный 23andMe v5 dataset.

## FILTER policy

Конвертер не вводит новый порог качества, которого не было в задаче. Called SNP с `FILTER != PASS` сохраняется и отдельно считается в QC report.

Если конкретному downstream workflow нужны только `PASS`, фильтрация должна быть отдельным, явно сформулированным analytical contract, а не скрытой частью преобразования формата.

## UI

Для generic single-sample VCF показывается минимальный профиль:

```text
VCF CONVERSION
 2. Convert VCF -> 23andMe v5-compatible SNP TXT
 0. Exit
```

Docker для этой операции не требуется. Нужен только Python 3; интернет нужен только для дополнительного rsID enrichment.

## Что намеренно не изменено

Release `9.1.0` не меняет:

- существующий Dante WGS analytical contract;
- существующий Genotek/genotyping-array workflow;
- Stage 24 finding rules;
- коммерческий/research split;
- правила CNV/SV;
- интерпретацию отсутствующего варианта как negative evidence.

Generic VCF converter является отдельным runtime adapter с fallback на прежние реализации.
