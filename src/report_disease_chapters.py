#!/usr/bin/env python3
"""Deterministic clinical chapter taxonomy for Stage-24 report publication.

No AI is used. Classification is based only on text already present in Stage-24
ClinVar conditions / disease relations (plus open-context rows for research overlays).
Every event receives exactly one primary chapter for rendering and may keep additional
secondary chapters for auditability. Source condition strings are never modified.
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Iterable

TAXONOMY_SCHEMA = "dante-disease-chapters-v1"

CHAPTERS: list[dict[str, object]] = [
    {"code": "D01_ONCOLOGY", "title": "Онкологические заболевания", "description": "Наследственная предрасположенность к опухолям, опухолевые синдромы и варианты, аннотированные в контексте злокачественных или потенциально опухолевых состояний. Наличие варианта в этом разделе не означает наличие опухоли у пациента."},
    {"code": "D02_CARDIOVASCULAR", "title": "Сердечно-сосудистые заболевания и липидный обмен", "description": "Кардиомиопатии, аритмии, сосудистые заболевания, атеросклеротические фенотипы и наследственные нарушения липопротеинов, имеющие прямое отношение к сердечно-сосудистому риску."},
    {"code": "D03_ENDOCRINE", "title": "Эндокринные заболевания и регуляция глюкозы", "description": "Сахарный диабет, гиперинсулинизм/гипогликемия, заболевания щитовидной, надпочечниковой и других эндокринных систем, а также связанные эндокринные количественные признаки."},
    {"code": "D04_INBORN_METABOLIC", "title": "Врожденные ошибки обмена и митохондриальные заболевания", "description": "Наследственные нарушения аминокислотного, органического, билирубинового, железного и иных метаболических путей, врожденные нарушения гликозилирования, митохондриальные и пероксисомные заболевания."},
    {"code": "D05_NEUROLOGY", "title": "Неврологические и нервно-мышечные заболевания", "description": "Эпилепсии, атаксии, нейропатии, спастические параплегии, заболевания двигательного нейрона, мигрень, каналопатии, миотонии и другие неврологические/нервно-мышечные состояния."},
    {"code": "D06_NEURODEVELOPMENTAL", "title": "Нейроразвитие, психиатрия и когнитивные фенотипы", "description": "Нарушения нейроразвития, интеллектуального развития, аутистического спектра и психиатрические/когнитивные фенотипы. Этот раздел не является психиатрическим диагнозом."},
    {"code": "D07_HEMATOLOGY", "title": "Гематологические заболевания и гемостаз", "description": "Анемии, тромбоцитарные нарушения, наследственные дефекты свертывания, гемолиз, гемоглобинопатии и другие заболевания крови и гемостаза."},
    {"code": "D08_IMMUNE", "title": "Иммунные, аутоиммунные и воспалительные заболевания", "description": "Первичные иммунодефициты, аутоиммунные, аутовоспалительные и другие иммунологические состояния. Органоспецифические воспалительные заболевания по возможности относятся к соответствующему органному разделу."},
    {"code": "D09_RENAL", "title": "Заболевания почек и мочевыводящих путей", "description": "Нефропатии, гломерулярные заболевания, нефротические синдромы, камнеобразование и другие наследственные/ассоциированные заболевания почек и мочевыводящих путей."},
    {"code": "D10_GASTRO", "title": "Желудочно-кишечные, печеночные и панкреатические заболевания", "description": "Заболевания желудочно-кишечного тракта, печени, желчных путей и поджелудочной железы, включая наследственные панкреатиты и воспалительные заболевания кишечника."},
    {"code": "D11_RESPIRATORY", "title": "Заболевания дыхательной системы", "description": "Наследственные и ассоциированные заболевания легких и дыхательных путей, включая обструктивные фенотипы, муковисцидоз и первичную цилиарную дискинезию."},
    {"code": "D12_OPHTHALMOLOGY", "title": "Заболевания органов зрения", "description": "Наследственные заболевания сетчатки, макулы, зрительного нерва, роговицы, хрусталика и другие офтальмологические фенотипы."},
    {"code": "D13_HEARING_ENT", "title": "Нарушения слуха и ЛОР-заболевания", "description": "Наследственная тугоухость/глухота и синдромы с ведущим нарушением слуха, а также другие ЛОР-фенотипы."},
    {"code": "D14_MUSCULOSKELETAL", "title": "Костно-мышечные, соединительнотканные и скелетные заболевания", "description": "Скелетные дисплазии, заболевания соединительной ткани, костей и суставов, врожденные деформации и другие преимущественно костно-мышечные фенотипы."},
    {"code": "D15_DERMATOLOGY", "title": "Кожа, волосы и пигментация", "description": "Наследственные и ассоциированные заболевания кожи, волос и пигментации, включая альбинизм, кератозы и псориатические фенотипы."},
    {"code": "D16_REPRODUCTIVE", "title": "Репродуктивное здоровье и половое развитие", "description": "Фертильность, сперматогенез, крипторхизм, нарушения полового развития и другие наследственные репродуктивные фенотипы."},
    {"code": "D17_CONGENITAL_MULTISYSTEM", "title": "Врожденные пороки и мультисистемные синдромы", "description": "Врожденные аномалии и синдромы, затрагивающие несколько систем и не имеющие более подходящего ведущего клинического раздела."},
    {"code": "D18_INFECTIOUS", "title": "Инфекционная восприимчивость и ответ на инфекцию", "description": "Генетические варианты, описанные в контексте восприимчивости, тяжести или ответа на инфекционные заболевания."},
    {"code": "D19_PHARMACOGENETICS", "title": "Фармакогенетика и ответ на лекарственные препараты", "description": "Drug-response/toxicity фенотипы из ClinVar. Этот раздел фиксирует конкретные variant-level аннотации, но не заменяет полноценный гаплотипный PGx-анализ (например CYP star alleles), если он отдельно не реализован."},
    {"code": "D20_TRAITS", "title": "Неболезненные и количественные признаки", "description": "Количественные признаки, защитные эффекты и другие фенотипы, не являющиеся самостоятельным заболеванием (например группа крови, запаховые или лабораторные количественные признаки)."},
    {"code": "D99_UNCLASSIFIED", "title": "Недостаточно данных для клинической группировки", "description": "Stage 24 не содержит достаточно специфичного названия состояния для надежного детерминированного отнесения к одному из клинических разделов. Исходные строки и ссылки сохраняются без изменений."},
]

CHAPTER_BY_CODE = {str(item["code"]): item for item in CHAPTERS}
CHAPTER_ORDER = [str(item["code"]) for item in CHAPTERS]
CHAPTER_RANK = {code: i for i, code in enumerate(CHAPTER_ORDER)}

_GENERIC_PATTERNS = [
    r"^not provided$", r"^not specified$", r"^unspecified$", r"^see cases$",
    r"^family history$", r"^inborn genetic diseases$", r"^genetic disease$",
    r"^reclassified variant of unknown significance$", r"^variant of unknown significance$",
    r"^[a-z0-9]+ related disorder$", r"^[a-z0-9]+ related condition$",
]

_RULES: list[tuple[str, tuple[str, ...]]] = [
    ("D19_PHARMACOGENETICS", (r"\bdrug response\b", r"\bresponse\b", r"\btoxicity\b", r"\bmetabolizer\b", r"\bclopidogrel\b", r"\birinotecan\b", r"\bfluorouracil\b", r"\bcapecitabine\b", r"\bwarfarin\b", r"\bstatin\b", r"\bthalidomide\b", r"\binhibitor response\b")),
    ("D01_ONCOLOGY", (r"\bcancer\b", r"\bcarcinoma\b", r"\btumou?r\b", r"\bneoplasm\b", r"\bneoplasia\b", r"\bleukemia\b", r"\blymphoma\b", r"\bmelanoma\b", r"\bsarcoma\b", r"\bglioma\b", r"\bosteosarcoma\b", r"\bneuroblastoma\b", r"\bparaganglioma\b", r"\bpheochromocytoma\b", r"\bmalignant\b", r"\blynch syndrome\b", r"\bcancer predisposition\b")),
    ("D02_CARDIOVASCULAR", (r"\bcardiomyopathy\b", r"\bcardiovascular\b", r"\bmyocardial infarction\b", r"\bcoronary\b", r"\batherosclero", r"\barrhythm", r"\blong qt\b", r"\bshort qt\b", r"\bbrugada\b", r"\baortic\b", r"\baorta\b", r"\barteriopathy\b", r"\bvascular\b", r"\bhypertension\b", r"\bhypercholester", r"\bhypobetalipoprotein", r"\bhyperlipid", r"\blipoprotein\b", r"\bldl\b", r"\bheart disease\b", r"\bcardiac\b", r"\bseptal defect\b", r"\bheart defect\b")),
    ("D03_ENDOCRINE", (r"\bdiabetes\b", r"\bdiabetic\b", r"\bhyperinsulin", r"\bhypoglyc", r"\binsulin", r"\bobesity\b", r"\bthyroid\b", r"\bthyroxine\b", r"\bparathyroid\b", r"\bpituitary\b", r"\bgrowth hormone\b", r"\badrenal\b", r"\bcortisol\b", r"\baldosteron", r"\bhyperglyc", r"\bglucose homeostasis\b")),
    ("D04_INBORN_METABOLIC", (r"\bmitochond", r"\bperoxisom", r"\bglycosylation\b", r"\bmaple syrup urine\b", r"\bhyperprolinemia\b", r"\bhemochromatosis\b", r"\bporphyria\b", r"\bcrigler najjar\b", r"\bgilbert syndrome\b", r"\burea cycle\b", r"\borganic acid", r"\bamino acid", r"\bmetabolic disease\b", r"\bmetabolic disorder\b", r"\bbiogenesis disorder\b", r"\b3 hydroxy 3 methylglutaryl coa synthase deficiency\b", r"\btransferrin.*trait\b", r"\bbilirubin\b", r"\bdihydropyrimidine dehydrogenase deficiency\b")),
    ("D07_HEMATOLOGY", (r"\banemia\b", r"\bhaemol", r"\bhemol", r"\bthrombocyt", r"\bthrombocythemia\b", r"\bplatelet\b", r"\bcoagulat", r"\bbleeding\b", r"\bhemophilia\b", r"\bfactor [ivx0-9]+\b", r"\bhemoglobin", r"\bspherocyt", r"\belliptocyt", r"\bpyropoikilocyt", r"\bagranulocyt", r"\bneutropenia\b", r"\bpolycythemia\b", r"\bthrombosis\b", r"\bamegakaryocytic\b", r"\bblood disorder\b")),
    ("D09_RENAL", (r"\bnephro", r"\bkidney\b", r"\bglomerul", r"\brenal(?! cell carcinoma)", r"\bnephrolith", r"\burolith", r"\burinary\b")),
    ("D10_GASTRO", (r"\binflammatory bowel\b", r"\bcrohn\b", r"\bulcerative colitis\b", r"\bpancreatitis\b", r"\bpancreatic\b", r"\bhepatic\b", r"\bliver\b", r"\bcholest", r"\bbile acid\b", r"\bgallbladder\b", r"\bmalabsorption\b", r"\bceliac\b", r"\bgastrointestinal\b", r"\bintestinal\b")),
    ("D11_RESPIRATORY", (r"\bpulmonary\b", r"\blung\b", r"\basthma\b", r"\bobstructive\b", r"\bcystic fibrosis\b", r"\bciliary dyskinesia\b", r"\brespiratory\b", r"\bbronch", r"\binterstitial lung\b")),
    ("D08_IMMUNE", (r"\brheumatoid\b", r"\blupus\b", r"\bautoimmune\b", r"\bautoinflamm", r"\bimmune\b", r"\bimmunodefic", r"\bagammaglobul", r"\bcomplement\b", r"\ballerg", r"\bmannose binding lectin deficiency\b", r"\binflammat", r"\bimmune dysregulation\b")),
    ("D06_NEURODEVELOPMENTAL", (r"\bintellectual dis", r"\bdevelopmental disorder\b", r"\bdevelopmental delay\b", r"\bautis", r"\badhd\b", r"\bschizophren", r"\bbipolar\b", r"\bdepression\b", r"\banxiety\b", r"\bpsychosis\b", r"\bpitt hopkins\b", r"\bkleefstra\b")),
    ("D05_NEUROLOGY", (r"\bepilep", r"\bseizure\b", r"\bataxia\b", r"\bneuropath", r"\bcharcot marie tooth\b", r"\bspastic paraplegia\b", r"\bamyotrophic lateral sclerosis\b", r"\bparkinson\b", r"\balzheimer\b", r"\bdementia\b", r"\bmigraine\b", r"\bmyotonia\b", r"\bmyasthen", r"\bparalysis\b", r"\bdystonia\b", r"\bencephal", r"\bcerebell", r"\bcerebral\b", r"\bneurolog", r"\bmuscular dystrophy\b", r"\bmyopathy\b", r"\bazorean disease\b")),
    ("D12_OPHTHALMOLOGY", (r"\bretina", r"\bretinitis\b", r"\bmacular\b", r"\boptic\b", r"\bcornea\b", r"\bcorneal\b", r"\bcataract\b", r"\bglaucoma\b", r"\bocular\b", r"\beye\b", r"\bnystagmus\b")),
    ("D13_HEARING_ENT", (r"\bhearing loss\b", r"\bdeaf", r"\busher syndrome\b", r"\bauditory\b", r"\botologic\b")),
    ("D14_MUSCULOSKELETAL", (r"\bdysplasia\b", r"\bskeletal\b", r"\bosteogenesis\b", r"\behlers danlos\b", r"\bconnective tissue\b", r"\bbone\b", r"\bscoliosis\b", r"\bsynost", r"\bchondrodysplasia\b", r"\bosteodysplastic\b", r"\bspondylo", r"\bjoint\b", r"\bcartilage\b", r"\bshort stature\b")),
    ("D15_DERMATOLOGY", (r"\bskin\b", r"\bhair\b", r"\bpigment", r"\balbinism\b", r"\bpsoriasis\b", r"\bpityriasis\b", r"\bkeratosis\b", r"\bpiebaldism\b", r"\bdermat")),
    ("D16_REPRODUCTIVE", (r"\bcryptorchid", r"\binfertil", r"\bfertility\b", r"\bspermat", r"\btesticular\b", r"\bovarian(?! cancer)", r"\bsex development\b", r"\bgonadal\b", r"\breproductive\b")),
    ("D18_INFECTIOUS", (r"\binfectious\b", r"\binfection\b", r"\bviral\b", r"\bbacterial\b", r"\bhiv\b", r"\bmalaria\b", r"\btuberculosis\b", r"\bcovid\b", r"\bhepatitis [abc]\b")),
    ("D17_CONGENITAL_MULTISYSTEM", (r"\bfetal anomal", r"\bcongenital anomal", r"\bmultiple congenital\b", r"\bnoonan\b", r"\b3mc syndrome\b", r"\bfrank ter haar\b", r"\bsaul wilson\b", r"\bprimordial dwarfism\b", r"\bsyndrome\b")),
    ("D20_TRAITS", (r"\bquantitative trait\b", r"\bpolymorphism\b", r"\bability to smell\b", r"\bblood group\b", r"\bprotection against\b", r"\bprotective\b", r"\bsusceptibility to\b", r"\bfamily history\b", r"\bserum level\b", r"\bphenotype\b")),
]

_COMPILED_GENERIC = [re.compile(p, re.I) for p in _GENERIC_PATTERNS]
_COMPILED_RULES = [(code, [re.compile(p, re.I) for p in patterns]) for code, patterns in _RULES]


def normalize_condition(text: str) -> str:
    text = str(text or "").strip().replace("_", " ")
    text = re.sub(r"([a-z])([A-Z])", r"\1 \2", text)
    text = text.lower()
    text = text.replace("/", " ").replace("-", " ").replace(":", " ").replace(",", " ")
    text = re.sub(r"[^a-z0-9+ ]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def split_conditions(value: object) -> list[str]:
    if value is None:
        return []
    result: list[str] = []
    seen: set[str] = set()
    for token in re.split(r"[|;]+", str(value)):
        raw = token.strip()
        if not raw:
            continue
        key = normalize_condition(raw)
        if not key or key in seen:
            continue
        seen.add(key)
        result.append(raw)
    return result


def is_generic_condition(text: str) -> bool:
    normalized = normalize_condition(text)
    if not normalized:
        return True
    return any(pattern.search(normalized) for pattern in _COMPILED_GENERIC)


def classify_condition(text: str) -> str:
    normalized = normalize_condition(text)
    if not normalized or is_generic_condition(text):
        return "D99_UNCLASSIFIED"
    for code, patterns in _COMPILED_RULES:
        if any(pattern.search(normalized) for pattern in patterns):
            return code
    return "D99_UNCLASSIFIED"


def informative_conditions(values: Iterable[object]) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        for raw in split_conditions(value):
            key = normalize_condition(raw)
            if not key or key in seen or is_generic_condition(raw):
                continue
            seen.add(key)
            result.append(raw.replace("_", " "))
    return result


def assign_conditions(values: Iterable[object]) -> dict[str, object]:
    conditions = informative_conditions(values)
    counts: Counter[str] = Counter()
    condition_map: list[dict[str, str]] = []
    for condition in conditions:
        chapter = classify_condition(condition)
        counts[chapter] += 1
        condition_map.append({"condition": condition, "chapter": chapter})

    non_unknown = [(code, n) for code, n in counts.items() if code != "D99_UNCLASSIFIED" and n > 0]
    if non_unknown:
        primary = min(non_unknown, key=lambda item: (-item[1], CHAPTER_RANK[item[0]]))[0]
    else:
        primary = "D99_UNCLASSIFIED"

    secondary = [code for code, n in sorted(counts.items(), key=lambda item: (-item[1], CHAPTER_RANK[item[0]])) if n > 0 and code not in {primary, "D99_UNCLASSIFIED"}]
    return {"schema": TAXONOMY_SCHEMA, "primary": primary, "secondary": secondary, "conditions": conditions, "condition_map": condition_map, "chapter_counts": dict(counts)}


def chapter_title(code: str) -> str:
    return str(CHAPTER_BY_CODE.get(code, CHAPTER_BY_CODE["D99_UNCLASSIFIED"])["title"])


def chapter_description(code: str) -> str:
    return str(CHAPTER_BY_CODE.get(code, CHAPTER_BY_CODE["D99_UNCLASSIFIED"])["description"])


def self_test() -> int:
    cases = [
        ("Breast_and_colorectal_cancer,_susceptibility_to", "D01_ONCOLOGY"),
        ("Familial_hypercholesterolemia", "D02_CARDIOVASCULAR"),
        ("Maturity-onset_diabetes_of_the_young", "D03_ENDOCRINE"),
        ("Maple_syrup_urine_disease", "D04_INBORN_METABOLIC"),
        ("Mitochondrial_complex_I_deficiency,_nuclear_type_1", "D04_INBORN_METABOLIC"),
        ("Amyotrophic_lateral_sclerosis", "D05_NEUROLOGY"),
        ("Intellectual_developmental_disorder,_X-linked_111", "D06_NEURODEVELOPMENTAL"),
        ("Hereditary_spherocytosis_type_3", "D07_HEMATOLOGY"),
        ("Rheumatoid_arthritis", "D08_IMMUNE"),
        ("Nephrolithiasis,_calcium_oxalate", "D09_RENAL"),
        ("Hereditary_pancreatitis", "D10_GASTRO"),
        ("Pulmonary_disease,_chronic_obstructive,_susceptibility_to", "D11_RESPIRATORY"),
        ("Retinitis_pigmentosa_58", "D12_OPHTHALMOLOGY"),
        ("Autosomal_recessive_nonsyndromic_hearing_loss_59", "D13_HEARING_ENT"),
        ("Diastrophic_dysplasia", "D14_MUSCULOSKELETAL"),
        ("Oculocutaneous_albinism_type_1A", "D15_DERMATOLOGY"),
        ("Bilateral_cryptorchidism", "D16_REPRODUCTIVE"),
        ("Noonan_syndrome_and_Noonan-related_syndrome", "D17_CONGENITAL_MULTISYSTEM"),
        ("COVID-19_susceptibility", "D18_INFECTIOUS"),
        ("Clopidogrel_response", "D19_PHARMACOGENETICS"),
        ("C3HEX,_ability_to_smell", "D20_TRAITS"),
        ("not_provided", "D99_UNCLASSIFIED"),
    ]
    for condition, expected in cases:
        actual = classify_condition(condition)
        assert actual == expected, (condition, expected, actual)
    mixed = assign_conditions(["Maturity-onset_diabetes_of_the_young|Type_2_diabetes_mellitus|Nonpapillary_renal_cell_carcinoma"])
    assert mixed["primary"] == "D03_ENDOCRINE", mixed
    assert "D01_ONCOLOGY" in mixed["secondary"], mixed
    ugt = assign_conditions(["Gilbert_syndrome|Crigler-Najjar_syndrome|Irinotecan_response"])
    assert ugt["primary"] == "D04_INBORN_METABOLIC", ugt
    assert "D19_PHARMACOGENETICS" in ugt["secondary"], ugt
    print("SELF-TEST PASS | deterministic disease chapter taxonomy")
    return 0


if __name__ == "__main__":
    raise SystemExit(self_test())
