# =============================================================================
# 02_codebooks.R — Codes ICD et item_id MIMIC-IV pour les covariables
# =============================================================================

# -----------------------------------------------------------------------------
# Codes ICD-9 et ICD-10 pour les comorbidités
# -----------------------------------------------------------------------------
get_icd_codes <- function() {
  list(

    # Insuffisance rénale chronique
    ckd = list(
      icd9  = c("585", "5851", "5852", "5853", "5854", "5855", "5856",
                "5859", "V420", "V451", "V4511", "V4512"),
      icd10 = c("N18", "N180", "N181", "N182", "N183", "N184", "N185",
                "N186", "N189", "Z940", "Z992", "Z9115")
    ),

    # ESRD / Dialyse chronique
    esrd = list(
      icd9  = c("5856", "V451", "V4511", "V4512", "V560", "V561",
                "V562", "V5631", "V5632"),
      icd10 = c("N186", "Z992", "Z9115", "Z9113", "Z9114")
    ),

    # Diabète (avec ou sans complications)
    diabetes = list(
      icd9  = c("250", "2500", "2501", "2502", "2503", "2504", "2505",
                "2506", "2507", "2508", "2509"),
      icd10 = c("E10", "E11", "E12", "E13", "E14")
    ),

    # Insuffisance cardiaque
    heart_failure = list(
      icd9  = c("428", "4280", "4281", "4282", "4289", "40201", "40211",
                "40291", "40401", "40411", "40491"),
      icd10 = c("I50", "I500", "I501", "I502", "I503", "I504", "I509")
    ),

    # Maladie hépatique (chronique / cirrhose)
    liver_disease = list(
      icd9  = c("571", "5710", "5711", "5712", "5713", "5714", "5715",
                "5716", "5718", "5719", "07022", "07023", "07032",
                "07033", "07044", "07054"),
      icd10 = c("K70", "K71", "K72", "K73", "K74", "K760", "K766",
                "K767", "K768", "K769", "B15", "B16", "B17", "B18", "B19")
    ),

    # Cancer
    cancer = list(
      icd9  = c("140", "141", "142", "143", "144", "145", "146", "147",
                "148", "149", "150", "151", "152", "153", "154", "155",
                "156", "157", "158", "159", "160", "161", "162", "163",
                "164", "165", "170", "171", "172", "173", "174", "175",
                "176", "179", "180", "181", "182", "183", "184", "185",
                "186", "187", "188", "189", "190", "191", "192", "193",
                "194", "195", "196", "197", "198", "199", "200", "201",
                "202", "203", "204", "205", "206", "207", "208"),
      icd10 = c("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8",
                "C9")
    ),

    # Immunosuppression / transplantation
    immunosuppression = list(
      icd9  = c("V420", "V421", "V422", "V424", "V425", "V426", "V427",
                "V428", "V429", "V4281", "V4282", "V4283", "V4284",
                "V4286", "V4287", "V4289",
                "279", "2790", "2791", "2792", "2793", "2794", "2798", "2799"),
      icd10 = c("Z94", "Z940", "Z941", "Z942", "Z943", "Z944", "Z945",
                "Z946", "Z947", "Z948", "Z949", "D80", "D81", "D82",
                "D83", "D84", "D85", "D86", "D87", "D88", "D89")
    ),

    # Sepsis / Choc septique
    sepsis = list(
      icd9  = c("99591", "99592", "0031", "0362", "03811", "03812",
                "03819", "0389", "7907"),
      icd10 = c("A41", "A410", "A411", "A412", "A413", "A414", "A415",
                "A418", "A419", "R6520", "R6521", "R6521", "R578")
    ),

    # Ventilation mécanique (procédure ICD)
    mech_vent_proc = list(
      icd9  = c("9670", "9671", "9672"),
      icd10 = c("5A1935Z", "5A1945Z", "5A1955Z")
    ),

    # RRT / Dialyse en aigu (procédure)
    rrt_proc = list(
      icd9  = c("3995", "5498"),
      icd10 = c("5A1D70Z", "5A1D80Z", "5A1D90Z", "3E1M39Z")
    )
  )
}

# -----------------------------------------------------------------------------
# Item IDs MIMIC-IV pour chartevents / labevents
# -----------------------------------------------------------------------------
get_item_ids <- function() {
  list(
    # Créatinine sérique (labevents d_labitems)
    creatinine_lab_label   = "creatinine",

    # Chartevents ICU - vasopresseurs (d_items)
    vasopressor_itemids = c(
      221662L, # Dopamine
      221653L, # Dobutamine
      221289L, # Epinephrine
      221906L, # Norepinephrine
      222315L, # Vasopressin
      229617L  # Phenylephrine
    ),

    # Chartevents ICU - ventilation mécanique
    vent_itemids = c(
      223849L, # Ventilator mode
      223848L, # Ventilator mode (simplifié)
      224688L, # Respiratory Rate (Set)
      224687L, # Minute Volume
      220339L, # Peep set
      224696L  # Plateau pressure
    ),

    # Chartevents - fréquence cardiaque
    hr_itemids = c(220045L),

    # Chartevents - pression artérielle moyenne
    map_itemids = c(220052L, 220181L, 225312L),

    # Chartevents - température
    temp_itemids = c(223761L, 223762L),

    # Chartevents - SpO2
    spo2_itemids = c(220277L),

    # Chartevents - GCS total
    gcs_itemids = c(220739L, 223900L, 223901L, 227011L, 227012L),

    # Chartevents - FiO2
    fio2_itemids = c(223835L, 3420L, 3422L),

    # Inputevents - Vasopresseurs item IDs (pour détection ICU)
    vasopressor_input_itemids = c(
      221662L, 221653L, 221289L, 221906L, 222315L, 229617L
    ),

    # Diurèse outputevents
    urine_output_itemids = c(226559L, 226560L, 226561L, 226584L,
                              226563L, 226564L, 226565L, 226567L,
                              226557L, 226558L)
  )
}

# -----------------------------------------------------------------------------
# Codes Charlson par catégorie (ICD-9 et ICD-10) — version simplifiée
# -----------------------------------------------------------------------------
get_charlson_codes <- function() {
  # Utiliser les codes ICD pour calculer le score de Charlson
  # Basé sur : Quan et al. Med Care 2005
  list(
    myocardial_infarction = list(
      icd9  = c("410","412"),
      icd10 = c("I21","I22","I252")
    ),
    congestive_heart_failure = list(
      icd9  = c("39891","40201","40211","40291","40401","40411","40491","4280","4281","4289"),
      icd10 = c("I099","I110","I130","I132","I255","I420","I425","I426","I427","I428",
                "I429","I43","I50","P290")
    ),
    peripheral_vascular = list(
      icd9  = c("0930","4373","440","4410","4411","4412","4413","4414","4415","4416",
                "4419","4431","4432","4438","4439","4471","5571","5579","V434"),
      icd10 = c("I70","I71","I731","I738","I739","I771","I790","I792","K551","K558",
                "K559","Z958","Z959")
    ),
    cerebrovascular = list(
      icd9  = c("36234","430","431","4320","4321","4329","433","434","435","436","437","438"),
      icd10 = c("G45","G46","H340","H341","H342","I60","I61","I62","I63","I64","I65",
                "I66","I67","I68","I69")
    ),
    dementia = list(
      icd9  = c("290","2941","3312"),
      icd10 = c("F00","F01","F02","F03","F051","G30","G311")
    ),
    chronic_pulmonary = list(
      icd9  = c("4168","4169","490","491","492","493","494","495","496","500","501",
                "502","503","504","505","5064","5081","5088"),
      icd10 = c("I278","I279","J40","J41","J42","J43","J44","J45","J46","J47","J60",
                "J61","J62","J63","J64","J65","J66","J67","J684","J701","J703")
    ),
    rheumatic_disease = list(
      icd9  = c("4465","7100","7101","7102","7103","7104","7108","7109","7140","7141",
                "7142","7148","725"),
      icd10 = c("M05","M06","M315","M32","M33","M34","M351","M353","M360")
    ),
    peptic_ulcer = list(
      icd9  = c("531","532","533","534"),
      icd10 = c("K25","K26","K27","K28")
    ),
    mild_liver_disease = list(
      icd9  = c("07022","07023","07032","07033","07044","07054","0706","0709","5710",
                "5712","5713","5714","5718","5719","5730","5731","5732","5733","5734",
                "5738","5739","V427"),
      icd10 = c("B18","K700","K701","K702","K703","K709","K713","K714","K715","K717",
                "K73","K74","K760","K762","K763","K764","K768","K769","Z944")
    ),
    diabetes_wo_complications = list(
      icd9  = c("2500","2501","2502","2503","2508","2509"),
      icd10 = c("E100","E101","E106","E108","E109","E110","E111","E116","E118","E119",
                "E120","E121","E126","E128","E129","E130","E131","E136","E138","E139",
                "E140","E141","E146","E148","E149")
    ),
    diabetes_w_complications = list(
      icd9  = c("2504","2505","2506","2507"),
      icd10 = c("E102","E103","E104","E105","E107","E112","E113","E114","E115","E117",
                "E122","E123","E124","E125","E127","E132","E133","E134","E135","E137",
                "E142","E143","E144","E145","E147")
    ),
    hemiplegia_paraplegia = list(
      icd9  = c("3341","342","343","3440","3441","3442","3443","3444","3445","3446","3449"),
      icd10 = c("G041","G114","G801","G802","G81","G82","G830","G831","G832","G833",
                "G834","G839")
    ),
    renal_disease = list(
      icd9  = c("40301","40311","40391","40402","40412","40492","582","5830","5831",
                "5832","5834","5836","5837","5838","5839","585","586","5880","V420","V451"),
      icd10 = c("I120","I131","N032","N033","N034","N035","N036","N037","N052","N053",
                "N054","N055","N056","N057","N18","N19","N250","Z490","Z491","Z492",
                "Z940","Z992")
    ),
    malignancy = list(
      icd9  = c("140","141","142","143","144","145","146","147","148","149","150","151",
                "152","153","154","155","156","157","158","159","160","161","162","163",
                "164","165","166","167","168","169","170","171","172","174","175","176",
                "179","180","181","182","183","184","185","186","187","188","189","190",
                "191","192","193","194","195","196","197","198","199","200","201","202",
                "203","204","205","206","207","208","2386"),
      icd10 = c("C0","C1","C2","C3","C4","C5","C6","C70","C71","C72","C73","C74","C75",
                "C76","C81","C82","C83","C84","C85","C86","C87","C88","C89","C90","C91",
                "C92","C93","C94","C95","C96","C97")
    ),
    severe_liver_disease = list(
      icd9  = c("4560","4561","4562","5722","5723","5724","5725","5726","5727","5728"),
      icd10 = c("I850","I859","I864","I982","K704","K711","K721","K729","K765","K766",
                "K767")
    ),
    metastatic_solid_tumor = list(
      icd9  = c("196","197","198","199"),
      icd10 = c("C77","C78","C79","C80")
    ),
    aids = list(
      icd9  = c("042","0440"),
      icd10 = c("B20","B21","B22","B24")
    )
  )
}

# Poids Charlson
CHARLSON_WEIGHTS <- c(
  myocardial_infarction = 1, congestive_heart_failure = 1,
  peripheral_vascular = 1, cerebrovascular = 1, dementia = 1,
  chronic_pulmonary = 1, rheumatic_disease = 1, peptic_ulcer = 1,
  mild_liver_disease = 1, diabetes_wo_complications = 1,
  diabetes_w_complications = 2, hemiplegia_paraplegia = 2,
  renal_disease = 2, malignancy = 2, severe_liver_disease = 3,
  metastatic_solid_tumor = 6, aids = 6
)

message("[02_codebooks] Codebooks ICD et item_ids chargés.")
