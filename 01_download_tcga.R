# 분석에 사용할 폴더 만들기
dir.create("data", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# TCGA 데이터를 다루기 위한 패키지 설치

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
packageVersion("TCGAbiolinks")

# TCGAbiolinks 불러오기
library(TCGAbiolinks)

# TCGA-LUAD 임상정보 가져오기
clinical_luad <- GDCquery_clinic(
  project = "TCGA-LUAD",
  type = "clinical"
)

# 데이터 크기 확인
dim(clinical_luad)

# 앞부분 확인
head(clinical_luad)

# 실제 고유 환자 수
length(unique(clinical_luad$bcr_patient_barcode))

# 같은 환자가 여러 번 들어있는지 확인
sum(duplicated(clinical_luad$bcr_patient_barcode))

# 생존 상태 확인
table(clinical_luad$vital_status, useNA = "ifany")

# 방사선치료 여부 확인
table(
  clinical_luad$treatments_radiation_treatment_or_therapy,
  useNA = "ifany"
)

# 생존분석용 데이터 만들기
clinical_surv <- clinical_luad

# 사망 여부: Dead = 1, Alive = 0
clinical_surv$OS_event <- ifelse(
  clinical_surv$vital_status == "Dead", 1,
  ifelse(clinical_surv$vital_status == "Alive", 0, NA)
)

# 생존기간:
# 사망자는 days_to_death
# 생존자는 days_to_last_follow_up
clinical_surv$OS_days <- ifelse(
  clinical_surv$OS_event == 1,
  clinical_surv$days_to_death,
  clinical_surv$days_to_last_follow_up
)

# 숫자형으로 변환
clinical_surv$OS_days <- as.numeric(clinical_surv$OS_days)

# 확인
table(clinical_surv$OS_event, useNA = "ifany")
summary(clinical_surv$OS_days)

clinical_surv$RT_status <- ifelse(
  clinical_surv$treatments_radiation_treatment_or_therapy == "yes", "RT",
  ifelse(
    clinical_surv$treatments_radiation_treatment_or_therapy == "no",
    "No_RT",
    NA
  )
)

table(clinical_surv$RT_status, useNA = "ifany")

# TCGA-LUAD somatic mutation 데이터 검색
query_mut <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access = "open"
)

query_mut

# mutation 데이터 다운로드
GDCdownload(
  query_mut,
  method = "api",
  directory = "data/GDC_mutation"
)

unlink("data/GDC_mutation", recursive = TRUE)

GDCdownload(
  query_mut,
  method = "api",
  directory = "data/GDC_mutation",
  files.per.chunk = 50
)

# 다운로드된 파일 수 확인
length(list.files("data/GDC_mutation", recursive = TRUE))

mut_luad <- GDCprepare(
  query_mut,
  directory = "data/GDC_mutation"
)

dim(mut_luad)

colnames(mut_luad)

head(mut_luad)

# 환자 ID 만들기
# TCGA-XX-YYYY-01A... → TCGA-XX-YYYY
mut_luad$patient_id <- substr(
  mut_luad$Tumor_Sample_Barcode,
  1, 12
)

# KEAP1 mutation만 추출
keap1_mut <- mut_luad[
  mut_luad$Hugo_Symbol == "KEAP1",
]

# KEAP1 mutation 기록 수
nrow(keap1_mut)

# KEAP1 mutation이 있는 고유 환자 수
length(unique(keap1_mut$patient_id))

# 어떤 종류의 KEAP1 mutation인지 확인
table(
  keap1_mut$Variant_Classification,
  useNA = "ifany"
)

# 단백질에 영향을 줄 가능성이 있는 KEAP1 mutation만 선택
functional_classes <- c(
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Splice_Site"
)

keap1_mut_functional <- keap1_mut[
  keap1_mut$Variant_Classification %in% functional_classes,
]

# 기능성 KEAP1 mutation 환자 수
length(unique(keap1_mut_functional$patient_id))

# Silent mutation 환자 ID
silent_patients <- unique(
  keap1_mut$patient_id[
    keap1_mut$Variant_Classification == "Silent"
  ]
)

silent_patients

# 이 환자가 다른 KEAP1 mutation도 가지고 있는지 확인
keap1_mut[
  keap1_mut$patient_id %in% silent_patients,
  c(
    "patient_id",
    "Variant_Classification",
    "HGVSp_Short"
  )
]

# query에서 실제 mutation 파일 정보 확인
mut_results <- getResults(query_mut)

# 결과 크기
dim(mut_results)

# 어떤 컬럼들이 있는지 확인
colnames(mut_results)

# 환자 / case 관련 컬럼 확인
mut_results[
  1:6,
  grep(
    "case|patient|barcode|sample",
    colnames(mut_results),
    ignore.case = TRUE
  ),
  drop = FALSE
]

# mutation sequencing이 수행된 환자 ID 추출
sequenced_patients <- unique(
  substr(mut_results$cases, 1, 12)
)

# 환자 수 확인
length(sequenced_patients)

# 앞부분 확인
head(sequenced_patients)

# 기능성/non-silent KEAP1 mutation 환자 목록
keap1_mut_patients <- unique(
  keap1_mut_functional$patient_id
)

# 혹시 KEAP1 MUT 환자가 모두 sequenced 환자 안에 있는지 확인
all(keap1_mut_patients %in% sequenced_patients)

# 환자별 KEAP1 상태 만들기
keap1_status <- data.frame(
  patient_id = sequenced_patients
)

keap1_status$KEAP1_status <- ifelse(
  keap1_status$patient_id %in% keap1_mut_patients,
  "MUT",
  "WT"
)

# WT / MUT 환자 수 확인
table(keap1_status$KEAP1_status)

# clinical 데이터에 patient_id 만들기
clinical_surv$patient_id <- clinical_surv$bcr_patient_barcode

# clinical + KEAP1 상태 합치기
analysis_df <- merge(
  clinical_surv,
  keap1_status,
  by = "patient_id"
)

# 합쳐진 환자 수
nrow(analysis_df)

# KEAP1 WT / MUT 환자 수
table(analysis_df$KEAP1_status)

# 생존정보까지 있는 환자 수 확인
table(
  analysis_df$KEAP1_status,
  is.na(analysis_df$OS_days)
)

# 생존정보가 있는 환자만 추출
surv_df <- analysis_df[
  !is.na(analysis_df$OS_days) &
    !is.na(analysis_df$OS_event),
]

# 최종 생존분석 환자 수
nrow(surv_df)

# KEAP1 WT / MUT 환자 수
table(surv_df$KEAP1_status)

# 각 그룹에서 Alive(0) / Dead(1) 환자 수
table(
  surv_df$KEAP1_status,
  surv_df$OS_event
)

# survival 패키지 설치/불러오기
if (!requireNamespace("survival", quietly = TRUE)) {
  install.packages("survival")
}

if (!requireNamespace("survminer", quietly = TRUE)) {
  install.packages("survminer")
}

library(survival)
library(survminer)

# Kaplan-Meier survival model
fit_keap1 <- survfit(
  Surv(OS_days, OS_event) ~ KEAP1_status,
  data = surv_df
)

# 결과 확인
fit_keap1

# Log-rank test
survdiff(
  Surv(OS_days, OS_event) ~ KEAP1_status,
  data = surv_df
)

ggsurvplot(
  fit_keap1,
  data = surv_df,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Days",
  ylab = "Overall survival probability",
  legend.title = "KEAP1",
  legend.labs = c("MUT", "WT")
)

# TCGA-LUAD miRNA-seq 데이터 검색
query_mirna <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  sample.type = "Primary Tumor"
)

# 검색 결과를 표로 저장
mirna_results <- getResults(query_mirna)

# 몇 개 파일이 검색됐는지 확인
dim(mirna_results)

# 앞부분 확인
head(mirna_results)

# 오늘 작업 상태 저장
save.image("data/workspace_checkpoint.RData")

file.exists("data/workspace_checkpoint.RData")

# miRNA 데이터 다운로드
GDCdownload(
  query_mirna,
  method = "api",
  directory = "data/GDC_miRNA",
  files.per.chunk = 50
)

# 다운로드된 파일 수 확인
length(list.files("data/GDC_miRNA", recursive = TRUE))

# 다운로드한 miRNA 파일들을 하나의 R 객체로 합치기
mirna_luad <- GDCprepare(
  query_mirna,
  directory = "data/GDC_miRNA"
)

# 어떤 형태의 데이터인지 확인
class(mirna_luad)

# 데이터 크기 확인
dim(mirna_luad)

# 컬럼 이름 확인
colnames(mirna_luad)

# 앞부분 확인
head(mirna_luad)

# miR-130b가 들어간 miRNA 이름 찾기
mirna_luad[
  grepl("130b", mirna_luad$miRNA_ID, ignore.case = TRUE),
  "miRNA_ID"
]

grep(
  "130b",
  mirna_luad$miRNA_ID,
  value = TRUE,
  ignore.case = TRUE
)

query_mirna_iso <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Isoform Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  sample.type = "Primary Tumor"
)

iso_results <- getResults(query_mirna_iso)

dim(iso_results)
head(iso_results)

# isoform 파일의 고유 환자 수
length(unique(iso_results$cases.submitter_id))

# 같은 환자가 여러 파일을 가지고 있는 경우가 몇 개인지
sum(duplicated(iso_results$cases.submitter_id))

# 앞에서 받은 miRNA quantification과
# isoform quantification의 샘플 구성이 동일한지 확인
setequal(
  mirna_results$sample.submitter_id,
  iso_results$sample.submitter_id
)

# 환자별 isoform 파일 개수 세기
patient_file_count <- table(
  iso_results$cases.submitter_id
)

# 2개 이상의 파일을 가진 환자만 확인
patient_file_count[
  patient_file_count > 1
]

duplicate_patients <- names(
  patient_file_count[
    patient_file_count > 1
  ]
)

iso_results[
  iso_results$cases.submitter_id %in% duplicate_patients,
  c(
    "cases.submitter_id",
    "sample.submitter_id",
    "cases",
    "file_name"
  )
]

GDCdownload(
  query_mirna_iso,
  method = "api",
  directory = "data/GDC_miRNA_isoform",
  files.per.chunk = 50
)

length(
  list.files(
    "data/GDC_miRNA_isoform",
    recursive = TRUE
  )
)

mirna_iso_luad <- GDCprepare(
  query_mirna_iso,
  directory = "data/GDC_miRNA_isoform"
)

class(mirna_iso_luad)
dim(mirna_iso_luad)
colnames(mirna_iso_luad)[1:15]
head(mirna_iso_luad[, 1:15])
head(mirna_iso_luad)
head(mirna_iso_luad[, 1:7])

# hsa-mir-130b isoform만 추출
mir130b_iso <- mirna_iso_luad[
  mirna_iso_luad$miRNA_ID == "hsa-mir-130b",
]

# 몇 행이 있는지 확인
dim(mir130b_iso)

# 앞부분 확인
head(mir130b_iso)

# miR-130b의 region 종류 확인
unique(mir130b_iso$miRNA_region)

table(
  mir130b_iso$miRNA_region,
  useNA = "ifany"
)

# miR-130b-5p에 해당하는 mature isoform만 추출
mir130b_5p_iso <- mir130b_iso[
  mir130b_iso$miRNA_region == "mature,MIMAT0004680",
]

# 크기 확인
dim(mir130b_5p_iso)

# 몇 개의 샘플이 들어있는지 확인
length(unique(mir130b_5p_iso$barcode))

# 앞부분 확인
head(mir130b_5p_iso)

# 같은 aliquot 안의 miR-130b-5p isoform RPM 합산
mir130b_5p_aliquot <- aggregate(
  reads_per_million_miRNA_mapped ~ barcode,
  data = mir130b_5p_iso,
  FUN = sum
)

# 컬럼 이름 변경
colnames(mir130b_5p_aliquot)[2] <- "miR130b_5p_RPM"

# 확인
dim(mir130b_5p_aliquot)
head(mir130b_5p_aliquot)

# sample ID: 예) TCGA-05-4244-01A
mir130b_5p_aliquot$sample_id <- substr(
  mir130b_5p_aliquot$barcode,
  1, 16
)

# patient ID: 예) TCGA-05-4244
mir130b_5p_aliquot$patient_id <- substr(
  mir130b_5p_aliquot$barcode,
  1, 12
)

# 확인
head(mir130b_5p_aliquot)

# 같은 tumor sample의 반복 aliquot를 평균
mir130b_5p_sample <- aggregate(
  miR130b_5p_RPM ~ sample_id + patient_id,
  data = mir130b_5p_aliquot,
  FUN = mean
)

dim(mir130b_5p_sample)
head(mir130b_5p_sample)

# 환자당 miR-130b-5p 값 하나로 정리
mir130b_5p_patient <- aggregate(
  miR130b_5p_RPM ~ patient_id,
  data = mir130b_5p_sample,
  FUN = mean
)

# 최종 환자 수 확인
dim(mir130b_5p_patient)
# 앞부분 확인
head(mir130b_5p_patient)

mir_keap1_df <- merge(
  mir130b_5p_patient,
  keap1_status,
  by = "patient_id"
)

# 합쳐진 환자 수
dim(mir_keap1_df)

# WT / MUT 각각 몇 명인지
table(mir_keap1_df$KEAP1_status)

# 앞부분 확인
head(mir_keap1_df)

# WT / MUT별 miR-130b-5p 발현량 요약
aggregate(
  miR130b_5p_RPM ~ KEAP1_status,
  data = mir_keap1_df,
  FUN = function(x) c(
    n = length(x),
    median = median(x),
    mean = mean(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
)

wilcox.test(
  miR130b_5p_RPM ~ KEAP1_status,
  data = mir_keap1_df,
  exact = FALSE
)

# hsa-mir-130b 행 찾기
mir130b_precursor <- mirna_luad[
  mirna_luad$miRNA_ID == "hsa-mir-130b",
]

dim(mir130b_precursor)

rpm_cols <- grep(
  "^reads_per_million_miRNA_mapped_",
  colnames(mirna_luad),
  value = TRUE
)

length(rpm_cols)

# 519개 RPM 값을 세로 형태로 변환
mir130b_precursor_aliquot <- data.frame(
  barcode = sub(
    "^reads_per_million_miRNA_mapped_",
    "",
    rpm_cols
  ),
  miR130b_precursor_RPM = as.numeric(
    unlist(
      mir130b_precursor[1, rpm_cols],
      use.names = FALSE
    )
  )
)

dim(mir130b_precursor_aliquot)

head(mir130b_precursor_aliquot)

# sample ID와 patient ID 만들기
mir130b_precursor_aliquot$sample_id <- substr(
  mir130b_precursor_aliquot$barcode,
  1, 16
)

mir130b_precursor_aliquot$patient_id <- substr(
  mir130b_precursor_aliquot$barcode,
  1, 12
)

# 같은 sample의 반복 aliquot 평균
mir130b_precursor_sample <- aggregate(
  miR130b_precursor_RPM ~ sample_id + patient_id,
  data = mir130b_precursor_aliquot,
  FUN = mean
)

# 환자당 하나의 값으로 평균
mir130b_precursor_patient <- aggregate(
  miR130b_precursor_RPM ~ patient_id,
  data = mir130b_precursor_sample,
  FUN = mean
)

dim(mir130b_precursor_sample)
dim(mir130b_precursor_patient)

precursor_keap1_df <- merge(
  mir130b_precursor_patient,
  keap1_status,
  by = "patient_id"
)

dim(precursor_keap1_df)
table(precursor_keap1_df$KEAP1_status)

# KEAP1 WT / MUT별 precursor hsa-mir-130b 발현량 요약
aggregate(
  miR130b_precursor_RPM ~ KEAP1_status,
  data = precursor_keap1_df,
  FUN = function(x) c(
    n = length(x),
    median = median(x),
    mean = mean(x),
    Q1 = quantile(x, 0.25),
    Q3 = quantile(x, 0.75)
  )
)

wilcox.test(
  miR130b_precursor_RPM ~ KEAP1_status,
  data = precursor_keap1_df,
  exact = FALSE
)

# mature miR-130b-5p + KEAP1 + survival 정보 합치기
mir_surv_df <- merge(
  mir_keap1_df,
  clinical_surv[
    , c(
      "patient_id",
      "OS_days",
      "OS_event"
    )
  ],
  by = "patient_id"
)

# 전체 크기
dim(mir_surv_df)

# KEAP1 MUT 환자만 추출
mut_mir_surv <- mir_surv_df[
  mir_surv_df$KEAP1_status == "MUT" &
    !is.na(mir_surv_df$OS_days) &
    !is.na(mir_surv_df$OS_event),
]

# 분석 가능한 MUT 환자 수
nrow(mut_mir_surv)

# Alive / Dead 수
table(mut_mir_surv$OS_event)

# miR-130b-5p 값 확인
summary(mut_mir_surv$miR130b_5p_RPM)

# KEAP1 MUT + miRNA 생존 cohort에서
# OS_days가 0인 환자가 있는지 확인
table(mut_mir_surv$OS_days == 0)

# 만약 있다면 어떤 환자인지 확인
mut_mir_surv[
  mut_mir_surv$OS_days == 0,
  c(
    "patient_id",
    "OS_days",
    "OS_event",
    "miR130b_5p_RPM"
  )
]

clinical_surv[
  clinical_surv$patient_id %in% c(
    "TCGA-05-4395",
    "TCGA-86-8281"
  ),
  c(
    "patient_id",
    "vital_status",
    "days_to_death",
    "days_to_last_follow_up",
    "OS_days",
    "OS_event"
  )
]

# 생존기간이 0보다 큰 KEAP1 MUT 환자만 사용
mut_mir_surv_clean <- mut_mir_surv[
  mut_mir_surv$OS_days > 0,
]

# 확인
nrow(mut_mir_surv_clean)
table(mut_mir_surv_clean$OS_event)

cut_mir <- surv_cutpoint(
  mut_mir_surv_clean,
  time = "OS_days",
  event = "OS_event",
  variables = "miR130b_5p_RPM",
  minprop = 0.1
)

cut_mir$cutpoint

# cutoff를 이용해서 Low / High 그룹 만들기
cut_mir_cat <- surv_categorize(cut_mir)

# 각 그룹 환자 수
table(cut_mir_cat$miR130b_5p_RPM)

# Low / High 그룹 붙이기
mut_mir_surv_clean$miR_group <- cut_mir_cat$miR130b_5p_RPM

# 순서 명확하게 지정
mut_mir_surv_clean$miR_group <- factor(
  mut_mir_surv_clean$miR_group,
  levels = c("low", "high")
)

# 확인
table(mut_mir_surv_clean$miR_group)

# cutoff를 직접 적용
mut_mir_surv_clean$miR_group <- ifelse(
  mut_mir_surv_clean$miR130b_5p_RPM <= 7.245196,
  "low",
  "high"
)

mut_mir_surv_clean$miR_group <- factor(
  mut_mir_surv_clean$miR_group,
  levels = c("low", "high")
)

table(mut_mir_surv_clean$miR_group)

fit_mir <- survfit(
  Surv(OS_days, OS_event) ~ miR_group,
  data = mut_mir_surv_clean
)

fit_mir

# Log-rank test
survdiff(
  Surv(OS_days, OS_event) ~ miR_group,
  data = mut_mir_surv_clean
)

ggsurvplot(
  fit_mir,
  data = mut_mir_surv_clean,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = TRUE,
  xlab = "Days",
  ylab = "Overall survival probability",
  legend.title = "miR-130b-5p",
  legend.labs = c("Low", "High")
)

# 발현량 분포가 한쪽으로 치우쳐 있으므로 log2 변환
mut_mir_surv_clean$log2_miR130b_5p <- log2(
  mut_mir_surv_clean$miR130b_5p_RPM + 1
)

# miR-130b-5p 발현량 자체와 생존의 관계
cox_mir <- coxph(
  Surv(OS_days, OS_event) ~ log2_miR130b_5p,
  data = mut_mir_surv_clean
)

summary(cox_mir)

# miRNA 데이터가 있는 KEAP1 MUT 환자 ID
current_mut_ids <- mir_keap1_df$patient_id[
  mir_keap1_df$KEAP1_status == "MUT"
]

# 이 환자들의 KEAP1 mutation 상세정보
keap1_mut_check <- keap1_mut[
  keap1_mut$patient_id %in% current_mut_ids,
  c(
    "patient_id",
    "Variant_Classification",
    "HGVSp_Short"
  )
]

# mutation 종류별 개수
table(keap1_mut_check$Variant_Classification)

# 환자 수 확인
length(unique(keap1_mut_check$patient_id))

table(
  mir_surv_df$KEAP1_status,
  complete.cases(
    mir_surv_df$OS_days,
    mir_surv_df$OS_event
  )
)
# 생존정보가 있고, 추적기간도 0보다 큰 환자만
mir_surv_positive <- mir_surv_df[
  !is.na(mir_surv_df$OS_days) &
    !is.na(mir_surv_df$OS_event) &
    mir_surv_df$OS_days > 0,
]
nrow(mir_surv_positive)
table(mir_surv_positive$KEAP1_status)

table(
  mir_surv_complete$KEAP1_status,
  mir_surv_complete$OS_days == 0
)

mir_surv_complete <- mir_surv_df[
  !is.na(mir_surv_df$OS_days) &
    !is.na(mir_surv_df$OS_event),
]

table(
  mir_surv_complete$KEAP1_status,
  mir_surv_complete$OS_days == 0
)

mir_keap1_df$KEAP1_status <- factor(
  mir_keap1_df$KEAP1_status,
  levels = c("WT", "MUT")
)

table(mir_keap1_df$KEAP1_status)
mir_keap1_df$miR130b_5p_log2 <- log2(
  mir_keap1_df$miR130b_5p_RPM + 1
)
summary(mir_keap1_df$miR130b_5p_log2)

aggregate(
  miR130b_5p_log2 ~ KEAP1_status,
  data = mir_keap1_df,
  FUN = function(x) c(
    n = length(x),
    median = median(x),
    mean = mean(x)
  )
)

wilcox_log2 <- wilcox.test(
  miR130b_5p_log2 ~ KEAP1_status,
  data = mir_keap1_df,
  exact = FALSE
)

wilcox_log2
library(ggplot2)
ggplot(
  mir_keap1_df,
  aes(
    x = KEAP1_status,
    y = miR130b_5p_log2,
    fill = KEAP1_status
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  labs(
    x = "KEAP1 status",
    y = "miR-130b-5p (log2 RPM + 1)",
    title = "miR-130b-5p expression by KEAP1 status"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

ggplot(
  mir_keap1_df,
  aes(
    x = KEAP1_status,
    y = miR130b_5p_log2,
    fill = KEAP1_status
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.3,
    size = 1
  ) +
  labs(
    x = "KEAP1 status",
    y = "miR-130b-5p (log2 RPM + 1)",
    title = "miR-130b-5p expression by KEAP1 status"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )
p_text <- paste0(
  "Wilcoxon p = ",
  signif(wilcox_log2$p.value, 3)
)

p_text

ggplot(
  mir_keap1_df,
  aes(
    x = KEAP1_status,
    y = miR130b_5p_log2,
    fill = KEAP1_status
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.3,
    size = 1
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 8.0,
    label = p_text,
    size = 5
  ) +
  labs(
    x = "KEAP1 status",
    y = "miR-130b-5p (log2 RPM + 1)",
    title = "miR-130b-5p expression by KEAP1 status"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

mir_sample_keap1_df <- merge(
  mir130b_5p_sample,
  keap1_status,
  by = "patient_id"
)

# 전체 sample 수
dim(mir_sample_keap1_df)

# WT / MUT sample 수
table(mir_sample_keap1_df$KEAP1_status)

# 앞부분 확인
head(mir_sample_keap1_df)

mir_sample_keap1_df$miR130b_5p_log2 <- log2(
  mir_sample_keap1_df$miR130b_5p_RPM + 1
)

aggregate(
  miR130b_5p_log2 ~ KEAP1_status,
  data = mir_sample_keap1_df,
  FUN = function(x) c(
    n = length(x),
    median = median(x),
    mean = mean(x)
  )
)


wilcox_sample <- wilcox.test(
  miR130b_5p_log2 ~ KEAP1_status,
  data = mir_sample_keap1_df,
  exact = FALSE
)

wilcox_sample


mir130b_5p_count_aliquot <- aggregate(
  read_count ~ barcode,
  data = mir130b_5p_iso,
  FUN = sum
)

colnames(mir130b_5p_count_aliquot)[2] <- "miR130b_5p_count"

dim(mir130b_5p_count_aliquot)

head(mir130b_5p_count_aliquot)


# sample ID와 patient ID 만들기
mir130b_5p_count_aliquot$sample_id <- substr(
  mir130b_5p_count_aliquot$barcode,
  1, 16
)

mir130b_5p_count_aliquot$patient_id <- substr(
  mir130b_5p_count_aliquot$barcode,
  1, 12
)

# 같은 tumor sample의 반복 aliquot를 평균
mir130b_5p_count_sample <- aggregate(
  miR130b_5p_count ~ sample_id + patient_id,
  data = mir130b_5p_count_aliquot,
  FUN = mean
)

dim(mir130b_5p_count_sample)

head(mir130b_5p_count_sample)


mir_count_sample_keap1 <- merge(
  mir130b_5p_count_sample,
  keap1_status,
  by = "patient_id"
)

dim(mir_count_sample_keap1)

table(mir_count_sample_keap1$KEAP1_status)

mir_count_sample_keap1$miR130b_5p_log2count <- log2(
  mir_count_sample_keap1$miR130b_5p_count + 1
)


aggregate(
  miR130b_5p_log2count ~ KEAP1_status,
  data = mir_count_sample_keap1,
  FUN = function(x) c(
    n = length(x),
    median = median(x),
    mean = mean(x)
  )
)

summary(
  mir_count_sample_keap1$miR130b_5p_log2count
)

wilcox_count <- wilcox.test(
  miR130b_5p_log2count ~ KEAP1_status,
  data = mir_count_sample_keap1,
  exact = FALSE
)

wilcox_count


ggplot(
  mir_count_sample_keap1,
  aes(
    x = KEAP1_status,
    y = miR130b_5p_log2count,
    fill = KEAP1_status
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.3,
    size = 1
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 9.7,
    label = "Wilcoxon p = 0.0429",
    size = 5
  ) +
  labs(
    x = "KEAP1 status",
    y = "miR-130b-5p (log2 count + 1)",
    title = "miR-130b-5p expression by KEAP1 status"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

mir_count_sample_keap1$KEAP1_status <- factor(
  mir_count_sample_keap1$KEAP1_status,
  levels = c("WT", "MUT")
)

ggplot(
  mir_count_sample_keap1,
  aes(
    x = KEAP1_status,
    y = miR130b_5p_log2count,
    fill = KEAP1_status
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.3,
    size = 1
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 9.7,
    label = "Wilcoxon p = 0.0429",
    size = 5
  ) +
  labs(
    x = "KEAP1 status",
    y = "miR-130b-5p (log2 count + 1)",
    title = "miR-130b-5p expression by KEAP1 status"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )


ggplot(
  mir_count_sample_keap1,
  aes(
    x = KEAP1_status,
    y = miR130b_5p_log2count,
    fill = KEAP1_status
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.3,
    size = 1
  ) +
  annotate(
    "text",
    x = 1.5,
    y = 9.7,
    label = "Wilcoxon p = 0.0429",
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "WT" = "#00BFC4",
      "MUT" = "#F28E2B"
    )
  ) +
  labs(
    x = "KEAP1 status",
    y = "miR-130b-5p (log2 count + 1)",
    title = "miR-130b-5p expression by KEAP1 status"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

length(unique(mir130b_5p_iso$isoform_coords))

sort(
  table(mir130b_5p_iso$isoform_coords),
  decreasing = TRUE
)[1:20]

mir130b_5p_exact <- mir130b_5p_iso[
  mir130b_5p_iso$isoform_coords ==
    "hg38:chr22:21653316-21653336:+",
]

dim(mir130b_5p_exact)

length(unique(mir130b_5p_exact$barcode))

head(mir130b_5p_exact)

# 전체 519개 tumor aliquot 목록
all_barcodes <- data.frame(
  barcode = unique(mir130b_5p_aliquot$barcode)
)

# exact mature miR-130b-5p count만 가져오기
mir130b_5p_exact_count <- as.data.frame(
  mir130b_5p_exact[, c("barcode", "read_count")]
)

# 전체 519개 샘플에 붙이기
mir130b_5p_exact_all <- merge(
  all_barcodes,
  mir130b_5p_exact_count,
  by = "barcode",
  all.x = TRUE
)

# exact isoform row가 없었던 샘플은 count = 0
mir130b_5p_exact_all$read_count[
  is.na(mir130b_5p_exact_all$read_count)
] <- 0

# 확인
dim(mir130b_5p_exact_all)

table(mir130b_5p_exact_all$read_count == 0)

head(mir130b_5p_exact_all)


# sample ID / patient ID 만들기
mir130b_5p_exact_all$sample_id <- substr(
  mir130b_5p_exact_all$barcode,
  1, 16
)

mir130b_5p_exact_all$patient_id <- substr(
  mir130b_5p_exact_all$barcode,
  1, 12
)

# 같은 tumor sample의 반복 aliquot를 평균
mir130b_5p_exact_sample <- aggregate(
  read_count ~ sample_id + patient_id,
  data = mir130b_5p_exact_all,
  FUN = mean
)

# 확인
dim(mir130b_5p_exact_sample)



mir_exact_sample_keap1 <- merge(
  mir130b_5p_exact_sample,
  keap1_status,
  by = "patient_id"
)

dim(mir_exact_sample_keap1)

table(mir_exact_sample_keap1$KEAP1_status)

head(mir_exact_sample_keap1)


mir_exact_sample_keap1$miR130b_5p_exact_log2 <- log2(
  mir_exact_sample_keap1$read_count + 1
)


aggregate(
  miR130b_5p_exact_log2 ~ KEAP1_status,
  data = mir_exact_sample_keap1,
  FUN = function(x) c(
    n = length(x),
    median = median(x),
    mean = mean(x)
  )
)


wilcox_exact <- wilcox.test(
  miR130b_5p_exact_log2 ~ KEAP1_status,
  data = mir_exact_sample_keap1,
  exact = FALSE
)

wilcox_exact


length(unique(mirna_iso_luad$miRNA_region))

length(unique(mirna_iso_luad$miRNA_ID))

mirna_mature_all <- mirna_iso_luad[
  grepl("^mature,", mirna_iso_luad$miRNA_region),
]

dim(mirna_mature_all)

length(unique(mirna_mature_all$miRNA_region))

length(unique(mirna_mature_all$barcode))

head(mirna_mature_all)



mirna_mature_sum <- aggregate(
  read_count ~ miRNA_region + barcode,
  data = mirna_mature_all,
  FUN = sum
)

dim(mirna_mature_sum)

head(mirna_mature_sum)


mirna_count_matrix <- xtabs(
  read_count ~ miRNA_region + barcode,
  data = mirna_mature_sum
)

mirna_count_matrix <- as.matrix(mirna_count_matrix)

dim(mirna_count_matrix)

mirna_count_matrix[1:5, 1:5]

"mature,MIMAT0004680" %in% rownames(mirna_count_matrix)
head(mir130b_5p_exact_sample)


library(DESeq2)

# 519개 sample 정보를 DESeq2 형식으로 만들기
mirna_coldata <- data.frame(
  row.names = colnames(mirna_count_matrix),
  group = rep("all", ncol(mirna_count_matrix))
)

# DESeq2 객체 생성
dds_mirna <- DESeqDataSetFromMatrix(
  countData = mirna_count_matrix,
  colData = mirna_coldata,
  design = ~ 1
)

# 샘플별 sequencing depth 보정값 계산
dds_mirna <- estimateSizeFactors(
  dds_mirna,
  type = "poscounts"
)

# size factor 확인
summary(sizeFactors(dds_mirna))


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("DESeq2")

library(DESeq2)



mirna_coldata <- data.frame(
  row.names = colnames(mirna_count_matrix),
  group = rep("all", ncol(mirna_count_matrix))
)

dds_mirna <- DESeqDataSetFromMatrix(
  countData = mirna_count_matrix,
  colData = mirna_coldata,
  design = ~ 1
)

dds_mirna <- estimateSizeFactors(
  dds_mirna,
  type = "poscounts"
)

summary(sizeFactors(dds_mirna))



# DESeq2 normalized count matrix
mirna_norm_counts <- counts(
  dds_mirna,
  normalized = TRUE
)

# miR-130b-5p만 추출
mir130b_5p_norm <- data.frame(
  barcode = colnames(mirna_norm_counts),
  miR130b_5p_norm = as.numeric(
    mirna_norm_counts["mature,MIMAT0004680", ]
  )
)

dim(mir130b_5p_norm)

head(mir130b_5p_norm)

summary(mir130b_5p_norm$miR130b_5p_norm)



# sample ID와 patient ID 만들기
mir130b_5p_norm$sample_id <- substr(
  mir130b_5p_norm$barcode,
  1, 16
)

mir130b_5p_norm$patient_id <- substr(
  mir130b_5p_norm$barcode,
  1, 12
)

# 같은 tumor sample의 반복 aliquot는 평균
mir130b_5p_norm_sample <- aggregate(
  miR130b_5p_norm ~ sample_id + patient_id,
  data = mir130b_5p_norm,
  FUN = mean
)

dim(mir130b_5p_norm_sample)

head(mir130b_5p_norm_sample)


mir_norm_sample_keap1 <- merge(
  mir130b_5p_norm_sample,
  keap1_status,
  by = "patient_id"
)

dim(mir_norm_sample_keap1)

table(mir_norm_sample_keap1$KEAP1_status)



mir_norm_sample_keap1$miR130b_5p_log2norm <- log2(
  mir_norm_sample_keap1$miR130b_5p_norm + 1
)


aggregate(
  miR130b_5p_log2norm ~ KEAP1_status,
  data = mir_norm_sample_keap1,
  FUN = function(x) c(
    n = length(x),
    median = median(x),
    mean = mean(x)
  )
)

wilcox_norm <- wilcox.test(
  miR130b_5p_log2norm ~ KEAP1_status,
  data = mir_norm_sample_keap1,
  exact = FALSE
)

wilcox_norm


library(ggplot2)

# 1) WT가 왼쪽, MUT가 오른쪽에 오도록 순서 고정
mir_norm_sample_keap1$KEAP1_status <- factor(
  mir_norm_sample_keap1$KEAP1_status,
  levels = c("WT", "MUT")
)

# 2) p-value 계산
wilcox_norm <- wilcox.test(
  miR130b_5p_log2norm ~ KEAP1_status,
  data = mir_norm_sample_keap1,
  exact = FALSE
)

p_text <- paste0(
  "Wilcoxon p = ",
  signif(wilcox_norm$p.value, 4)
)

# 3) 최종 그림
p_final <- ggplot(
  mir_norm_sample_keap1,
  aes(
    x = KEAP1_status,
    y = miR130b_5p_log2norm,
    fill = KEAP1_status
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.3,
    size = 1
  ) +
  annotate(
    "text",
    x = 1.5,
    y = max(mir_norm_sample_keap1$miR130b_5p_log2norm, na.rm = TRUE) + 0.2,
    label = p_text,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "WT" = "#00BFC4",
      "MUT" = "#F28E2B"
    )
  ) +
  labs(
    x = "KEAP1 status",
    y = "miR-130b-5p (log2 normalized count + 1)",
    title = "miR-130b-5p expression by KEAP1 status"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# 4) 출력
p_final


mir130b_5p_norm_patient <- aggregate(
  miR130b_5p_norm ~ patient_id,
  data = mir130b_5p_norm_sample,
  FUN = mean
)

dim(mir130b_5p_norm_patient)

head(mir130b_5p_norm_patient)


mir130b_5p_norm_patient$miR130b_5p_log2norm <- log2(
  mir130b_5p_norm_patient$miR130b_5p_norm + 1
)

summary(
  mir130b_5p_norm_patient$miR130b_5p_log2norm
)


# miR-130b-5p + KEAP1 상태 합치기
mir_surv_norm <- merge(
  mir130b_5p_norm_patient,
  keap1_status,
  by = "patient_id"
)

# 생존정보(OS) 합치기
mir_surv_norm <- merge(
  mir_surv_norm,
  clinical_surv[, c("patient_id", "OS_days", "OS_event")],
  by = "patient_id"
)

# 전체 확인
dim(mir_surv_norm)

table(mir_surv_norm$KEAP1_status)



mut_mir_surv_norm <- mir_surv_norm[
  mir_surv_norm$KEAP1_status == "MUT" &
    !is.na(mir_surv_norm$OS_days) &
    !is.na(mir_surv_norm$OS_event) &
    mir_surv_norm$OS_days > 0,
]

dim(mut_mir_surv_norm)

table(mut_mir_surv_norm$OS_event)

summary(mut_mir_surv_norm$miR130b_5p_log2norm)



library(survminer)

cut_mir_norm <- surv_cutpoint(
  mut_mir_surv_norm,
  time = "OS_days",
  event = "OS_event",
  variables = "miR130b_5p_log2norm",
  minprop = 0.1
)

cut_mir_norm

cut_mir_norm$cutpoint



cutoff_mir <- 4.479162

mut_mir_surv_norm$miR_group <- ifelse(
  mut_mir_surv_norm$miR130b_5p_log2norm <= cutoff_mir,
  "Low",
  "High"
)

# 순서 고정
mut_mir_surv_norm$miR_group <- factor(
  mut_mir_surv_norm$miR_group,
  levels = c("Low", "High")
)

# 각 그룹 환자 수
table(mut_mir_surv_norm$miR_group)


table(
  mut_mir_surv_norm$miR_group,
  mut_mir_surv_norm$OS_event
)


mut_mir_surv_cat <- surv_categorize(
  cut_mir_norm
)

table(mut_mir_surv_cat$miR130b_5p_log2norm)

head(mut_mir_surv_cat)

table(
  mut_mir_surv_cat$miR130b_5p_log2norm,
  mut_mir_surv_cat$OS_event
)

library(survival)

# 그룹 이름을 보기 좋게 정리
mut_mir_surv_cat$miR_group <- factor(
  mut_mir_surv_cat$miR130b_5p_log2norm,
  levels = c("low", "high"),
  labels = c("Low", "High")
)

# Kaplan-Meier 모델
fit_mir_norm <- survfit(
  Surv(OS_days, OS_event) ~ miR_group,
  data = mut_mir_surv_cat
)

fit_mir_norm

logrank_mir_norm <- survdiff(
  Surv(OS_days, OS_event) ~ miR_group,
  data = mut_mir_surv_cat
)

logrank_mir_norm


p_logrank <- 1 - pchisq(
  logrank_mir_norm$chisq,
  df = 1
)

p_logrank

cut_mir_norm_40 <- surv_cutpoint(
  mut_mir_surv_norm,
  time = "OS_days",
  event = "OS_event",
  variables = "miR130b_5p_log2norm",
  minprop = 0.4
)

cut_mir_norm_40


mut_mir_surv_cat_40 <- surv_categorize(
  cut_mir_norm_40
)

table(
  mut_mir_surv_cat_40$miR130b_5p_log2norm
)


# 보기 좋은 그룹 이름 만들기
mut_mir_surv_cat_40$miR_group <- factor(
  mut_mir_surv_cat_40$miR130b_5p_log2norm,
  levels = c("low", "high"),
  labels = c("Low", "High")
)

# Kaplan-Meier
fit_mir_norm_40 <- survfit(
  Surv(OS_days, OS_event) ~ miR_group,
  data = mut_mir_surv_cat_40
)

fit_mir_norm_40



logrank_mir_norm_40 <- survdiff(
  Surv(OS_days, OS_event) ~ miR_group,
  data = mut_mir_surv_cat_40
)

logrank_mir_norm_40


p_logrank_40 <- 1 - pchisq(
  logrank_mir_norm_40$chisq,
  df = 1
)

p_logrank_40