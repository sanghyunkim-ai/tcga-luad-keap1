# ============================================================
# TCGA-LUAD KEAP1 / miR-130b-5p / PIK3CB reproduction analysis
# ============================================================

# 현재 프로젝트 위치 확인
getwd()

# 프로젝트 안의 파일/폴더 확인
list.files

# ============================================================
# 1. 저장된 TCGA 데이터 확인
# ============================================================

list.files("data")

list.files(
  "data",
  recursive = FALSE,
  full.names = TRUE
)

# ============================================================
# 2. 필요한 패키지 설치 여부 확인
# ============================================================

packages_needed <- c(
  "TCGAbiolinks",
  "DESeq2",
  "edgeR",
  "survival",
  "survminer",
  "ggplot2"
)

sapply(
  packages_needed,
  requireNamespace,
  quietly = TRUE
)

library(edgeR)
requireNamespace("edgeR", quietly = TRUE)

# ============================================================
# 3. TCGA-LUAD clinical data
# ============================================================

library(TCGAbiolinks)

clinical_luad <- GDCquery_clinic(
  project = "TCGA-LUAD",
  type = "clinical"
)

dim(clinical_luad)

head(
  clinical_luad[, c(
    "bcr_patient_barcode",
    "vital_status",
    "days_to_death",
    "days_to_last_follow_up"
  )]
)

length(unique(clinical_luad$bcr_patient_barcode))

# ============================================================
# 4. Overall survival (OS) 만들기
# ============================================================

clinical_surv <- clinical_luad

# patient ID
clinical_surv$patient_id <- clinical_surv$bcr_patient_barcode

# 사망 여부
clinical_surv$OS_event <- ifelse(
  clinical_surv$vital_status == "Dead", 1,
  ifelse(
    clinical_surv$vital_status == "Alive", 0,
    NA
  )
)

# 전체생존기간
clinical_surv$OS_days <- ifelse(
  clinical_surv$OS_event == 1,
  clinical_surv$days_to_death,
  clinical_surv$days_to_last_follow_up
)

clinical_surv$OS_days <- as.numeric(
  clinical_surv$OS_days
)



table(
  clinical_surv$OS_event,
  useNA = "ifany"
)

summary(
  clinical_surv$OS_days
)


# ============================================================
# 5. Radiotherapy status 만들기
# ============================================================

clinical_surv$RT_status <- ifelse(
  clinical_surv$treatments_radiation_treatment_or_therapy == "yes",
  "RT",
  ifelse(
    clinical_surv$treatments_radiation_treatment_or_therapy == "no",
    "No_RT",
    NA
  )
)

table

# ============================================================
# 6. KEAP1 mutation status 만들기
# ============================================================

query_mut <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access = "open"
)

mut_luad <- GDCprepare(
  query_mut,
  directory = "data/GDC_mutation"
)

dim(mut_luad)

keap1_mut <- mut_luad[
  mut_luad$Hugo_Symbol == "KEAP1",
]

keap1_mut$patient_id <- substr(
  keap1_mut$Tumor_Sample_Barcode,
  1,
  12
)

dim(keap1_mut)

length(unique(keap1_mut$patient_id))

table(keap1_mut$Variant_Classification)


# ============================================================
# 7. KEAP1 WT / MUT 정의
# ============================================================

keap1_classes <- c(
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Splice_Site"
)

keap1_mut_nonsilent <- keap1_mut[
  keap1_mut$Variant_Classification %in% keap1_classes,
]

# non-silent KEAP1 mutation 환자
keap1_mut_patients <- unique(
  keap1_mut_nonsilent$patient_id
)

length(keap1_mut_patients)


mut_results <- getResults(query_mut)

sequenced_patients <- unique(
  substr(mut_results$cases, 1, 12)
)

length(sequenced_patients)


keap1_status <- data.frame(
  patient_id = sequenced_patients
)

keap1_status$KEAP1_status <- ifelse(
  keap1_status$patient_id %in% keap1_mut_patients,
  "MUT",
  "WT"
)

table(keap1_status$KEAP1_status)


# ============================================================
# 8. TCGA-LUAD miRNA isoform data 불러오기
# ============================================================

query_mirna_iso <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Isoform Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling",
  sample.type = "Primary Tumor"
)

mirna_iso_luad <- GDCprepare(
  query_mirna_iso,
  directory = "data/GDC_miRNA_isoform"
)

dim(mirna_iso_luad)
length(unique(mirna_iso_luad$barcode))
colnames(mirna_iso_luad)



# ============================================================
# 9. Mature miRNA raw count matrix 만들기
# ============================================================

# mature miRNA만 선택
mirna_mature <- mirna_iso_luad[
  grepl("^mature,", mirna_iso_luad$miRNA_region),
]

dim(mirna_mature)

length(unique(mirna_mature$miRNA_region))


mirna_mature_sum <- aggregate(
  read_count ~ miRNA_region + barcode,
  data = mirna_mature,
  FUN = sum
)

dim(mirna_mature_sum)


mirna_count_matrix <- stats::xtabs(
  read_count ~ miRNA_region + barcode,
  data = mirna_mature_sum
)

mirna_count_matrix <- as.matrix(
  mirna_count_matrix
)

dim(mirna_count_matrix)

"mature,MIMAT0004680" %in%
  rownames(mirna_count_matrix)

# ============================================================
# 10. DESeq2 normalization
# ============================================================

library(DESeq2)

# 519개 aliquot 정보
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

# miRNA 데이터는 0 count가 많기 때문에 poscounts 사용
dds_mirna <- estimateSizeFactors(
  dds_mirna,
  type = "poscounts"
)

# 각 sample의 보정계수 확인
summary(sizeFactors(dds_mirna))


mirna_norm_counts <- counts(
  dds_mirna,
  normalized = TRUE
)

dim(mirna_norm_counts)

summary(
  mirna_norm_counts["mature,MIMAT0004680", ]
)

# ============================================================
# 11. miR-130b-5p normalized expression table 만들기
# ============================================================

mir130b_sample <- data.frame(
  barcode = colnames(mirna_norm_counts),
  miR130b_5p_norm = as.numeric(
    mirna_norm_counts["mature,MIMAT0004680", ]
  )
)

# TCGA barcode에서 patient / sample ID 추출
mir130b_sample$patient_id <- substr(
  mir130b_sample$barcode,
  1, 12
)

mir130b_sample$sample_id <- substr(
  mir130b_sample$barcode,
  1, 16
)

# log2 normalized expression
mir130b_sample$miR130b_5p_log2norm <- log2(
  mir130b_sample$miR130b_5p_norm + 1
)

dim(mir130b_sample)

length(unique(mir130b_sample$patient_id))

length(unique(mir130b_sample$sample_id))

head(mir130b_sample)

# ============================================================
# 12. miR-130b-5p + KEAP1 status 결합
# ============================================================

mir130b_keap1 <- merge(
  mir130b_sample,
  keap1_status,
  by = "patient_id"
)

dim(mir130b_keap1)

table(mir130b_keap1$KEAP1_status)

length(unique(mir130b_keap1$patient_id))

length(unique(mir130b_keap1$sample_id))


# ============================================================
# 13. Overall survival 정보 결합
# ============================================================

mir130b_keap1_os <- merge(
  mir130b_keap1,
  clinical_surv[, c(
    "patient_id",
    "OS_days",
    "OS_event"
  )],
  by = "patient_id"
)

dim(mir130b_keap1_os)

table(
  mir130b_keap1_os$KEAP1_status
)


table(
  is.na(mir130b_keap1_os$OS_days),
  mir130b_keap1_os$KEAP1_status
)

table(
  is.na(mir130b_keap1_os$OS_event),
  mir130b_keap1_os$KEAP1_status
)


tcga_mir_primary <- mir130b_keap1_os[
  !is.na(mir130b_keap1_os$OS_days) &
    !is.na(mir130b_keap1_os$OS_event),
]

dim(tcga_mir_primary)

table(
  tcga_mir_primary$KEAP1_status
)

table(
  tcga_mir_primary$OS_days == 0,
  tcga_mir_primary$KEAP1_status
)


# ============================================================
# 14. 현재 후보 cohort의 sample / patient 구조 확인
# ============================================================

# aliquot 수
nrow(tcga_mir_primary)

# unique tumor sample 수
length(unique(tcga_mir_primary$sample_id))

# unique patient 수
length(unique(tcga_mir_primary$patient_id))

# KEAP1별 unique patient 수
tapply(
  tcga_mir_primary$patient_id,
  tcga_mir_primary$KEAP1_status,
  function(x) length(unique(x))
)

tcga_mir_os_positive <- tcga_mir_primary[
  tcga_mir_primary$OS_days > 0,
]

dim(tcga_mir_os_positive)

table(tcga_mir_os_positive$KEAP1_status)

length(unique(tcga_mir_os_positive$patient_id))


# ============================================================
# 15. TCGA-LUAD miRNA 전체 sample type 확인
# ============================================================

query_mirna_iso_all <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Isoform Expression Quantification",
  workflow.type = "BCGSC miRNA Profiling"
)

mirna_iso_results_all <- getResults(
  query_mirna_iso_all
)

dim(mirna_iso_results_all)

table(
  mirna_iso_results_all$sample_type,
  useNA = "ifany"
)

length(
  unique(mirna_iso_results_all$cases)
)

table(
  substr(mirna_iso_results_all$sample.submitter_id, 14, 15)
)

# ============================================================
# 16. Tumor miRNA sample의 환자 중복 구조 확인
# ============================================================

mirna_tumor_results <- mirna_iso_results_all[
  mirna_iso_results_all$sample_type %in%
    c("Primary Tumor", "Recurrent Tumor"),
]

# patient ID 만들기
mirna_tumor_results$patient_id <- substr(
  mirna_tumor_results$sample.submitter_id,
  1, 12
)

# tumor sample 수
nrow(mirna_tumor_results)

# unique patient 수
length(unique(mirna_tumor_results$patient_id))

# 환자당 tumor sample 수
sort(
  table(mirna_tumor_results$patient_id),
  decreasing = TRUE
)[1:20]

mirna_primary_results <- mirna_iso_results_all[
  mirna_iso_results_all$sample_type == "Primary Tumor",
]

mirna_primary_results$patient_id <- substr(
  mirna_primary_results$sample.submitter_id,
  1, 12
)

nrow(mirna_primary_results)

length(unique(mirna_primary_results$patient_id))

sort(
  table(mirna_primary_results$patient_id),
  decreasing = TRUE
)


# ============================================================
# 17. TCGA-LUAD RNA-seq metadata 확인
# ============================================================

query_rna <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

rna_results <- getResults(query_rna)

dim(rna_results)

table(
  rna_results$sample_type,
  useNA = "ifany"
)


rna_primary_results <- rna_results[
  rna_results$sample_type == "Primary Tumor",
]

nrow(rna_primary_results)


rna_primary_results$patient_id <- substr(
  rna_primary_results$sample.submitter_id,
  1, 12
)

length(
  unique(rna_primary_results$patient_id)
)


sort(
  table(rna_primary_results$patient_id),
  decreasing = TRUE
)[1:20]


# ============================================================
# 18. Multi-omics + OS patient intersection
# ============================================================

# miRNA-seq가 있는 환자
mirna_patients <- unique(
  mirna_primary_results$patient_id
)

# RNA-seq가 있는 환자
rna_patients <- unique(
  rna_primary_results$patient_id
)

# mutation profiling이 된 환자
mutation_patients <- unique(
  sequenced_patients
)

# OS 정보가 있는 환자
os_patients <- unique(
  clinical_surv$patient_id[
    !is.na(clinical_surv$OS_days) &
      !is.na(clinical_surv$OS_event)
  ]
)

length(mirna_patients)
length(rna_patients)
length(mutation_patients)
length(os_patients)


common_patients <- Reduce(
  intersect,
  list(
    mirna_patients,
    rna_patients,
    mutation_patients,
    os_patients
  )
)

length(common_patients)


common_keap1 <- keap1_status[
  keap1_status$patient_id %in% common_patients,
]

table(common_keap1$KEAP1_status)


# ============================================================
# 19. RNA-seq / miRNA-seq 동일 sample 교집합 확인
# ============================================================

# miRNA primary tumor sample ID
mirna_primary_results$sample_id <- substr(
  mirna_primary_results$sample.submitter_id,
  1, 16
)

# RNA primary tumor sample ID
rna_primary_results$sample_id <- substr(
  rna_primary_results$sample.submitter_id,
  1, 16
)

# unique sample 수
length(unique(mirna_primary_results$sample_id))

length(unique(rna_primary_results$sample_id))

common_sample_ids <- intersect(
  unique(mirna_primary_results$sample_id),
  unique(rna_primary_results$sample_id)
)

length(common_sample_ids)


common_sample_patients <- unique(
  substr(common_sample_ids, 1, 12)
)

common_sample_patients_complete <- Reduce(
  intersect,
  list(
    common_sample_patients,
    mutation_patients,
    os_patients
  )
)

length(common_sample_patients_complete)


common_sample_keap1 <- keap1_status[
  keap1_status$patient_id %in% common_sample_patients_complete,
]

table(common_sample_keap1$KEAP1_status)


common_complete_sample_ids <- common_sample_ids[
  substr(common_sample_ids, 1, 12) %in%
    common_sample_patients_complete
]

length(common_complete_sample_ids)

# ============================================================
# 20. OS_days > 0 적용한 multi-omics cohort 확인
# ============================================================

os_positive_patients <- unique(
  clinical_surv$patient_id[
    !is.na(clinical_surv$OS_days) &
      !is.na(clinical_surv$OS_event) &
      clinical_surv$OS_days > 0
  ]
)

common_patients_os_positive <- Reduce(
  intersect,
  list(
    mirna_patients,
    rna_patients,
    mutation_patients,
    os_positive_patients
  )
)

length(common_patients_os_positive)

common_keap1_os_positive <- keap1_status[
  keap1_status$patient_id %in% common_patients_os_positive,
]

table(common_keap1_os_positive$KEAP1_status)


# ============================================================
# 21. TCGA-LUAD RNA-seq 다운로드
# ============================================================

dir.create(
  "data/GDC_RNA",
  showWarnings = FALSE
)


exists("common_patients_os_positive")

length(common_patients_os_positive)

table(common_keap1_os_positive$KEAP1_status)


# ============================================================
# 21. TCGA-LUAD RNA-seq 다운로드
# ============================================================

dir.create(
  "data/GDC_RNA",
  showWarnings = FALSE
)

GDCdownload(
  query_rna,
  directory = "data/GDC_RNA",
  files.per.chunk = 50
)

length(
  list.files(
    "data/GDC_RNA",
    recursive = TRUE
  )
)

# ============================================================
# 22. RNA-seq 데이터 준비
# ============================================================

rna_luad <- GDCprepare(
  query_rna,
  directory = "data/GDC_RNA"
)

class(rna_luad)

dim(rna_luad)

assayNames(rna_luad)


# ============================================================
# 23. PIK3CB gene 확인
# ============================================================

colnames(rowData(rna_luad))

head(
  as.data.frame(rowData(rna_luad))
)

# ============================================================
# 24. PIK3CB gene row 찾기
# ============================================================

pik3cb_idx <- which(
  rowData(rna_luad)$gene_name == "PIK3CB"
)

pik3cb_idx

length(pik3cb_idx)

as.data.frame(
  rowData(rna_luad)[pik3cb_idx, ]
)

assay(
  rna_luad,
  "unstranded"
)[pik3cb_idx, 1:10]


# ============================================================
# 25. Primary Tumor PIK3CB expression 추출
# ============================================================

table(
  colData(rna_luad)$sample_type,
  useNA = "ifany"
)


rna_primary_idx <- which(
  colData(rna_luad)$sample_type == "Primary Tumor"
)

length(rna_primary_idx)




pik3cb_sample <- data.frame(
  barcode = colnames(rna_luad)[rna_primary_idx],
  PIK3CB_raw_count = as.numeric(
    assay(rna_luad, "unstranded")[
      pik3cb_idx,
      rna_primary_idx
    ]
  )
)

# patient / sample ID 생성
pik3cb_sample$patient_id <- substr(
  pik3cb_sample$barcode,
  1, 12
)

pik3cb_sample$sample_id <- substr(
  pik3cb_sample$barcode,
  1, 16
)

dim(pik3cb_sample)

length(unique(pik3cb_sample$patient_id))

length(unique(pik3cb_sample$sample_id))

head(pik3cb_sample)



# ============================================================
# 26. RNA-seq DESeq2 normalization
# ============================================================

# Primary Tumor의 전체 gene raw count matrix
rna_count_primary <- assay(
  rna_luad,
  "unstranded"
)[, rna_primary_idx]

dim(rna_count_primary)


# sample 정보
rna_coldata <- data.frame(
  row.names = colnames(rna_count_primary),
  group = rep("all", ncol(rna_count_primary))
)

# DESeq2 객체 생성
dds_rna <- DESeqDataSetFromMatrix(
  countData = rna_count_primary,
  colData = rna_coldata,
  design = ~ 1
)

# 표준 median-of-ratios 방식으로 size factor 계산
dds_rna <- estimateSizeFactors(dds_rna)

summary(sizeFactors(dds_rna))


rna_norm_counts <- counts(
  dds_rna,
  normalized = TRUE
)

dim(rna_norm_counts)

summary(
  rna_norm_counts[pik3cb_idx, ]
)

# ============================================================
# 27. PIK3CB normalized expression + miRNA sample matching
# ============================================================

pik3cb_norm_aliquot <- data.frame(
  barcode = colnames(rna_norm_counts),
  PIK3CB_norm = as.numeric(
    rna_norm_counts[pik3cb_idx, ]
  )
)

pik3cb_norm_aliquot$patient_id <- substr(
  pik3cb_norm_aliquot$barcode,
  1, 12
)

pik3cb_norm_aliquot$sample_id <- substr(
  pik3cb_norm_aliquot$barcode,
  1, 16
)

dim(pik3cb_norm_aliquot)


pik3cb_norm_sample <- aggregate(
  PIK3CB_norm ~ patient_id + sample_id,
  data = pik3cb_norm_aliquot,
  FUN = mean
)

pik3cb_norm_sample$PIK3CB_log2norm <- log2(
  pik3cb_norm_sample$PIK3CB_norm + 1
)

dim(pik3cb_norm_sample)


mir130b_norm_sample2 <- aggregate(
  miR130b_5p_norm ~ patient_id + sample_id,
  data = mir130b_sample,
  FUN = mean
)

mir130b_norm_sample2$miR130b_5p_log2norm <- log2(
  mir130b_norm_sample2$miR130b_5p_norm + 1
)

dim(mir130b_norm_sample2)


multi_expr_sample <- merge(
  mir130b_norm_sample2,
  pik3cb_norm_sample,
  by = c("patient_id", "sample_id")
)

dim(multi_expr_sample)

length(unique(multi_expr_sample$patient_id))

head(multi_expr_sample)

# ============================================================
# 28. RNA + miRNA + KEAP1 + OS 최종 cohort 확인
# ============================================================

# KEAP1 상태 붙이기
multi_expr_keap1 <- merge(
  multi_expr_sample,
  keap1_status,
  by = "patient_id"
)

dim(multi_expr_keap1)

table(multi_expr_keap1$KEAP1_status)


multi_expr_keap1_os <- merge(
  multi_expr_keap1,
  clinical_surv[, c(
    "patient_id",
    "OS_days",
    "OS_event"
  )],
  by = "patient_id"
)

dim(multi_expr_keap1_os)

table(multi_expr_keap1_os$KEAP1_status)


multi_expr_complete <- multi_expr_keap1_os[
  !is.na(multi_expr_keap1_os$OS_days) &
    !is.na(multi_expr_keap1_os$OS_event),
]

dim(multi_expr_complete)

table(multi_expr_complete$KEAP1_status)

length(unique(multi_expr_complete$patient_id))


table(
  multi_expr_complete$OS_days == 0,
  multi_expr_complete$KEAP1_status
)


# ============================================================
# 29. Current-GDC master cohort
# ============================================================

master_cohort <- multi_expr_complete[
  multi_expr_complete$OS_days > 0,
]

dim(master_cohort)

table(master_cohort$KEAP1_status)

length(unique(master_cohort$patient_id))


master_cohort <- merge(
  master_cohort,
  clinical_surv[, c(
    "patient_id",
    "RT_status"
  )],
  by = "patient_id",
  all.x = TRUE
)

table(
  master_cohort$RT_status,
  useNA = "ifany"
)

# ============================================================
# 30. miR-130b-5p expression: KEAP1 WT vs MUT
# ============================================================

master_cohort$KEAP1_status <- factor(
  master_cohort$KEAP1_status,
  levels = c("WT", "MUT")
)

aggregate(
  miR130b_5p_log2norm ~ KEAP1_status,
  data = master_cohort,
  FUN = function(x) c(
    n = length(x),
    median = median(x),
    mean = mean(x)
  )
)

wilcox_mir_master <- wilcox.test(
  miR130b_5p_log2norm ~ KEAP1_status,
  data = master_cohort,
  exact = FALSE
)

wilcox_mir_master



# ============================================================
# 31. Final violin plot
# ============================================================

library(ggplot2)

p_text <- paste0(
  "Wilcoxon p = ",
  signif(wilcox_mir_master$p.value, 4)
)

p_mir_final <- ggplot(
  master_cohort,
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
    y = max(master_cohort$miR130b_5p_log2norm, na.rm = TRUE) + 0.2,
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
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = "black"
    )
  )

p_mir_final


# ============================================================
# 32. KEAP1 MUT subgroup for miR-130b survival analysis
# ============================================================

mir_mut_surv <- master_cohort[
  master_cohort$KEAP1_status == "MUT",
]

dim(mir_mut_surv)

table(mir_mut_surv$OS_event)

summary(
  mir_mut_surv$miR130b_5p_log2norm
)

# KEAP1 MUT sample 수와 실제 환자 수 비교
nrow(mir_mut_surv)

length(unique(mir_mut_surv$patient_id))


# ============================================================
# 33. Maximally selected rank statistics
# ============================================================

library(survminer)

cut_mir_master <- surv_cutpoint(
  mir_mut_surv,
  time = "OS_days",
  event = "OS_event",
  variables = "miR130b_5p_log2norm",
  minprop = 0.1
)

cut_mir_master

cut_mir_master$cutpoint



# ============================================================
# 34. Low / High group 만들기
# ============================================================

mir_mut_cat <- surv_categorize(
  cut_mir_master
)

table(
  mir_mut_cat$miR130b_5p_log2norm
)


table(
  mir_mut_cat$miR130b_5p_log2norm,
  mir_mut_cat$OS_event
)


# ============================================================
# 35. Kaplan-Meier survival analysis
# ============================================================

library(survival)

mir_mut_cat$miR_group <- factor(
  mir_mut_cat$miR130b_5p_log2norm,
  levels = c("low", "high"),
  labels = c("Low", "High")
)

fit_mir_mut <- survfit(
  Surv(OS_days, OS_event) ~ miR_group,
  data = mir_mut_cat
)

fit_mir_mut


logrank_mir_mut <- survdiff(
  Surv(OS_days, OS_event) ~ miR_group,
  data = mir_mut_cat
)

logrank_mir_mut

p_mir_mut <- 1 - pchisq(
  logrank_mir_mut$chisq,
  df = 1
)

p_mir_mut


# ============================================================
# 36. KEAP1 MUT miR-130b-5p Kaplan-Meier plot
# ============================================================

library(survminer)

km_mir_mut <- ggsurvplot(
  fit_mir_mut,
  data = mir_mut_cat,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  censor = TRUE,
  
  palette = c(
    "Low" = "#00A6D6",
    "High" = "#E84A5F"
  ),
  
  xlab = "Days",
  ylab = "Survival probability",
  
  legend.title = "miR-130b-5p",
  legend.labs = c("Low", "High"),
  
  title = "miR-130b-5p in KEAP1 MUT",
  
  risk.table.height = 0.25,
  ggtheme = theme_bw()
)

# ============================================================
# miR-130b-5p median split
# ============================================================

mir_cutoff_median <- median(
  mir_mut_surv$miR130b_5p_log2norm,
  na.rm = TRUE
)

mir_cutoff_median


mir_mut_surv$miR_group_median <- ifelse(
  mir_mut_surv$miR130b_5p_log2norm <= mir_cutoff_median,
  "Low",
  "High"
)

mir_mut_surv$miR_group_median <- factor(
  mir_mut_surv$miR_group_median,
  levels = c("Low", "High")
)

table(mir_mut_surv$miR_group_median)

# ============================================================
# Median split Kaplan-Meier analysis
# ============================================================

fit_mir_median <- survfit(
  Surv(OS_days, OS_event) ~ miR_group_median,
  data = mir_mut_surv
)

fit_mir_median

logrank_mir_median <- survdiff(
  Surv(OS_days, OS_event) ~ miR_group_median,
  data = mir_mut_surv
)

logrank_mir_median

p_mir_median <- 1 - pchisq(
  logrank_mir_median$chisq,
  df = 1
)

p_mir_median

# ============================================================
# 37. PIK3CB survival analysis in KEAP1 WT
# ============================================================

pik3cb_wt_surv <- master_cohort[
  master_cohort$KEAP1_status == "WT",
]

dim(pik3cb_wt_surv)

table(pik3cb_wt_surv$OS_event)

summary(
  pik3cb_wt_surv$PIK3CB_log2norm
)
km_mir_mut


# ============================================================
# 38. PIK3CB cutoff in KEAP1 WT
# ============================================================

cut_pik3cb_wt <- surv_cutpoint(
  pik3cb_wt_surv,
  time = "OS_days",
  event = "OS_event",
  variables = "PIK3CB_log2norm",
  minprop = 0.1
)

cut_pik3cb_wt

cut_pik3cb_wt$cutpoint

pik3cb_wt_cat <- surv_categorize(
  cut_pik3cb_wt
)

table(
  pik3cb_wt_cat$PIK3CB_log2norm
)

table(
  pik3cb_wt_cat$PIK3CB_log2norm,
  pik3cb_wt_cat$OS_event
)


# ============================================================
# 39. Kaplan-Meier analysis: PIK3CB in KEAP1 WT
# ============================================================

pik3cb_wt_cat$PIK3CB_group <- factor(
  pik3cb_wt_cat$PIK3CB_log2norm,
  levels = c("low", "high"),
  labels = c("Low", "High")
)

fit_pik3cb_wt <- survfit(
  Surv(OS_days, OS_event) ~ PIK3CB_group,
  data = pik3cb_wt_cat
)

fit_pik3cb_wt

logrank_pik3cb_wt <- survdiff(
  Surv(OS_days, OS_event) ~ PIK3CB_group,
  data = pik3cb_wt_cat
)

logrank_pik3cb_wt

p_pik3cb_wt <- 1 - pchisq(
  logrank_pik3cb_wt$chisq,
  df = 1
)

p_pik3cb_wt

# ============================================================
# 40. Kaplan-Meier plot: PIK3CB in KEAP1 WT
# ============================================================

km_pik3cb_wt <- ggsurvplot(
  fit_pik3cb_wt,
  data = pik3cb_wt_cat,
  
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  censor = TRUE,
  
  palette = c(
    "Low" = "#00A6D6",
    "High" = "#E84A5F"
  ),
  
  xlab = "Days",
  ylab = "Survival probability",
  
  legend.title = "PIK3CB",
  legend.labs = c("Low", "High"),
  
  title = "PIK3CB in KEAP1 WT",
  
  risk.table.height = 0.25,
  ggtheme = theme_bw()
)

# ============================================================
# 41. PIK3CB survival analysis in KEAP1 MUT
# ============================================================

pik3cb_mut_surv <- master_cohort[
  master_cohort$KEAP1_status == "MUT",
]

dim(pik3cb_mut_surv)

table(pik3cb_mut_surv$OS_event)

summary(
  pik3cb_mut_surv$PIK3CB_log2norm
)

# ============================================================
# 42. PIK3CB cutoff in KEAP1 MUT
# ============================================================

cut_pik3cb_mut <- surv_cutpoint(
  pik3cb_mut_surv,
  time = "OS_days",
  event = "OS_event",
  variables = "PIK3CB_log2norm",
  minprop = 0.1
)

cut_pik3cb_mut

cut_pik3cb_mut$cutpoint

pik3cb_mut_cat <- surv_categorize(
  cut_pik3cb_mut
)

table(
  pik3cb_mut_cat$PIK3CB_log2norm
)

table(
  pik3cb_mut_cat$PIK3CB_log2norm,
  pik3cb_mut_cat$OS_event
)
km_pik3cb_wt

# ============================================================
# 43. Kaplan-Meier analysis: PIK3CB in KEAP1 MUT
# ============================================================

pik3cb_mut_cat$PIK3CB_group <- factor(
  pik3cb_mut_cat$PIK3CB_log2norm,
  levels = c("low", "high"),
  labels = c("Low", "High")
)

fit_pik3cb_mut <- survfit(
  Surv(OS_days, OS_event) ~ PIK3CB_group,
  data = pik3cb_mut_cat
)

fit_pik3cb_mut

logrank_pik3cb_mut <- survdiff(
  Surv(OS_days, OS_event) ~ PIK3CB_group,
  data = pik3cb_mut_cat
)

logrank_pik3cb_mut

p_pik3cb_mut <- 1 - pchisq(
  logrank_pik3cb_mut$chisq,
  df = 1
)

p_pik3cb_mut

# ============================================================
# 44. Kaplan-Meier plot: PIK3CB in KEAP1 MUT
# ============================================================

km_pik3cb_mut <- ggsurvplot(
  fit_pik3cb_mut,
  data = pik3cb_mut_cat,
  
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  censor = TRUE,
  
  palette = c(
    "Low" = "#00A6D6",
    "High" = "#E84A5F"
  ),
  
  xlab = "Days",
  ylab = "Survival probability",
  
  legend.title = "PIK3CB",
  legend.labs = c("Low", "High"),
  
  title = "PIK3CB in KEAP1 MUT",
  
  risk.table.height = 0.25,
  ggtheme = theme_bw()
)

km_pik3cb_mut

# ============================================================
# 45. KEAP1 WT + Radiotherapy subgroup
# ============================================================

wt_rt <- master_cohort[
  master_cohort$KEAP1_status == "WT" &
    master_cohort$RT_status == "RT",
]

dim(wt_rt)

table(wt_rt$OS_event)

length(unique(wt_rt$patient_id))

summary(wt_rt$miR130b_5p_log2norm)

summary(wt_rt$PIK3CB_log2norm)


# ============================================================
# 45. KEAP1 WT + Radiotherapy subgroup (수정)
# ============================================================

wt_rt <- master_cohort[
  !is.na(master_cohort$RT_status) &
    master_cohort$KEAP1_status == "WT" &
    master_cohort$RT_status == "RT",
]

dim(wt_rt)

table(wt_rt$OS_event)

length(unique(wt_rt$patient_id))

summary(wt_rt$miR130b_5p_log2norm)

summary(wt_rt$PIK3CB_log2norm)


# ============================================================
# 46. PIK3CB + miR-130b-5p composite risk score
# ============================================================

# PIK3CB Z-score
wt_rt$PIK3CB_Z <- as.numeric(
  scale(wt_rt$PIK3CB_log2norm)
)

# miR-130b-5p Z-score
wt_rt$miR130b_Z <- as.numeric(
  scale(wt_rt$miR130b_5p_log2norm)
)

# Composite risk score
wt_rt$composite_score <- (
  wt_rt$PIK3CB_Z -
    wt_rt$miR130b_Z
)

summary(wt_rt$PIK3CB_Z)

summary(wt_rt$miR130b_Z)

summary(wt_rt$composite_score)

head(
  wt_rt[, c(
    "patient_id",
    "PIK3CB_log2norm",
    "miR130b_5p_log2norm",
    "PIK3CB_Z",
    "miR130b_Z",
    "composite_score"
  )]
)

# ============================================================
# 47. Composite risk cutoff
# ============================================================

cut_composite <- surv_cutpoint(
  wt_rt,
  time = "OS_days",
  event = "OS_event",
  variables = "composite_score",
  minprop = 0.1
)

cut_composite

cut_composite$cutpoint

composite_cat <- surv_categorize(
  cut_composite
)

table(
  composite_cat$composite_score
)

table(
  composite_cat$composite_score,
  composite_cat$OS_event
)


# ============================================================
# 48. Kaplan-Meier analysis: Composite risk
# ============================================================

composite_cat$risk_group <- factor(
  composite_cat$composite_score,
  levels = c("low", "high"),
  labels = c("Low risk", "High risk")
)

fit_composite <- survfit(
  Surv(OS_days, OS_event) ~ risk_group,
  data = composite_cat
)

fit_composite


logrank_composite <- survdiff(
  Surv(OS_days, OS_event) ~ risk_group,
  data = composite_cat
)

logrank_composite


p_composite <- 1 - pchisq(
  logrank_composite$chisq,
  df = 1
)

p_composite


# ============================================================
# 49. Kaplan-Meier plot: Composite risk
# ============================================================

km_composite <- ggsurvplot(
  fit_composite,
  data = composite_cat,
  
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  censor = TRUE,
  
  palette = c(
    "Low risk" = "#00A6D6",
    "High risk" = "#E84A5F"
  ),
  
  xlab = "Days",
  ylab = "Survival probability",
  
  legend.title = "Composite risk",
  legend.labs = c("Low risk", "High risk"),
  
  title = "PIK3CB / miR-130b-5p composite risk in KEAP1 WT + RT",
  
  risk.table.height = 0.25,
  ggtheme = theme_bw()
)

km_composite

# ============================================================
# 50. ssGSEA 준비 - RNA-seq 데이터 확인
# ============================================================

# 사용 가능한 expression assay 확인
assayNames(rna_luad)

# gene symbol 확인
head(rowData(rna_luad)$gene_name)

# gene symbol이 없는 유전자 수
sum(
  is.na(rowData(rna_luad)$gene_name) |
    rowData(rna_luad)$gene_name == ""
)

# 같은 gene symbol이 중복되어 있는지 확인
sum(
  duplicated(rowData(rna_luad)$gene_name)
)

# Primary Tumor 개수 다시 확인
sum(colData(rna_luad)$sample_type == "Primary Tumor")


# ============================================================
# 51. ssGSEA용 expression matrix 만들기
# ============================================================

# Primary Tumor만 선택
rna_primary_idx <- which(
  colData(rna_luad)$sample_type == "Primary Tumor"
)

# TPM expression
expr_tpm <- assay(
  rna_luad,
  "tpm_unstrand"
)[, rna_primary_idx]

# log2(TPM + 1)
expr_log2tpm <- log2(expr_tpm + 1)

# gene symbol
gene_symbol <- rowData(rna_luad)$gene_name

# 각 gene row의 평균 발현
gene_mean <- rowMeans(
  expr_log2tpm,
  na.rm = TRUE
)

# 같은 gene symbol이 여러 개라면
# 평균 발현이 가장 높은 row를 먼저 오도록 정렬
ord <- order(
  gene_symbol,
  -gene_mean
)

# gene symbol당 하나만 유지
keep <- !duplicated(
  gene_symbol[ord]
)

expr_ssgsea <- expr_log2tpm[
  ord[keep],
  ,
  drop = FALSE
]

# row name을 gene symbol로 변경
rownames(expr_ssgsea) <- gene_symbol[
  ord[keep]
]

# 확인
dim(expr_ssgsea)

sum(
  duplicated(rownames(expr_ssgsea))
)

head(rownames(expr_ssgsea))

summary(
  as.numeric(expr_ssgsea["PIK3CB", ])
)

# ============================================================
# 52. PI3K pathway gene set 확인
# ============================================================

if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}

library(msigdbr)


msig_human <- msigdbr(
  species = "Homo sapiens"
)

dim(msig_human)

colnames(msig_human)


pi3k_sets <- unique(
  msig_human$gs_name[
    grepl(
      "PI3K",
      msig_human$gs_name,
      ignore.case = TRUE
    )
  ]
)

pi3k_sets

# ============================================================
# 53. Hallmark PI3K-AKT-mTOR gene set 추출
# ============================================================

pi3k_genes <- unique(
  msig_human$gene_symbol[
    msig_human$gs_name ==
      "HALLMARK_PI3K_AKT_MTOR_SIGNALING"
  ]
)

# gene set 크기
length(pi3k_genes)

# 우리 expression matrix에 존재하는 gene 수
sum(
  pi3k_genes %in% rownames(expr_ssgsea)
)

# 실제 겹치는 gene들
pi3k_genes_present <- intersect(
  pi3k_genes,
  rownames(expr_ssgsea)
)

head(pi3k_genes_present, 20)

# ============================================================
# 54. ssGSEA score 계산
# ============================================================

if (!requireNamespace("GSVA", quietly = TRUE)) {
  BiocManager::install("GSVA")
}

library(GSVA)

pi3k_gene_set <- list(
  HALLMARK_PI3K_AKT_MTOR_SIGNALING = pi3k_genes_present
)

ssgsea_param <- ssgseaParam(
  expr_ssgsea,
  pi3k_gene_set
)

pi3k_ssgsea <- gsva(
  ssgsea_param
)


# ============================================================
# 54-1. GSVA 불러오기
# ============================================================

library(GSVA)

packageVersion("GSVA")


ssgsea_param <- GSVA::ssgseaParam(
  expr_ssgsea,
  pi3k_gene_set
)

pi3k_ssgsea <- GSVA::gsva(
  ssgsea_param
)


dim(pi3k_ssgsea)

summary(
  as.numeric(pi3k_ssgsea[1, ])
)
# ============================================================
# 55. PIK3CB와 PI3K ssGSEA score 준비
# ============================================================

# PIK3CB가 pathway gene set 안에 포함되어 있는지 확인
"PIK3CB" %in% pi3k_genes_present


# ============================================================
# 55-1. PIK3CB expression + PI3K ssGSEA score 합치기
# ============================================================

pi3k_score_df <- data.frame(
  barcode = colnames(pi3k_ssgsea),
  
  PI3K_ssGSEA = as.numeric(
    pi3k_ssgsea[1, ]
  ),
  
  PIK3CB_log2TPM = as.numeric(
    expr_ssgsea[
      "PIK3CB",
      colnames(pi3k_ssgsea)
    ]
  )
)

pi3k_score_df$patient_id <- substr(
  pi3k_score_df$barcode,
  1, 12
)

pi3k_score_df$sample_id <- substr(
  pi3k_score_df$barcode,
  1, 16
)

dim(pi3k_score_df)

head(pi3k_score_df)


# ============================================================
# 56. sample 단위로 정리 + KEAP1 status 추가
# ============================================================

pi3k_score_sample <- aggregate(
  cbind(
    PI3K_ssGSEA,
    PIK3CB_log2TPM
  ) ~ patient_id + sample_id,
  data = pi3k_score_df,
  FUN = mean
)

# KEAP1 WT / MUT 붙이기
pi3k_score_sample <- merge(
  pi3k_score_sample,
  keap1_status,
  by = "patient_id"
)

dim(pi3k_score_sample)

table(pi3k_score_sample$KEAP1_status)

head(pi3k_score_sample)


# 전체 PIK3CB vs PI3K pathway correlation
cor_all <- cor.test(
  pi3k_score_sample$PIK3CB_log2TPM,
  pi3k_score_sample$PI3K_ssGSEA,
  method = "spearman",
  exact = FALSE
)

cor_all

# ============================================================
# 57. KEAP1 WT / MUT 각각 상관분석
# ============================================================

# KEAP1 WT
cor_wt <- cor.test(
  pi3k_score_sample$PIK3CB_log2TPM[
    pi3k_score_sample$KEAP1_status == "WT"
  ],
  pi3k_score_sample$PI3K_ssGSEA[
    pi3k_score_sample$KEAP1_status == "WT"
  ],
  method = "spearman",
  exact = FALSE
)

cor_wt


# KEAP1 MUT
cor_mut <- cor.test(
  pi3k_score_sample$PIK3CB_log2TPM[
    pi3k_score_sample$KEAP1_status == "MUT"
  ],
  pi3k_score_sample$PI3K_ssGSEA[
    pi3k_score_sample$KEAP1_status == "MUT"
  ],
  method = "spearman",
  exact = FALSE
)

cor_mut


# ============================================================
# 58. PIK3CB vs PI3K ssGSEA scatter plot
# ============================================================

library(ggplot2)

pi3k_score_sample$KEAP1_status <- factor(
  pi3k_score_sample$KEAP1_status,
  levels = c("WT", "MUT")
)

corr_label <- data.frame(
  KEAP1_status = factor(
    c("WT", "MUT"),
    levels = c("WT", "MUT")
  ),
  label = c(
    "Spearman rho = 0.027\np = 0.577",
    "Spearman rho = 0.392\np = 0.00011"
  )
)

p_pi3k_corr <- ggplot(
  pi3k_score_sample,
  aes(
    x = PIK3CB_log2TPM,
    y = PI3K_ssGSEA
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  facet_wrap(
    ~ KEAP1_status
  ) +
  geom_text(
    data = corr_label,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.2,
    size = 4
  ) +
  labs(
    x = "PIK3CB expression (log2 TPM + 1)",
    y = "PI3K-AKT-mTOR ssGSEA score",
    title = "PIK3CB expression vs PI3K pathway activity"
  ) +
  theme_bw()



# ============================================================
# 58-1. Scatter plot 수정
# ============================================================

p_pi3k_corr <- ggplot(
  pi3k_score_sample,
  aes(
    x = PIK3CB_log2TPM,
    y = PI3K_ssGSEA
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.8
  ) +
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  facet_wrap(
    ~ KEAP1_status
  ) +
  geom_text(
    data = corr_label,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.2,
    size = 4
  ) +
  labs(
    x = "PIK3CB expression (log2 TPM + 1)",
    y = "PI3K-AKT-mTOR ssGSEA score",
    title = "PIK3CB expression vs PI3K pathway activity"
  ) +
  theme_bw()

p_pi3k_corr
p_pi3k_corr


# ============================================================
# 58-2. 최종 권장 scatter plot
# ============================================================

p_pi3k_corr_final <- ggplot(
  pi3k_score_sample,
  aes(
    x = PIK3CB_log2TPM,
    y = PI3K_ssGSEA
  )
) +
  geom_point(
    alpha = 0.6,
    size = 1.8
  ) +
  facet_wrap(
    ~ KEAP1_status
  ) +
  geom_text(
    data = corr_label,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.2,
    size = 4
  ) +
  labs(
    x = "PIK3CB expression (log2 TPM + 1)",
    y = "PI3K-AKT-mTOR ssGSEA score",
    title = "PIK3CB expression vs PI3K pathway activity"
  ) +
  theme_bw()

p_pi3k_corr_final


# ============================================================
# 59. WT vs MUT Spearman correlation 차이 검정
# ============================================================

# 환자당 하나의 값으로 정리
pi3k_patient <- aggregate(
  cbind(
    PI3K_ssGSEA,
    PIK3CB_log2TPM
  ) ~ patient_id + KEAP1_status,
  data = pi3k_score_sample,
  FUN = mean
)

# 환자 수 확인
dim(pi3k_patient)

table(pi3k_patient$KEAP1_status)



# WT correlation
rho_wt_patient <- cor(
  pi3k_patient$PIK3CB_log2TPM[
    pi3k_patient$KEAP1_status == "WT"
  ],
  pi3k_patient$PI3K_ssGSEA[
    pi3k_patient$KEAP1_status == "WT"
  ],
  method = "spearman"
)

# MUT correlation
rho_mut_patient <- cor(
  pi3k_patient$PIK3CB_log2TPM[
    pi3k_patient$KEAP1_status == "MUT"
  ],
  pi3k_patient$PI3K_ssGSEA[
    pi3k_patient$KEAP1_status == "MUT"
  ],
  method = "spearman"
)

rho_wt_patient
rho_mut_patient

# 관찰된 두 correlation의 차이
observed_diff <- rho_mut_patient - rho_wt_patient

observed_diff



# ============================================================
# 60. Permutation test
# ============================================================

set.seed(123)

B <- 10000

perm_diff <- replicate(B, {
  
  perm_status <- sample(
    pi3k_patient$KEAP1_status
  )
  
  rho_wt_perm <- cor(
    pi3k_patient$PIK3CB_log2TPM[
      perm_status == "WT"
    ],
    pi3k_patient$PI3K_ssGSEA[
      perm_status == "WT"
    ],
    method = "spearman"
  )
  
  rho_mut_perm <- cor(
    pi3k_patient$PIK3CB_log2TPM[
      perm_status == "MUT"
    ],
    pi3k_patient$PI3K_ssGSEA[
      perm_status == "MUT"
    ],
    method = "spearman"
  )
  
  rho_mut_perm - rho_wt_perm
})

p_diff <- (
  sum(abs(perm_diff) >= abs(observed_diff)) + 1
) / (B + 1)

p_diff


# ============================================================
# 61. Patient-level WT / MUT correlation test
# ============================================================

cor_wt_patient <- cor.test(
  pi3k_patient$PIK3CB_log2TPM[
    pi3k_patient$KEAP1_status == "WT"
  ],
  pi3k_patient$PI3K_ssGSEA[
    pi3k_patient$KEAP1_status == "WT"
  ],
  method = "spearman",
  exact = FALSE
)

cor_mut_patient <- cor.test(
  pi3k_patient$PIK3CB_log2TPM[
    pi3k_patient$KEAP1_status == "MUT"
  ],
  pi3k_patient$PI3K_ssGSEA[
    pi3k_patient$KEAP1_status == "MUT"
  ],
  method = "spearman",
  exact = FALSE
)

cor_wt_patient
cor_mut_patient


# ============================================================
# 62. Final scatter plot
# ============================================================

library(ggplot2)

# factor 순서 고정
pi3k_patient$KEAP1_status <- factor(
  pi3k_patient$KEAP1_status,
  levels = c("WT", "MUT")
)

# 패널에 넣을 라벨 만들기
corr_label_patient <- data.frame(
  KEAP1_status = factor(
    c("WT", "MUT"),
    levels = c("WT", "MUT")
  ),
  label = c(
    paste0(
      "n = ", sum(pi3k_patient$KEAP1_status == "WT"),
      "\nSpearman rho = ", round(cor_wt_patient$estimate, 3),
      "\np = ", signif(cor_wt_patient$p.value, 3)
    ),
    paste0(
      "n = ", sum(pi3k_patient$KEAP1_status == "MUT"),
      "\nSpearman rho = ", round(cor_mut_patient$estimate, 3),
      "\np = ", signif(cor_mut_patient$p.value, 3)
    )
  )
)

p_pi3k_final <- ggplot(
  pi3k_patient,
  aes(
    x = PIK3CB_log2TPM,
    y = PI3K_ssGSEA,
    color = KEAP1_status
  )
) +
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8
  ) +
  facet_wrap(
    ~ KEAP1_status,
    nrow = 1
  ) +
  geom_text(
    data = corr_label_patient,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.1,
    vjust = 1.2,
    size = 4
  ) +
  scale_color_manual(
    values = c("WT" = "#2A9D8F", "MUT" = "#F4A261")
  ) +
  labs(
    x = "PIK3CB expression (log2 TPM + 1)",
    y = "PI3K-AKT-mTOR ssGSEA score",
    title = "Association between PIK3CB expression and PI3K pathway activity",
    subtitle = paste0(
      "Difference in Spearman correlation between WT and MUT: permutation p = ",
      signif(p_diff, 3)
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

p_pi3k_final


ggsave(
  filename = "figures/PIK3CB_PI3K_ssGSEA_WT_MUT_final.png",
  plot = p_pi3k_final,
  width = 10,
  height = 5,
  dpi = 300
)


# ============================================================
# 63. PI3K pathway heatmap용 patient-level matrix 만들기
# ============================================================

# 105개 PI3K pathway gene만 추출
# sample x gene 형태로 전환
expr_pi3k_sample <- t(
  expr_ssgsea[
    pi3k_genes_present,
    ,
    drop = FALSE
  ]
)

# 각 RNA sample의 patient ID
patient_id_rna <- substr(
  rownames(expr_pi3k_sample),
  1, 12
)

# 같은 환자에게 여러 sample이 있으면 합산
expr_pi3k_sum <- rowsum(
  expr_pi3k_sample,
  group = patient_id_rna,
  reorder = TRUE
)

# 환자별 sample 개수
n_per_patient <- table(patient_id_rna)

# 환자별 평균 expression
expr_pi3k_patient <- expr_pi3k_sum /
  as.numeric(
    n_per_patient[rownames(expr_pi3k_sum)]
  )

dim(expr_pi3k_patient)



# ============================================================
# 64. KEAP1 정보와 환자 맞추기
# ============================================================

common_patients <- intersect(
  rownames(expr_pi3k_patient),
  pi3k_patient$patient_id
)

# heatmap expression:
# gene x patient
expr_heatmap <- t(
  expr_pi3k_patient[
    common_patients,
    ,
    drop = FALSE
  ]
)

# 환자 annotation
heatmap_meta <- pi3k_patient[
  match(
    common_patients,
    pi3k_patient$patient_id
  ),
  c(
    "patient_id",
    "KEAP1_status",
    "PIK3CB_log2TPM",
    "PI3K_ssGSEA"
  )
]

rownames(heatmap_meta) <- heatmap_meta$patient_id

# 확인
dim(expr_heatmap)

table(heatmap_meta$KEAP1_status)

all(
  colnames(expr_heatmap) ==
    rownames(heatmap_meta)
)


# ============================================================
# 65. Heatmap용 Z-score + 환자 순서 정렬
# ============================================================

# 각 gene별 Z-score
# gene마다 환자들 사이의 상대적인 발현 높고/낮음을 나타냄
expr_heatmap_z <- t(
  scale(
    t(expr_heatmap)
  )
)

# NA가 생겼는지 확인
sum(is.na(expr_heatmap_z))

# 극단적인 값 때문에 색이 너무 찌그러지는 것을 방지
expr_heatmap_z[
  expr_heatmap_z > 2
] <- 2

expr_heatmap_z[
  expr_heatmap_z < -2
] <- -2

# WT -> MUT 순서,
# 각 그룹에서는 PI3K ssGSEA 낮음 -> 높음 순서
patient_order <- order(
  heatmap_meta$KEAP1_status,
  heatmap_meta$PI3K_ssGSEA
)

# expression과 annotation 모두 같은 순서로 정렬
expr_heatmap_z <- expr_heatmap_z[
  ,
  patient_order,
  drop = FALSE
]

heatmap_meta <- heatmap_meta[
  patient_order,
  ,
  drop = FALSE
]

# 확인
dim(expr_heatmap_z)

table(
  heatmap_meta$KEAP1_status
)

head(
  heatmap_meta[
    ,
    c(
      "KEAP1_status",
      "PIK3CB_log2TPM",
      "PI3K_ssGSEA"
    )
  ]
)

# ============================================================
# 66. PI3K pathway heatmap
# ============================================================

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)

# 위쪽 annotation
annotation_col <- data.frame(
  KEAP1_status = heatmap_meta$KEAP1_status,
  PIK3CB = heatmap_meta$PIK3CB_log2TPM,
  PI3K_score = heatmap_meta$PI3K_ssGSEA
)

rownames(annotation_col) <- rownames(heatmap_meta)

# WT가 끝나는 위치
wt_end <- sum(
  heatmap_meta$KEAP1_status == "WT"
)

wt_end


# Heatmap
pheatmap(
  expr_heatmap_z,
  
  # gene은 비슷한 패턴끼리 묶기
  cluster_rows = TRUE,
  
  # 환자 순서는 우리가 정한 순서 유지
  cluster_cols = FALSE,
  
  # 507명 이름은 표시하지 않음
  show_colnames = FALSE,
  
  # gene 이름은 표시
  show_rownames = TRUE,
  
  # 위쪽 annotation
  annotation_col = annotation_col,
  
  # WT와 MUT 사이에 간격
  gaps_col = wt_end,
  
  # Z-score 색 범위
  color = colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(100),
  
  breaks = seq(
    -2, 2,
    length.out = 101
  ),
  
  border_color = NA,
  
  fontsize_row = 6,
  
  main = "PI3K-AKT-mTOR pathway expression in TCGA-LUAD"
)


# ============================================================
# 67. KEAP1 WT vs MUT PI3K pathway activity 비교
# ============================================================

wilcox_pi3k <- wilcox.test(
  PI3K_ssGSEA ~ KEAP1_status,
  data = pi3k_patient,
  exact = FALSE
)

wilcox_pi3k



aggregate(
  PI3K_ssGSEA ~ KEAP1_status,
  data = pi3k_patient,
  FUN = function(x) c(
    n = length(x),
    median = median(x),
    mean = mean(x)
  )
)

# ============================================================
# 68. PI3K pathway에서 변동성이 큰 25개 gene 선택
# ============================================================

gene_sd <- apply(
  expr_heatmap_z,
  1,
  sd
)

top25_genes <- names(
  sort(
    gene_sd,
    decreasing = TRUE
  )
)[1:25]

top25_genes


expr_heatmap_top25 <- expr_heatmap_z[
  top25_genes,
  ,
  drop = FALSE
]

dim(expr_heatmap_top25)


# ============================================================
# 69. Top 25 PI3K pathway heatmap
# ============================================================

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}

library(pheatmap)


# ============================================================
# 70. Heatmap annotation 만들기
# ============================================================

annotation_col <- data.frame(
  KEAP1_status = heatmap_meta$KEAP1_status,
  PIK3CB_log2TPM = heatmap_meta$PIK3CB_log2TPM,
  PI3K_ssGSEA = heatmap_meta$PI3K_ssGSEA
)

rownames(annotation_col) <- rownames(heatmap_meta)

head(annotation_col)


all(colnames(expr_heatmap_top25) == rownames(annotation_col))

# ============================================================
# 71. annotation 색 설정
# ============================================================

ann_colors <- list(
  KEAP1_status = c(
    WT = "#4DB6AC",   # 청록색
    MUT = "#F4A259"   # 오렌지색
  )
)


# ============================================================
# 72. 최종 heatmap 그리기
# ============================================================

pheatmap(
  expr_heatmap_top25,
  
  color = colorRampPalette(
    c("#3B4CC0", "white", "#B40426")
  )(100),
  
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  
  show_colnames = FALSE,
  show_rownames = TRUE,
  
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  
  scale = "none",
  
  fontsize_row = 10,
  fontsize_col = 6,
  
  border_color = NA,
  
  main = "Top variable PI3K-AKT-mTOR pathway genes in TCGA-LUAD"
)

# ============================================================
# 73. Heatmap 저장
# ============================================================

pdf("figures/heatmap_top25_pi3k_pathway.pdf", width = 12, height = 8)

pheatmap(
  expr_heatmap_top25,
  color = colorRampPalette(
    c("#3B4CC0", "white", "#B40426")
  )(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  scale = "none",
  fontsize_row = 10,
  fontsize_col = 6,
  border_color = NA,
  main = "Top variable PI3K-AKT-mTOR pathway genes in TCGA-LUAD"
)

dev.off()

# ============================================================
# 74. Final polished PI3K heatmap
# ============================================================

library(pheatmap)

# WT 환자가 끝나는 위치
wt_end <- sum(
  heatmap_meta$KEAP1_status == "WT"
)

# annotation
annotation_col_final <- data.frame(
  KEAP1 = heatmap_meta$KEAP1_status,
  PIK3CB = heatmap_meta$PIK3CB_log2TPM,
  PI3K_score = heatmap_meta$PI3K_ssGSEA
)

rownames(annotation_col_final) <- rownames(heatmap_meta)

# KEAP1 색
annotation_colors_final <- list(
  KEAP1 = c(
    WT = "#4DB6AC",
    MUT = "#F4A261"
  )
)

# 최종 heatmap
pheatmap(
  expr_heatmap_top25,
  
  # gene은 유사한 발현 패턴끼리 clustering
  cluster_rows = TRUE,
  
  # 환자 순서는 WT -> MUT,
  # 각 군 안에서는 PI3K score 순서 유지
  cluster_cols = FALSE,
  
  # WT / MUT 사이 간격
  gaps_col = wt_end,
  
  show_colnames = FALSE,
  show_rownames = TRUE,
  
  annotation_col = annotation_col_final,
  annotation_colors = annotation_colors_final,
  
  color = colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(100),
  
  breaks = seq(
    -2, 2,
    length.out = 101
  ),
  
  border_color = NA,
  
  fontsize = 10,
  fontsize_row = 10,
  
  main = "PI3K-AKT-mTOR pathway expression in TCGA-LUAD"
)


# annotation 만들기
annotation_col_final <- data.frame(
  KEAP1 = heatmap_meta$KEAP1_status,
  PIK3CB = heatmap_meta$PIK3CB_log2TPM,
  PI3K_score = heatmap_meta$PI3K_ssGSEA
)

rownames(annotation_col_final) <- rownames(heatmap_meta)

# annotation 색
annotation_colors_final <- list(
  KEAP1 = c(
    WT = "#4DB6AC",
    MUT = "#F4A261"
  )
)

wt_end <- sum(
  heatmap_meta$KEAP1_status == "WT"
)

wt_end



pheatmap(
  expr_heatmap_top25,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = wt_end,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annotation_col_final,
  annotation_colors = annotation_colors_final,
  color = colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(100),
  breaks = seq(-2, 2, length.out = 101),
  border_color = NA,
  fontsize = 10,
  fontsize_row = 10,
  main = "PI3K-AKT-mTOR pathway expression in TCGA-LUAD"
)

dev.cur()
dev.list()


dev.off()


pheatmap(
  expr_heatmap_top25,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = wt_end,
  show_colnames = FALSE,
  show_rownames = TRUE,
  annotation_col = annotation_col_final,
  annotation_colors = annotation_colors_final,
  color = colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(100),
  breaks = seq(-2, 2, length.out = 101),
  border_color = NA,
  fontsize = 10,
  fontsize_row = 10,
  main = "PI3K-AKT-mTOR pathway expression in TCGA-LUAD"
)
  
  # ============================================================
  # 75. 105 PI3K genes: KEAP1 WT vs MUT gene-wise test
  # ============================================================
  
  # 각 patient의 KEAP1 상태를 expression matrix 순서에 맞춤
  keap1_for_expr <- pi3k_patient$KEAP1_status[
    match(
      colnames(expr_heatmap),
      pi3k_patient$patient_id
    )
  ]
  
  # 제대로 맞았는지 확인
  table(keap1_for_expr)
  
  sum(is.na(keap1_for_expr))
  
  
  # gene-wise Wilcoxon test
  gene_test_results <- lapply(
    rownames(expr_heatmap),
    function(gene) {
      
      x_wt <- expr_heatmap[
        gene,
        keap1_for_expr == "WT"
      ]
      
      x_mut <- expr_heatmap[
        gene,
        keap1_for_expr == "MUT"
      ]
      
      test <- wilcox.test(
        x_mut,
        x_wt,
        exact = FALSE
      )
      
      data.frame(
        gene = gene,
        
        WT_median = median(
          x_wt,
          na.rm = TRUE
        ),
        
        MUT_median = median(
          x_mut,
          na.rm = TRUE
        ),
        
        delta_median = median(
          x_mut,
          na.rm = TRUE
        ) - median(
          x_wt,
          na.rm = TRUE
        ),
        
        p_value = test$p.value
      )
    }
  )
  
  gene_test_results <- do.call(
    rbind,
    gene_test_results
  )
  
  # BH-FDR 보정
  gene_test_results$FDR <- p.adjust(
    gene_test_results$p_value,
    method = "BH"
  )
  
  # FDR가 작은 순서로 정렬
  gene_test_results <- gene_test_results[
    order(gene_test_results$FDR),
  ]
  
  rownames(gene_test_results) <- NULL

  
  # 상위 20개 결과
  head(gene_test_results, 20)
  
  # FDR < 0.05인 유전자 개수
  sum(gene_test_results$FDR < 0.05)
  
  # 유의한 유전자들의 방향
  table(
    significant = gene_test_results$FDR < 0.05,
    direction = ifelse(
      gene_test_results$delta_median > 0,
      "Higher_in_MUT",
      "Lower_in_MUT"
    )
  )
  
  # ============================================================
  # 76. FDR < 0.05인 PI3K pathway gene 추출
  # ============================================================
  
  sig_pi3k_genes <- gene_test_results$gene[
    gene_test_results$FDR < 0.05
  ]
  
  length(sig_pi3k_genes)
  
  sig_pi3k_genes
  
  expr_heatmap_sig <- expr_heatmap_z[
    sig_pi3k_genes,
    ,
    drop = FALSE
  ]
  
  dim(expr_heatmap_sig)
  
  
  # ============================================================
  # 77. 유의 유전자 방향 정리
  # ============================================================
  
  sig_info <- gene_test_results[
    gene_test_results$FDR < 0.05,
  ]
  
  sig_info$direction <- ifelse(
    sig_info$delta_median > 0,
    "Higher in MUT",
    "Lower in MUT"
  )
  
  table(sig_info$direction)
  
  
  # ============================================================
  # 78. Heatmap gene 순서 정리
  # ============================================================
  
  sig_info <- sig_info[
    order(
      sig_info$direction,
      -abs(sig_info$delta_median)
    ),
  ]
  
  gene_order <- sig_info$gene
  
  expr_heatmap_sig_ordered <- expr_heatmap_z[
    gene_order,
    ,
    drop = FALSE
  ]
  
  annotation_row_sig <- data.frame(
    Direction = sig_info$direction
  )
  
  rownames(annotation_row_sig) <- sig_info$gene
  
  
  # ============================================================
  # 79. Significant PI3K pathway genes heatmap
  # ============================================================
  
  annotation_col_final <- data.frame(
    KEAP1 = heatmap_meta$KEAP1_status,
    PIK3CB = heatmap_meta$PIK3CB_log2TPM,
    PI3K_score = heatmap_meta$PI3K_ssGSEA
  )
  
  rownames(annotation_col_final) <- rownames(heatmap_meta)
  
  annotation_colors_sig <- list(
    KEAP1 = c(
      WT = "#4DB6AC",
      MUT = "#F4A261"
    ),
    Direction = c(
      "Higher in MUT" = "#E76F51",
      "Lower in MUT" = "#457B9D"
    )
  )
  
  wt_end <- sum(
    heatmap_meta$KEAP1_status == "WT"
  )
  
  pheatmap(
    expr_heatmap_sig_ordered,
    
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    
    gaps_col = wt_end,
    
    show_colnames = FALSE,
    show_rownames = TRUE,
    
    annotation_col = annotation_col_final,
    annotation_row = annotation_row_sig,
    
    annotation_colors = annotation_colors_sig,
    
    color = colorRampPalette(
      c("#2166AC", "white", "#B2182B")
    )(100),
    
    breaks = seq(
      -2, 2,
      length.out = 101
    ),
    
    border_color = NA,
    
    fontsize_row = 7,
    
    main = "Differentially expressed PI3K-AKT-mTOR pathway genes"
  )
  
  # ============================================================
  # 지금까지 만든 Figure 전체 저장
  # ============================================================
  
  dir.create("figures", showWarnings = FALSE)
  library(ggplot2)
  
  master_cohort$KEAP1_status <- factor(
    master_cohort$KEAP1_status,
    levels = c("WT", "MUT")
  )
  
  p_mir_expression <- ggplot(
    master_cohort,
    aes(
      x = KEAP1_status,
      y = miR130b_5p_log2norm,
      fill = KEAP1_status
    )
  ) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.08,
      alpha = 0.35,
      size = 1
    ) +
    scale_fill_manual(
      values = c(
        "WT" = "#00BFC4",
        "MUT" = "#F28E2B"
      )
    ) +
    annotate(
      "text",
      x = 1.5,
      y = max(master_cohort$miR130b_5p_log2norm) * 0.98,
      label = "p = 0.0153",
      size = 5
    ) +
    labs(
      x = "KEAP1 status",
      y = "miR-130b-5p (log2 normalized count + 1)",
      title = "miR-130b-5p expression by KEAP1 status"
    ) +
    theme_bw() +
    theme(
      legend.position = "none"
    )
  
  ggsave(
    "figures/01_miR130b_KEAP1_expression.png",
    p_mir_expression,
    width = 6,
    height = 5,
    dpi = 300
  )
  library(survival)
  library(survminer)
  
  mir_mut_cat <- surv_categorize(
    cut_mir_master
  )
  
  mir_mut_cat$miR_group <- factor(
    mir_mut_cat$miR130b_5p_log2norm,
    levels = c("low", "high"),
    labels = c("Low", "High")
  )
  
  fit_mir_mut <- survfit(
    Surv(OS_days, OS_event) ~ miR_group,
    data = mir_mut_cat
  )
  
  km_mir_mut <- ggsurvplot(
    fit_mir_mut,
    data = mir_mut_cat,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    censor = TRUE,
    palette = c(
      "Low" = "#00A6D6",
      "High" = "#E84A5F"
    ),
    xlab = "Days",
    ylab = "Survival probability",
    legend.title = "miR-130b-5p",
    legend.labs = c("Low", "High"),
    title = "miR-130b-5p in KEAP1 MUT",
    risk.table.height = 0.25,
    ggtheme = theme_bw()
  )
  
  png(
    "figures/02_miR130b_survival_KEAP1_MUT.png",
    width = 1800,
    height = 1500,
    res = 200
  )
  
  print(km_mir_mut)
  
  dev.off()
  
  km_pik3cb_wt <- ggsurvplot(
    fit_pik3cb_wt,
    data = pik3cb_wt_cat,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    censor = TRUE,
    palette = c(
      "Low" = "#00A6D6",
      "High" = "#E84A5F"
    ),
    xlab = "Days",
    ylab = "Survival probability",
    legend.title = "PIK3CB",
    legend.labs = c("Low", "High"),
    title = "PIK3CB in KEAP1 WT",
    risk.table.height = 0.25,
    ggtheme = theme_bw()
  )
  
  png(
    "figures/03_PIK3CB_survival_KEAP1_WT.png",
    width = 1800,
    height = 1500,
    res = 200
  )
  
  print(km_pik3cb_wt)
  
  dev.off()
  
  
  km_pik3cb_mut <- ggsurvplot(
    fit_pik3cb_mut,
    data = pik3cb_mut_cat,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    censor = TRUE,
    palette = c(
      "Low" = "#00A6D6",
      "High" = "#E84A5F"
    ),
    xlab = "Days",
    ylab = "Survival probability",
    legend.title = "PIK3CB",
    legend.labs = c("Low", "High"),
    title = "PIK3CB in KEAP1 MUT",
    risk.table.height = 0.25,
    ggtheme = theme_bw()
  )
  
  png(
    "figures/04_PIK3CB_survival_KEAP1_MUT.png",
    width = 1800,
    height = 1500,
    res = 200
  )
  
  print(km_pik3cb_mut)
  
  dev.off()
  
  
  km_composite <- ggsurvplot(
    fit_composite,
    data = composite_cat,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    censor = TRUE,
    palette = c(
      "Low risk" = "#00A6D6",
      "High risk" = "#E84A5F"
    ),
    xlab = "Days",
    ylab = "Survival probability",
    legend.title = "Composite risk",
    legend.labs = c(
      "Low risk",
      "High risk"
    ),
    title = "PIK3CB / miR-130b-5p composite risk in KEAP1 WT + RT",
    risk.table.height = 0.25,
    ggtheme = theme_bw()
  )
  
  png(
    "figures/05_composite_risk_KEAP1_WT_RT.png",
    width = 1800,
    height = 1500,
    res = 200
  )
  
  print(km_composite)
  
  dev.off()
  
  
  ggsave(
    "figures/06_PIK3CB_PI3K_ssGSEA_WT_MUT.png",
    p_pi3k_final,
    width = 10,
    height = 5,
    dpi = 300
  )
  
  
  library(pheatmap)
  
  png(
    "figures/07_PI3K_heatmap_105_genes.png",
    width = 2600,
    height = 2200,
    res = 200
  )
  
  pheatmap(
    expr_heatmap_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    annotation_col = annotation_col_final,
    annotation_colors = annotation_colors_final,
    gaps_col = wt_end,
    color = colorRampPalette(
      c("#2166AC", "white", "#B2182B")
    )(100),
    breaks = seq(-2, 2, length.out = 101),
    border_color = NA,
    fontsize_row = 6,
    main = "PI3K-AKT-mTOR pathway genes in TCGA-LUAD"
  )
  
  dev.off()
  
  
  png(
    "figures/08_PI3K_heatmap_top25_variable.png",
    width = 2600,
    height = 1800,
    res = 200
  )
  
  pheatmap(
    expr_heatmap_top25,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    gaps_col = wt_end,
    show_colnames = FALSE,
    show_rownames = TRUE,
    annotation_col = annotation_col_final,
    annotation_colors = annotation_colors_final,
    color = colorRampPalette(
      c("#2166AC", "white", "#B2182B")
    )(100),
    breaks = seq(-2, 2, length.out = 101),
    border_color = NA,
    fontsize_row = 10,
    main = "Top variable PI3K-AKT-mTOR pathway genes"
  )
  
  dev.off()
  
  
  
  png(
    "figures/09_PI3K_heatmap_54_significant_genes.png",
    width = 2800,
    height = 2200,
    res = 200
  )
  
  pheatmap(
    expr_heatmap_sig_ordered,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_col = wt_end,
    show_colnames = FALSE,
    show_rownames = TRUE,
    annotation_col = annotation_col_final,
    annotation_row = annotation_row_sig,
    annotation_colors = annotation_colors_sig,
    color = colorRampPalette(
      c("#2166AC", "white", "#B2182B")
    )(100),
    breaks = seq(-2, 2, length.out = 101),
    border_color = NA,
    fontsize_row = 7,
    main = "Significant PI3K-AKT-mTOR pathway genes"
  )
  
  dev.off()
  
  list.files(
    "figures",
    full.names = TRUE
  )
  
  
  # ============================================================
  # 09. 유의한 54개 PI3K pathway gene heatmap 저장
  # ============================================================
  
  png(
    "figures/09_PI3K_heatmap_54_significant_genes.png",
    width = 2800,
    height = 2200,
    res = 200
  )
  
  pheatmap(
    expr_heatmap_sig_ordered,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_col = wt_end,
    show_colnames = FALSE,
    show_rownames = TRUE,
    annotation_col = annotation_col_final,
    annotation_row = annotation_row_sig,
    annotation_colors = annotation_colors_sig,
    color = colorRampPalette(
      c("#2166AC", "white", "#B2182B")
    )(100),
    breaks = seq(-2, 2, length.out = 101),
    border_color = NA,
    fontsize_row = 7,
    main = "Significant PI3K-AKT-mTOR pathway genes"
  )
  
  dev.off()
  
  
  
  list.files("figures")
  
  
  # ============================================================
  # 80. Genome-wide DEG 분석 시작
  # 필요한 객체 확인
  # ============================================================
  
  getwd()
  
  exists("rna_luad")
  exists("keap1_status")
  
  dim(rna_luad)
  table(keap1_status$KEAP1_status)
  
  
  # RNA-seq sample metadata 확인
  colnames(colData(rna_luad))
  
  # sample type 확인
  table(colData(rna_luad)$sample_type, useNA = "ifany")
  
  # 사용 가능한 expression assay 확인
  assayNames(rna_luad)
  
  
  # ============================================================
  # 81. Primary Tumor 환자 중복 확인
  # ============================================================
  
  rna_meta <- as.data.frame(colData(rna_luad))
  
  primary_meta <- rna_meta[
    rna_meta$sample_type == "Primary Tumor",
    c(
      "barcode",
      "patient",
      "sample",
      "sample_submitter_id",
      "sample_type"
    )
  ]
  
  # Primary Tumor sample 수
  nrow(primary_meta)
  
  # Primary Tumor 환자 수
  length(unique(primary_meta$patient))
  
  # 환자당 sample 개수 분포
  table(table(primary_meta$patient))
  
  # 2개 이상의 Primary Tumor sample을 가진 환자
  dup_patients <- names(
    which(table(primary_meta$patient) > 1)
  )
  
  length(dup_patients)
  
  # 중복 환자의 실제 sample 정보
  primary_meta[
    primary_meta$patient %in% dup_patients,
  ][order(
    primary_meta[
      primary_meta$patient %in% dup_patients,
      "patient"
    ]
  ), ]
  
  # ============================================================
  # 82. 중복 Primary Tumor sample의 library size 확인
  # ============================================================
  
  # Primary Tumor 열 위치
  primary_idx <- which(
    colData(rna_luad)$sample_type == "Primary Tumor"
  )
  
  # raw count
  primary_counts <- assay(
    rna_luad,
    "unstranded"
  )[, primary_idx]
  
  # metadata 순서 확인
  primary_meta2 <- as.data.frame(
    colData(rna_luad)[primary_idx, ]
  )
  
  # 각 sample의 전체 sequencing count
  primary_meta2$library_size <- colSums(primary_counts)
  
  # 보기 편하게 million 단위
  primary_meta2$library_M <- primary_meta2$library_size / 1e6
  
  # 중복 환자만 확인
  dup_qc <- primary_meta2[
    primary_meta2$patient %in% dup_patients,
    c(
      "patient",
      "sample",
      "barcode",
      "library_M"
    )
  ]
  
  dup_qc <- dup_qc[
    order(dup_qc$patient, -dup_qc$library_M),
  ]
  
  dup_qc
  
  
  
  # ============================================================
  # 83. 환자당 Primary Tumor 1개 선택
  # 규칙:
  # 1) 01A 우선
  # 2) 같은 vial이면 library size가 큰 sample 우선
  # ============================================================
  
  primary_meta2$vial_priority <- ifelse(
    grepl("-01A$", primary_meta2$sample), 1,
    ifelse(
      grepl("-01B$", primary_meta2$sample), 2,
      ifelse(grepl("-01C$", primary_meta2$sample), 3, 9)
    )
  )
  
  # 환자 → vial 우선순위 → library size 순으로 정렬
  primary_meta2 <- primary_meta2[
    order(
      primary_meta2$patient,
      primary_meta2$vial_priority,
      -primary_meta2$library_size
    ),
  ]
  
  # 환자당 첫 번째 sample만 선택
  primary_unique <- primary_meta2[
    !duplicated(primary_meta2$patient),
  ]
  
  # 확인
  nrow(primary_unique)
  length(unique(primary_unique$patient))
  table(primary_unique$vial_priority)
  
  # 중복 환자 12명에서 최종 선택된 sample 확인
  primary_unique[
    primary_unique$patient %in% dup_patients,
    c("patient", "sample", "barcode", "library_M")
  ]
  
  
  # ============================================================
  # 84. Primary Tumor cohort에 KEAP1 mutation status 매칭
  # ============================================================
  
  # KEAP1 status 데이터 구조 확인
  colnames(keap1_status)
  head(keap1_status)
  
  
  # ============================================================
  # 85. KEAP1 status 매칭 및 최종 DEG cohort 확인
  # ============================================================
  
  primary_unique$KEAP1_status <- keap1_status$KEAP1_status[
    match(
      primary_unique$patient,
      keap1_status$patient_id
    )
  ]
  
  # 517명 중 KEAP1 상태 매칭 결과
  table(primary_unique$KEAP1_status, useNA = "ifany")
  
  # mutation profiling 정보가 없는 환자 수
  sum(is.na(primary_unique$KEAP1_status))
  
  # 최종 DEG cohort
  deg_meta <- primary_unique[
    !is.na(primary_unique$KEAP1_status),
  ]
  
  # 확인
  nrow(deg_meta)
  length(unique(deg_meta$patient))
  table(deg_meta$KEAP1_status)
  
  # ============================================================
  # 86. 최종 507명 genome-wide raw count matrix 생성
  # ============================================================
  
  # rna_luad에서 최종 환자들의 sample 위치 찾기
  deg_idx <- match(
    deg_meta$barcode,
    colnames(rna_luad)
  )
  
  # 매칭 실패가 있는지 확인
  sum(is.na(deg_idx))
  
  # raw count matrix 추출
  deg_counts <- assay(
    rna_luad,
    "unstranded"
  )[, deg_idx]
  
  # sample 이름을 barcode로 유지
  colnames(deg_counts) <- deg_meta$barcode
  
  # 확인
  dim(deg_counts)
  
  # metadata와 count matrix 순서가 정확히 같은지 확인
  all(colnames(deg_counts) == deg_meta$barcode)
  
  # 그룹 확인
  table(deg_meta$KEAP1_status)
  
  
  # ============================================================
  # 87. Low-count gene filtering
  # 기준: count >= 10인 sample이 최소 10개 이상
  # ============================================================
  
  keep_gene <- rowSums(deg_counts >= 10) >= 10
  
  # filtering 전/후 gene 수
  c(
    before = nrow(deg_counts),
    after  = sum(keep_gene),
    removed = sum(!keep_gene)
  )
  
  # filtering
  deg_counts_filt <- deg_counts[keep_gene, ]
  
  # 최종 크기 확인
  dim(deg_counts_filt)
  
  # count matrix 기본 확인
  summary(rowSums(deg_counts_filt))
  
  
  
  # ============================================================
  # 88. DESeq2 객체 생성 + VST 준비
  # ============================================================
  
  library(DESeq2)
  
  # WT를 기준(reference)으로 설정
  deg_meta$KEAP1_status <- factor(
    deg_meta$KEAP1_status,
    levels = c("WT", "MUT")
  )
  
  # DESeq2용 metadata
  deg_coldata <- data.frame(
    KEAP1_status = deg_meta$KEAP1_status,
    row.names = deg_meta$barcode
  )
  
  # 순서가 정확히 맞는지 확인
  all(colnames(deg_counts_filt) == rownames(deg_coldata))
  
  # DESeq2 객체 생성
  dds_deg <- DESeqDataSetFromMatrix(
    countData = round(deg_counts_filt),
    colData = deg_coldata,
    design = ~ KEAP1_status
  )
  
  # normalization
  dds_deg <- estimateSizeFactors(dds_deg)
  
  # PCA용 variance stabilizing transformation
  vsd_deg <- vst(
    dds_deg,
    blind = TRUE
  )
  
  # 확인
  dim(dds_deg)
  dim(vsd_deg)
  
  summary(sizeFactors(dds_deg))
  
  table(colData(dds_deg)$KEAP1_status)
  
  
  # ============================================================
  # 89. PCA QC - KEAP1 WT vs MUT
  # ============================================================
  
  library(ggplot2)
  
  # PCA 데이터 추출
  pca_data <- plotPCA(
    vsd_deg,
    intgroup = "KEAP1_status",
    returnData = TRUE
  )
  
  # 각 PC가 설명하는 분산 비율
  percent_var <- round(
    100 * attr(pca_data, "percentVar"),
    1
  )
  
  # 확인
  head(pca_data)
  
  percent_var
  
  # PCA plot
  p_pca <- ggplot(
    pca_data,
    aes(
      x = PC1,
      y = PC2,
      color = KEAP1_status
    )
  ) +
    geom_point(
      size = 2,
      alpha = 0.7
    ) +
    xlab(
      paste0("PC1: ", percent_var[1], "% variance")
    ) +
    ylab(
      paste0("PC2: ", percent_var[2], "% variance")
    ) +
    ggtitle(
      "TCGA-LUAD RNA-seq PCA: KEAP1 WT vs MUT"
    ) +
    theme_classic(base_size = 13)
  
  p_pca
  
  print(p_pca)
  ggsave(
    "figures/10_DEG_PCA_KEAP1_WT_MUT.png",
    plot = p_pca,
    width = 7,
    height = 5.5,
    dpi = 300
  )
  
  
  # ============================================================
  # 90. Genome-wide DESeq2: KEAP1 MUT vs WT
  # ============================================================
  
  # DESeq2 모델 fitting
  dds_deg <- DESeq(dds_deg)
  
  # MUT vs WT 결과
  # log2FoldChange > 0 : MUT에서 높음
  # log2FoldChange < 0 : MUT에서 낮음
  res_deg <- results(
    dds_deg,
    contrast = c("KEAP1_status", "MUT", "WT"),
    alpha = 0.05
  )
  
  # 기본 결과 확인
  summary(res_deg)
  
  # data.frame으로 변환
  res_deg_df <- as.data.frame(res_deg)
  
  # 전체 gene 수
  nrow(res_deg_df)
  
  # FDR < 0.05
  sum(res_deg_df$padj < 0.05, na.rm = TRUE)
  
  # FDR < 0.05 + |log2FC| >= 1
  sum(
    res_deg_df$padj < 0.05 &
      abs(res_deg_df$log2FoldChange) >= 1,
    na.rm = TRUE
  )
  
  
  # ============================================================
  # 91. RNA-seq gene annotation 구조 확인
  # ============================================================
  
  colnames(rowData(rna_luad))
  
  head(
    as.data.frame(rowData(rna_luad))
  )
  
  
  # ============================================================
  # 92. DEG 결과에 gene annotation 추가
  # ============================================================
  
  # Ensembl ID 저장
  res_deg_df$gene_id <- rownames(res_deg_df)
  
  # RNA annotation
  gene_anno <- as.data.frame(rowData(rna_luad))[
    ,
    c("gene_id", "gene_name", "gene_type")
  ]
  
  # DEG 결과와 annotation 매칭
  anno_idx <- match(
    res_deg_df$gene_id,
    gene_anno$gene_id
  )
  
  res_deg_df$gene_name <- gene_anno$gene_name[anno_idx]
  res_deg_df$gene_type <- gene_anno$gene_type[anno_idx]
  
  # annotation 확인
  head(
    res_deg_df[
      ,
      c(
        "gene_id",
        "gene_name",
        "gene_type",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj"
      )
    ]
  )
  
  # gene symbol 매칭 실패 확인
  sum(is.na(res_deg_df$gene_name))
  
  
  # ============================================================
  # 93. Significant DEG 정의
  # 기준: padj < 0.05 & |log2FC| >= 1
  # ============================================================
  
  sig_deg <- res_deg_df[
    !is.na(res_deg_df$padj) &
      res_deg_df$padj < 0.05 &
      abs(res_deg_df$log2FoldChange) >= 1,
  ]
  
  sig_deg$direction <- ifelse(
    sig_deg$log2FoldChange > 0,
    "Up_in_MUT",
    "Down_in_MUT"
  )
  
  # 전체 / 방향별 개수
  nrow(sig_deg)
  table(sig_deg$direction)
  
  
  
  # MUT에서 증가한 유전자 Top 20
  top_up <- sig_deg[
    order(-sig_deg$log2FoldChange),
  ]
  
  head(
    top_up[
      ,
      c(
        "gene_name",
        "gene_type",
        "baseMean",
        "log2FoldChange",
        "padj"
      )
    ],
    20
  )
  
  # MUT에서 감소한 유전자 Top 20
  top_down <- sig_deg[
    order(sig_deg$log2FoldChange),
  ]
  
  head(
    top_down[
      ,
      c(
        "gene_name",
        "gene_type",
        "baseMean",
        "log2FoldChange",
        "padj"
      )
    ],
    20
  )
  
  
  # ============================================================
  # 94. MA plot
  # ============================================================
  
  png(
    "figures/11_DEG_MA_KEAP1_WT_MUT.png",
    width = 1800,
    height = 1500,
    res = 250
  )
  
  plotMA(
    res_deg,
    alpha = 0.05,
    ylim = c(-6, 6),
    main = "KEAP1 MUT vs WT - MA plot"
  )
  
  dev.off()
  
  # RStudio 화면에도 표시
  plotMA(
    res_deg,
    alpha = 0.05,
    ylim = c(-6, 6),
    main = "KEAP1 MUT vs WT - MA plot"
  )
  
  plotMA(
    res_deg,
    alpha = 0.05,
    ylim = c(-6, 6),
    main = "KEAP1 MUT vs WT - MA plot"
  )

  graphics.off()
  
  plotMA(
    res_deg,
    alpha = 0.05,
    ylim = c(-6, 6),
    main = "KEAP1 MUT vs WT - MA plot"
  )
  
  
  # ============================================================
  # 95. Volcano plot
  # padj < 0.05 & |log2FC| >= 1
  # ============================================================
  
  volcano_df <- res_deg_df[
    !is.na(res_deg_df$padj) &
      !is.na(res_deg_df$log2FoldChange),
  ]
  
  # -log10(FDR)
  volcano_df$neglog10_padj <- -log10(
    pmax(volcano_df$padj, .Machine$double.xmin)
  )
  
  # DEG 분류
  volcano_df$DEG <- "Not significant"
  
  volcano_df$DEG[
    volcano_df$padj < 0.05 &
      volcano_df$log2FoldChange >= 1
  ] <- "Up in MUT"
  
  volcano_df$DEG[
    volcano_df$padj < 0.05 &
      volcano_df$log2FoldChange <= -1
  ] <- "Down in MUT"
  
  # 개수 확인
  table(volcano_df$DEG)
  
  dev.cur()
  
  
  graphics.off()
  
  
  library(ggplot2)
  
  print(p_volcano)
  
  
  # Volcano plot용 데이터 만들기
  volcano_df <- res_deg_df[
    !is.na(res_deg_df$padj) &
      !is.na(res_deg_df$log2FoldChange),
  ]
  
  volcano_df$neglog10_padj <- -log10(
    pmax(volcano_df$padj, .Machine$double.xmin)
  )
  
  volcano_df$DEG <- "Not significant"
  
  volcano_df$DEG[
    volcano_df$padj < 0.05 &
      volcano_df$log2FoldChange >= 1
  ] <- "Up in MUT"
  
  volcano_df$DEG[
    volcano_df$padj < 0.05 &
      volcano_df$log2FoldChange <= -1
  ] <- "Down in MUT"
  
  # 개수 확인
  table(volcano_df$DEG)
  
  
  library(ggplot2)
  
  p_volcano <- ggplot(
    volcano_df,
    aes(
      x = log2FoldChange,
      y = neglog10_padj,
      color = DEG
    )
  ) +
    geom_point(
      alpha = 0.5,
      size = 1.3
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    xlab("log2 Fold Change (MUT vs WT)") +
    ylab("-log10 adjusted p-value") +
    ggtitle("KEAP1 MUT vs WT - Volcano plot") +
    theme_classic(base_size = 13)
  
  
  # label을 더 정확하게 수정
  volcano_df$DEG <- "Other"
  
  volcano_df$DEG[
    volcano_df$padj < 0.05 &
      volcano_df$log2FoldChange >= 1
  ] <- "Up in MUT"
  
  volcano_df$DEG[
    volcano_df$padj < 0.05 &
      volcano_df$log2FoldChange <= -1
  ] <- "Down in MUT"
  
  table(volcano_df$DEG)
  
  p_volcano <- ggplot(
    volcano_df,
    aes(
      x = log2FoldChange,
      y = neglog10_padj,
      color = DEG
    )
  ) +
    geom_point(alpha = 0.5, size = 1.3) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    xlab("log2 Fold Change (MUT vs WT)") +
    ylab("-log10 adjusted p-value") +
    ggtitle("KEAP1 MUT vs WT - Volcano plot") +
    theme_classic(base_size = 13)
  
  print(p_volcano)
  
  ggsave(
    "figures/12_DEG_Volcano_KEAP1_WT_MUT.png",
    plot = p_volcano,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  
  # ============================================================
  # 96. 해석용 robust protein-coding DEG 확인
  # DEG 정의 자체는 그대로 유지
  # ============================================================
  
  robust_deg <- sig_deg[
    sig_deg$gene_type == "protein_coding" &
      sig_deg$baseMean >= 50,
  ]
  
  # 개수
  nrow(robust_deg)
  
  table(robust_deg$direction)
  
  
  # MUT에서 증가
  robust_up <- robust_deg[
    robust_deg$direction == "Up_in_MUT",
  ]
  
  robust_up <- robust_up[
    order(-robust_up$log2FoldChange),
  ]
  
  head(
    robust_up[
      ,
      c(
        "gene_name",
        "baseMean",
        "log2FoldChange",
        "padj"
      )
    ],
    20
  )
  
  # MUT에서 감소
  robust_down <- robust_deg[
    robust_deg$direction == "Down_in_MUT",
  ]
  
  robust_down <- robust_down[
    order(robust_down$log2FoldChange),
  ]
  
  head(
    robust_down[
      ,
      c(
        "gene_name",
        "baseMean",
        "log2FoldChange",
        "padj"
      )
    ],
    20
  )
  
  
  # ============================================================
  # 97. Top 40 robust DEG heatmap
  # Up 20 + Down 20
  # ============================================================
  
  library(pheatmap)
  
  # Top 20 Up / Top 20 Down
  heatmap_genes <- c(
    head(robust_up$gene_id, 20),
    head(robust_down$gene_id, 20)
  )
  
  # VST expression matrix
  vst_mat <- assay(vsd_deg)
  
  # 선택한 gene만 추출
  heatmap_mat <- vst_mat[
    rownames(vst_mat) %in% heatmap_genes,
    ,
    drop = FALSE
  ]
  
  # gene symbol 붙이기
  heatmap_symbols <- res_deg_df$gene_name[
    match(rownames(heatmap_mat), res_deg_df$gene_id)
  ]
  
  rownames(heatmap_mat) <- heatmap_symbols
  
  # 각 gene별 Z-score
  heatmap_z <- t(
    scale(
      t(heatmap_mat)
    )
  )
  
  # 샘플 annotation
  heatmap_anno <- data.frame(
    KEAP1 = deg_meta$KEAP1_status
  )
  
  rownames(heatmap_anno) <- deg_meta$barcode
  
  # sample을 WT -> MUT 순서로 정렬
  sample_order <- order(deg_meta$KEAP1_status)
  
  heatmap_z <- heatmap_z[, sample_order]
  
  heatmap_anno <- heatmap_anno[
    sample_order,
    ,
    drop = FALSE
  ]
  
  # 확인
  dim(heatmap_z)
  
  
  pheatmap(
    heatmap_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    annotation_col = heatmap_anno,
    fontsize_row = 8,
    main = "Top robust DEGs: KEAP1 MUT vs WT"
  )
  
  
  png(
    "figures/13_DEG_Top40_heatmap_KEAP1_WT_MUT.png",
    width = 2400,
    height = 1500,
    res = 250
  )
  
  pheatmap(
    heatmap_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    annotation_col = heatmap_anno,
    fontsize_row = 8,
    main = "Top robust DEGs: KEAP1 MUT vs WT"
  )
  
  dev.off()
  
  
  requireNamespace("clusterProfiler", quietly = TRUE)
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
  requireNamespace("enrichplot", quietly = TRUE)
  
  
  
  c(
    clusterProfiler = requireNamespace("clusterProfiler", quietly = TRUE),
    org.Hs.eg.db = requireNamespace("org.Hs.eg.db", quietly = TRUE),
    enrichplot = requireNamespace("enrichplot", quietly = TRUE)
  )
  
  
  
  # Bioconductor 설치 도구가 없으면 먼저 설치
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # GO/KEGG 분석 패키지 설치
  BiocManager::install(
    c(
      "clusterProfiler",
      "org.Hs.eg.db",
      "enrichplot"
    ),
    ask = FALSE,
    update = FALSE
  )
  
  c(
    clusterProfiler = requireNamespace("clusterProfiler", quietly = TRUE),
    org.Hs.eg.db = requireNamespace("org.Hs.eg.db", quietly = TRUE),
    enrichplot = requireNamespace("enrichplot", quietly = TRUE)
  )
  
  
  # ============================================================
  # 98. GO enrichment용 Gene ID 매핑
  # ============================================================
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  
  # Ensembl 버전 번호 제거
  # 예: ENSG00000125851.10 -> ENSG00000125851
  res_deg_df$ensembl_clean <- sub(
    "\\..*$",
    "",
    res_deg_df$gene_id
  )
  
  sig_deg$ensembl_clean <- sub(
    "\\..*$",
    "",
    sig_deg$gene_id
  )
  
  # Ensembl -> Entrez ID 변환
  res_deg_df$ENTREZID <- mapIds(
    org.Hs.eg.db,
    keys = res_deg_df$ensembl_clean,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  sig_deg$ENTREZID <- mapIds(
    org.Hs.eg.db,
    keys = sig_deg$ensembl_clean,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  # 전체 background
  go_universe <- unique(
    na.omit(res_deg_df$ENTREZID)
  )
  
  # MUT에서 증가한 DEG
  go_up <- unique(
    na.omit(
      sig_deg$ENTREZID[
        sig_deg$direction == "Up_in_MUT"
      ]
    )
  )
  
  # MUT에서 감소한 DEG
  go_down <- unique(
    na.omit(
      sig_deg$ENTREZID[
        sig_deg$direction == "Down_in_MUT"
      ]
    )
  )
  
  # 매핑 결과 확인
  c(
    tested_genes = nrow(res_deg_df),
    universe_mapped = length(go_universe),
    up_DEG_original = sum(sig_deg$direction == "Up_in_MUT"),
    up_DEG_mapped = length(go_up),
    down_DEG_original = sum(sig_deg$direction == "Down_in_MUT"),
    down_DEG_mapped = length(go_down)
  )
  
  
  # ============================================================
  # 98. GO enrichment용 Gene ID 매핑
  # ============================================================
  
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  
  # Ensembl 버전 번호 제거
  # 예: ENSG00000125851.10 -> ENSG00000125851
  res_deg_df$ensembl_clean <- sub(
    "\\..*$",
    "",
    res_deg_df$gene_id
  )
  
  sig_deg$ensembl_clean <- sub(
    "\\..*$",
    "",
    sig_deg$gene_id
  )
  
  # Ensembl -> Entrez ID 변환
  res_deg_df$ENTREZID <- mapIds(
    org.Hs.eg.db,
    keys = res_deg_df$ensembl_clean,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  sig_deg$ENTREZID <- mapIds(
    org.Hs.eg.db,
    keys = sig_deg$ensembl_clean,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )
  
  # 전체 background
  go_universe <- unique(
    na.omit(res_deg_df$ENTREZID)
  )
  
  # MUT에서 증가한 DEG
  go_up <- unique(
    na.omit(
      sig_deg$ENTREZID[
        sig_deg$direction == "Up_in_MUT"
      ]
    )
  )
  
  # MUT에서 감소한 DEG
  go_down <- unique(
    na.omit(
      sig_deg$ENTREZID[
        sig_deg$direction == "Down_in_MUT"
      ]
    )
  )
  
  # 매핑 결과 확인
  c(
    tested_genes = nrow(res_deg_df),
    universe_mapped = length(go_universe),
    up_DEG_original = sum(sig_deg$direction == "Up_in_MUT"),
    up_DEG_mapped = length(go_up),
    down_DEG_original = sum(sig_deg$direction == "Down_in_MUT"),
    down_DEG_mapped = length(go_down)
  )
  
  
  # ============================================================
  # 99. GO Biological Process enrichment
  # Up / Down DEG separately
  # ============================================================
  
  go_bp_up <- enrichGO(
    gene          = go_up,
    universe      = go_universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  go_bp_down <- enrichGO(
    gene          = go_down,
    universe      = go_universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  
  
  # 유의한 GO BP term 개수
  nrow(as.data.frame(go_bp_up))
  nrow(as.data.frame(go_bp_down))
  
  
  # MUT에서 증가한 DEG의 GO
  head(
    as.data.frame(go_bp_up)[
      ,
      c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "Count")
    ],
    15
  )
  
  # MUT에서 감소한 DEG의 GO
  head(
    as.data.frame(go_bp_down)[
      ,
      c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "Count")
    ],
    15
  )
  
  
  # ============================================================
  # 100. GO BP dotplot
  # ============================================================
  
  library(enrichplot)
  library(ggplot2)
  
  p_go_up <- dotplot(
    go_bp_up,
    showCategory = 15
  ) +
    ggtitle("GO Biological Process - Up in KEAP1 MUT")
  
  p_go_down <- dotplot(
    go_bp_down,
    showCategory = 15
  ) +
    ggtitle("GO Biological Process - Down in KEAP1 MUT")
  
  # Plots 창에서 보기
  print(p_go_up)
  
  
  # ============================================================
  # 100. GO Biological Process dotplot
  # KEAP1 MUT에서 증가한 DEG
  # ============================================================
  
  library(enrichplot)
  library(ggplot2)
  
  p_go_up <- dotplot(
    go_bp_up,
    showCategory = 15
  ) +
    ggtitle("GO Biological Process - Up in KEAP1 MUT")
  
  print(p_go_up)
  
  
  # 그래픽 장치 초기화
  graphics.off()
  
  library(enrichplot)
  library(ggplot2)
  
  # Up GO dotplot 다시 생성
  p_go_up <- dotplot(
    go_bp_up,
    showCategory = 15
  ) +
    ggtitle("GO Biological Process - Up in KEAP1 MUT")
  
  # Plots 창에 출력
  print(p_go_up)
  
  
  # ============================================================
  # 101. GO Biological Process dotplot
  # KEAP1 MUT에서 감소한 DEG
  # ============================================================
  
  graphics.off()
  
  p_go_down <- dotplot(
    go_bp_down,
    showCategory = 15
  ) +
    ggtitle("GO Biological Process - Down in KEAP1 MUT")
  
  print(p_go_down)
  
  
  ggsave(
    "figures/15_GO_BP_Down_KEAP1_MUT.png",
    plot = p_go_down,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  ggsave(
    "figures/14_GO_BP_Up_KEAP1_MUT.png",
    plot = p_go_up,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  # ============================================================
  # 102. KEGG enrichment
  # Up / Down DEG separately
  # ============================================================
  
  kegg_up <- enrichKEGG(
    gene          = go_up,
    universe      = go_universe,
    organism      = "hsa",
    keyType       = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  kegg_down <- enrichKEGG(
    gene          = go_down,
    universe      = go_universe,
    organism      = "hsa",
    keyType       = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  # 유의 pathway 수
  nrow(as.data.frame(kegg_up))
  nrow(as.data.frame(kegg_down))
  
  # ============================================================
  # 103. KEGG enrichment 결과 확인
  # ============================================================
  
  # KEAP1 MUT에서 증가한 DEG
  head(
    as.data.frame(kegg_up)[
      ,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ],
    15
  )
  
  
  # KEAP1 MUT에서 감소한 DEG
  head(
    as.data.frame(kegg_down)[
      ,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "p.adjust",
        "Count"
      )
    ],
    15
  )
  
  
  # ============================================================
  # 104. KEGG enrichment dotplot
  # ============================================================
  
  library(enrichplot)
  library(ggplot2)
  
  graphics.off()
  
  # MUT에서 증가
  p_kegg_up <- dotplot(
    kegg_up,
    showCategory = 15
  ) +
    ggtitle("KEGG Pathways - Up in KEAP1 MUT")
  
  print(p_kegg_up)
  
  p_kegg_down <- dotplot(
    kegg_down,
    showCategory = 8
  ) +
    ggtitle("KEGG Pathways - Down in KEAP1 MUT")
  
  print(p_kegg_down)
  
  ggsave(
    "figures/16_KEGG_Up_KEAP1_MUT.png",
    plot = p_kegg_up,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    "figures/17_KEGG_Down_KEAP1_MUT.png",
    plot = p_kegg_down,
    width = 8,
    height = 5,
    dpi = 300
  )
  
  
  # ============================================================
  # 105. GSEA용 ranked gene list 생성
  # 전체 유전자 사용
  # ============================================================
  
  gsea_df <- res_deg_df[
    !is.na(res_deg_df$stat) &
      !is.na(res_deg_df$ENTREZID),
    c(
      "gene_name",
      "ENTREZID",
      "log2FoldChange",
      "stat",
      "padj"
    )
  ]
  
  # 같은 Entrez ID가 여러 번 매핑된 경우
  # |stat|가 가장 큰 gene 하나만 유지
  gsea_df <- gsea_df[
    order(-abs(gsea_df$stat)),
  ]
  
  gsea_df <- gsea_df[
    !duplicated(gsea_df$ENTREZID),
  ]
  
  # GSEA용 named vector 생성
  gene_rank <- gsea_df$stat
  
  names(gene_rank) <- gsea_df$ENTREZID
  
  # GSEA에서는 반드시 큰 값 -> 작은 값 순으로 정렬
  gene_rank <- sort(
    gene_rank,
    decreasing = TRUE
  )
  
  # 확인
  length(gene_rank)
  
  head(gene_rank, 10)
  
  tail(gene_rank, 10)
  
  
  # ============================================================
  # 106. Hallmark gene set 준비
  # ============================================================
  
  library(msigdbr)
  
  hallmark_df <- msigdbr(
    species = "Homo sapiens",
    collection = "H"
  )
  
  # 구조 확인
  dim(hallmark_df)
  
  colnames(hallmark_df)
  
  # Hallmark pathway 개수
  length(unique(hallmark_df$gs_name))
  
  head(hallmark_df)
  
  
  dim(hallmark_df)
  
  colnames(hallmark_df)
  
  length(unique(hallmark_df$gs_name))
  
  
  
  # ============================================================
  # 107. Hallmark GSEA
  # ============================================================
  
  # GSEA용 Hallmark TERM2GENE 만들기
  hallmark_t2g <- hallmark_df[
    !is.na(hallmark_df$ncbi_gene),
    c("gs_name", "ncbi_gene")
  ]
  
  # Entrez ID를 문자형으로 통일
  hallmark_t2g$ncbi_gene <- as.character(
    hallmark_t2g$ncbi_gene
  )
  
  # 중복 제거
  hallmark_t2g <- unique(hallmark_t2g)
  
  # gene_rank의 이름도 문자형 확인
  names(gene_rank) <- as.character(names(gene_rank))
  
  # Hallmark GSEA 실행
  gsea_hallmark <- GSEA(
    geneList      = gene_rank,
    TERM2GENE     = hallmark_t2g,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    minGSSize     = 10,
    maxGSSize     = 500,
    verbose       = FALSE,
    seed          = TRUE
  )
  
  # ============================================================
  # 107-1. Hallmark GSEA 재실행
  # 매우 작은 p-value까지 정확하게 계산
  # ============================================================
  
  gsea_hallmark <- GSEA(
    geneList      = gene_rank,
    TERM2GENE     = hallmark_t2g,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    minGSSize     = 10,
    maxGSSize     = 500,
    eps           = 0,
    verbose       = FALSE,
    seed          = TRUE
  )
  
  
  gsea_hallmark_df <- as.data.frame(gsea_hallmark)
  
  # 분석된 pathway 수
  nrow(gsea_hallmark_df)
  
  # FDR < 0.05 pathway 수
  sum(
    gsea_hallmark_df$p.adjust < 0.05,
    na.rm = TRUE
  )
  
  # 결과 컬럼
  colnames(gsea_hallmark_df)
  
  
  # MUT에서 가장 강하게 증가한 pathway
  head(
    gsea_hallmark_df[
      order(-gsea_hallmark_df$NES),
      c("ID", "NES", "p.adjust", "core_enrichment")
    ],
    10
  )
  
  # MUT에서 가장 강하게 감소한 pathway
  head(
    gsea_hallmark_df[
      order(gsea_hallmark_df$NES),
      c("ID", "NES", "p.adjust", "core_enrichment")
    ],
    10
  )
  
  
  # ============================================================
  # 108. Hallmark GSEA Top pathways NES plot
  # ============================================================
  
  library(ggplot2)
  
  # FDR < 0.05만 선택
  gsea_sig <- gsea_hallmark_df[
    !is.na(gsea_hallmark_df$p.adjust) &
      gsea_hallmark_df$p.adjust < 0.05,
  ]
  
  # MUT 쪽 enrichment: NES > 0
  gsea_pos <- gsea_sig[
    gsea_sig$NES > 0,
  ]
  
  gsea_pos <- gsea_pos[
    order(-gsea_pos$NES),
  ]
  
  gsea_pos_top <- head(gsea_pos, 10)
  
  # WT 쪽 enrichment: NES < 0
  gsea_neg <- gsea_sig[
    gsea_sig$NES < 0,
  ]
  
  gsea_neg <- gsea_neg[
    order(gsea_neg$NES),
  ]
  
  gsea_neg_top <- head(gsea_neg, 10)
  
  # 두 그룹 합치기
  gsea_top20 <- rbind(
    gsea_neg_top,
    gsea_pos_top
  )
  
  # pathway 이름 보기 좋게 정리
  gsea_top20$Pathway <- gsub(
    "^HALLMARK_",
    "",
    gsea_top20$ID
  )
  
  gsea_top20$Pathway <- gsub(
    "_",
    " ",
    gsea_top20$Pathway
  )
  
  # 방향 표시
  gsea_top20$Direction <- ifelse(
    gsea_top20$NES > 0,
    "Enriched in MUT",
    "Enriched in WT"
  )
  
  # NES 순서대로 factor 설정
  gsea_top20$Pathway <- factor(
    gsea_top20$Pathway,
    levels = gsea_top20$Pathway[
      order(gsea_top20$NES)
    ]
  )
  
  # 확인
  gsea_top20[
    ,
    c("Pathway", "NES", "p.adjust", "Direction")
  ]
  
  
  
  # ============================================================
  # 109. GSEA NES plot
  # ============================================================
  
  p_gsea_nes <- ggplot(
    gsea_top20,
    aes(
      x = NES,
      y = Pathway,
      fill = Direction
    )
  ) +
    geom_col(width = 0.75) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed"
    ) +
    xlab("Normalized Enrichment Score (NES)") +
    ylab(NULL) +
    ggtitle(
      "Hallmark GSEA: KEAP1 MUT vs WT"
    ) +
    theme_classic(base_size = 12)
  
  print(p_gsea_nes)
  
  
  
  # ============================================================
  # 110. 대표 GSEA enrichment curve
  # ROS pathway
  # ============================================================
  
  library(enrichplot)
  
  graphics.off()
  
  p_gsea_ros <- gseaplot2(
    gsea_hallmark,
    geneSetID = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    title = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
  )
  
  print(p_gsea_ros)
  
  
  ggsave(
    "figures/19_GSEA_ROS_KEAP1_MUT_WT.png",
    plot = p_gsea_ros,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  # ============================================================
  # 111. Oxidative phosphorylation GSEA curve
  # ============================================================
  
  graphics.off()
  
  p_gsea_oxphos <- gseaplot2(
    gsea_hallmark,
    geneSetID = "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    title = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
  )
  
  print(p_gsea_oxphos)
  
  
  ggsave(
    "figures/20_GSEA_OXPHOS_KEAP1_MUT_WT.png",
    plot = p_gsea_oxphos,
    width = 8,
    height = 6,
    dpi = 300
  )
    
    # ============================================================
    # 112. DEG 결과 저장
    # ============================================================
    
    # 전체 DEG 결과
    write.csv(
      res_deg_df,
      "results/DEG_all_KEAP1_MUT_vs_WT.csv",
      row.names = FALSE
    )
    
    # significant DEG
    write.csv(
      sig_deg,
      "results/DEG_significant_KEAP1_MUT_vs_WT.csv",
      row.names = FALSE
    )
    
    # robust protein-coding DEG
    write.csv(
      robust_deg,
      "results/DEG_robust_protein_coding_KEAP1_MUT_vs_WT.csv",
      row.names = FALSE
    )
    
    
    # ============================================================
    # 113. Pathway 결과 저장
    # ============================================================
    
    write.csv(
      as.data.frame(go_bp_up),
      "results/GO_BP_Up_KEAP1_MUT.csv",
      row.names = FALSE
    )
    
    write.csv(
      as.data.frame(go_bp_down),
      "results/GO_BP_Down_KEAP1_MUT.csv",
      row.names = FALSE
    )
    
    write.csv(
      as.data.frame(kegg_up),
      "results/KEGG_Up_KEAP1_MUT.csv",
      row.names = FALSE
    )
    
    write.csv(
      as.data.frame(kegg_down),
      "results/KEGG_Down_KEAP1_MUT.csv",
      row.names = FALSE
    )
    
    write.csv(
      gsea_hallmark_df,
      "results/GSEA_Hallmark_KEAP1_MUT_vs_WT.csv",
      row.names = FALSE
    )
  
  
  list.files("results")
  
  
  list.files("figures")
  
  
  # ============================================================
  # 114. Hallmark GSEA NES plot 저장
  # ============================================================
  
  ggsave(
    "figures/18_GSEA_Hallmark_NES_KEAP1_MUT_WT.png",
    plot = p_gsea_nes,
    width = 9,
    height = 7,
    dpi = 300
  )
  
  file.exists(
    "figures/18_GSEA_Hallmark_NES_KEAP1_MUT_WT.png"
  )
  
  
  # ============================================================
  # 115. Clinical confounder 구조 확인
  # ============================================================
  
  # Age
  summary(deg_meta$age_at_diagnosis)
  
  # Sex
  table(
    deg_meta$KEAP1_status,
    deg_meta$sex_at_birth,
    useNA = "ifany"
  )
  
  # Pathologic stage
  table(
    deg_meta$KEAP1_status,
    deg_meta$ajcc_pathologic_stage,
    useNA = "ifany"
  )
  
  # Smoking status
  table(
    deg_meta$KEAP1_status,
    deg_meta$tobacco_smoking_status,
    useNA = "ifany"
  )
  
  
  # ============================================================
  # 116. Clinical confounder statistical tests
  # ============================================================
  
  # ------------------------------------------------------------
  # 1. Age
  # ------------------------------------------------------------
  
  deg_meta$age_years <- deg_meta$age_at_diagnosis / 365.25
  
  aggregate(
    age_years ~ KEAP1_status,
    data = deg_meta,
    FUN = function(x) c(
      n = sum(!is.na(x)),
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE)
    )
  )
  
  wilcox.test(
    age_years ~ KEAP1_status,
    data = deg_meta
  )
  
  
  # ------------------------------------------------------------
  # 2. Sex
  # ------------------------------------------------------------
  
  sex_tab <- table(
    deg_meta$KEAP1_status,
    deg_meta$sex_at_birth
  )
  
  sex_tab
  
  chisq.test(sex_tab)
  
  
  # ------------------------------------------------------------
  # 3. Pathologic stage
  # ------------------------------------------------------------
  
  stage_raw <- as.character(
    deg_meta$ajcc_pathologic_stage
  )
  
  deg_meta$stage_group <- ifelse(
    grepl("^Stage IV", stage_raw), "IV",
    ifelse(
      grepl("^Stage III", stage_raw), "III",
      ifelse(
        grepl("^Stage II", stage_raw), "II",
        ifelse(
          grepl("^Stage I", stage_raw), "I",
          NA
        )
      )
    )
  )
  
  stage_tab <- table(
    deg_meta$KEAP1_status,
    deg_meta$stage_group
  )
  
  stage_tab
  
  chisq.test(stage_tab)
  
  
  # ------------------------------------------------------------
  # 4. Smoking
  # ------------------------------------------------------------
  
  smoke_raw <- as.character(
    deg_meta$tobacco_smoking_status
  )
  
  deg_meta$smoking_group <- ifelse(
    smoke_raw == "Lifelong Non-Smoker",
    "Never",
    ifelse(
      smoke_raw %in% c(
        "Current Smoker",
        "Current Reformed Smoker for < or = 15 yrs",
        "Current Reformed Smoker for > 15 yrs",
        "Current Reformed Smoker, Duration Not Specified"
      ),
      "Ever",
      NA
    )
  )
  
  smoking_tab <- table(
    deg_meta$KEAP1_status,
    deg_meta$smoking_group
  )
  
  smoking_tab
  
  chisq.test(smoking_tab)
  
  
  # ============================================================
  # 117. Stage 정확한 검정 + adjusted analysis cohort 확인
  # ============================================================
  
  # Stage: 작은 셀이 있으므로 Fisher exact test
  fisher.test(
    stage_tab,
    simulate.p.value = TRUE,
    B = 10000
  )
  
  # adjusted DEG에 사용할 변수 정리
  deg_meta$sex_group <- factor(
    deg_meta$sex_at_birth
  )
  
  deg_meta$stage_group <- factor(
    deg_meta$stage_group,
    levels = c("I", "II", "III", "IV")
  )
  
  deg_meta$smoking_group <- factor(
    deg_meta$smoking_group,
    levels = c("Never", "Ever")
  )
  
  # Sex + Stage + Smoking이 모두 있는 환자만 선택
  adjust_keep <- complete.cases(
    deg_meta[
      ,
      c(
        "KEAP1_status",
        "sex_group",
        "stage_group",
        "smoking_group"
      )
    ]
  )
  
  # 최종 보정 분석 환자 수
  sum(adjust_keep)
  
  # WT / MUT 수
  table(
    deg_meta$KEAP1_status[adjust_keep]
  )
  
  # 각 변수 분포
  table(
    deg_meta$sex_group[adjust_keep],
    deg_meta$KEAP1_status[adjust_keep]
  )
  
  table(
    deg_meta$stage_group[adjust_keep],
    deg_meta$KEAP1_status[adjust_keep]
  )
  
  table(
    deg_meta$smoking_group[adjust_keep],
    deg_meta$KEAP1_status[adjust_keep]
  )
  
  
  # ============================================================
  # 118. Adjusted DESeq2 cohort 준비
  # Covariates: Sex + Stage + Smoking
  # ============================================================
  
  # complete-case 환자만 count matrix에서 선택
  adj_counts <- deg_counts_filt[, adjust_keep]
  
  # metadata도 같은 환자만 선택
  adj_meta <- deg_meta[adjust_keep, ]
  
  # row name을 barcode로 설정
  rownames(adj_meta) <- adj_meta$barcode
  
  # factor와 reference level 명확히 설정
  adj_meta$KEAP1_status <- factor(
    adj_meta$KEAP1_status,
    levels = c("WT", "MUT")
  )
  
  adj_meta$sex_group <- factor(
    adj_meta$sex_group,
    levels = c("female", "male")
  )
  
  adj_meta$stage_group <- factor(
    adj_meta$stage_group,
    levels = c("I", "II", "III", "IV")
  )
  
  adj_meta$smoking_group <- factor(
    adj_meta$smoking_group,
    levels = c("Never", "Ever")
  )
  
  # count와 metadata 순서 확인
  all(
    colnames(adj_counts) == rownames(adj_meta)
  )
  
  # 크기 확인
  dim(adj_counts)
  
  table(adj_meta$KEAP1_status)
  
  
  
  # design matrix 확인
  design_mat <- model.matrix(
    ~ sex_group + stage_group + smoking_group + KEAP1_status,
    data = adj_meta
  )
  
  # rank와 column 수가 같아야 함
  c(
    rank = qr(design_mat)$rank,
    columns = ncol(design_mat)
  )
  
  colnames(design_mat)
  
  
  
  # ============================================================
  # 119. Covariate-adjusted DESeq2
  # Sex + Stage + Smoking 보정
  # ============================================================
  
  library(DESeq2)
  
  # DESeq2에 필요한 metadata만 구성
  adj_coldata <- data.frame(
    sex_group     = adj_meta$sex_group,
    stage_group   = adj_meta$stage_group,
    smoking_group = adj_meta$smoking_group,
    KEAP1_status  = adj_meta$KEAP1_status,
    row.names     = adj_meta$barcode
  )
  
  # 순서 확인
  all(
    colnames(adj_counts) == rownames(adj_coldata)
  )
  
  # DESeq2 객체 생성
  dds_adj <- DESeqDataSetFromMatrix(
    countData = round(adj_counts),
    colData   = adj_coldata,
    design    = ~ sex_group +
      stage_group +
      smoking_group +
      KEAP1_status
  )
  
  # DESeq2 실행
  dds_adj <- DESeq(dds_adj)
  
  # ============================================================
  # 120. Adjusted KEAP1 MUT vs WT 결과
  # ============================================================
  
  res_adj <- results(
    dds_adj,
    contrast = c(
      "KEAP1_status",
      "MUT",
      "WT"
    ),
    alpha = 0.05
  )
  
  summary(res_adj)
  
  res_adj_df <- as.data.frame(res_adj)
  
  # FDR < 0.05
  sum(
    res_adj_df$padj < 0.05,
    na.rm = TRUE
  )
  
  # FDR < 0.05 + |log2FC| >= 1
  sum(
    res_adj_df$padj < 0.05 &
      abs(res_adj_df$log2FoldChange) >= 1,
    na.rm = TRUE
  )
  
  
  # ============================================================
  # 121. Unadjusted vs Adjusted DEG robustness 비교
  # ============================================================
  
  # adjusted 결과에 gene ID 추가
  res_adj_df$gene_id <- rownames(res_adj_df)
  
  # 같은 gene끼리 매칭
  adj_idx <- match(
    res_deg_df$gene_id,
    res_adj_df$gene_id
  )
  
  deg_compare <- data.frame(
    gene_id = res_deg_df$gene_id,
    gene_name = res_deg_df$gene_name,
    
    unadj_log2FC = res_deg_df$log2FoldChange,
    unadj_padj   = res_deg_df$padj,
    
    adj_log2FC = res_adj_df$log2FoldChange[adj_idx],
    adj_padj   = res_adj_df$padj[adj_idx]
  )
  
  # 둘 다 log2FC가 존재하는 gene
  compare_ok <- complete.cases(
    deg_compare[, c("unadj_log2FC", "adj_log2FC")]
  )
  
  # 전체 gene의 effect-size correlation
  cor(
    deg_compare$unadj_log2FC[compare_ok],
    deg_compare$adj_log2FC[compare_ok],
    method = "spearman"
  )
  
  cor(
    deg_compare$unadj_log2FC[compare_ok],
    deg_compare$adj_log2FC[compare_ok],
    method = "pearson"
  )
  
  
  # 기존 stringent DEG
  deg_compare$unadj_sig <- (
    !is.na(deg_compare$unadj_padj) &
      deg_compare$unadj_padj < 0.05 &
      abs(deg_compare$unadj_log2FC) >= 1
  )
  
  # adjusted stringent DEG
  deg_compare$adj_sig <- (
    !is.na(deg_compare$adj_padj) &
      deg_compare$adj_padj < 0.05 &
      abs(deg_compare$adj_log2FC) >= 1
  )
  
  # 개수
  c(
    unadjusted = sum(deg_compare$unadj_sig),
    adjusted   = sum(deg_compare$adj_sig),
    overlap    = sum(
      deg_compare$unadj_sig &
        deg_compare$adj_sig
    )
  )
  
  # overlap gene 중 방향까지 같은지
  overlap_idx <- (
    deg_compare$unadj_sig &
      deg_compare$adj_sig
  )
  
  same_direction <- sign(
    deg_compare$unadj_log2FC[overlap_idx]
  ) == sign(
    deg_compare$adj_log2FC[overlap_idx]
  )
  
  c(
    overlap_genes = sum(overlap_idx),
    same_direction = sum(same_direction),
    concordance_percent =
      mean(same_direction) * 100
  )
  
  
  # ============================================================
  # 122. Adjusted Hallmark GSEA 준비
  # ============================================================
  
  # adjusted 결과에 기존 Entrez ID 매칭
  adj_match <- match(
    res_adj_df$gene_id,
    res_deg_df$gene_id
  )
  
  res_adj_df$ENTREZID <- res_deg_df$ENTREZID[adj_match]
  
  # GSEA에 사용할 gene
  gsea_adj_df <- res_adj_df[
    !is.na(res_adj_df$stat) &
      !is.na(res_adj_df$ENTREZID),
    c(
      "gene_id",
      "ENTREZID",
      "log2FoldChange",
      "stat",
      "padj"
    )
  ]
  
  # 같은 Entrez ID가 여러 개면 |stat|가 가장 큰 것 유지
  gsea_adj_df <- gsea_adj_df[
    order(-abs(gsea_adj_df$stat)),
  ]
  
  gsea_adj_df <- gsea_adj_df[
    !duplicated(gsea_adj_df$ENTREZID),
  ]
  
  # ranked list
  gene_rank_adj <- gsea_adj_df$stat
  names(gene_rank_adj) <- as.character(gsea_adj_df$ENTREZID)
  
  gene_rank_adj <- sort(
    gene_rank_adj,
    decreasing = TRUE
  )
  
  length(gene_rank_adj)
  
  
  # ============================================================
  # 123. Covariate-adjusted Hallmark GSEA
  # ============================================================
  
  gsea_hallmark_adj <- GSEA(
    geneList      = gene_rank_adj,
    TERM2GENE     = hallmark_t2g,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    minGSSize     = 10,
    maxGSSize     = 500,
    eps           = 0,
    verbose       = FALSE,
    seed          = TRUE
  )
  
  gsea_hallmark_adj_df <- as.data.frame(
    gsea_hallmark_adj
  )
  
  # 분석 pathway 수
  nrow(gsea_hallmark_adj_df)
  
  # FDR < 0.05
  sum(
    gsea_hallmark_adj_df$p.adjust < 0.05,
    na.rm = TRUE
  )
  
  
  
  # Adjusted: MUT 쪽 Top 10
  head(
    gsea_hallmark_adj_df[
      order(-gsea_hallmark_adj_df$NES),
      c("ID", "NES", "p.adjust")
    ],
    10
  )
  
  # Adjusted: WT 쪽 Top 10
  head(
    gsea_hallmark_adj_df[
      order(gsea_hallmark_adj_df$NES),
      c("ID", "NES", "p.adjust")
    ],
    10
  )
  
  
  # ============================================================
  # 124. Unadjusted vs Adjusted Hallmark GSEA 비교
  # ============================================================
  
  gsea_compare <- merge(
    gsea_hallmark_df[
      ,
      c("ID", "NES", "p.adjust")
    ],
    gsea_hallmark_adj_df[
      ,
      c("ID", "NES", "p.adjust")
    ],
    by = "ID",
    suffixes = c("_unadj", "_adj")
  )
  
  # pathway 수
  nrow(gsea_compare)
  
  # NES correlation
  cor(
    gsea_compare$NES_unadj,
    gsea_compare$NES_adj,
    method = "spearman"
  )
  
  cor(
    gsea_compare$NES_unadj,
    gsea_compare$NES_adj,
    method = "pearson"
  )
  
  # 방향 일치율
  same_direction_gsea <- sign(
    gsea_compare$NES_unadj
  ) == sign(
    gsea_compare$NES_adj
  )
  
  c(
    pathways = nrow(gsea_compare),
    same_direction = sum(same_direction_gsea),
    concordance_percent =
      mean(same_direction_gsea) * 100
  )
  
  
  # ============================================================
  # 125. 보정 전후 방향이 달라진 Hallmark pathway 확인
  # ============================================================
  
  gsea_compare[
    sign(gsea_compare$NES_unadj) !=
      sign(gsea_compare$NES_adj),
    c(
      "ID",
      "NES_unadj",
      "p.adjust_unadj",
      "NES_adj",
      "p.adjust_adj"
    )
  ]
  
  
  # ============================================================
  # 126. Adjusted sensitivity analysis 결과 저장
  # ============================================================
  
  # Adjusted DEG 전체
  write.csv(
    res_adj_df,
    "results/DEG_adjusted_KEAP1_MUT_vs_WT.csv",
    row.names = FALSE
  )
  
  # Unadjusted vs adjusted DEG 비교
  write.csv(
    deg_compare,
    "results/DEG_unadjusted_vs_adjusted_comparison.csv",
    row.names = FALSE
  )
  
  # Adjusted Hallmark GSEA
  write.csv(
    gsea_hallmark_adj_df,
    "results/GSEA_Hallmark_adjusted_KEAP1_MUT_vs_WT.csv",
    row.names = FALSE
  )
  
  # Hallmark GSEA 보정 전후 비교
  write.csv(
    gsea_compare,
    "results/GSEA_Hallmark_unadjusted_vs_adjusted.csv",
    row.names = FALSE
  )
  
  # 확인
  list.files("results")
  
  
  x <- readLines("02_tcga_reproduction.R")
  
  cat(
    sprintf(
      "%4d: %s\n",
      2865:2885,
      x[2865:2885]
    ),
    sep = ""
  )
  
  source("~/Documents/tcga-luad-keap1/02_tcga_reproduction.R")
  
  
  x <- readLines("02_tcga_reproduction.R")
  
  cat(
    sprintf(
      "%4d: %s\n",
      5025:5045,
      x[5025:5045]
    ),
    sep = ""
  )
  
  cat(
    sprintf(
      "%4d: %s\n",
      5015:5030,
      x[5015:5030]
    ),
    sep = ""
  )
  parse(file = "~/Documents/tcga-luad-keap1/02_tcga_reproduction.R")
  
  x <- readLines("~/Documents/tcga-luad-keap1/02_tcga_reproduction.R")
  
  cat(
    sprintf(
      "%4d: %s\n",
      5075:5095,
      x[5075:5095]
    ),
    sep = ""
  )
  
  parse(file = "~/Documents/tcga-luad-keap1/02_tcga_reproduction.R")
  