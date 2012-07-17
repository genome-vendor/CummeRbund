-- Creator:       MySQL Workbench 5.2.33/ExportSQLite plugin 2009.12.02
-- Author:        Loyal Goff
-- Caption:       New Model
-- Project:       Name of the project
-- Changed:       2012-05-01 15:20
-- Created:       2011-05-02 12:52
PRAGMA foreign_keys = OFF;

-- Schema: cuffData
BEGIN;
CREATE TABLE "genes"(
  "gene_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT
);
CREATE INDEX "genes.gsn_index" ON "genes"("gene_short_name");
CREATE INDEX "genes.cc_index" ON "genes"("class_code");
CREATE TABLE "biasData"(
  "biasData_id" INTEGER PRIMARY KEY NOT NULL
);
CREATE TABLE "samples"(
  "sample_index" INTEGER NOT NULL,
  "sample_name" VARCHAR(45) PRIMARY KEY NOT NULL
);
CREATE TABLE "TSS"(
  "TSS_group_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_id" VARCHAR(45) NOT NULL,
  "gene_short_name" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_TSS_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "TSS.fk_TSS_genes1" ON "TSS"("gene_id");
CREATE TABLE "TSSData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_TSSData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "TSSData.fk_TSSData_TSS1" ON "TSSData"("TSS_group_id");
CREATE INDEX "TSSData.fk_TSSData_samples1" ON "TSSData"("sample_name");
CREATE TABLE "CDS"(
  "CDS_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "gene_id" VARCHAR(45),
  "TSS_group_id" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  "gene_short_name" VARCHAR(45),
  CONSTRAINT "fk_CDS_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_CDS_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
CREATE INDEX "CDS.fk_CDS_genes1" ON "CDS"("gene_id");
CREATE INDEX "CDS.fk_CDS_TSS1" ON "CDS"("TSS_group_id");
CREATE TABLE "CDSData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_CDSData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "CDSData.fk_CDSData_CDS1" ON "CDSData"("CDS_id");
CREATE INDEX "CDSData.fk_CDSData_samples1" ON "CDSData"("sample_name");
CREATE TABLE "splicingDiffData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_splicingDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_splicingDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_splicingDiffData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_splicingDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples1" ON "splicingDiffData"("sample_1");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_samples2" ON "splicingDiffData"("sample_2");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_TSS1" ON "splicingDiffData"("TSS_group_id");
CREATE INDEX "splicingDiffData.fk_splicingDiffData_genes1" ON "splicingDiffData"("gene_id");
CREATE TABLE "TSSExpDiffData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "ln_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_TSSExpDiffData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_TSSExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_TSS1" ON "TSSExpDiffData"("TSS_group_id");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples1" ON "TSSExpDiffData"("sample_1");
CREATE INDEX "TSSExpDiffData.fk_TSSExpDiffData_samples2" ON "TSSExpDiffData"("sample_2");
CREATE INDEX "TSSExpDiffData.TSSExpDiffData_sig_index" ON "TSSExpDiffData"("test_stat","p_value","q_value","significant");
CREATE TABLE "CDSDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_CDSDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples1" ON "CDSDiffData"("sample_1");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_samples2" ON "CDSDiffData"("sample_2");
CREATE INDEX "CDSDiffData.fk_CDSDiffData_genes1" ON "CDSDiffData"("gene_id");
CREATE TABLE "CDSExpDiffData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "ln_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_CDSExpDiffData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_CDSExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_CDS1" ON "CDSExpDiffData"("CDS_id");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples1" ON "CDSExpDiffData"("sample_1");
CREATE INDEX "CDSExpDiffData.fk_CDSExpDiffData_samples2" ON "CDSExpDiffData"("sample_2");
CREATE INDEX "CDSExpDiffData.CDSExpDiffData_sig_index" ON "CDSExpDiffData"("test_stat","p_value","q_value","significant");
CREATE TABLE "promoterDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "JS_dist" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_promoterDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_promoterDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_promoterDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "promoterDiffData.fk_promoterDiffData_genes1" ON "promoterDiffData"("gene_id");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples1" ON "promoterDiffData"("sample_1");
CREATE INDEX "promoterDiffData.fk_promoterDiffData_samples2" ON "promoterDiffData"("sample_2");
CREATE TABLE "geneFeatures"(
  "gene_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_geneFeatures_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "geneFeatures.fk_geneFeatures_genes1" ON "geneFeatures"("gene_id");
CREATE TABLE "TSSFeatures"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_TSSFeatures_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id")
);
CREATE INDEX "TSSFeatures.fk_TSSFeatures_TSS1" ON "TSSFeatures"("TSS_group_id");
CREATE TABLE "CDSFeatures"(
  "CDS_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_CDSFeatures_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id")
);
CREATE INDEX "CDSFeatures.fk_CDSFeatures_CDS1" ON "CDSFeatures"("CDS_id");
CREATE TABLE "model_transcripts"(
  "model_transcript_id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL
);
CREATE TABLE "geneCount"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "uncertainty" FLOAT,
  "dispersion" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_geneCount_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "geneCount.fk_geneCount_samples1" ON "geneCount"("sample_name");
CREATE INDEX "geneCount.fk_geneCount_genes1" ON "geneCount"("gene_id");
CREATE TABLE "CDSCount"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "uncertainty" FLOAT,
  "dispersion" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_CDSCount_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "CDSCount.fk_CDSCount_CDS1" ON "CDSCount"("CDS_id");
CREATE INDEX "CDSCount.fk_CDSCount_samples1" ON "CDSCount"("sample_name");
CREATE TABLE "TSSCount"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "uncertainty" FLOAT,
  "dispersion" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_TSSCount_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "TSSCount.fk_TSSCount_TSS1" ON "TSSCount"("TSS_group_id");
CREATE INDEX "TSSCount.fk_TSSCount_samples1" ON "TSSCount"("sample_name");
CREATE TABLE "replicates"(
  "file" INTEGER NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "replicate" VARCHAR(45),
  "rep_name" VARCHAR(45) PRIMARY KEY NOT NULL,
  "total_mass" FLOAT,
  "norm_mass" FLOAT,
  "internal_scale" FLOAT,
  "external_scale" FLOAT,
  CONSTRAINT "fk_replicates_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "replicates.fk_replicates_samples1" ON "replicates"("sample_name");
CREATE TABLE "geneReplicateData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "replicate" INTEGER,
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneData_genes10"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneReplicateData_replicates1"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name"),
  CONSTRAINT "fk_geneReplicateData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "geneReplicateData.fk_geneReplicateData_genes1" ON "geneReplicateData"("gene_id");
CREATE INDEX "geneReplicateData.fk_geneReplicateData_replicates1" ON "geneReplicateData"("rep_name");
CREATE INDEX "geneReplicateData.fk_geneReplicateData_samples1" ON "geneReplicateData"("sample_name");
CREATE TABLE "CDSReplicateData"(
  "CDS_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "replicate" INTEGER,
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneReplicateData_replicates100"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name"),
  CONSTRAINT "fk_CDSReplicateData_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_CDSReplicateData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "CDSReplicateData.fk_CDSReplicateData_replicates1" ON "CDSReplicateData"("rep_name");
CREATE INDEX "CDSReplicateData.fk_CDSReplicateData_CDS1" ON "CDSReplicateData"("CDS_id");
CREATE INDEX "CDSReplicateData.fk_CDSReplicateData_samples1" ON "CDSReplicateData"("sample_name");
CREATE TABLE "TSSReplicateData"(
  "TSS_group_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "replicate" VARCHAR(45),
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneReplicateData_replicates10000"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name"),
  CONSTRAINT "fk_TSSReplicateData_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_TSSReplicateData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "TSSReplicateData.fk_TSSReplicateData_replicates1" ON "TSSReplicateData"("rep_name");
CREATE INDEX "TSSReplicateData.fk_TSSReplicateData_TSS1" ON "TSSReplicateData"("TSS_group_id");
CREATE INDEX "TSSReplicateData.fk_TSSReplicateData_samples1" ON "TSSReplicateData"("sample_name");
CREATE TABLE "runInfo"(
  "param" VARCHAR(45),
  "value" TEXT
);
CREATE TABLE "geneData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_geneData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "geneData.fk_geneData_genes1" ON "geneData"("gene_id");
CREATE INDEX "geneData.fk_geneData_samples1" ON "geneData"("sample_name");
CREATE TABLE "phenoData"(
  "sample_name" VARCHAR(45) NOT NULL,
  "parameter" VARCHAR(45) NOT NULL,
  "value" VARCHAR(45),
  CONSTRAINT "fk_phenoData_samples"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "phenoData.fk_phenoData_samples" ON "phenoData"("sample_name");
CREATE TABLE "geneExpDiffData"(
  "gene_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "ln_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_geneExpDiffData_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_geneExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_geneExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_genes1" ON "geneExpDiffData"("gene_id");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples1" ON "geneExpDiffData"("sample_1");
CREATE INDEX "geneExpDiffData.fk_geneExpDiffData_samples2" ON "geneExpDiffData"("sample_2");
CREATE INDEX "geneExpDiffData.geneExpDiff_status_index" ON "geneExpDiffData"("status");
CREATE INDEX "geneExpDiffData.geneExpDiff_sig_index" ON "geneExpDiffData"("significant","p_value","q_value","test_stat");
CREATE TABLE "isoforms"(
  "isoform_id" VARCHAR(45) PRIMARY KEY NOT NULL,
  "gene_id" VARCHAR(45),
  "CDS_id" VARCHAR(45),
  "gene_short_name" VARCHAR(45),
  "TSS_group_id" VARCHAR(45),
  "class_code" VARCHAR(45),
  "nearest_ref_id" VARCHAR(45),
  "locus" VARCHAR(45),
  "length" INTEGER,
  "coverage" FLOAT,
  CONSTRAINT "fk_isoforms_TSS1"
    FOREIGN KEY("TSS_group_id")
    REFERENCES "TSS"("TSS_group_id"),
  CONSTRAINT "fk_isoforms_CDS1"
    FOREIGN KEY("CDS_id")
    REFERENCES "CDS"("CDS_id"),
  CONSTRAINT "fk_isoforms_genes1"
    FOREIGN KEY("gene_id")
    REFERENCES "genes"("gene_id")
);
CREATE INDEX "isoforms.fk_isoforms_TSS1" ON "isoforms"("TSS_group_id");
CREATE INDEX "isoforms.fk_isoforms_CDS1" ON "isoforms"("CDS_id");
CREATE INDEX "isoforms.fk_isoforms_genes1" ON "isoforms"("gene_id");
CREATE TABLE "isoformData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "fpkm" FLOAT NOT NULL,
  "conf_hi" FLOAT,
  "conf_lo" FLOAT,
  "quant_status" VARCHAR(45),
  CONSTRAINT "fk_isoformData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_isoformData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
CREATE INDEX "isoformData.fk_isoformData_samples1" ON "isoformData"("sample_name");
CREATE INDEX "isoformData.fk_isoformData_isoforms1" ON "isoformData"("isoform_id");
CREATE TABLE "isoformExpDiffData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_1" VARCHAR(45) NOT NULL,
  "sample_2" VARCHAR(45) NOT NULL,
  "status" VARCHAR(45),
  "value_1" FLOAT,
  "value_2" FLOAT,
  "ln_fold_change" FLOAT,
  "test_stat" FLOAT,
  "p_value" FLOAT,
  "q_value" FLOAT,
  "significant" VARCHAR(45),
  CONSTRAINT "fk_isoformExpDiffData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id"),
  CONSTRAINT "fk_isoformExpDiffData_samples1"
    FOREIGN KEY("sample_1")
    REFERENCES "samples"("sample_name"),
  CONSTRAINT "fk_isoformExpDiffData_samples2"
    FOREIGN KEY("sample_2")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_isoforms1" ON "isoformExpDiffData"("isoform_id");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples1" ON "isoformExpDiffData"("sample_1");
CREATE INDEX "isoformExpDiffData.fk_isoformExpDiffData_samples2" ON "isoformExpDiffData"("sample_2");
CREATE INDEX "isoformExpDiffData.isoformExpDiffData_sig_index" ON "isoformExpDiffData"("test_stat","p_value","q_value","significant");
CREATE TABLE "isoformFeatures"(
  "isoform_id" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_isoformFeatures_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
CREATE INDEX "isoformFeatures.fk_isoformFeatures_isoforms1" ON "isoformFeatures"("isoform_id");
CREATE TABLE "features"(
--   GTF Features (all lines/records from reference .gtf file)
  "feature_id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  "genes_gene_id" VARCHAR(45) NOT NULL,
  "isoforms_isoform_id" VARCHAR(45) NOT NULL,
  "seqname" VARCHAR(45) NOT NULL,
  "source" VARCHAR(45) NOT NULL,
  "type_id" INTEGER,
  "start" INTEGER,
  "end" INTEGER,
  "score" FLOAT,
  "strand" VARCHAR(45),
  "frame" VARCHAR(45),
  CONSTRAINT "fk_features_genes1"
    FOREIGN KEY("genes_gene_id")
    REFERENCES "genes"("gene_id"),
  CONSTRAINT "fk_features_isoforms1"
    FOREIGN KEY("isoforms_isoform_id")
    REFERENCES "isoforms"("isoform_id")
);
CREATE INDEX "features.features_seqname_index" ON "features"("seqname");
CREATE INDEX "features.features_type_index" ON "features"("type_id");
CREATE INDEX "features.features_strand_index" ON "features"("strand");
CREATE INDEX "features.features_start_end_index" ON "features"("start","end");
CREATE INDEX "features.fk_features_genes1" ON "features"("genes_gene_id");
CREATE INDEX "features.fk_features_isoforms1" ON "features"("isoforms_isoform_id");
CREATE TABLE "attributes"(
  "attribute_lookup_id" INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
  "feature_id" INTEGER NOT NULL,
  "attribute" VARCHAR(45) NOT NULL,
  "value" VARCHAR(45) NOT NULL,
  CONSTRAINT "fk_attribute_lookup_features1"
    FOREIGN KEY("feature_id")
    REFERENCES "features"("feature_id")
);
CREATE INDEX "attributes.fk_attributes_feature_id" ON "attributes"("feature_id");
CREATE INDEX "attributes.attributes_attribute_index" ON "attributes"("attribute");
CREATE INDEX "attributes.attributes_value_index" ON "attributes"("value");
CREATE TABLE "isoformCount"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "count" FLOAT,
  "variance" FLOAT,
  "uncertainty" FLOAT,
  "dispersion" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_isoformCount_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id"),
  CONSTRAINT "fk_isoformCount_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "isoformCount.fk_isoformCount_isoforms1" ON "isoformCount"("isoform_id");
CREATE INDEX "isoformCount.fk_isoformCount_samples1" ON "isoformCount"("sample_name");
CREATE TABLE "isoformReplicateData"(
  "isoform_id" VARCHAR(45) NOT NULL,
  "sample_name" VARCHAR(45) NOT NULL,
  "replicate" INTEGER,
  "rep_name" VARCHAR(45) NOT NULL,
  "raw_frags" FLOAT,
  "internal_scaled_frags" FLOAT,
  "external_scaled_frags" FLOAT,
  "fpkm" FLOAT,
  "effective_length" FLOAT,
  "status" VARCHAR(45),
  CONSTRAINT "fk_geneReplicateData_replicates10"
    FOREIGN KEY("rep_name")
    REFERENCES "replicates"("rep_name"),
  CONSTRAINT "fk_isoformReplicateData_isoforms1"
    FOREIGN KEY("isoform_id")
    REFERENCES "isoforms"("isoform_id"),
  CONSTRAINT "fk_isoformReplicateData_samples1"
    FOREIGN KEY("sample_name")
    REFERENCES "samples"("sample_name")
);
CREATE INDEX "isoformReplicateData.fk_isoformReplicateData_replicates1" ON "isoformReplicateData"("rep_name");
CREATE INDEX "isoformReplicateData.fk_isoformReplicateData_isoforms1" ON "isoformReplicateData"("isoform_id");
CREATE INDEX "isoformReplicateData.fk_isoformReplicateData_samples1" ON "isoformReplicateData"("sample_name");
COMMIT;
