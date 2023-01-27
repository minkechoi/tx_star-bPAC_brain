
###ZFIN_DATA
zfin_1to1 = read.table("../../../ZFIN_data.v.1.0.2.4/ensembl_1_to_1.txt",sep = "\t")
zfin_uniprot = read.table("../../../ZFIN_data.v.1.0.2.4/uniprot.txt",sep = "\t")
zfin_GenBank = read.table("../../../ZFIN_data.v.1.0.2.4/genbank.txt",sep = "\t")
zfin_RefSeq = read.table("../../../ZFIN_data.v.1.0.2.4/refseq.txt",sep = "\t")

c_1to1 = c("ZFIN ID","Gene_SO_ID","Symbol","Ensembl_ID")
c_z_uniprot = c("ZFIN ID","Gene_SO_ID","Symbol","UniProt_ID")
c_z_GenBank = c("ZFIN ID","Gene_SO_ID","Symbol","GenBank_ID")
c_z_RefSeq = c("ZFIN ID","Gene_SO_ID","Symbol","RefSeq_ID")

colnames(zfin_1to1)=c_1to1
colnames(zfin_uniprot)=c_z_uniprot
colnames(zfin_GenBank)=c_z_GenBank
colnames(zfin_RefSeq)=c_z_RefSeq
  
#
zfin_phenotypic_zh = read.table("../../../ZFIN_data.v.1.0.2.4/ortho.txt",sep = "\t",fill = TRUE )
zfin_human_ortho = read.table("../../../ZFIN_data.v.1.0.2.4/human_orthos.txt",sep = "\t",fill = TRUE )
zfin_mouse_ortho = read.table("../../../ZFIN_data.v.1.0.2.4/mouse_orthos.txt",sep = "\t",fill = TRUE )

c_human_ortho = c("ZFIN_ID","ZFIN_Symbol","ZFIN_Name","Human_Symbol","Human_Name","OMIM_ID","Gene_ID","HGNC_ID","Evidence","Pub_ID")
c_mouse_ortho = c("ZFIN_ID","ZFIN_Symbol","ZFIN_Name","Mouse_Symbol","Mouse_Name","MGI_ID","Gene_ID","Evidence","Pub_ID")
c_phenotypic_zh = c("ZFIN_ID","ZFIN_Symbol","Entrez_Zebrafish_Gene_ID","Human_Gene_Symbol","Entrez_Human_Gene_ID")

colnames(zfin_phenotypic_zh)=c_phenotypic_zh
colnames(zfin_human_ortho)=c_human_ortho
colnames(zfin_mouse_ortho)=c_mouse_ortho

#
zfin_gene_disease = read.table("../../../ZFIN_data.v.1.0.2.4/gene2DiseaseViaOrthology.txt",sep = "\t",fill = TRUE )
zfin_disease_model = read.table("../../../ZFIN_data.v.1.0.2.4/fish_model_disease.txt",sep = "\t",fill = TRUE )
zfin_exp_pheno = read.table("../../../ZFIN_data.v.1.0.2.4/gene_expression_phenotype.txt",sep = "\t",fill = TRUE )

c_Human_Disease_Models= c("Fish_ZDB_ID","Environment_ZDB_ID","is_a_model","DO_Term_ID","DO_Term_Name","Publication_ZDB_ID","PubMed_ID","Evidence_Code")
c_gene_disease= c("Zebrafish_Gene_ID","Zebrafish_Gene_Symbol","Human_Ortholog_Entrez_Gene_Id","Human_Ortholog_Symbol","DO_Term_Name","DO_Term_ID","OMIM_Term_Name","OMIM_ID","Evidence_Code","Publication")
c_exp_pheno = c("Gene_Symbol","Gene_ID","RO_Term","RO_ID","Superstructure_Name","Superstructure_ID","Substructure_Name","Substructure_ID","Phenotype_Keyword_Name","Phenotype_Keyword_ID","Phenotype_Tag","Start_Stage_Name","Start_Stage_ID","End_Stage_Name","End_Stage_ID","Assay","Assay_MMO_ID","Probe_ID","Antibody_Name","Antibody_ID","Fish_ID","Environment_ID","Figure_ID","Publication_ID","PubMed_ID")

colnames(zfin_gene_disease)=c_gene_disease
colnames(zfin_disease_model)=c_Human_Disease_Models
colnames(zfin_exp_pheno)=c_exp_pheno

#
zfin_phenotype_zh = read.table("../../../ZFIN_data.v.1.0.2.4/pheno_fish.txt",sep = "\t",fill = TRUE )
zfin_phenotype_z = read.table("../../../ZFIN_data.v.1.0.2.4/phenoGeneCleanData_fish.txt",sep = "\t",fill = TRUE )

c_pheno_zh = c("Gene_ID","Entrez_Zebrafish_Gene_ID","Entrez_Human_Gene_ID","ZFIN_Gene_Symbol","Affected_Structure_or_Process_1_subterm_OBO_ID","Affected_Structure_or_Process_1_subterm_name","Post-Composed_Relationship_ID","Post-Composed_Relationship_Name","Affected_Structure_or_Process_1_superterm_OBO_ID","Affected_Structure_or_Process_1_superterm_name","Phenotype_Keyword_OBO_ID","Phenotype_Quality","Phenotype_Tag","Affected_Structure_or_Process_2_subterm_OBO_ID","Affected_Structure_or_Process_2_subterm_name","Post-Composed_Relationship_ID","Post-Composed_Relationship_Name","Affected_Structure_or_Process_2_superterm_OBO_ID","Affected_Structure_or_Process_2_superterm_name")
c_pheno_z = c("ID","Gene_Symbol","Gene_ID","Affected_Structure_or_Process_1_subterm_ID","Affected_Structure_or_Process_1_subterm_Name","Post-composed_Relationship_ID","Post-composed_Relationship_Name","Affected_Structure_or_Process_1_superterm_ID","Affected_Structure_or_Process_1_superterm_Name","Phenotype_Keyword_ID","Phenotype_Keyword_Name","Phenotype_Tag","Affected_Structure_or_Process_2_subterm_ID","Affected_Structure_or_Process_2_subterm_name","Post-composed_Relationship_(rel)_ID","Post-composed_Relationship_(rel)_Name","Affected_Structure_or_Process_2_superterm_ID","Affected_Structure_or_Process_2_superterm_name","Fish_ID","Fish_Display_Name","Start_Stage_ID","End_Stage_ID","Fish_Environment_ID","Publication_ID","Figure_ID")

colnames(zfin_phenotype_zh)=c_pheno_zh
colnames(zfin_phenotype_z)=c_pheno_z




