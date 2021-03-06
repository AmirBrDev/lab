ALTER TABLE signif_chimeras_of_iron_limitation_cl CHANGE poly_u_of_rna1 rna1_poly_u integer;
ALTER TABLE signif_chimeras_of_iron_limitation_cl CHANGE poly_u_of_rna2 rna2_poly_u integer;
ALTER TABLE signif_chimeras_of_iron_limitation_cl MODIFY COLUMN start_1 integer;
ALTER TABLE signif_chimeras_of_iron_limitation_cl MODIFY COLUMN start_2 integer;
ALTER TABLE signif_chimeras_of_iron_limitation_cl MODIFY COLUMN end_1 integer;
ALTER TABLE signif_chimeras_of_iron_limitation_cl MODIFY COLUMN end_2 integer;
ALTER TABLE signif_chimeras_of_iron_limitation_cl MODIFY COLUMN interactions integer;

ALTER TABLE signif_chimeras_of_stat_limitation_cl CHANGE poly_u_of_rna1 rna1_poly_u integer;
ALTER TABLE signif_chimeras_of_stat_limitation_cl CHANGE poly_u_of_rna2 rna2_poly_u integer;
ALTER TABLE signif_chimeras_of_stat_limitation_cl MODIFY COLUMN start_1 integer;
ALTER TABLE signif_chimeras_of_stat_limitation_cl MODIFY COLUMN start_2 integer;
ALTER TABLE signif_chimeras_of_stat_limitation_cl MODIFY COLUMN end_1 integer;
ALTER TABLE signif_chimeras_of_stat_limitation_cl MODIFY COLUMN end_2 integer;
ALTER TABLE signif_chimeras_of_stat_limitation_cl MODIFY COLUMN interactions integer;

ALTER TABLE signif_chimeras_of_log_limitation_cl CHANGE poly_u_of_rna1 rna1_poly_u integer;
ALTER TABLE signif_chimeras_of_log_limitation_cl CHANGE poly_u_of_rna2 rna2_poly_u integer;
ALTER TABLE signif_chimeras_of_log_limitation_cl MODIFY COLUMN start_1 integer;
ALTER TABLE signif_chimeras_of_log_limitation_cl MODIFY COLUMN start_2 integer;
ALTER TABLE signif_chimeras_of_log_limitation_cl MODIFY COLUMN end_1 integer;
ALTER TABLE signif_chimeras_of_log_limitation_cl MODIFY COLUMN end_2 integer;
ALTER TABLE signif_chimeras_of_log_limitation_cl MODIFY COLUMN interactions integer;