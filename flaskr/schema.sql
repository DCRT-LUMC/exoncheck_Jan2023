-- Initialize the database.
-- Drop any existing data and create empty tables.

DROP TABLE IF EXISTS user;
-- post refers to the variant
DROP TABLE IF EXISTS variant;
DROP TABLE IF EXISTS exon;

-- Each element is a column in the table, each entry will be a row that behaves as a dict in Python.
-- user is a table to store users and their login details.
CREATE TABLE user (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  username TEXT UNIQUE NOT NULL,
  password TEXT NOT NULL
);

-- Contains data collected for the input variants by the user
CREATE TABLE variant (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  author_id INTEGER NOT NULL,
  created TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  input_variant TEXT,

  --if input is a gene (i.e. title1)
  NM_id TEXT,
  NC_id TEXT,
  in_frame_exons TEXT,
  out_of_frame_exons TEXT,

  --If variant is NP_, create a column to store json object of alternative NM_variants
  NM_variants TEXT,

  -- hg-version
  hg TEXT,

  --gene
  gene_symbol TEXT,
  ENSG_gene_id TEXT,
  elig TEXT,

  --protein
  prot_name TEXT,
  short_name TEXT,

  --variant
  NC_variant TEXT,
  strand TEXT,
  hg_variant TEXT,
  MANE_select_NM_variant TEXT,
  MANE_select_ENST_variant TEXT,
  consequence_variant TEXT,

  --exon
  exon_number TEXT,
  total_exons TEXT,
  exon_number_interpretation TEXT,
  coding_exons TEXT,
  NC_exon TEXT,
  exon_length TEXT,
  total_protein_length TEXT,
  percentage_length TEXT,
  length_condition TEXT,
  nearest_splice_distant TEXT,
  nearest_end TEXT,
  frame TEXT,
  splice_dist_interpretation TEXT,
  MANE_select_NM_exon TEXT,
  r_exon_skip TEXT,
  MANE_select_ENST_exon TEXT,
  consequence_skipping TEXT,

  --domain
  uniprot_link TEXT,
  domain_info TEXT,

  --expression
  gtex_link TEXT,
  expression_brain TEXT,
  expression_fibroblasts TEXT,
  expression_tibial_nerve TEXT,
  expression_blood TEXT,
  expression_transformed_lymphocytes TEXT,
  expression_periphery_retina TEXT,
  expression_center_retina TEXT,

  --lovd
  lovd_output TEXT,

  --identifiers
  omim_id TEXT,

  --links
  omim_link TEXT,
  gnomAD_link TEXT,
  decipher_link TEXT,
  clinvar_link TEXT,

  FOREIGN KEY (author_id) REFERENCES user (id)
);

-- Contains data collected for the selected exon corresponding to the gene input by the user
CREATE TABLE exon (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  author_id INTEGER NOT NULL,
  created TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  input_gene TEXT,

  --gene
  NM_id TEXT,
  NC_id TEXT,
  elig TEXT,

  --exon
  exon_number TEXT,
  total_exons TEXT,
  exon_number_interpretation TEXT,
  coding_exons TEXT,
  NC_exon TEXT,
  exon_length TEXT,
  total_protein_length TEXT,
  percentage_length TEXT,
  length_condition TEXT,
  frame TEXT,
  MANE_select_NM_exon TEXT,
  r_exon_skip TEXT,
  consequence_skipping TEXT,

  FOREIGN KEY (author_id) REFERENCES user (id)
);
