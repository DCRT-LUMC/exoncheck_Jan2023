-- Initialize the database.
-- Drop any existing data and create empty tables.

DROP TABLE IF EXISTS user;
DROP TABLE IF EXISTS post;

CREATE TABLE user (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  username TEXT UNIQUE NOT NULL,
  password TEXT NOT NULL
);

CREATE TABLE post (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  author_id INTEGER NOT NULL,
  created TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  title TEXT,

  --gene
  gene_symbol TEXT,
  ENSG_gene_id TEXT,

  --variant
  NC_variant TEXT,
  strand TEXT,
  hg38_variant TEXT,
  MANE_select_NM_variant TEXT,
  MANE_select_ENST_variant TEXT,
  consequence_variant TEXT,

  --exon
  exon_number TEXT,
  total_exons TEXT,
  exon_number_interpretation TEXT,
  NC_exon TEXT,
  exon_length TEXT,
  total_protein_length TEXT,
  percentage_length TEXT,
  nearest_splice_distant TEXT,
  frame TEXT,
  MANE_select_NM_exon TEXT,
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



