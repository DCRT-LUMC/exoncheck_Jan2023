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
  NC_exon TEXT,
  exon_length TEXT,
  total_protein_length TEXT,
  percentage_length TEXT,
  frame TEXT,
  MANE_select_NM_exon TEXT,
  MANE_select_ENST_exon TEXT,
  consequence_skipping TEXT,

  --domain
  uniprot_link TEXT,
  domain_info TEXT,

  --expression
  gtex_link TEXT,
  expression_eye TEXT,
  expression_brain TEXT,
  expression_fibroblasts TEXT,
  expression_tibial_nerve TEXT,
  expression_blood TEXT,
  expression_transformed_lymphocytes TEXT,

  --lovd
  no_exact_lovd_matches TEXT,
  no_partial_lovd_matches TEXT,
  exact_lovd_match_link TEXT,
  lovd_link TEXT,
  partial_lovd_match1 TEXT,
  partial_lovd_match2 TEXT,
  partial_lovd_match3 TEXT,
  partial_lovd_match4 TEXT,
  partial_lovd_match5 TEXT,
  partial_lovd_match6 TEXT,
  partial_lovd_match7 TEXT,
  partial_lovd_match8 TEXT,
  partial_lovd_match9 TEXT,
  partial_lovd_match10 TEXT,

    --identifiers
  omim_id TEXT,

    --links
  omim_link TEXT,
  gnomAD_link TEXT,
  decipher_link TEXT,
  clinvar_link TEXT,

  FOREIGN KEY (author_id) REFERENCES user (id)
);



