# UNITE and SILVA as BLAST databases

Download the reference databases UNITE (ITS of fungi and other
eukaryotes) and SILVA (SSU rRNA, 16S/18S) and format them as BLAST
databases that work like the ones of the NCBI in BLASTr: each sequence
has an NCBI taxid (so
[`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md),
[`blast_lca()`](https://heronoh.github.io/BLASTr/reference/blast_lca.md)
and the summaries work), its accession as the subject id, and the
classification of UNITE or SILVA in the subject header.

- `make_unite_db()`: the UNITE "general FASTA release" (representative
  and reference sequences of the species hypotheses, SH), for all
  eukaryotes or only fungi. UNITE has no NCBI taxids: each sequence gets
  the taxid of the deepest name of its UNITE classification found in the
  NCBI taxonomy (the species, else the genus, the family...). Names in
  open nomenclature (`"Fungi_sp"`, `"..._gen_Incertae_sedis"`) are
  skipped, and homonyms are resolved by the higher ranks of the
  classification.

- `make_silva_db()`: SILVA SSU Ref NR99, with RNA converted to DNA. The
  taxid is the NCBI taxid SILVA gives each sequence (from the organism
  of the INSDC record), or, when missing, the one found from the names
  as for UNITE. The SILVA classification, which may differ from the
  organism of the record, is in the header.

The table `<name>_taxids.tsv`, next to the database, records for each
sequence the taxid, how it was found and the name it was found by.

Both databases are licensed under CC BY 4.0: cite the release (the DOI
is printed and kept in `<name>_source.txt`).

## Usage

``` r
make_unite_db(
  dir,
  set = c("eukaryotes", "fungi"),
  doi = NULL,
  singletons = FALSE,
  taxonomy_db = NULL
)

make_silva_db(dir, release = "current", taxonomy_db = NULL)
```

## Arguments

- dir:

  Directory where the database directory is created.

- set:

  UNITE: `"eukaryotes"` (all eukaryotes) or `"fungi"`.

- doi:

  UNITE: the DOI of a general FASTA release (e.g.
  `"10.15156/BIO/3671843"`); `NULL` takes the newest one of `set` listed
  at <https://unite.ut.ee/repository.php>.

- singletons:

  UNITE: use the file that also has the singletons (the "\_dev" file of
  the release).

- taxonomy_db:

  Taxonomy database to find the taxids from names (see
  [`get_tax()`](https://heronoh.github.io/BLASTr/reference/get_tax.md));
  `NULL` uses (or builds) the current one.

- release:

  SILVA: release number (e.g. `"138.2"`), or `"current"`.

## Value

The path of the database (as `db_path` of
[`parallel_blast()`](https://heronoh.github.io/BLASTr/reference/parallel_blast.md)),
invisibly.

## See also

Other database functions:
[`make_blast_db()`](https://heronoh.github.io/BLASTr/reference/make_blast_db.md),
[`ncbi_blast_dbs`](https://heronoh.github.io/BLASTr/reference/ncbi_blast_dbs.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# downloads UNITE and SILVA and builds the BLAST databases
its <- make_unite_db("~/blast_databases")
ssu <- make_silva_db("~/blast_databases")
parallel_blast(seqs, db_path = its, total_cores = 8)
} # }
```
