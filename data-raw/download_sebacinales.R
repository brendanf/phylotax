# download reannotated RDP training set sequences from Sebacinales
# the ones that are annotated to genus get kept as references.
# the rest have their names removed and are used for query sequences.
# also grab Dacrymycetes as an outgroup.

library(magrittr)

rdpfile <- tempfile(fileext <- ".zip")
utils::download.file(
  paste0('https://sourceforge.net/projects/rdp-classifier/files/',
         'RDP_Classifier_TrainingData/',
         'RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata.zip'),
  destfile = rdpfile)
utils::unzip(rdpfile, exdir = tempdir())
rdp <- Biostrings::readBStringSet(file.path(
  tempdir(),
  "RDPClassifier_fungiLSU_trainsetNo11_rawtrainingdata",
  "fungiLSU_train_012014.fa"
))

rdp <- chartr("Uu", "Tt", rdp)
rdp <- Biostrings::DNAStringSet(rdp)

sebacinales <- rdp[grepl("Sebacinales", names(rdp))]
sebacinales_accno <-
  vapply(strsplit(names(sebacinales), "\t", fixed = TRUE), dplyr::first, "")
sebacinales_oldnames <- gsub("\\t.+;", "=", names(sebacinales))
names(sebacinales_oldnames) <- sebacinales_accno
names(sebacinales) <- sebacinales_accno

# an ENTREZ key is probably necessary to avoid ketting kicked off the NCBI server.
if (file.exists("ENTREZ_KEY")) {
  Sys.setenv(ENTREZ_KEY = readLines("ENTREZ_KEY"))
}
seb_taxa <-
  taxa::lookup_tax_data(sebacinales_accno, type = "seq_id", database = "ncbi")

seb_classifications <-
  taxa::get_data(seb_taxa, "query_data")[[1]] %>%
  tibble::enframe(name = "taxon_id", value = "accno") %>%
  dplyr::left_join(
    taxa::taxonomy_table(seb_taxa, add_id_col = TRUE),
    by = "taxon_id"
  ) %>%
  # the taxonomy table does not have ranks, so we need to figure them out.
  dplyr::rename_if(~ "Fungi" %in% ., ~"kingdom") %>%
  dplyr::rename_if(~ "Basidiomycota" %in% ., ~"phylum") %>%
  dplyr::rename_if(~ "Agaricomycetes" %in% ., ~"class") %>%
  dplyr::rename_if(~ "Sebacinales" %in% ., ~"order") %>%
  dplyr::rename_if(~ "Sebacinaceae" %in% ., ~"family") %>%
  dplyr::rename_if(~ "Sebacina" %in% ., ~"genus") %>%
  dplyr::rename_if(~ "Serendipita vermifera" %in% ., ~"species")

# reference sequences must actually belong to a known species!
# this means that NCBI will update the taxonomy if the species is moved
# For instance, sequences annotated as "Sebacina vermifera" are now called
# Serendipita vermifera, and located in Serendipitaceae.
# sequences annotated as "Sebacina cf. vermifera" are still called Sebacina
# cf. vermifera, and are still in Sebacinaceae!
seb_refs <-
  dplyr::filter(
    seb_classifications,
    !grepl("(unclassified|uncultured)", genus),
    !grepl("(unclassified|unidentified|uncultured|sp\\.|aff\\.|cf\\.)", species)
  ) %>%
  dplyr::transmute(
    accno = accno,
    sintax = sprintf("tax=k:%s,p:%s,c:%s,o:%s,f:%s,g:%s",
                     kingdom, phylum, class, order, family, genus),
    dada2 = paste(kingdom, phylum, class, order, family, genus, sep = ";")
  )

# write the reference databases
refseqs <- sebacinales[seb_refs$accno]
Biostrings::writeXStringSet(
  magrittr::set_names(refseqs, seb_refs$sintax),
  here::here("inst", "extdata", "sebacinales.sintax.fasta.gz"),
  compress = TRUE
)
Biostrings::writeXStringSet(
  magrittr::set_names(refseqs, seb_refs$dada2),
  here::here("inst", "extdata", "sebacinales.dada2.fasta.gz"),
  compress = TRUE
)

seb_unknowns <-
  dplyr::anti_join(seb_classifications, seb_refs, by = "accno")

# write the unknowns
Biostrings::writeXStringSet(
  sebacinales[seb_unknowns$accno],
)

# Now pull out the Dacrymycetes to use as an outgroup
# just give accno, family, and genus
dacrymycetes <- rdp[grepl("Dacrymycetes", names(rdp))]
names(dacrymycetes) <- gsub("\\t.+Dacrymycetales;", "=", names(dacrymycetes))
Biostrings::writeXStringSet(
  dacrymycetes,
  here::here("inst", "extdata", "dacrymycetes.fasta.gz"),
  compress = TRUE
)

usethis::use_data(sebacinales_oldnames)
