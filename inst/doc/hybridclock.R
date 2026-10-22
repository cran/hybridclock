## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 5)

## -----------------------------------------------------------------------------
library(hybridclock)
ex <- hybridclock_example()
table(ex$groups)
ex$truth

## -----------------------------------------------------------------------------
mk <- diagnostic_markers(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
sh <- shared_loci(ex$genotypes, ex$groups, "species1", "species2", locus = ex$locus)
c(diagnostic = nrow(mk), shared = nrow(sh))
sc <- hybrid_scores(ex$genotypes, mk)
plot_triangle(sc, ex$groups)

## -----------------------------------------------------------------------------
cl <- classify_hybrids(sc)
table(ex$groups, cl$class)

## -----------------------------------------------------------------------------
ref <- build_reference(n_sims = 300, generations = 20, breeders = c(20, 150),
                       n_markers = nrow(mk), shared_freq = sh$freq, sample_size = 20)
clock <- train_clock(ref, committees = 5)
clock

## -----------------------------------------------------------------------------
hyb <- ex$groups %in% c("site_A", "site_B")
est <- estimate_history(clock, sc[hyb, ], ex$groups[hyb],
                        diag = orient_genotypes(ex$genotypes, mk)[hyb, ],
                        shared = ex$genotypes[hyb, sh$column])
est[, c("population", "method", "hybrid_age_both", "hybrid_age_both_lower", "hybrid_age_both_upper",
        "hybrid_age_rare", "backcrosses")]

## -----------------------------------------------------------------------------
markers_needed(1:8)

