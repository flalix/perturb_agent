#!/usr/bin/env Rscript
# =============================================================================
# combat_seq.R — Correção de lote em nível de CONTAGEM (TCGA vs CPTAC3),
#                compatível com BayesPrism.
#
#   Rscript combat_seq.R --counts counts_brutos.tsv --design design.csv \
#                        --out counts_combatseq.tsv
#
# Instalação:  BiocManager::install("sva")
# =============================================================================

suppressPackageStartupMessages({ library(sva); library(optparse) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts", help = "TSV genes x amostras, contagens BRUTAS inteiras"),
  make_option("--design", help = "CSV com colunas: sample, study, group"),
  make_option("--out",    default = "counts_combatseq.tsv"),
  make_option("--covars", default = "",
              help = "Colunas extras a preservar, separadas por vírgula (ex: sex,stage)")
)))

cts <- as.matrix(read.delim(opt$counts, row.names = 1, check.names = FALSE))
des <- read.csv(opt$design, row.names = 1)
des <- des[colnames(cts), , drop = FALSE]

stopifnot(all(cts == floor(cts)), min(cts) >= 0)   # ComBat_seq exige contagens

# --- Confundimento: ComBat_seq não avisa, então checamos aqui ----------------
tab <- table(des$study, des$group)
print(tab)
if (any(tab == 0))
  stop("Estudo e grupo estão confundidos (célula vazia na tabela cruzada). ",
       "Corrigir aqui remove o efeito biológico junto. Use ~ study + group ",
       "no modelo de DE em vez disto.")
if (min(tab) < 5)
  warning("Desbalanço forte (célula mínima = ", min(tab), "). ",
          "ComBat_seq vai encolher o efeito de grupo.")

# --- Filtros exigidos pelo modelo -------------------------------------------
# Gene com contagem zero em TODAS as amostras de um lote não tem como ter
# dispersão estimada naquele lote; ComBat_seq falha ou devolve NA.
ok <- Reduce(`&`, lapply(split(seq_len(ncol(cts)), des$study),
                         function(i) rowSums(cts[, i, drop = FALSE]) > 0))
message(sprintf("Genes: %d -> %d (expressos em todos os lotes)", nrow(cts), sum(ok)))
cts <- cts[ok, ]

# --- Correção ----------------------------------------------------------------
# `group` é o que impede a remoção do sinal biológico: o modelo estima o
# efeito de lote CONDICIONAL ao grupo. Omitir esse argumento é o erro mais
# comum e o mais caro — a saída fica visualmente limpa e biologicamente vazia.
#
# `covar_mod` preserva covariáveis adicionais pelo mesmo mecanismo.
covar_mod <- NULL
if (nzchar(opt$covars)) {
  vars <- strsplit(opt$covars, ",")[[1]]
  covar_mod <- model.matrix(
    as.formula(paste("~", paste(vars, collapse = "+"))), data = des
  )
}

message("Rodando ComBat_seq...")
adj <- ComBat_seq(counts    = cts,
                  batch     = as.factor(des$study),
                  group     = as.factor(des$group),
                  covar_mod = covar_mod,
                  full_mod  = TRUE)

stopifnot(all(adj == floor(adj)), min(adj) >= 0)
storage.mode(adj) <- "integer"

write.table(adj, opt$out, sep = "\t", quote = FALSE, col.names = NA)
message("Escrito: ", opt$out, sprintf("  (%d genes x %d amostras)", nrow(adj), ncol(adj)))

# --- Diagnóstico -------------------------------------------------------------
# Separação por estudo em PCA, antes e depois. Serve para inspeção, não como
# prova: PCA limpa é fácil de obter removendo sinal demais. A prova real é o
# validate_theta() depois da deconvolução.
lcpm <- function(m) log2(t(t(m) / colSums(m)) * 1e6 + 1)
pdf(sub("\\.tsv$", "_pca.pdf", opt$out), width = 11, height = 5)
par(mfrow = c(1, 2))
for (nm in c("antes", "depois")) {
  m <- if (nm == "antes") cts else adj
  p <- prcomp(t(lcpm(m)), scale. = TRUE)
  ve <- round(100 * summary(p)$importance[2, 1:2])
  plot(p$x[, 1], p$x[, 2],
       col = as.integer(as.factor(des$study)),
       pch = c(1, 19)[as.integer(as.factor(des$group))],
       xlab = sprintf("PC1 (%d%%)", ve[1]), ylab = sprintf("PC2 (%d%%)", ve[2]),
       main = paste("PCA", nm))
  legend("topright", bty = "n", cex = 0.8,
         legend = c(levels(as.factor(des$study)), levels(as.factor(des$group))),
         col = c(seq_along(levels(as.factor(des$study))), 1, 1),
         pch = c(rep(19, nlevels(as.factor(des$study))), 1, 19))
}
dev.off()
