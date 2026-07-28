#!/usr/bin/env Rscript
# =============================================================================
# run_bayesprism.R — deconvolução BayesPrism do bulk PAAD
#
# Instalação:
#   install.packages(c("devtools","Matrix","snowfall","NMF","BiocManager"))
#   BiocManager::install("scran")
#   devtools::install_github("Danko-Lab/BayesPrism/BayesPrism")
#
# Uso:
#   Rscript run_bayesprism.R --indir bp_input --outdir bp_output \
#                            --key malignant --cores 32
# =============================================================================

suppressPackageStartupMessages({
  library(BayesPrism)
  library(Matrix)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--indir",  default = "bp_input"),
  make_option("--outdir", default = "bp_output"),
  make_option("--key",    default = "malignant",
              help = "Rótulo em cell.type.labels das células malignas."),
  make_option("--cores",  default = 16L, type = "integer"),
  make_option("--update-gibbs", default = TRUE, type = "logical")
)))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Carregar
# -----------------------------------------------------------------------------
# ATENÇÃO À ORIENTAÇÃO. BayesPrism (implementação em R) espera:
#   bk.dat : amostras (linhas) x genes (colunas)
#   sc.dat : células  (linhas) x genes (colunas)
# A documentação do gateway web usa a convenção transposta. Se você inverter,
# a função roda e devolve lixo silenciosamente.
message("[1/6] Carregando matrizes...")

bk.dat <- as.matrix(read.delim(file.path(opt$indir, "bulk.tsv"),
                               row.names = 1, check.names = FALSE))
storage.mode(bk.dat) <- "numeric"

sc.dat <- readMM(file.path(opt$indir, "sc_matrix.mtx"))          # células x genes
sc.dat <- as(sc.dat, "CsparseMatrix")                            # dgCMatrix (v2.2+)
rownames(sc.dat) <- readLines(file.path(opt$indir, "sc_cells.txt"))
colnames(sc.dat) <- readLines(file.path(opt$indir, "sc_genes.txt"))

cell.type.labels  <- readLines(file.path(opt$indir, "cell_type_labels.txt"))
cell.state.labels <- readLines(file.path(opt$indir, "cell_state_labels.txt"))

stopifnot(
  nrow(sc.dat) == length(cell.type.labels),
  nrow(sc.dat) == length(cell.state.labels),
  opt$key %in% cell.type.labels
)

message(sprintf("      bulk: %d amostras x %d genes", nrow(bk.dat), ncol(bk.dat)))
message(sprintf("      ref : %d células x %d genes | %d tipos, %d estados",
                nrow(sc.dat), ncol(sc.dat),
                length(unique(cell.type.labels)),
                length(unique(cell.state.labels))))

# -----------------------------------------------------------------------------
# 2. QC dos rótulos
# -----------------------------------------------------------------------------
# Correlação de Pearson entre os perfis. Estados com r > ~0.95 são
# praticamente colineares — o Gibbs não os separa e o theta deles fica
# instável (CV alto). Funda-os antes de rodar.
message("[2/6] QC de correlação entre perfis...")
pdf(file.path(opt$outdir, "qc_cor_phi.pdf"), width = 10, height = 9)
plot.cor.phi(input = sc.dat, input.labels = cell.state.labels,
             title = "Correlação entre cell states",
             cexRow = 0.35, cexCol = 0.35, margins = c(10, 10))
plot.cor.phi(input = sc.dat, input.labels = cell.type.labels,
             title = "Correlação entre cell types",
             cexRow = 0.8, cexCol = 0.8, margins = c(10, 10))
dev.off()

# -----------------------------------------------------------------------------
# 3. Filtragem de genes sensíveis a efeito de lote
# -----------------------------------------------------------------------------
# Esta etapa é o que mais importa para o SEU caso. O eixo técnico que
# apareceu na clusterização (lncRNAs de clone, antisenso, pseudogenes,
# KCNQ1OT1, MALAT1) é exatamente a classe de transcritos que difere entre
# plataformas de bulk e de single-cell. select.gene.type("protein_coding")
# remove esse eixo do problema por construção.
message("[3/6] Limpando genes...")
sc.stat <- plot.scRNA.outlier(
  input = sc.dat, cell.type.labels = cell.type.labels,
  species = "hs", return.raw = TRUE, pdf.prefix = file.path(opt$outdir, "qc_sc")
)
bk.stat <- plot.bulk.outlier(
  bulk.input = bk.dat, sc.input = sc.dat, cell.type.labels = cell.type.labels,
  species = "hs", return.raw = TRUE, pdf.prefix = file.path(opt$outdir, "qc_bk")
)
write.csv(sc.stat, file.path(opt$outdir, "gene_stat_sc.csv"))
write.csv(bk.stat, file.path(opt$outdir, "gene_stat_bulk.csv"))

sc.dat.filtered <- cleanup.genes(
  input          = sc.dat,
  input.type     = "count.matrix",
  species        = "hs",
  gene.group     = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
  exp.cells      = 5
)

# Restringir a protein-coding: remove ~40% dos genes e quase todo o ruído
# de anotação, com ganho consistente de acurácia em theta.
sc.dat.pc <- select.gene.type(sc.dat.filtered, gene.type = "protein_coding")
message(sprintf("      genes: %d -> %d -> %d (pc)",
                ncol(sc.dat), ncol(sc.dat.filtered), ncol(sc.dat.pc)))

# Diagnóstico da concordância entre plataformas.
pdf(file.path(opt$outdir, "qc_bulk_vs_sc.pdf"), width = 6, height = 6)
plot.bulk.vs.sc(sc.input = sc.dat.pc, bulk.input = bk.dat,
                pdf.prefix = NULL)
dev.off()

# -----------------------------------------------------------------------------
# 4. Construir o objeto prism
# -----------------------------------------------------------------------------
# `key` identifica o compartimento maligno. BayesPrism trata esse tipo de
# forma especial: o prior dele NÃO é compartilhado entre amostras (cada
# tumor tem seu próprio perfil), enquanto os tipos não-malignos usam um
# prior comum. É isso que permite recuperar Z_tumor específico da amostra.
# Se você passar key = NA, tudo vira prior compartilhado e você perde
# a heterogeneidade tumoral.
#
# outlier.cut/outlier.fraction: remove da MISTURA genes que representam
# >1% do total em >10% das amostras. Em PAAD isso normalmente pega
# hemoglobinas e algum tripsinogênio residual de ácino.
message("[4/6] Construindo prism...")
myPrism <- new.prism(
  reference         = sc.dat.pc,
  mixture           = bk.dat,
  input.type        = "count.matrix",
  cell.type.labels  = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key               = opt$key,
  outlier.cut       = 0.01,
  outlier.fraction  = 0.10
)
saveRDS(myPrism, file.path(opt$outdir, "prism_object.rds"))

# -----------------------------------------------------------------------------
# 5. Rodar
# -----------------------------------------------------------------------------
# ~180 amostras TCGA x ~15k genes x ~40 estados: da ordem de 1-3 h com 32 cores.
# Se for proibitivo, InstaPrism (github.com/humengying0907/InstaPrism)
# reimplementa o mesmo modelo em forma fechada e roda em minutos; use-o para
# varrer parâmetros e o BayesPrism para o resultado final.
message("[5/6] Rodando Gibbs (isso demora)...")
t0 <- Sys.time()
bp.res <- run.prism(prism = myPrism, n.cores = opt$cores,
                    update.gibbs = opt$`update-gibbs`)
message(sprintf("      concluído em %.1f min",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

saveRDS(bp.res, file.path(opt$outdir, "bp_res.rds"))

# -----------------------------------------------------------------------------
# 6. Extrair
# -----------------------------------------------------------------------------
message("[6/6] Extraindo posteriores...")

# theta: fração por tipo celular. which.theta="final" usa o theta atualizado
# pelo Gibbs (mais acurado); "first" é o do passo inicial.
theta <- get.fraction(bp = bp.res, which.theta = "final", state.or.type = "type")
write.csv(theta, file.path(opt$outdir, "fractions_type.csv"))

theta.state <- get.fraction(bp = bp.res, which.theta = "final",
                            state.or.type = "state")
write.csv(theta.state, file.path(opt$outdir, "fractions_state.csv"))

# CV do theta: incerteza por amostra/tipo. Frações com CV > ~0.5 não devem
# entrar em teste estatístico downstream.
theta.cv <- bp.res@posterior.theta_f@theta.cv
write.csv(theta.cv, file.path(opt$outdir, "fractions_type_cv.csv"))

# Z: expressão específica de tipo celular, escala de contagens,
# amostras x genes, uma matriz por tipo. Isto é o produto principal.
for (ct in colnames(theta)) {
  Z <- get.exp(bp = bp.res, state.or.type = "type", cell.name = ct)
  write.csv(Z, file.path(opt$outdir,
                         sprintf("Z_%s.csv", gsub("[^A-Za-z0-9]+", "_", ct))))
}

message("Pronto. Resultados em ", normalizePath(opt$outdir))
