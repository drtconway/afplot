suppressMessages(library(docopt))
"Usage: afplot (locus | density | histogram | boxplot | violin) [options] <vcf-file> [<range>...]
        afplot help

Options:
    -o FILE, --output FILE      file to write plots to [default: afplot.pdf]
    -i FIELD, --info FIELD      get the allele frequency data from the named INFO field
    -g, --genotype              compute the allele frequency data from the AD field in the sample genotypes.
    -G GENOME, --genome GENOME  genome version to use [default: hg38]
" -> doc;

opts <- docopt(doc);
print(opts);

if (opts$help) {
    cat(doc);
    quit(save="no");
}

if (is.null(opts$info) && !opts$genotype) {
    stop("please give a method for getting allele frequency information: either --info FIELD, or --genotype");
}

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(VariantAnnotation))

sanitize.locus <- function(vcf, locus) {
    orig <- locus;
  if (str_detect(locus, regex("^[^:]+$"))) {
    n <- seqlengths(seqinfo(vcf)[locus]);
    if (is.na(n)) {
      stop(sprintf("could not find length for '%s'", orig));
    }
    locus <- sprintf("%s:1-%d", locus, n);
  }
  if (!str_detect(locus, regex("^[^:]+:[0-9]+-[0-9]+"))) {
    stop(sprintf("locus must be <chrom> or <chrom>:start-end; not '%s'", orig));
  }
  return(locus);
}

select.loci <- function(vcf, loci) {
  return(subsetByOverlaps(vcf, GRanges(loci, seqinfo=seqinfo(vcf))));
}

cons.table.from.info <- function(vcf, sample.name, field.name) {
  info.fields <- colnames(info(vcf));
  if (!(field.name %in% info.fields)) {
    stop(sprintf("could not find field '%s' in fields:\n\t%s", field.name, paste(info.fields, collapse=", ")));
  }
  res <- data.table(chrom=as.character(seqnames(vcf)), pos=start(vcf), sample=sample.name, af=info(vcf)[[field.name]]);
  res <- res[!is.na(freq)];
  return(res);
}

ad.to.af <- function(ad.entry) {
  x <- unlist(ad.entry);
  tot <- sum(x);
  return((tot - x[1])/tot);
}

compute.alt.frac.from.ad <- function(vcf) {
  res <- data.table(chrom=as.character(seqnames(vcf)), pos=start(vcf))
  sample.names <- samples(header(vcf));
  for (idx in (1:length(sample.names))) {
    sample.name <- sample.names[idx];
    ad <- geno(vcf)$AD[, idx];
    res[, (sample.name) := unlist(lapply(ad, ad.to.af))];
  }
  return(res);
}

prepare.af.data <- function(tbl, sample.names) {
  res <- melt(tbl, id.vars=c("chrom", "pos"), measure.vars=sample.names, variable.name="sample", value.name="af");
  res <- res[!is.na(af)];
  return(res);
}

make.plot.by.locus <- function(tbl, locus) {
  nchrom <- length(unique(tbl$chrom));
  if (nchrom > 1) {
    stop(sprintf("locus based plots must be from a single chromosome (%d given).", nchrom));
  }
  g <- ggplot(tbl, aes(pos, af, colour=sample));
  g <- g + geom_point();
  g <- g + facet_wrap(~sample, ncol=1);
  g <- g + theme_minimal();
  g <- g + labs(title=locus, x="position", y="AF");
  g <- g + theme(legend.position = "none");
  return(g);
}

make.density.plot <- function(tbl) {
  g <- ggplot(tbl, aes(af, fill=sample));
  g <- g + geom_density();
  g <- g + facet_wrap(~sample, ncol=1);
  g <- g + theme_minimal();
  g <- g + labs(title="alt allele frequency distribution", x="AF", y="density");
  g <- g + theme(legend.position = "none");
  return(g);

}

make.histogram.plot <- function(tbl) {
  g <- ggplot(tbl, aes(af, fill=sample));
  g <- g + geom_histogram(bins=100);
  g <- g + facet_wrap(~sample, ncol=1);
  g <- g + theme_minimal();
  g <- g + labs(title="alt allele frequency distribution", x="AF", y="count");
  g <- g + theme(legend.position = "none");
  return(g);

}

make.violin.plot <- function(tbl) {
  g <- ggplot(tbl, aes(sample, af, fill=sample));
  g <- g + geom_violin();
  g <- g + theme_minimal();
  g <- g + labs(title="alt allele frequency distribution", x="sample", y="AF");
  g <- g + theme(legend.position = "none");
  return(g);

}

make.box.plot <- function(tbl, split = FALSE) {
  g <- ggplot(tbl, aes(sample, af, fill=sample));
  g <- g + geom_boxplot();
  if (split) {
    g <- g + facet_wrap(~chrom);
  }
  g <- g + theme_minimal();
  g <- g + labs(title="alt allele frequency distribution", x="sample", y="AF");
  g <- g + theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1));
  return(g);
}

vcf <- readVcf(opts$vcf_file, opts$genome);
sample.names <- samples(header(vcf));
loci <- sapply(unlist(opts$range), function(loc) sanitize.locus(vcf, loc));
if (is.list(loci) && length(loci) > 0) {
    vcf <- select.loci(vcf, loci);
}

if (opts$genotype) {
    tbl <- compute.alt.frac.from.ad(vcf);
    tbl <- prepare.af.data(tbl, sample.names);
} else {
    tbl <- cons.table.from.info(vcf, sample.names[1], opts$info);
}

if (opts$locus) {
    g <- make.plot.by.locus(tbl, "chr1");
} else if (opts$density) {
    g <- make.density.plot(tbl);
} else if (opts$histogram) {
    g <- make.histogram.plot(tbl);
} else if (opts$boxplot) {
    g <- make.box.plot(tbl);
} else if (opts$violin) {
    g <- make.violin.plot(tbl);
} else {
    stop("internal error: plot type not handled!");
}

ggsave(opts$output, g, width=8, height=5, units="in", dpi=300);
