library(ape)
library(xml2)
library(data.table)

fn <- list.files("sifter/hundred", full.names = TRUE, pattern = "pli$")
fn

dat <- parallel::mclapply(fn, function(f) {
  xml2::as_list(xml2::read_html(f))
}, mc.cores = 4L)


dat2 <- parallel::mclapply(dat, function(d) {
  res <- lapply(
    d$html$body$family[-1],
    function(b) {
      name   <- b$proteinname
      number <- b$proteinnumber
      go     <- b$gonumber
      moc    <- b$moc

      nann <- max(1, length(go))

      data.table(
	name   = rep(unname(name), nann),
	number = rep(unname(number), nann),
	go     = if (is.null(go)) NA else unname(as.vector(go)),
	moc    = if (is.null(moc)) NA else unname(as.vector(moc))
      )
      
    }
    )
  do.call(rbind, res)
})

famids <- sapply(dat, function(d) d$html$body$family$familyid)
famids <- unlist(famids, recursive = TRUE)
famids <- rep(famids, sapply(dat2, nrow))

dat2 <- do.call(rbind, dat2)
dat2$famid <- famids

fwrite(dat2, "sifter/hundredfamilies_annotations.csv")


