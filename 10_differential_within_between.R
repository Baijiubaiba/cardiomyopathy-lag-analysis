suppressPackageStartupMessages(library(data.table))

paired_model_fit <- function(data, terms) {
  X <- cbind(`(Intercept)`=1,as.matrix(data[,..terms]))
  fit <- lm.fit(X,data$response)
  if(fit$rank != ncol(X) || any(!is.finite(fit$coefficients))) stop("Rank-deficient or nonfinite paired design")
  list(coefficients=setNames(fit$coefficients,colnames(X)),rank=fit$rank,condition=kappa(X,exact=TRUE))
}

paired_lag_curve <- function(beta, names, prepared) {
  as.numeric(prepared$lag_spline %*% solve(prepared$lag_transform,unname(beta[names])))
}

# Supply the in-memory object returned by prepare_gbd_dlm_data().
# The original standardized lag features and their transformation are retained.
# This post-review analysis does not fit a new Bayesian model or evaluate prediction.
analyse_differential_within_between <- function(prepared, replicates=300L,
                                               seed=20260906L, minimum_success=285L) {
  ac <- copy(prepared$d[cause_ac==1L]); oc <- copy(prepared$d[cause_ac==0L])
  setorder(ac,location_id,year); setorder(oc,location_id,year)
  if(!identical(ac[,.(location_id,year)],oc[,.(location_id,year)])) stop("Unmatched outcome pairs")
  pair <- copy(ac); pair[,response:=ac$y_log-oc$y_log]
  if(nrow(pair)!=4896L || uniqueN(pair$location_id)!=204L || any(pair[,.N,by=location_id]$N!=24L)) stop("Unexpected analysis coverage")
  z <- prepared$lag_names; calendar <- prepared$time_names
  terms <- list(baseline=z,cause_sdi=c(z,"sdi_conditional"),
                cause_calendar=c(z,calendar),cause_both=c(z,"sdi_conditional",calendar))
  for(term in c(z,"sdi_conditional")) {
    pair[,(paste0("between_",term)):=mean(get(term)),by=location_id]
    pair[,(paste0("within_",term)):=get(term)-get(paste0("between_",term))]
  }
  wn <- paste0("within_",z); bn <- paste0("between_",z)
  terms$within_between <- c(wn,bn,"within_sdi_conditional","between_sdi_conditional",calendar)
  summarize_fit <- function(data,id) {
    fit <- paired_model_fit(data,terms[[id]])
    if(id=="within_between") {
      within <- mean(paired_lag_curve(fit$coefficients,wn,prepared))
      between <- mean(paired_lag_curve(fit$coefficients,bn,prepared))
      data.table(model=id,component=c("within","between","within_minus_between"),
                 cumulative=c(within,between,within-between),rank=fit$rank,condition=fit$condition)
    } else {
      data.table(model=id,component="pooled",cumulative=mean(paired_lag_curve(fit$coefficients,z,prepared)),rank=fit$rank,condition=fit$condition)
    }
  }
  full <- rbindlist(lapply(names(terms),function(id)summarize_fit(pair,id)))
  ids <- sort(unique(pair$location_id)); set.seed(seed)
  draws <- lapply(seq_len(replicates),function(b)sample(ids,length(ids),replace=TRUE))
  bootstrap <- rbindlist(lapply(seq_len(replicates),function(b) {
    data <- rbindlist(lapply(draws[[b]],function(id)pair[location_id==id]))
    rbindlist(lapply(names(terms),function(id) {
      tryCatch({
        row <- summarize_fit(data,id)
        row[,`:=`(replicate=b,status="success",message="")]
        row
      },error=function(e) {
        data.table(model=id,component=if(id=="within_between")c("within","between","within_minus_between") else "pooled",cumulative=NA_real_,rank=NA_integer_,condition=NA_real_,replicate=b,status="failed",message=conditionMessage(e))
      })
    }))
  }))
  intervals <- bootstrap[status=="success",.(bootstrap_success=.N,lower=quantile(cumulative,.025),upper=quantile(cumulative,.975)),by=.(model,component)]
  summary <- merge(full,intervals,by=c("model","component"),all.x=TRUE)
  if(anyNA(summary$bootstrap_success) || any(summary$bootstrap_success<minimum_success)) stop("Insufficient successful location resamples")
  list(summary=summary,bootstrap=bootstrap,
       boundary="Within and between coefficients are separate estimands, not additive attribution shares. All intervals condition on the supplied point estimates and original preprocessing; no causality or independent validation is established.")
}
