
suppressPackageStartupMessages({
  library(relaimpo)
  library(R.matlab)
  library(foreach)
  library(doParallel)
})

args = commandArgs(trailingOnly=FALSE)
script = sub('^--file=', '', grep('^--file=', args, value=TRUE)[1])
project = normalizePath(file.path(dirname(script), '..'), mustWork=TRUE)
indir = file.path(project, 'results', 'multiplecorr')

num_cores = 32
nperms = 1000L

# ===== Select tasks =====
run_tasks = c('nm', 'dv_1')
# run_tasks = c('nm')
# run_tasks = c('nm', 'dv_1')
# ========================

stopifnot(
  length(run_tasks) > 0L,
  all(run_tasks %in% c('nm', 'dv_1')),
  dir.exists(indir)
)

cat('Input/output:', indir, '\n')
cat('Selected tasks:', paste(run_tasks, collapse=', '), '\n')

cl = parallel::makeCluster(max(1L, as.integer(num_cores)))
registerDoParallel(cl)

run_one = function(input, output) {
  dat = readMat(input)
  X = dat$x
  y = as.numeric(dat$y)
  yp = dat$yperm

  stopifnot(
    nrow(X) == length(y),
    nrow(yp) == length(y),
    ncol(yp) == nperms
  )

  xdf = as.data.frame(X)
  ic_data = cbind(y, yp)

  model = lm(y ~ ., data=xdf)
  predicted_y = as.numeric(predict(model, newdata=xdf))
  observed_y = y

  result = foreach(
    p = seq_len(nperms + 1L),
    .combine = 'cbind',
    .packages = 'relaimpo',
    .export = c('xdf', 'ic_data')
  ) %dopar% {
    current_y = ic_data[, p]
    m = lm(current_y ~ ., data=xdf)
    importance = calc.relimp(m, type='genizi', rela=TRUE)
    c(summary(m)$adj.r.squared, as.numeric(importance@genizi))
  }

  adj_r2 = as.numeric(result[1, ])
  weis = t(result[-1, , drop=FALSE])

  writeMat(
    output,
    adj_r2=adj_r2,
    weis=weis,
    predicted_y=predicted_y,
    observed_y=observed_y
  )

  cat(
    'Saved ', basename(output),
    ' adjusted R2=', format(adj_r2[1], digits=9),
    '\n', sep=''
  )
}

tryCatch({
  for(task in run_tasks) {
    nk = if(task == 'nm') 75L else 74L

    for(md in c('ct', 'nt', 'mc')) {
      for(k in seq_len(nk)) {

        infile = file.path(
          indir, sprintf('%s_%s_%d.mat', task, md, k)
        )
        outfile = file.path(
          indir, sprintf('results_%s_%s_%d.mat', task, md, k)
        )

        if(!file.exists(infile)) {
          stop('Missing S2 input ', infile)
        }

        cat('Running ', task, ' ', md, ' ', k, '\n')
        run_one(infile, outfile)
      }
    }
  }
}, finally={
  parallel::stopCluster(cl)
  registerDoSEQ()
})

cat('Completed:', paste(run_tasks, collapse=', '), '\n')

