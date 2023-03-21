build: NAMESPACE
		R CMD build . && R CMD check --no-manual *_*.tar.gz
  
NAMESPACE: R/*
		Rscript --vanilla -e 'roxygen2::roxygenize()'
  
.PHONY: build
