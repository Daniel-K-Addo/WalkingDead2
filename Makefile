build: NAMESPACE
	R CMD build . 

check: install
	R CMD check --no-manual *_*.tar.gz

install: build
	R CMD INSTALL *_*.tar.gz
	
  
NAMESPACE: R/*
	Rscript --vanilla -e 'roxygen2::roxygenize()'

docs: NAMESPACE
  
.PHONY: build check install docs
