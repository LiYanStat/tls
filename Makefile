objects := $(wildcard R/*.R) DESCRIPTION
version := $(shell grep "Version" DESCRIPTION | sed "s/Version: //")
pkg := $(shell grep "Package" DESCRIPTION | sed "s/Package: //")
tar := $(pkg)_$(version).tar.gz
checkLog := $(pkg).Rcheck/00check.log
yr := $(shell date +"%Y")
dt := $(shell date +"%Y-%m-%d")

.PHONY: check
check: $(checkLog)

.PHONY: build
build: $(tar)

$(tar): $(objects)
	@make -s updateDate
	Rscript -e "library(methods); devtools::document();";
	R CMD build .

$(checkLog): $(tar)
	R CMD check --as-cran $(tar)

.PHONY: install
install: $(tar)
	R CMD INSTALL $(tar)

## update copyright year in COPYRIGHT, R script and date in DESCRIPTION
.PHONY: updateDate
updateDate:
	@for Rfile in R/*.R; do \
	if ! grep -q 'Copyright (C)' $$Rfile;\
	then cat COPYRIGHT $$Rfile  > tmp;\
	mv tmp $$Rfile;\
	fi;\
	sed -i '' 's/Copyright (C) 2018-[0-9]*/Copyright (C) 2018-$(yr)/' $$Rfile;\
	done;
	@echo "updating Date: $(dt)"
	@sed -i '' 's/Date: [0-9]\{4\}-[0-9]\{1,2\}-[0-9]\{1,2\}/Date: $(dt)/' DESCRIPTION
	@sed -i '' 's/Copyright (C) 2018-[0-9]\{4\}/Copyright (C) 2018-$(yr)/' COPYRIGHT

## update version number
.PHONY: newVersion
newVersion:
	@read -p "new version number: " NEWVER;\
	sed -i '' 's/^Version: [0-9]\.[0-9]\.[0-9]\.*[0-9]*[0-9]*[0-9]*[0-9]*/Version: '$$NEWVER'/' DESCRIPTION;\

.PHONY: clean
clean:
	rm -rf *~ */*~ *.Rhistroy *.tar.gz *.Rcheck/ .\#*
