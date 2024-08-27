.PHONY: clean docs start setwd check install load render
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys
from urllib.request import pathname2url

# The input is expected to be the full HTML filename
filename = sys.argv[1]
filepath = os.path.abspath(os.path.join("./vignettes/", filename))
webbrowser.open("file://" + pathname2url(filepath))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python3 -c "$$BROWSER_PYSCRIPT"

help:
	@python3 -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

clean: ## remove all build, test, coverage and Python artifacts
	rm -f .Rbuildignore
	rm -f .Rhistory
	rm -f *.RData
	rm -f *.Rproj
	rm -rf .Rproj.user

start: ## start or restart R session
	Rscript -e "system('R')"

setwd: ## set working directory
	Rscript -e "setwd(getwd())"

docs: clean setwd ## generate docs		
	Rscript -e "devtools::document('.')"

check: clean setwd ## check package
	Rscript -e "devtools::check('.')"

install: clean setwd ## install package
	Rscript -e "devtools::install('.')"	

load: clean setwd ## load all (when developing the package)
	Rscript -e "devtools::load_all('.')"

render: ## run markdown file in /vignettes
	@read -p "Enter the name of the Rmd file (without extension): " filename; \
	Rscript -e "rmarkdown::render(paste0('./vignettes/', '$$filename', '.Rmd'))"; \
	python3 -c "$$BROWSER_PYSCRIPT" "$$filename.html"