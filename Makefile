data: explore-variables.Rmd function_library.R
	Rscript -e "knitr::knit('explore-variables.Rmd')" > explore-out.Rout 2>&1

analysis:  tmle.Rmd function_library.R explore-variables.Rmd
	nice Rscript -e "knitr::knit('tmle.Rmd')" > tmle-out.Rout 2>&1
	Rscript -e "markdown::markdownToHTML('tmle.md', 'tmle.html')"

bootstrap:  tmle-boot.Rmd function_library.R explore-variables.Rmd
	nice Rscript -e "knitr::knit('tmle-boot.Rmd')" > tmle-boot.Rout 2>&1
	Rscript -e "markdown::markdownToHTML('tmle-boot.md', 'tmle-boot.html')"


clean:
	rm -f data/* cache/*
