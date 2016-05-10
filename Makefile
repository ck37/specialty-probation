data: explore-variables.Rmd function_library.R
	Rscript -e "knitr::knit('explore-variables.Rmd')" > explore-out.Rout 2>&1

tmle: data tmle.Rmd function_library.R explore-variables.Rmd
	nice Rscript -e "knitr::knit('tmle.Rmd')" > tmle-out.Rout 2>&1
	Rscript -e "markdown::markdownToHTML('tmle.md', 'tmle.html')"

bootstrap: data tmle-boot.Rmd function_library.R explore-variables.Rmd
	nice Rscript -e "knitr::knit('tmle-boot.Rmd')" > tmle-boot.Rout 2>&1
	Rscript -e "markdown::markdownToHTML('tmle-boot.md', 'tmle-boot.html')"

robust:  data lib-robustness.Rmd function_library.R explore-variables.Rmd
	nice Rscript -e "knitr::knit('lib-robustness.Rmd')" > lib-robustness.Rout 2>&1
	Rscript -e "markdown::markdownToHTML('lib-robustness.md', 'lib-robustness.html')"

clean:
	rm -rf data/* cache/* tmle_cache/*
