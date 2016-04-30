data: explore-variables.Rmd function_library.R
	Rscript -e "knitr::knit('explore-variables.Rmd')" > explore-out.Rout 2>&1

analysis:  tmle.Rmd function_library.R
	Rscript -e "knitr::knit('tmle.Rmd')" > tmle-out.Rout 2>&1
	Rscript -e "markdown::markdownToHTML('tmle.md', 'tmle.html')"

clean:
	rm -f data/*
