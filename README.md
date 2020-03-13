# drugSimDB
Drug Similarity information

## Introduction

R is a powerful programming language for statistical computing with many packages and tools.

The goal of this article is to outline some ways to combine Shiny & R Markdown.

We'll cover the main approaches below:

1. **Interactive documents**

2. **Including Markdown content in a Shiny app**

3. **Generating downloadable reports from Shiny apps**

4. **Embedding Shiny apps via their URLs in R Markdown output**

5. **Shiny components in rtutor**

## Interactive Documents - Overview
- R Markdown leverages Shiny at its core to make this possible. Add runtime: shiny to the document’s YAML header.
- [http://rmarkdown.rstudio.com/authoring_shiny.html](http://rmarkdown.rstudio.com/authoring_shiny.html)
- [http://shiny.rstudio.com/articles/interactive-docs.html](http://shiny.rstudio.com/articles/interactive-docs.html)
- [https://blog.rstudio.org/2014/06/19/interactive-documents-an-incredibly-easy-way-to-use-shiny/](https://blog.rstudio.org/2014/06/19/interactive-documents-an-incredibly-easy-way-to-use-shiny/)
- [http://rmarkdown.rstudio.com/flexdashboard/shiny.html](http://rmarkdown.rstudio.com/flexdashboard/shiny.html)
- As of Spring 17 - Shiny apps can’t currently be embedded into R Notebooks.

## Interactive Documents - an R Markdown file that contains Shiny components and can be an app or used within a report

- Types of Interactive Documents:
- 1. RMD Interactive Document/Report with Shiny components 
- Example: Demo by using File > New > RMarkdown > Shiny > Shiny Document
- 2. HTML presentations with Embed Shiny components
- Example: Demo by using File > New > RMarkdown > Shiny > Shiny Presentation
- *1/2 Note*. It is also possible to embed an entire shiny application within an R Markdown document using the `shinyAppDir` function or `shinyAppFile` or source
- http://rmarkdown.rstudio.com/authoring_embedded_shiny.html
- 3. R Markdown report as Shiny app
- Example [https://jjallaire.shinyapps.io/shiny-ggplot2-brushing/](https://jjallaire.shinyapps.io/shiny-ggplot2-brushing/)

## Include Markdown content in a Shiny app - Overview
- It's possible to include files with HTML, text, or Markdown content in a Shiny app.
- [http://shiny.rstudio.com/gallery/including-html-text-and-markdown-files.html](http://shiny.rstudio.com/gallery/including-html-text-and-markdown-files.html)
- The includeMarkdown function requires the markdown package.
- [https://shiny.rstudio.com/reference/shiny/latest/include.html](https://shiny.rstudio.com/reference/shiny/latest/include.html)

## Include Markdown content in a Shiny app - Example
- [https://github.com/abhik1368/openfda-dashboard](https://github.com/abhik1368/openfda-dashboard)
- [https://raw.githubusercontent.com/abhik1368/openfda-dashboard/master/ui.R](https://raw.githubusercontent.com/abhik1368/openfda-dashboard/master/ui.R)
- [https://diabetescheck.shinyapps.io/openfda-dashboard/](https://diabetescheck.shinyapps.io/openfda-dashboard/)

## Generating downloadable reports from Shiny app
- It is often required after interacting with a Shiny application to download a report in HTML or PDF format. This can be done with knitr and rmarkdown. The app takes one input value, and passes it as a parameter to an .Rmd document, which is rendered (or knitted) into an HTML output file.
- [https://shiny.rstudio.com/articles/generating-reports.html](https://shiny.rstudio.com/articles/generating-reports.html)
- [http://shiny.rstudio.com/gallery/download-knitr-reports.html](http://shiny.rstudio.com/gallery/download-knitr-reports.html)
- Another way to generate a PDF is to pass in a value to output_format argument of rmarkdown::render(). It will also require pdflatex to be installed on your system.
- http://shiny.rstudio.com/gallery/download-knitr-reports.html

## Generating downloadable reports from Shiny app - Example

- [https://github.com/aaronjfisher/interAdapt/tree/master/r_package/inst/interAdapt](https://github.com/aaronjfisher/interAdapt/tree/master/r_package/inst/interAdapt)
- [https://michaelrosenblum.shinyapps.io/interAdapt/](https://michaelrosenblum.shinyapps.io/interAdapt/)

## Embedding Shiny apps via their URLs in R Markdown output - Overview

- With the knitr package, the functions knitr::include_url() and include_app() can be used by inserting an iframe or a screenshot/image in the output. iframes are basically a box on one web page to embed another web page. knitr::include_app() was designed for embedding Shiny apps via their URLs in RMD output.
- [https://bookdown.org/yihui/bookdown/web-pages-and-shiny-apps.html](https://bookdown.org/yihui/bookdown/web-pages-and-shiny-apps.html)


## Embedding Shiny apps via their URLs in R Markdown output - Examples

- RMD HTML Output - http://r4ds.had.co.nz/r-markdown-formats.html
- https://bookdown.org/
- http://rmarkdown.rstudio.com/rmarkdown_websites.html
- https://github.com/rstudio/blogdown
  
## Shiny Components in rtutor - Overview & Examples
- rtutor package makes it easy to turn any R Markdown document into an interactive tutorial
- rtutor package uses runtime: shiny_prerendered to turn regular R Markdown documents into live tutorials
- https://rstudio.github.io/rtutor/index.html
- https://rstudio.github.io/rtutor/index.html#shiny_components
- https://rstudio.github.io/rtutor/examples.html

## Further Reading

- http://rmarkdown.rstudio.com/authoring_shiny_prerendered.html
- http://colorado.rstudio.com:3939/content/188/ by Edgar Ruiz
- [http://rmarkdown.rstudio.com/authoring_shiny_advanced.html](http://rmarkdown.rstudio.com/authoring_shiny_advanced.html)
