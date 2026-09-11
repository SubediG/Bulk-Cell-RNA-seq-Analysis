FROM bioconductor/bioconductor_docker:3.22-R-4.5.2

WORKDIR /project

COPY . .

RUN R -e "renv::restore(prompt = FALSE)"

CMD ["Rscript", "DESeq2_mTNBC_Analysis.R"]
