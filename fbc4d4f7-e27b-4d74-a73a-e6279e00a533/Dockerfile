FROM rocker/r-ver:4.3.0

WORKDIR /workspace

RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('ggplot2', 'dplyr'), repos='https://cloud.r-project.org/')"

COPY plot.R /workspace/
COPY extracted_solution.csv /workspace/

CMD ["Rscript", "plot.R"]
