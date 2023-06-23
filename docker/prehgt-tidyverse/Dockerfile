FROM arcadiascience/rbase_4.2.2:0.0.1

# Install the required dependencies and R
RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
        dirmngr \
        gpg-agent \
        libcurl4-openssl-dev \
        libicu-dev \
        libssl-dev \
        libtiff-dev \
        libxml2-dev \
        make \
        pandoc \
        zlib1g-dev

# Install R packages required to run R modules
RUN Rscript -e "install.packages('tidyverse', version = '2.0.0', repos='https://cloud.r-project.org/')"
RUN Rscript -e "install.packages('janitor', version = '2.2.0', repos='https://cloud.r-project.org/')"
RUN Rscript -e "install.packages('fastcluster', version = '1.2.3', repos='https://cloud.r-project.org/')"
RUN Rscript -e "install.packages('dbplyr', version = '2.3.2', repos='https://cloud.r-project.org/')"
RUN Rscript -e "install.packages('dbi', version = '1.1.3', repos='https://cloud.r-project.org/')"
RUN Rscript -e "install.packages('rsqlite', version = '2.3.0', repos='https://cloud.r-project.org/')"

# Clean up
RUN apt-get -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set the default command to launch an interactive Ubuntu shell
CMD ["/bin/bash", "-l"]
