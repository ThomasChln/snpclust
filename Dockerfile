from thomaschln/r-publish

# Install plink, copied from
# https://github.com/GELOG/plink/blob/v1.07/plink-1.07-bin/docker/Dockerfile 
# Changed to http://zzz.bwh.harvard.edu

# Environment variables
env PLINK_VERSION       1.07
env PLINK_HOME          /usr/local/plink
env PATH                $PLINK_HOME:$PATH

run wget http://zzz.bwh.harvard.edu/plink/dist/plink-$PLINK_VERSION-x86_64.zip && \
    unzip plink-$PLINK_VERSION-x86_64.zip -d /usr/local/ && \
    rm plink-$PLINK_VERSION-x86_64.zip && \
    cd /usr/local && \
    ln -s plink-$PLINK_VERSION-x86_64 $PLINK_HOME && \
    rm -rf /var/lib/apt/lists/*

# Install SHAPEIT
# https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html
run wget -O shapeit.tgz https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.17.linux.tar.gz \
  && tar -zxvf shapeit.tgz \
  && mv /shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit /bin/


# Install SNPClust 

run apt-get update && apt-get install -y libnetcdf-dev
run R -e "install.packages('BiocManager');BiocManager::install('GWASTools')"
run R -e "devtools::install_github('zhengxwen/SNPRelate')"

add ./ snpclust/
run R -e "devtools::install('snpclust', dependencies = TRUE)"
