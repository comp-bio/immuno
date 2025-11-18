# check=skip=FromPlatformFlagConstDisallowed
FROM --platform=linux/amd64 python:3.8-slim

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    tcsh gawk build-essential nano \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

COPY IEDB_MHC_I-3.1.6.tar.gz /IEDB_MHC_I-3.1.6.tar.gz
COPY IEDB_NetChop-3.0.tar.gz /IEDB_NetChop-3.0.tar.gz

WORKDIR /app

RUN cd / && tar -xvzf IEDB_NetChop-3.0.tar.gz
RUN cd / && tar -xvzf IEDB_MHC_I-3.1.6.tar.gz

RUN ln -s /usr/local/bin/python3.8 /usr/bin/python || true
RUN chmod +x /netchop/build.sh && cd /netchop && ./build.sh || echo "netchop build failed"
RUN chmod +x /mhc_i/configure && cd /mhc_i && ./configure || echo "configure failed or was already run"

SHELL ["/bin/tcsh", "-c"]
CMD ["/bin/tcsh"]
