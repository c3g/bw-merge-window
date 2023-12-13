FROM python:3.11-bookworm

RUN apt-get update -y && \
    apt-get install -y \
      build-essential \
      libcurl4-openssl-dev \
      libz-dev && \
    rm -rf /var/lib/apt/lists/*

RUN wget http://nz2.archive.ubuntu.com/ubuntu/pool/main/o/openssl/libssl1.1_1.1.1f-1ubuntu2_amd64.deb && \
    dpkg -i libssl1.1_1.1.1f-1ubuntu2_amd64.deb

WORKDIR /test

COPY run_tests.bash .

RUN pip install --no-cache-dir poetry

COPY pyproject.toml .
COPY poetry.lock .

RUN poetry install --no-root

COPY README.md .
COPY bw_merge_window bw_merge_window
COPY tests tests

RUN poetry install

CMD ["/bin/bash", "/test/run_tests.bash"]
