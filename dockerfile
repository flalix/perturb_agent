FROM rocker/r-ver:4.3.2

ENV PYTHONUNBUFFERED=1
ENV UV_SYSTEM_PYTHON=1
ENV STORAGE_DIR=/opt/render/data
ENV RENV_PATHS_CACHE=/opt/renv/cache

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    curl \
    git \
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
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install uv

WORKDIR /app

COPY pyproject.toml uv.lock ./
RUN uv sync --frozen

COPY scripts/setup_renv.R scripts/setup_renv.R
RUN Rscript scripts/setup_renv.R

COPY . .

RUN mkdir -p ${STORAGE_DIR}

EXPOSE 10000

CMD uv run streamlit run app.py \
    --server.address 0.0.0.0 \
    --server.port ${PORT:-10000}