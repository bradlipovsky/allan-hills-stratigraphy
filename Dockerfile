FROM firedrakeproject/firedrake-vanilla:2025-01

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && rm -rf /var/lib/apt/lists/*

USER firedrake
RUN bash -lc "source /home/firedrake/firedrake/bin/activate \
    && python -m pip install --no-cache-dir \
        git+https://github.com/icepack/icepack.git \
        nbclient \
        nbconvert \
        openpyxl"

WORKDIR /workspace
