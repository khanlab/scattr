# Stage: build
# NOTE: g++ and libdatrie are required for poetry install
FROM python:3.9-slim-bullseye AS build
COPY ./poetry.lock ./pyproject.toml /
RUN mkdir -p /opt \
    && apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
        g++=4:10.2.1-1 \
        libdatrie1=0.2.13-1 \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && pip install --prefer-binary --no-cache-dir \
        poetry==1.2.2 \
    && poetry config virtualenvs.create false \
    && poetry install --only main \
    && apt-get purge -y -q g++ \
    && apt-get --purge -y -qq autoremove

# Stage: runtime
FROM build AS runtime
COPY ./dbsc /opt/dbsc
# ENTRYPOINT ["/opt/dbsc/run.py"]