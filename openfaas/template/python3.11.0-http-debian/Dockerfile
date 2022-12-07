FROM ghcr.io/openfaas/of-watchdog:0.9.10 AS watchdog

# Build image stage
FROM python:3.11.0-slim-bullseye AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG ADDITIONAL_BUILD_PACKAGE

# Upgrade pip and setuptools
RUN pip install --upgrade pip setuptools

RUN apt-get -qy update && \
	apt-get -qy install ${ADDITIONAL_BUILD_PACKAGE} && \
	rm -fr /var/lib/apt/lists/*

# Add non root user
RUN addgroup --system app && adduser app --system --ingroup app
USER app
WORKDIR /home/app

ENV PATH=/home/app/.local/bin:$PATH

COPY --chown=app index.py .
COPY --chown=app requirements.txt .
RUN pip install --no-cache-dir --user -r requirements.txt

RUN mkdir function
WORKDIR /home/app/function

RUN touch __init__.py
COPY --chown=app function/requirements.txt .
RUN pip install --no-cache-dir --user -r requirements.txt

COPY --chown=app function/ .

RUN rm -fr ../.cache

# Test image stage
FROM builder AS tester

ARG TEST_COMMAND=tox
ARG TEST_ENABLED=true
RUN [ "$TEST_ENABLED" = "false" ] && echo "skipping tests" || \
	(pip install tox==3.*; eval "$TEST_COMMAND" && rm -fr .tox)

# Production image stage
FROM python:3.11.0-slim-bullseye

COPY --from=watchdog /fwatchdog /usr/bin/fwatchdog
RUN chmod +x /usr/bin/fwatchdog

ARG DEBIAN_FRONTEND
ARG ADDITIONAL_PACKAGE
RUN apt-get -qy update && \
	apt-get -qy install ${ADDITIONAL_PACKAGE} && \
	rm -rf /var/lib/apt/lists/*

# Add non root user
RUN addgroup --system app && adduser app --system --ingroup app
USER app
WORKDIR /home/app

COPY --from=builder --chown=app /home/app /home/app
ENV PATH=/home/app/.local/bin:$PATH

# Set up of-watchdog for HTTP mode
ENV fprocess="python index.py"
ENV cgi_headers="true"
ENV mode="http"
ENV upstream_url="http://127.0.0.1:5000"

HEALTHCHECK --interval=5s CMD [ -e /tmp/.lock ] || exit 1

CMD ["fwatchdog"]
