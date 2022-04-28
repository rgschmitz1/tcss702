FROM openfaas/of-watchdog:0.7.7 as watchdog

# Build image stage
FROM python:3.7-slim-buster as build

COPY --from=watchdog /fwatchdog /usr/bin/fwatchdog
RUN chmod +x /usr/bin/fwatchdog

ARG ADDITIONAL_BUILD_PACKAGE
# Alternatively use ADD https:// (which will not be cached by Docker builder)

RUN apt-get -qy update && apt-get -qy install ${ADDITIONAL_BUILD_PACKAGE}

# Add non root user
RUN addgroup --system app && adduser app --system --ingroup app
USER app
WORKDIR /home/app

RUN python -m venv venv
ENV PATH=/home/app/venv/bin:$PATH

COPY --chown=app index.py .
COPY --chown=app requirements.txt .
RUN pip install -r requirements.txt

RUN mkdir -p function
RUN touch ./function/__init__.py
WORKDIR /home/app/function
COPY --chown=app function/requirements.txt .
RUN pip install -r requirements.txt

COPY --chown=app function/ .

ARG TEST_COMMAND=tox
ARG TEST_ENABLED=true
RUN [ "$TEST_ENABLED" = "false" ] && echo "skipping tests" || eval "$TEST_COMMAND"

# Production image stage
FROM python:3.7-slim-buster

COPY --from=watchdog /fwatchdog /usr/bin/fwatchdog
RUN chmod +x /usr/bin/fwatchdog

ARG ADDITIONAL_PACKAGE
RUN apt-get -qy update && \
	apt-get -qy install ${ADDITIONAL_PACKAGE}; \
	rm -rf /var/lib/apt/lists/*

# Add non root user
RUN addgroup --system app && adduser app --system --ingroup app
USER app
WORKDIR /home/app

COPY --from=build --chown=app /home/app /home/app
ENV PATH=/home/app/venv/bin:$PATH

# Set up of-watchdog for HTTP mode
ENV fprocess="python index.py"
ENV cgi_headers="true"
ENV mode="http"
ENV upstream_url="http://127.0.0.1:5000"

HEALTHCHECK --interval=5s CMD [ -e /tmp/.lock ] || exit 1

CMD ["fwatchdog"]
