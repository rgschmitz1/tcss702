ARG FUNCTION_DIR=/home/app/.local

# There is an issue in bullseye and qemu which prevents compilation
#FROM python:slim-bullseye as build-image
FROM python:slim-buster as build-image

# Install aws-lambda-cpp build dependencies
RUN apt-get update && \
	apt-get install -y \
	g++ \
	make \
	cmake \
	unzip \
	libcurl4-openssl-dev

# Add non root user
RUN addgroup --system app && adduser app --system --ingroup app
RUN chown app /home/app

USER app  
ENV PATH=$PATH:/home/app/.local/bin
WORKDIR /home/app

# Install the runtime interface client and other app requirements
COPY requirements.txt .
RUN pip install --user -r requirements.txt
# END OF build-image
#-------------------

FROM openfaas/of-watchdog:0.7.7 as watchdog
FROM python:slim-buster

COPY --from=watchdog /fwatchdog /usr/bin/fwatchdog
RUN chmod +x /usr/bin/fwatchdog

ARG FUNCTION_DIR
ARG ADDITIONAL_PACKAGE
# Alternatively use ADD https:// (which will not be cached by Docker builder)

RUN apt-get -qy update && \
	apt-get -qy install gcc make ${ADDITIONAL_PACKAGE} && \
	rm -rf /var/lib/apt/lists/*

# Add non root user
RUN addgroup --system app && adduser app --system --ingroup app
RUN chown app /home/app

USER app

ENV PATH=$PATH:/home/app/.local/bin

WORKDIR /home/app/

COPY index.py .
COPY requirements.txt .
USER root
RUN pip install -r requirements.txt
USER app

#COPY function/requirements.txt	.
#RUN pip install --user -r requirements.txt

COPY --from=build-image ${FUNCTION_DIR} ${FUNCTION_DIR}

# Create function directory
RUN mkdir -p /home/app/function
WORKDIR /home/app/function/

# Set NLTK_DATA environment variable
ENV NLTK_DATA=/home/app/function/nltk_data

USER root
COPY function/ .
RUN chown -R app:app ../

ARG TEST_COMMAND=tox
ARG TEST_ENABLED=true
RUN if [ "$TEST_ENABLED" = "false" ]; then \
    echo "skipping tests"; \
    else \
    eval "$TEST_COMMAND"; \
    fi

WORKDIR /home/app/

USER app

# Set up of-watchdog for HTTP mode
ENV fprocess="python index.py"
ENV cgi_headers="true"
ENV mode="http"
ENV upstream_url="http://127.0.0.1:5000"

HEALTHCHECK --interval=5s CMD [ -e /tmp/.lock ] || exit 1

CMD ["fwatchdog"]
