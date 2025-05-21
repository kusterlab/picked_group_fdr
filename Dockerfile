FROM python:3.10.17

MAINTAINER Matthew The "matthew.the@tum.de"

LABEL website=https://gitlab.lrz.de/proteomics/picked_group_fdr
LABEL description="Picked group FDR"
LABEL tags="protein grouping"
LABEL documentation=https://gitlab.lrz.de/proteomics/picked_group_fdr

# https://github.com/pyenv/pyenv/wiki/common-build-problems#prerequisites
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get --no-install-recommends -y install \
      git && \
    rm -rf \
      /var/lib/apt/lists/* \
      /usr/share/doc \
      /usr/share/doc-base \
      /usr/share/man \
      /usr/share/locale \
      /usr/share/zoneinfo

ENV HOME=/root
WORKDIR /root/

# see https://github.com/python-poetry/poetry/issues/1427
ENV LANG=C.UTF-8

RUN pip install poetry==1.8.3
# poetry uses virtualenvs by default -> we want global installation
RUN poetry config virtualenvs.create false

# urllib3 v2.0 only supports OpenSSL 1.1.1+
RUN pip install urllib3==1.26.6

ADD pyproject.toml /root/pyproject.toml
ADD poetry.lock /root/poetry.lock
ADD README.md /root/
ADD ./picked_group_fdr/ /root/picked_group_fdr
ADD ./tests/unit_tests/ /root/tests/unit_tests
RUN poetry install --only main

# for bootstrapping, see the "bootstrap" rule in the Makefile
ADD Makefile* /root/
ADD config.py /root/

RUN cd /root/
