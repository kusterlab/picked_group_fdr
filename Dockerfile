FROM ubuntu:16.04

MAINTAINER Matthew The "matthew.the@tum.de"

LABEL website=https://gitlab.lrz.de/proteomics/picked_group_fdr
LABEL description="Picked group FDR"
LABEL tags="protein grouping"
LABEL documentation=https://gitlab.lrz.de/proteomics/picked_group_fdr

# https://github.com/pyenv/pyenv/wiki/common-build-problems#prerequisites
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get --no-install-recommends -y install \
      git ca-certificates build-essential libssl-dev zlib1g-dev libbz2-dev \
      libreadline-dev libsqlite3-dev curl llvm libncurses5-dev libncursesw5-dev \
      xz-utils tk-dev libffi-dev liblzma-dev python-openssl && \
    rm -rf \
      /var/lib/apt/lists/* \
      /usr/share/doc \
      /usr/share/doc-base \
      /usr/share/man \
      /usr/share/locale \
      /usr/share/zoneinfo

## PYENV FOR MODERN PYTHON
# https://github.com/jprjr/docker-pyenv/blob/master/Dockerfile
# This allows an easy installation
ENV HOME  /root
ENV PYENV_ROOT $HOME/.pyenv
ENV PATH $PYENV_ROOT/shims:$PYENV_ROOT/bin:$PATH

RUN git clone https://github.com/pyenv/pyenv.git ~/.pyenv
##### THIS IS JUST FOR BEING INSIDE THE CONTAINER - NOT FOR INSTALLATION
RUN echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
RUN echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
RUN echo 'eval "$(pyenv init -)"' >> ~/.bashrc

WORKDIR /root/
RUN pyenv install 3.8.12
RUN pyenv local 3.8.12
RUN pyenv rehash

# see https://github.com/python-poetry/poetry/issues/1427
ENV LANG C.UTF-8

RUN pip install -U pip
RUN pip install poetry

ADD pyproject.toml /root/pyproject.toml
ADD poetry.lock /root/poetry.lock

# poetry useses virtualenvs by default -> we want global installation
RUN poetry config virtualenvs.create false
RUN poetry install --no-dev

ADD ./picked_group_fdr/ /root/picked_group_fdr

# RUN cd /root/picked_group_fdr && python setup.py build_ext --inplace

RUN cd /root/
