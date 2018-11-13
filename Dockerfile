FROM continuumio/miniconda

RUN apt-get update && \
    apt-get install -y git vim gcc python-tk software-properties-common && \
    apt-get install -y heimdal-dev && \
    apt-get install -y libkrb5-dev && \
    apt-get install -y build-essential && \
    apt-get install -y libopenblas-dev libeigen3-dev sqlite3 libsqlite3-dev libboost-dev libboost-system-dev libboost-thread-dev libboost-serialization-dev libboost-python-dev libboost-regex-dev libcairo2 libcairo2-dev libjpeg-dev libgif-dev && \
    rm -rf /var/lib/apt/lists/*

RUN useradd -ms /bin/bash askcos
RUN echo "alias ltr='ls -lhtr'" >> /home/askcos/.bashrc

RUN conda update conda && \
    conda create -n askcos python=2.7 && \
    conda install --yes -n askcos mkl mkl-service && \
    conda install --yes -n askcos numpy scipy scikit-learn numexpr && \
    conda install --yes -n askcos tensorflow=1.4.1 && \
    conda install --yes -n askcos -c conda-forge theano=0.9.0 keras pycairo && \
    conda install --yes -n askcos -c rdkit rdkit=2017.03.3 && \
    conda install --yes -n askcos pip

RUN echo "source activate askcos" >> /home/askcos/.bashrc

COPY requirements.txt requirements.txt
RUN bash -c "source activate askcos && pip install -r requirements.txt && rm requirements.txt"

WORKDIR /home/askcos/
RUN mkdir /home/askcos/ASKCOS

COPY . /home/askcos/ASKCOS/


RUN chown -R askcos:askcos /home/askcos/ASKCOS

USER askcos

ENV PYTHONPATH="/home/askcos/ASKCOS:/home/askcos/ASKCOS/askcos:${PYTHONPATH}"
ENV KERAS_BACKEND=theano
