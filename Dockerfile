FROM python:3.5-stretch

COPY --from=mefortunato/rdkit:2017.03-py35 /usr/local/rdkit-2017-03/rdkit /usr/local/rdkit-2017-03/rdkit
COPY --from=mefortunato/rdkit:2017.03-py35 /usr/local/rdkit-2017-03/lib /usr/local/rdkit-2017-03/lib
COPY requirements.txt requirements.txt

RUN apt-get update && \
    apt-get install -y libboost-thread-dev libboost-python-dev python-tk libopenblas-dev libeigen3-dev libcairo2-dev pkg-config python-dev python-mysqldb && \
    pip install Cython 'numpy>=1.16.4' && \
    pip install -r requirements.txt && rm requirements.txt && \
    useradd -ms /bin/bash askcos

COPY --chown=askcos:askcos . /usr/local/ASKCOS


WORKDIR /home/askcos
USER askcos

ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"/usr/local/rdkit-2017-03/lib"
ENV PYTHONPATH=${PYTHONPATH}:/usr/local/rdkit-2017-03:/usr/local/ASKCOS:/usr/local/ASKCOS/askcos/
ENV KERAS_BACKEND=theano

RUN cd /usr/local/ASKCOS/docs && python -m sphinx.ext.apidoc -f -o ./ .. && make html && cd -
