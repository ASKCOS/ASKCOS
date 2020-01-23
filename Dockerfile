FROM python:3.5-stretch

COPY --from=registry.gitlab.com/mlpds_mit/askcos/askcos/rdkit:2019.03-py35 /usr/local/rdkit-2019-03/rdkit /usr/local/rdkit-2019-03/rdkit
COPY --from=registry.gitlab.com/mlpds_mit/askcos/askcos/rdkit:2019.03-py35 /usr/local/rdkit-2019-03/lib /usr/local/rdkit-2019-03/lib
COPY requirements.txt requirements.txt

RUN apt-get update && \
    apt-get install -y libboost-thread-dev libboost-python-dev libboost-iostreams-dev python-tk libopenblas-dev libeigen3-dev libcairo2-dev pkg-config python-dev python-mysqldb && \
    pip install Cython 'numpy>=1.16.4' && \
    pip install -r requirements.txt && rm requirements.txt && \
    useradd -ms /bin/bash askcos

COPY --from=registry.gitlab.com/mlpds_mit/askcos/makeit-data:0.4.1 /data /usr/local/ASKCOS/makeit/data

COPY --chown=askcos:askcos . /usr/local/ASKCOS


WORKDIR /home/askcos
USER askcos

ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"/usr/local/rdkit-2019-03/lib"
ENV PYTHONPATH=${PYTHONPATH}:/usr/local/rdkit-2019-03:/usr/local/ASKCOS:/usr/local/ASKCOS/askcos/

RUN python /usr/local/ASKCOS/askcos/manage.py collectstatic --noinput
