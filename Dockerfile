FROM debian:stretch as rdkit

RUN apt-get update
RUN apt-get install -y git gcc cmake software-properties-common build-essential python-dev libopenblas-dev libeigen3-dev sqlite3 libsqlite3-dev libboost-dev libboost-system-dev libboost-thread-dev libboost-serialization-dev libboost-python-dev libboost-regex-dev libcairo2 libcairo2-dev libjpeg-dev libgif-dev
RUN apt-get install -y python-pip
RUN pip install numpy
RUN export RDBASE=/usr/local/rdkit-2017-03 && \
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$RDBASE/lib && \
    export PYTHONPATH=$PYTHONPATH:$RDBASE && \
    git clone -b Release_2017_03 https://github.com/rdkit/rdkit.git $RDBASE && \
    cd $RDBASE && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j8 && \
    make install

FROM debian:stretch

COPY --from=rdkit /usr/local/rdkit-2017-03 /usr/local/rdkit-2017-03
COPY requirements.txt requirements.txt

RUN apt-get update && \
    apt-get install -y python python-pip libboost-thread-dev libboost-python-dev python-tk libopenblas-dev libeigen3-dev libcairo2-dev pkg-config && \
    pip install -r requirements.txt && rm requirements.txt && \
    useradd -ms /bin/bash askcos

COPY --chown=askcos:askcos . /usr/local/ASKCOS

WORKDIR /home/askcos
USER askcos

ENV LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:"/usr/local/rdkit-2017-03/lib"
ENV PYTHONPATH=${PYTHONPATH}:/usr/local/rdkit-2017-03:/usr/local/ASKCOS:/usr/local/ASKCOS/askcos/
ENV KERAS_BACKEND=theano
