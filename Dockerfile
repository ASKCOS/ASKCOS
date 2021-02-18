FROM askcos/askcos-base:2019.03.4-gh2855-py35

RUN apt-get update && \
    apt-get install -y libboost-thread-dev libboost-python-dev libboost-iostreams-dev python-tk libopenblas-dev libeigen3-dev libcairo2-dev pkg-config python-dev python-mysqldb && \
    useradd -ms /bin/bash askcos

COPY requirements.txt requirements.txt
RUN pip install --no-cache-dir -r requirements.txt && rm requirements.txt

COPY --chown=askcos:askcos --from=askcos/askcos-data:0.4.1 /data /usr/local/ASKCOS/makeit/data
COPY --chown=askcos:askcos . /usr/local/ASKCOS

WORKDIR /home/askcos
USER askcos

ENV PYTHONPATH=/usr/local/ASKCOS/askcos:/usr/local/ASKCOS${PYTHONPATH:+:${PYTHONPATH}}

RUN python /usr/local/ASKCOS/askcos/manage.py collectstatic --noinput
