FROM mongo:4.4

COPY *.json.gz /init/
COPY seed.sh /init/
COPY *.js /init/
CMD bash /init/seed.sh
