FROM mongo

COPY *.json.gz /init/
COPY seed.sh /init/
COPY *.js /init/
CMD bash /init/seed.sh
