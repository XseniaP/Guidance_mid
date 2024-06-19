# syntax=docker/dockerfile:1
FROM --platform=linux/amd64 ubuntu:24.04
WORKDIR /

RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN apt-get -y update
RUN apt install python3-pip -y
RUN apt install python3.12-venv -y
RUN apt install mafft -y
RUN apt install prank -y
COPY guidance_Linux guidance_Linux
COPY test test
COPY start.sh start.sh
RUN chmod +x ./start.sh
RUN python3 -m venv .venv
RUN cd /guidance_Linux/script/programs/msa_set_score/ && chmod +x ./msa_set_score && cd /
RUN source .venv/bin/activate && pip install setuptools && python3 -m pip install -U wheel setuptools && python3 -m pip install -r ./guidance_Linux/requirements.txt
ENTRYPOINT ["./start.sh"]

