FROM ubuntu:14.04
MAINTAINER Mike Halagan <mhalagan@nmdp.org>

RUN apt-get update -q \
    && apt-get dist-upgrade -qy \
    && apt-get install -qyy wget curl build-essential cpp git \
    && apt-get -qyy install python3 python3-pip python3-dev python3-setuptools uwsgi-plugin-python3 \
    && cd opt/ && git clone https://github.com/nmdp-bioinformatics/service-act && cd service-act \
    && sudo pip3 install --upgrade pip \
    && sudo pip3 install -r requirements.txt \
    && pip install -e 'git+https://github.com/nmdp-bioinformatics/service-gfe-submission.git#egg=1.0.0&subdirectory=client-python' \
    && sudo python3 setup.py install 

CMD python3 /opt/service-act/app.py

