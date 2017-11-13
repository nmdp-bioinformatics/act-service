FROM ubuntu:16.10
MAINTAINER Mike Halagan <mhalagan@nmdp.org>

RUN apt-get update -q \
    && apt-get dist-upgrade -qy \
    && apt-get install -qyy wget curl build-essential cpp git \
    && apt-get -qyy install python3.6 python3-pip python3-dev python3-setuptools uwsgi-plugin-python3 \
    && pip install -e 'git+https://github.com/nmdp-bioinformatics/service-gfe-submission.git#egg=1.0.0&subdirectory=client-python' \
    && cd opt/ && git clone https://github.com/nmdp-bioinformatics/service-act && cd service-act \
    && pip3 install --upgrade pip \
    && pip3 install -r requirements.txt \
    && python3 setup.py install 

CMD uwsgi --http :8080 --plugin python --wsgi-file opt/service-act/app.py --callable app -p 10

