FROM ubuntu:16.10
MAINTAINER Mike Halagan <mhalagan@nmdp.org>

RUN apt-get update -q \
    && apt-get dist-upgrade -qy \
    && apt-get install -qyy wget curl build-essential cpp git \
    && apt-get -qyy install python3.6 python3-pip python3-dev python3-setuptools uwsgi-plugin-python3

RUN apt-get install python3.6-dev -qy

RUN cd opt/ && git clone https://github.com/nmdp-bioinformatics/service-act && cd service-act \
    && curl https://bootstrap.pypa.io/get-pip.py | python3.6 \
    && pip install --upgrade pip \
	&& pip install -e 'git+https://github.com/nmdp-bioinformatics/service-gfe-submission.git#egg=1.0.0&subdirectory=client-python'

RUN cd opt/service-act && pip install -r requirements.txt \
    && python3.6 setup.py install 

RUN apt-get install clustalo -y
RUN apt-get install ncbi-blast+ -y

CMD uwsgi --http :8080 --plugin python3.6 --wsgi-file opt/service-act/app.py --callable app -p 10

