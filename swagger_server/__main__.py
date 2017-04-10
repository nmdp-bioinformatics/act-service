#!/usr/bin/env python3

import connexion
from .encoder import JSONEncoder


if __name__ == '__main__':
    app = connexion.App(__name__, specification_dir='./swagger/')
    app.app.json_encoder = JSONEncoder
    app.add_api('swagger.yaml', arguments={'title': 'The Allele Calling  service provides an API for converting raw sequence data to GFE. It provides both a RESTful API and a simple user interface for converting raw sequence data to GFE results. Sequences can be submitted one at a time or as a fasta file. This service uses &lt;a href&#x3D;\&quot;https://github.com/nmdp-bioinformatics/service-feature\&quot;&gt;nmdp-bioinformatics/service-feature&lt;/a&gt; for encoding the raw sequence data and &lt;a href&#x3D;\&quot;https://github.com/nmdp-bioinformatics/HSA\&quot;&gt;nmdp-bioinformatics/HSA&lt;/a&gt; for aligning the raw sequence data. The code is open source, and available on &lt;a href&#x3D;\&quot;https://github.com/nmdp-bioinformatics/service-act\&quot;&gt;GitHub&lt;/a&gt;.&lt;br&gt;&lt;br&gt;Go to &lt;a href&#x3D;\&quot;http://service-act.readthedocs.io\&quot;&gt;service-act.readthedocs.io&lt;/a&gt; for more information'})
    app.run(port=8080)
