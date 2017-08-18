#!/usr/bin/env python3

import connexion
from flask import redirect
from swagger_server.encoder import JSONEncoder

app = connexion.App(__name__, specification_dir='./swagger_server/swagger/')
app.app.json_encoder = JSONEncoder
app.add_api('swagger.yaml', arguments={'title': 'The Allele Calling  service provides an API for converting raw sequence data to GFE and HLA allele calls.'})


@app.route('/')
def found():
    return redirect("/ui")

app.run(host="0.0.0.0")
application = app.app

