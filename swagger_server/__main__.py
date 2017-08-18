#!/usr/bin/env python3

import connexion
from flask import redirect
from .encoder import JSONEncoder



if __name__ == '__main__':
    app = connexion.App(__name__, specification_dir='./swagger/')
    app.app.json_encoder = JSONEncoder
    app.add_api('swagger.yaml', arguments={'title': 'The Allele Calling  service provides an API for converting raw sequence data to GFE and HLA allele calls.'})
    @app.route('/')
    def found():
        return redirect("/ui")
    app.run(port=5000, host="0.0.0.0")
