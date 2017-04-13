#!/usr/bin/env python3

import connexion
from swagger_server.encoder import JSONEncoder


if __name__ == '__main__':
    app = connexion.App(__name__, specification_dir='./swagger_server/swagger/')
    app.app.json_encoder = JSONEncoder
    app.add_api('swagger.yaml', arguments={'title': 'The Allele Calling  service provides an API for converting raw sequence data to GFE and HLA allele calls.'})
    app.run(host="0.0.0.0")
