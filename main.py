#!/usr/bin/env python3

import connexion
from flask import redirect
from swagger_server import encoder

app = connexion.App(__name__, specification_dir='./swagger_server/swagger/')
app.app.json_encoder = encoder.JSONEncoder
app.add_api('swagger.yaml', arguments={'title': 'Allele Calling Service'})
app.app.config['MAX_CONTENT_LENGTH'] = 1600 * 1024 * 1024


@app.route("/")
def basepath():
    return redirect("/ui")


if __name__ == '__main__':
    app.run(port=9000)
