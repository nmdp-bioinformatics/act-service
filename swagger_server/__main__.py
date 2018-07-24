#!/usr/bin/env python3

import connexion
from flask import redirect
from swagger_server import encoder

app = connexion.App(__name__, specification_dir='./swagger/')
app.app.json_encoder = encoder.JSONEncoder
app.add_api('swagger.yaml', arguments={'title': 'Allele Calling Service'})


@app.route("/")
def basepath():
    return redirect("/ui")


if __name__ == '__main__':
    app.run(port=8080)
