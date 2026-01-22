from blueprints.smarts_filter import smarts_filter
from blueprints.load_smarts import load_smarts_bp
from flask import Blueprint

def register_routes(app, in_production: bool, version_url_prefix: str):
    version = Blueprint("version", __name__, url_prefix=version_url_prefix)
    version.register_blueprint(smarts_filter)
    version.register_blueprint(load_smarts_bp)
    app.register_blueprint(version)
