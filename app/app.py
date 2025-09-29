import os
import yaml
from flask import Flask
from flask_cors import CORS
from flasgger import LazyJSONEncoder, Swagger
from dotenv import load_dotenv
from blueprints.version import register_routes


def _load_api_spec() -> dict:
    """Load Swagger YAML spec."""
    with open("swagger_template.yml", "r") as file:
        return yaml.safe_load(file)


def _get_updated_paths(paths_dict: dict, version_prefix: str, in_production: bool, url_prefix: str = ""):
    updated_paths = {}
    for path, definition in paths_dict.items():
        if in_production:
            # Production: URL_PREFIX + version prefix + path
            new_path = f"/{url_prefix.strip('/')}{version_prefix}{path}"
        else:
            # Dev: version prefix + path
            new_path = f"{version_prefix}{path}"
        updated_paths[new_path] = definition
    return updated_paths


def create_app():
    app = Flask(__name__)
    app.json_encoder = LazyJSONEncoder
    load_dotenv(".env")

    app.config.from_pyfile("config.py")
    CORS(app, resources={r"/*": {"origins": "*"}})

    # Load Swagger template
    swagger_template = _load_api_spec()
    VERSION = swagger_template["info"]["version"]
    VERSION_URL_PREFIX = f"/api/v{VERSION}"

    URL_PREFIX = app.config.get("URL_PREFIX", "")
    IN_PROD = app.config.get("FLASK_ENV", "") == "production"
    print("IN_PROD value:", IN_PROD)

    # Configure Flasgger
    swagger_config = {
        "headers": [],
        "specs": [
            {
                "endpoint": "apispec_1",
                "route": "/apidocs/apispec_1.json",
            }
        ],
        "static_url_path": f"/{URL_PREFIX}/flasgger_static",
        "swagger_ui": True,
        "specs_route": "/apidocs/",
        "ui_params": {
            "url_prefix": (f"/{URL_PREFIX}") if IN_PROD else "",
        },
        "auth": {},
    }

    # Update Swagger template paths so docs match real routes
    swagger_template["swaggerUiPrefix"] = f"/{URL_PREFIX}" if IN_PROD else ""
    swagger_template["paths"] = _get_updated_paths(
        swagger_template["paths"], VERSION_URL_PREFIX, IN_PROD, URL_PREFIX
    )

    # Setup Swagger
    swagger = Swagger(app, config=swagger_config, template=swagger_template)

    # Register routes
    register_routes(app, IN_PROD, VERSION_URL_PREFIX)

    return app


app = create_app()

if __name__ == "__main__":
    app.run(
        host="0.0.0.0",
        port=app.config.get("APP_PORT"),
        debug=True
    )
