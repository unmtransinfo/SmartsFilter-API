from flask import Blueprint, request, jsonify
import requests
from functools import lru_cache

load_smarts_bp = Blueprint("load_smarts", __name__, url_prefix="/load_smarts")

GITHUB_RAW_BASE = (
    "https://raw.githubusercontent.com/unmtransinfo/unm_biocomp/master/"
    "biocomp_war/src/main/webapp/data/smarts/"
)

# Cache up to 32 files to avoid repeated network calls
@lru_cache(maxsize=32)
def fetch_sma(filename: str) -> str:
    url = GITHUB_RAW_BASE + filename
    resp = requests.get(url, timeout=10)
    resp.raise_for_status()
    return resp.text


@load_smarts_bp.route("/load", methods=["GET"])
def load_smarts():
    """
    Load a SMARTS (.sma) file from GitHub
    ---
    tags:
      - SMARTS FILTER
    summary: Load SMARTS patterns from external source.
    description: Fetches SMARTS patterns from a predefined external URL.
    parameters:
      - name: file_name
        in: query
        type: string
        required: true
        description: Name of the .sma file to fetch from GitHub.
        example: unm_reactive.sma
    responses:
      200:
        description: Raw SMARTS file content in JSON
        content:
          application/json:
            schema:
              type: object
              properties:
                content:
                  type: string
                  example: "O~N(=O)-c(:*):* aromatic NO2\n[2H] deuterium"
      400:
        description: Invalid file_name parameter
      500:
        description: Failed to fetch file from GitHub
    """
    filename = request.args.get("file_name")
    if not filename or not filename.endswith(".sma"):
        return jsonify({"error": "Valid .sma filename required"}), 400

    try:
        text = fetch_sma(filename)
        # Return as JSON
        return jsonify({"content": text}), 200
    except requests.HTTPError as e:
        return jsonify({"error": f"Failed to fetch file: {str(e)}"}), 500
    except Exception as e:
        return jsonify({"error": str(e)}), 500
