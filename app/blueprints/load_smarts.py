from flask import Blueprint, request, jsonify
import requests
from functools import lru_cache
from utils.smarts_parser import SmartsFile

load_smarts_bp = Blueprint("load_smarts", __name__, url_prefix="/load_smarts")

GITHUB_RAW_BASE = (
    "https://raw.githubusercontent.com/unmtransinfo/unm_biocomp/master/"
    "biocomp_war/src/main/webapp/data/smarts/"
)


@lru_cache(maxsize=32)
def fetch_sma(filename: str) -> str:
    url = GITHUB_RAW_BASE + filename
    resp = requests.get(url, timeout=10)
    resp.raise_for_status()
    return resp.text


def preprocess_and_expand(raw_text: str) -> list:
    """Preprocess raw .sma text using SmartsFile parser.

    Handles comments, define macros, inline comments, RDKit normalization,
    recovery of edge-case patterns, and name sanitization.

    Returns:
        List of strings in the format 'SMARTS_PATTERN name'
    """
    sf = SmartsFile()
    sf.parse_file(raw_text, strict=False)
    return [str(s) for s in sf.smartses]


def parse_smarts_text(raw_text: str):
    """Parse raw SMARTS text into validated patterns and names.

    Uses the full SmartsFile pipeline: comment stripping, macro expansion,
    RDKit normalization, and pattern recovery.
    """
    sf = SmartsFile()
    sf.parse_file(raw_text, strict=False)

    smarts = []
    names = []
    for s in sf.smartses:
        smarts.append(s.smarts)
        names.append(s.name)

    invalid = []
    for f in sf.failed_smarts:
        invalid.append({
            "line": f["line"],
            "pattern": f["raw"],
            "name": f.get("name", ""),
            "error": f.get("error", "Unknown error"),
        })

    return smarts, names, invalid


# ---------------------------------------------------------------------------
# Flask routes
# ---------------------------------------------------------------------------

@load_smarts_bp.route("/load", methods=["GET", "POST"])
def load_smarts():
    """
    Load and parse SMARTS patterns.
    GET: fetch a .sma file from GitHub.
    POST: parse raw SMARTS text submitted in the request body.
    ---
    tags:
      - SMARTS FILTER
    summary: Load and parse SMARTS patterns.
    description: |
      GET  - Fetches a .sma file from GitHub and returns its raw content.
      POST - Accepts raw SMARTS text, preprocesses (skips comments, expands
             define macros, normalizes for RDKit), validates, and returns
             structured arrays.
    parameters:
      - name: file_name
        in: query
        type: string
        required: false
        description: (GET only) Name of the .sma file to fetch from GitHub.
        example: unm_reactive.sma
    requestBody:
      content:
        application/json:
          schema:
            type: object
            required:
              - content
            properties:
              content:
                type: string
                description: Raw SMARTS text (one pattern per line, optional name after whitespace).
                example: "O~N(=O)-c(:*):* aromatic_NO2\\n[2H] deuterium"
    responses:
      200:
        description: Parsed SMARTS patterns
        content:
          application/json:
            schema:
              type: object
              properties:
                smarts:
                  type: array
                  items:
                    type: string
                names:
                  type: array
                  items:
                    type: string
                invalid:
                  type: array
                  items:
                    type: object
      400:
        description: Bad request
      500:
        description: Server error
    """
    if request.method == "GET":
        filename = request.args.get("file_name")
        if not filename or not filename.endswith(".sma"):
            return jsonify({"error": "Valid .sma filename required"}), 400
        try:
            text = fetch_sma(filename)
            return jsonify({"content": text}), 200
        except requests.HTTPError as e:
            return jsonify({"error": f"Failed to fetch file: {str(e)}"}), 500
        except Exception as e:
            return jsonify({"error": str(e)}), 500

    # POST: accept raw SMARTS text and return parsed/validated arrays
    json_data = request.get_json()
    if not json_data or "content" not in json_data:
        return jsonify({"error": "JSON body with 'content' field required"}), 400

    raw_text = json_data["content"]
    if not isinstance(raw_text, str) or not raw_text.strip():
        return jsonify({"error": "Content must be a non-empty string"}), 400

    smarts, names, invalid = parse_smarts_text(raw_text)

    if not smarts:
        return jsonify({
            "error": "No valid SMARTS patterns found",
            "invalid": invalid
        }), 400

    return jsonify({
        "smarts": smarts,
        "names": names,
        "invalid": invalid
    }), 200


@load_smarts_bp.route("/preprocess_expand", methods=["POST"])
def preprocess_expand_endpoint():
    """
    Preprocess and expand macro SMARTS text.
    ---
    tags:
      - SMARTS FILTER
    summary: Expand define macros in SMARTS text.
    description: |
      Accepts raw .sma text with 'define' macros, comments, etc.
      Returns expanded SMARTS lines (macros resolved, comments stripped).
    requestBody:
      content:
        application/json:
          schema:
            type: object
            required:
              - smarts_text
            properties:
              smarts_text:
                type: string
    responses:
      200:
        description: Expanded SMARTS lines
      400:
        description: Bad request
    """
    json_data = request.get_json()
    if not json_data:
        return jsonify({"error": "JSON body required"}), 400

    raw_text = json_data.get("smarts_text") or json_data.get("content", "")
    if not raw_text or not raw_text.strip():
        return jsonify({"error": "Non-empty smarts_text required"}), 400

    expanded = preprocess_and_expand(raw_text)
    return jsonify({"expanded_smarts": expanded}), 200
