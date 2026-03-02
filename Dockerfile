FROM python:3.12-slim

WORKDIR /app

# Install system dependencies needed for git and potentially rdkit-tools
RUN apt-get update && apt-get install -y git curl && rm -rf /var/lib/apt/lists/*

COPY app/requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install git+https://github.com/jeremyjyang/rdkit-tools 
COPY app/ .

RUN mkdir -p flasgger_static && \
    cp -r /usr/local/lib/python3.12/site-packages/flasgger/ui3/static/* flasgger_static

ENV APP_PORT=8000
CMD gunicorn --bind "0.0.0.0:${APP_PORT}" app:app
