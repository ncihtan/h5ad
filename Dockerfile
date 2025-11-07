# The image is used by DCQC: https://github.com/Sage-Bionetworks-Workflows/py-dcqc
# The images is hosted at ghcr.io/sage-bionetworks-workflows/htan-h5ad-validator
FROM python:3.13

RUN pip install --no-cache-dir HTAN-h5ad-validate

COPY exemplars/* /usr/local/bin/
COPY h5ad.py /usr/local/bin/