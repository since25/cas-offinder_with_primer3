# Stage 1: Runtime dependencies
FROM python:3.9-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    bedtools \
    ocl-icd-opencl-dev \
    && rm -rf /var/lib/apt/lists/*

# Set up working directory
WORKDIR /app

# Copy application source and python dependencies
COPY pyproject.toml README.md ./
COPY src/ ./src/
COPY app/ ./app/

# Copy pre-compiled Linux cas-offinder binary
COPY build_linux/cas-offinder /usr/local/bin/cas-offinder
RUN chmod +x /usr/local/bin/cas-offinder

# Install python dependencies
RUN pip install --no-cache-dir .

# Set Python Path
ENV PYTHONPATH=/app/src
ENV STREAMLIT_SERVER_PORT=8501
ENV STREAMLIT_SERVER_ADDRESS=0.0.0.0

EXPOSE 8501

ENTRYPOINT ["streamlit", "run", "app/streamlit_app.py"]
