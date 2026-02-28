FROM python:3.10-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    wget \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone and build Cas-OFFinder globally
RUN git clone https://github.com/snugel/cas-offinder.git /opt/cas-offinder && \
    cd /opt/cas-offinder && \
    make

ENV PATH="/opt/cas-offinder:${PATH}"

WORKDIR /app

# Install python dependencies first for caching purposes
COPY pyproject.toml ./
# Because we only have pyproject.toml, let's install directly 
# We simulate a package install
RUN pip install --upgrade pip setuptools wheel
RUN pip install pandas primer3-py pyfaidx pyranges streamlit matplotlib openpyxl joblib pytest

# Copy application source
COPY . .

# Install the app itself in editable mode without isolation to ensure all hooks register correctly
RUN pip install -e .

EXPOSE 8501

CMD ["streamlit", "run", "app/streamlit_app.py", "--server.port=8501", "--server.address=0.0.0.0"]
