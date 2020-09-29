# Start from a core stack version
FROM jupyter/datascience-notebook

# Install from requirements.txt file
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt