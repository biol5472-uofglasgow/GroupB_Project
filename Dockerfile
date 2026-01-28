# Use Python 3.13 slim image as base
FROM python:3.13-slim

# Set working directory in container
WORKDIR /app

# Copy only dependency files first (for layer caching)
COPY pyproject.toml README.md LICENSE ./

# Install dependencies (this layer is cached if pyproject.toml doesn't change)
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir .

# Copy source code (changes here don't trigger dependency reinstall)
COPY src/ ./src/

# Reinstall to pick up code changes
RUN pip install --no-cache-dir --no-deps .

# Set the entry point to the CLI tool
ENTRYPOINT ["GroupB-tool"]

# Default command (show help if no arguments provided)
CMD ["--help"]
