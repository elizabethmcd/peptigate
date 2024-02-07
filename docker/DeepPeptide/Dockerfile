# Use an official NVIDIA CUDA base image with Python 3
FROM nvidia/cuda:11.3.1-cudnn8-runtime-ubuntu20.04

# Set a working directory
WORKDIR /DeepPeptide

# Clone the DeepPeptide repo
RUN apt-get update && apt-get install -y git
RUN git clone https://github.com/fteufel/DeepPeptide.git /app/DeepPeptide

# Install Python and pip
RUN apt-get install -y python3 python3-pip

# Install Python dependencies
RUN pip3 install -r /app/DeepPeptide/predictor/requirements.txt

# Set the working directory to the predictor directory inside the app
WORKDIR /app/DeepPeptide/predictor

# Set bash as the default entry point
ENTRYPOINT ["/bin/bash"]
