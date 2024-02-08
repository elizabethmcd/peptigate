We used the following command to build the docker container
```
docker build --progress=plain --platform linux/x86_64 --no-cache -t deeppeptide .
```

And we pushed the docker container to Docker Hub using:
```
docker login
docker images
docker tag deeppeptide arcadiascience/deeppeptide:f9906d9
docker push arcadiascience/deeppeptide:f9906d9
```

If you run this docker container on a computer with Nvidia GPUs, you can access the GPU with the `--gpus` flag.
The following command will execute the DeepPeptide `predict.py` command.
The `-v` flag specificies which directory to mount in the docker container, allowing the docker container to see the files in that directory and for you to see files the docker container writes in that directory after the `docker run` command finishes.
It assumes your FASTA file is called `fasta_file.fasta`, and will write results to the `testrun_outputs` directory in teh mounted directory.
```
docker run --gpus all -v /path/to/directory/with/your/fasta:/app/DeepPeptide/predictor -it deeppeptide -c "python3 predict.py -ff fasta_file.fasta -od testrun_outputs/"
```

If you run this docker container inside of a snakemake workflow using the `singularity` rule directive, you can use the following command to access the GPU.
```
snakemake -j all --use-singularity --singularity-args '\\-\\-nv'
```
