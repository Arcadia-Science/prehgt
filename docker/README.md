# Example instructions for how to build a docker container from a Dockerfile

The both images needed to be run with an explicit linux platform if built on M1 Mac.
```
docker build --platform linux/x86_64 --no-cache -t ncbi-genome-download ncbi-genome-download
docker login
docker images
docker tag ncbi-genome-download arcadiascience/ncbi-genome-download-patch:4c5c24e
docker push arcadiascience/ncbi-genome-download-patch:4c5c24e
```

```
docker build --progress=plain --no-cache --platform linux/x86_64 -t tidy-prehgt tidy-prehgt
docker tag tidy-prehgt arcadiascience/tidy-prehgt:2.0.0
docker push arcadiascience/tidy-prehgt:2.0.0
```
