# Tumors-produce-glucocorticoids-by-metabolite-recycling

This code accompanies the paper entitled:

<br>Taves MD, Otsuka S, Taylor MA, Donahue KM, Meyer TJ, Cam MC, Ashwell JD. Tumors produce glucocorticoids by metabolite recycling, not synthesis, and activate Tregs to promote growth. J Clin Invest. 2023 Sep 15;133(18):e164599. doi: 10.1172/JCI164599. PMID: 37471141; PMCID: PMC10503810.  [Link](https://pubmed.ncbi.nlm.nih.gov/37471141/)


To reproduce these results, follow these steps:

1.  Clone this GitHub repo (i.e. the page you are on):
    * ```git clone https://github.com/NIDAP-Community/Tumors-produce-glucocorticoids-by-metabolite-recycling.git```

2.  The input files for this pipeline will be available upon request. Please reach out to the authors before continue to following steps

3.  Install docker and build the docker container:
    * Navigate to the cloned repository directory. 
    * Move to the ./Docker_file/ directory of this repo

4.  Build the container:
    * ```docker build --tag Tumors-produce-glucocorticoids-by-metabolite-recycling .```

5.  Navidate to the cloned repository directory, Run the conainer by mounting the ./src/ directory of the repo to /tmp/ in the container:
    * ```docker run -ti -v $(pwd)/src:/mnt Tumors-produce-glucocorticoids-by-metabolite-recycling```
    
6.  Run the following code.
    * ```cd /mnt```
    * ```bash run_pipeline.sh```

