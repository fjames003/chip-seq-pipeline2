# Tutorial for general UNIX computers with docker

1. Download [cromwell](https://github.com/broadinstitute/cromwell).
    ```bash
    $ wget https://github.com/broadinstitute/cromwell/releases/download/34/cromwell-34.jar
    $ chmod +rx cromwell-34.jar
    ```

2. Git clone this pipeline and move into it.
    ```bash
    $ git clone https://github.com/ENCODE-DCC/chip-seq-pipeline2
    $ cd chip-seq-pipeline2
    ```

3. Download a SUBSAMPLED paired-end sample of [ENCSR936XTK](https://www.encodeproject.org/experiments/ENCSR936XTK/).
    ```bash
    $ wget https://storage.googleapis.com/encode-pipeline-test-samples/encode-chip-seq-pipeline/ENCSR936XTK/ENCSR936XTK_fastq_subsampled.tar
    $ tar xvf ENCSR936XTK_fastq_subsampled.tar
    ```

4. Download pre-built genome database for hg38.
    ```bash
    $ wget https://storage.googleapis.com/encode-pipeline-genome-data/test_genome_database_hg38_chip.tar
    $ tar xvf test_genome_database_hg38_chip.tar
    ```
    
5. Run a pipeline for the test sample.
    ```bash
    $ INPUT=examples/local/ENCSR936XTK_subsampled.json
    $ java -jar -Dconfig.file=backends/backend.conf cromwell-34.jar run chip.wdl -i ${INPUT} -o workflow_opts/docker.json
    ```

6. It will take about 6 hours. You will be able to find all outputs on `cromwell-executions/chip/[RANDOM_HASH_STRING]/`. See [output directory structure](output.md) for details.

7. See full specification for [input JSON file](input.md).
