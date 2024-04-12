# nf-sequenoscope

The Nextflow wrapper for this pipeline serves as a critical tool for bioinformatics workflows, especially in the context of ONT adaptive-sampling sequencing. The adaptive sampling technique, where sequencing efforts are dynamically focused on target DNA or RNA molecules, generates unique challenges in data management and interpretation. This pipeline addresses these through its three specialized modules:

## Analyze Module

- **Functionality**: Acts as the core component, taking input from FASTQ files, a reference FASTA file, and an optional sequencing summary.
- **Operations**: Performs filtering, mapping, and generates comprehensive sequencing statistics, aiding in the initial processing and quality assessment of the sequencing data.

## Plot Module

- **Functionality**: Provides advanced visual analytics capabilities, turning complex data into interpretable visual formats.
- **Importance**: This module is essential for comparing and visualizing data from adaptive sampling experiments, using input from the analyze module to generate plots that highlight key differences and trends in the data.

## Filter_ONT Module

- **Focus**: This module focuses on the precise filtering of raw reads based on custom criteria such as channel number, duration, and quality scores, among others.
- **Utility**: It’s tailored for researchers needing detailed control over data inclusion in analyses, particularly useful in experiments with complex sampling strategies.

Together, these modules form a pipeline that not only simplifies the data processing workflow but also enhances the ability of researchers to conduct detailed and meaningful analyses of ONT sequencing data. The Nextflow wrapper ensures the pipeline is scalable, reproducible, and easily integrable into larger bioinformatics operations, making it an invaluable resource for genomic research.

