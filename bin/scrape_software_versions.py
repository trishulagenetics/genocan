#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'trishulagenetics/genocan': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'fastP': ['v_fastp.txt', r"([\d\.]+)"],
    'Qualimap': ['v_qualimap.txt', r"QualiMap v.(\S+)"],
    'minimap2': ['v_minimap2.txt', r"(\S+)"],
    'samtools': ['v_samtools.txt', r"(\S+)"],
}

results = OrderedDict()
results['trishulagenetics/genocan'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['fastP'] = '<span style="color:#999999;\">N/A</span>'
results['Qualimap'] = '<span style="color:#999999;\">N/A</span>'
results['minimap2'] = '<span style="color:#999999;\">N/A</span>'
results['samtools'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'trishulagenetics/genocan-software-versions'
section_name: 'trishulagenetics/genocan Software Versions'
section_href: 'https://github.com/trishulagenetics/genocan'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")