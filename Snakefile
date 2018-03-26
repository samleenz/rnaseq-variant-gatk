# 2018-03-26
# Sam Lee

configfile: "config.yaml"
from os.path import join
import re
#


rule all:
    input:
      ...
    output:
      ...
    conda:
      "env/AAA.yaml"
    log:
      "logs/AAA/{BBB}.log"
    shell:
      """

      """