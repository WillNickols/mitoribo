# To run: $ python download_files.py -i to_download.txt -o OUTPUT_FOLDER
# to_download.txt should be one SRA per line with no header

from anadama2 import Workflow
import pandas as pd
import glob as glob
import os

workflow = Workflow()
args = workflow.parse_args()

output = os.path.abspath(args.output.rstrip("/")) + "/"
if not os.path.isdir(output):
	os.makedirs(output)

df = pd.read_csv(str(args.input), header = None)

files = (df.iloc[:, 0]).tolist()

for file in files:
    paths = glob.glob(output + file + '*.fastq*')
    if len(paths) == 0:
        workflow.add_task("fastq-dump " + file + " --gzip -O " + output)

workflow.go()