### Running Instructions
#### Directory
`testdata` is the directory for storing **small-scale test data**.

#### Script Running Commands
##### 1. Regenerate Results: Benchmark Test for All Methods
This command runs a benchmark test on all methods based on the task configuration file, and it is mainly used to regenerate test results.
```bash
# Run benchmark on all methods (based on the task configuration file, for regenerating results)
bash run.sh testdata/task_config/t00_TOY_s128_g1000.config
```

##### 2. For End Users: Run Axolotl Only
This command is a simplified version for end users, which only executes the logic related to Axolotl.
```bash
# For end users: Run axolotl only
bash axo.sh testdata/toy_small.tsv output
```

### Summary
1.  The `testdata` directory is used to store small-scale test data, which serves as the data source basis for script operation.
2.  The two types of script commands have different purposes: `run.sh` is used for the benchmark test of all methods and result regeneration, while `axo.sh` is a command for end users to run Axolotl independently.
3.  The Markdown format improves the hierarchy and readability of the command instructions significantly through hierarchical headings, code blocks, and lists.
