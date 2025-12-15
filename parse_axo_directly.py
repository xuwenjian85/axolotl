import sys
import os
import pandas as pd
import numpy as np
import yaml
import subprocess
from general import postprocess

def main():
    # Receive input parameters
    ctsfile = sys.argv[1]
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)
    
    # Read expression matrix
    cts = pd.read_csv(ctsfile, sep='\t', index_col='Gene')
    print(f'Input expression matrix shape: {cts.shape}')
    
    # Prepare output file paths
    base_name = os.path.splitext(os.path.basename(ctsfile))[0]
    outrider_out = os.path.join(output_dir, f'{base_name}_outrider.txt.gz')
    outsingle_out = os.path.join(output_dir, f'{base_name}_outsingle.txt.gz')
    mymethod_out = os.path.join(output_dir, f'{base_name}_axo.txt.gz')
    
    # Run OUTRIDER
    run_outrider(cts, outrider_out, output_dir)
    
    # Run OUTSINGLE
    run_outsingle(cts, outsingle_out, output_dir)
    
    # Run custom method
    run_mymethod(outsingle_out, outrider_out, mymethod_out)
    
    print(f'All results have been saved to: {output_dir}')

def run_outrider(cts, outrider_out, work_dir):
    print("Starting to run OUTRIDER...")
    
    # Calculate dim_q0 parameter
    Nsample = cts.shape[1]
    div_q0 = 3
    dim_q0 = int(Nsample / div_q0)
    print(f'OUTRIDER parameter: dim_q0 = {dim_q0}')
    
    # Save temporary expression file to working directory
    cts_workfile = os.path.join(work_dir, 'temp_cts_outrider.txt')
    cts.to_csv(cts_workfile, sep='\t')
    
    # Read configuration file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(current_dir, 'config.yaml'), 'r') as f:
        config = yaml.safe_load(f)
    
    # Prepare commands
    conda_path = config['paths']['conda']
    script_path = config['paths']['script1']
    script_name = config['paths']['outrider_script'] #'outrider_pred_dim.R'
    script = os.path.join(script_path, script_name)
    
    # Output file paths
    norm_cts_file = os.path.join(work_dir, 'temp_otrd_norm_cts.txt')
    result_file = os.path.join(work_dir, 'temp_otrd.txt')
    
    conda_cmds = [
        f'source {conda_path}',
        f'conda activate {config["paths"]["outrider_env"]}',
        f'Rscript {script} {cts_workfile} {dim_q0} {norm_cts_file} {result_file}',
        f'gzip -c {result_file} > {outrider_out}',
    ]
    
    # Execute commands
    subprocess.run(['bash', '-c', '; '.join(conda_cmds)], check=True)
    
    # Post-processing
    postprocess('outrider', outrider_out, cts_workfile)
    
    # Clean up intermediate files
    for f in [cts_workfile, norm_cts_file, result_file]:
        if os.path.exists(f):
            os.remove(f)
    
    print("OUTRIDER run completed")


def run_outsingle(cts, outsingle_out, work_dir):
    print("Starting to run OUTSINGLE...")
    
    # 1. Prepare temporary files and paths (corresponding to the original outsingle.sh logic)
    cts_workfile = os.path.join(work_dir, 'temp_cts_outsingle.txt')
    cts.to_csv(cts_workfile, sep='\t')  # Save expression file
    
    # 2. Parse paths and parameters from the original shell script
    out_osg = f"{cts_workfile}.osg.txt"  # Intermediate output file (uncompressed)
    
    # 3. Get user ID (replace the id command in the original shell script)
    # 4. Read configuration file to get environment and script paths
    current_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(current_dir, 'config.yaml'), 'r') as f:
        config = yaml.safe_load(f)
    
    # 5. Build execution commands (integrate Python calls and compression steps in the original shell script)
    conda_path = config['paths']['conda']
    script_path = config['paths']['script1']
    script_name = config['paths']['outsingle_script'] 
    script = os.path.join(script_path, script_name)
        
    outsingle_cmds = [
        f'source {conda_path}',
        f'conda activate {config["paths"]["outsingle_env"]}',  # Use configured environment
        # Call outsingle.py (path from configuration file)
        f'python {script} '
        f'{cts_workfile} '
        f'{out_osg} '
        f'{config["paths"]["outsingle_dir"]}',  # Third parameter: outsingle directory
        # Compress output file (corresponding to the original gzip command)
        f'gzip -c {out_osg} > {outsingle_out}'
    ]
    
    # 6. Execute commands
    subprocess.run(['bash', '-c', '; '.join(outsingle_cmds)], check=True)
    
    # 7. Post-processing
    postprocess('outsingle', outsingle_out, cts_workfile)
    
    # 8. Clean up intermediate files (including temporary expression file and intermediate output file)
    for f in [cts_workfile, out_osg]:
        if os.path.exists(f):
            os.remove(f)
    
    print("OUTSINGLE run completed")
    
def run_mymethod(outsingle_out, outrider_out, mymethod_out):
    print("Starting to run my method...")
    
    # Read configuration file to get environment information
    current_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(current_dir, 'config.yaml'), 'r') as f:
        config = yaml.safe_load(f)
    
    # Directly execute the logic in the original mymethod.sh (example command, need to adjust according to actual script content)
    conda_path = config['paths']['conda']
    script_path = config['paths']['script1']
    script_name = config['paths']['mymethod_script'] 
    script = os.path.join(script_path, script_name)    
    mymethod_cmds = [
        f'source {conda_path}',
        f'conda activate {config["paths"]["mymethod_env"]}',  # Assume the environment is configured
        # The following is the core command in the original mymethod.sh, modify according to actual content
        f'python {script}  {outsingle_out} {outrider_out} {mymethod_out}',  # Example Python script call
    ]
    
    subprocess.run(['bash', '-c', '; '.join(mymethod_cmds)], check=True)
    print("my method run completed")

if __name__ == "__main__":
    main()