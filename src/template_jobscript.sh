#!/bin/sh
# properties = {{properties}}
cd {metaLAFFA_directory}
export PYTHONPATH=${{PYTHONPATH}}:$({python} -c "from config import env; print(env['PYTHONPATH'])")
{{exec_job}}
