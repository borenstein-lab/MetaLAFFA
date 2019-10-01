#!/bin/sh
# properties = {{properties}}
cd {MetaLAFFA_directory}
export PYTHONPATH=${PYTHONPATH}:$({python} -c "from config import env; print(env['PYTHONPATH'])")
{{exec_job}}
