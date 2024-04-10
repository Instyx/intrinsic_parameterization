import subprocess
import sys

if len(sys.argv) > 1:
    kernelname = sys.argv[1]
else:
    kernelname = "kernel-manager.json"

    p = subprocess.Popen([sys.executable, "-c", f"from ipykernel.kernelapp import launch_new_instance;launch_new_instance(extra_arguments=[], connection_file='{kernelname}')"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

import jupyter_client, time

connection_file = jupyter_client.find_connection_file(f'{kernelname}')
client = jupyter_client.BlockingKernelClient()
client.load_connection_file(connection_file)
client.start_channels()
# print("lol")
x = client.execute_interactive(code="import taskMgr")
x = client.execute_interactive(code="taskMgr.start()")# , output_hook=msg_handler)

while True:
    try:
        time.sleep(1)
    except KeyboardInterrupt:
        print('\rShout down')
        break
        # Check for running Tasks

p.terminate()
