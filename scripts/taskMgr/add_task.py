import jupyter_client

connection_file = jupyter_client.find_connection_file('kernel-manager.json')
client = jupyter_client.BlockingKernelClient()
client.load_connection_file(connection_file)

import sys

if len(sys.argv) > 2:
    res = client.execute_interactive(code=f"taskMgr.tasks.insert(0,['{sys.argv[1]}','{sys.argv[2]}'])")
else:
    print("No Task defined")
