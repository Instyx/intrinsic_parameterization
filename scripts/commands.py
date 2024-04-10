import jupyter_client
from pathlib import Path
import sys


if not Path("kernel.txt").exists():
    print("No running Kernel")
    sys.exit(1)

open("kernel.txt", "r")
connection_file = jupyter_client.find_connection_file('/run/user/1000/jupyter/kernel-630afc98-481e-45b4-8eb1-6536827d28f2.json')

client = jupyter_client.BlockingKernelClient()
client.load_connection_file(connection_file)
client.start_channels()

x = client.execute_interactive("test_var")
