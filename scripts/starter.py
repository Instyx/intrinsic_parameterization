# import IPython
# IPython.start_ipython(config={'IPKernelApp': {'connection_file': '/run/user/1000/jupyter/kernel-630afc98-481e-45b4-8eb1-6536827d28f2.json'}})
%connect_info

import jupyter_client
connection_file = jupyter_client.find_connection_file('/run/user/1000/jupyter/kernel-630afc98-481e-45b4-8eb1-6536827d28f2.json')

client = jupyter_client.BlockingKernelClient()
client.load_connection_file(connection_file)
# client.start_channels()

 client.execute_interactive(code="test_var", output_hook=lambda msg: print(msg))
# print(x)
# reply = client.get_shell_msg(x["msg_id"])
# print(reply)
msg_id = client.execute("test_var")
# print(msg_id)
reply = client.get_shell_msg(msg_id)
# print(reply)

connection_file = jupyter_client.find_connection_file('kernel-manager.json')
from IPython import embed
embed()

dir(get_ipython().kernel).connection_file

import IPython
IPython.start_ipython()
