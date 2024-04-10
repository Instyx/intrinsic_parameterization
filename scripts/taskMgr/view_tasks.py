import jupyter_client
import os


connection_file = jupyter_client.find_connection_file('kernel-manager.json')
client = jupyter_client.BlockingKernelClient()
client.load_connection_file(connection_file)

import sys, time, json

LINE_UP = '\033[1A'
LINE_CLEAR = '\x1b[2K'

def callback(msg):
    if (msg["msg_type"] == "execute_result"):
        res = json.loads(eval(msg['content']["data"]['text/plain']))
        s = ""
        for r in res:
            key = list(r)[0]
            s += f'{key}: {r[key]["result"][0]}, {r[key]["result"][1]}\n'
        size = os.get_terminal_size()[0]
        print(s)
        for i in range(s.count("\n")+1):
            print(LINE_UP, end=LINE_CLEAR)

while True:
    try:
        res = client.execute_interactive(code=f"import json; json.dumps(taskMgr.done)", output_hook=callback)
        time.sleep(0.5)
    except KeyboardInterrupt:
        print('\rShout down')
        break
        # Check for running Tasks
