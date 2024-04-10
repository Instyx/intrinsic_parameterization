import os, glob, sys, math
import asyncio

tasks = []
done = []
MAX_JOBS = 5
running = []
stopped = False

async def run(cmd):
    # print("Run command:", cmd)
    proc = await asyncio.create_subprocess_shell(cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE)
    stdout, stderr = await proc.communicate()
    # print("Done command:", cmd)
    return stdout.decode("utf-8") if stdout is not None else None, stderr.decode("utf-8") if stderr is not None else None

def start_job(name, cmd):
    async def job():
        global done, running
        r = run(cmd)
        running.append(r)
        done.append({name: { "cmd": cmd, "result": await r }})
        running.remove(r)
    # print("Start job:", name)
    asyncio.ensure_future(job())

async def manager():
    global stopped
    stopped = False
    while not stopped:
        if len(tasks) > 0 and len(running) < MAX_JOBS:
            # print("Task gefunden:")
            start_job(*tasks.pop())
        await asyncio.sleep(0.1)
    print("Manager Stopped")

def stop():
    global stopped
    stopped = True

def start():
    asyncio.ensure_future(manager())

# %%

# import glob
# from pathlib import Path
# import numpy as np
# import time
#
# update_display("Test", display_id=10)
#
# for i in range(10):
#     time.sleep(1)
#     print("lol")
#
# len(done)
# len(tasks)
#
# done = []
#
# output = "../output"
# methods = ["dirichlet", "asap", "symdirichlet", "arap"]
# for i in filtered:
#     p = f"../meshes/cutted/{i}"
#     if not i.endswith("cutted.off"):
#         p = f"../meshes/thingi_data/{i}"
#
#     for m in methods:
#         cmd = f"../build/apps/cli {m} {p} {output} nodelaunay noprio"
#         tasks.insert(0,[f"{Path(i).name}", cmd])


# # for i in sorted(glob.glob("../meshes/thingi_data/*.off"), key = os.path.getsize):
# #     tasks.insert(0,[f"{Path(i).name}", f"../build/apps/filter {i}"])
# #
# # for i in sorted(glob.glob("../meshes/cutted/*.off"), key = os.path.getsize):
# #     tasks.insert(0,[f"{Path(i).name}", f"../build/apps/filter {i}"])
# #
# #
# # done
# #
# # MAX_JOBS = 5
# #
# # tasks = []
# # filtered = []
# # for i in done:
# #     if (i[list(i)[0]]["result"][0].endswith("Passed\n")):
# #         filtered.append(list(i)[0])
# # len(filtered)
# # open("filtered.txt", "w").write("\n".join(filtered))
