import os, glob, sys, math
import asyncio

tasks = []
done = []
MAX_JOBS = 5
running = []
stopped = False

async def run(cmd):
    print("Run command:", cmd)
    proc = await asyncio.create_subprocess_shell(cmd, stdout=asyncio.subprocess.PIPE, stderr=asyncio.subprocess.PIPE)
    stdout, stderr = await proc.communicate()
    print("Done command:", cmd)
    return stdout.decode("utf-8") if stdout is not None else None, stderr.decode("utf-8") if stderr is not None else None

def start_job(name, cmd):
    async def job():
        global done, running
        r = run(cmd)
        running.append(r)
        done.append({name: { "cmd": cmd, "result": await r }})
        running.remove(r)
    print("Start job:", name)
    asyncio.ensure_future(job())

async def manager():
    global stopped
    stopped = False
    while not stopped:
        if len(tasks) > 0 and len(running) < MAX_JOBS:
            print("Task gefunden:")
            start_job(*tasks.pop())
        await asyncio.sleep(0.1)
    print("Manager Stopped")

def stop():
    global stopped
    stopped = True

asyncio.ensure_future(manager())

for i in range(20):
    tasks.insert(0,[f"{i}.off", "lf"])

    tasks = [] + tasks

tasks = []
done
stop()

running
