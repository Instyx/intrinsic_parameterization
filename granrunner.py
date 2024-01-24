import os, glob, sys, math

methods = ["symdirichlet", "arap"]
grans = [5,10,50]
start = int(sys.argv[1]) if len(sys.argv) > 1 else 0
step = int(sys.argv[2]) if len(sys.argv) > 2 else 1
end = int(sys.argv[3]) if len(sys.argv) > 3 else None
folder = "thingi_data"
output = "gran_output"

num = math.floor((len(glob.glob(f"{folder}/*.off")[:end])-1-start)/step)

for k,i in enumerate(sorted(glob.glob(f"{folder}/*.off"), key = os.path.getsize)[start:end:step]):
    print("------------------------------------------------------------------------------------")
    print("------------------------------------------------------------------------------------")
    print(k,"of", num)
    print("------------------------------------------------------------------------------------")
    print("------------------------------------------------------------------------------------")
    for m in methods:
        for flip_gran in grans:
            cmd = f"build/cli {m} {i} {output} {flip_gran}"
            print(cmd)
            try:
                sts = os.system(cmd)
                if(sts == 2):
                    sys.exit(0)
            except Exception as e:
                print("Failed method", m, "for mesh", i)
            print()
        print()
        print()
